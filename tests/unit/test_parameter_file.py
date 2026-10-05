import csv
import json
import re

import pytest
from conftest import DATA_DIR, DEFAULTS, make_row

from pacybara_workflow.make_parameters import main, write_parameter_files
from pacybara_workflow.parameter_file import render_parameter_file
from pacybara_workflow.sample_sheet import SampleSheetError, read_sample_sheet

GOLDEN_DIR = DATA_DIR / "golden"


def test_matches_pacybara_prep_output(tmp_path, defaults):
    # golden/ is pacybara_prep's make_parameter.py run on
    # pacybara_prep_setup.csv, dated 20260101. Only INFASTQ (a "fastq.gz"
    # placeholder there) should differ.
    lines = (
        (DATA_DIR / "pacybara_prep_setup.csv").read_text(encoding="utf-8").splitlines()
    )
    lines = [lines[0] + ",fastq"] + [
        f"{line},{i}.fastq.gz" for i, line in enumerate(lines[1:])
    ]
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("\n".join(lines), encoding="utf-8")

    libraries = read_sample_sheet(sheet, defaults, "/data")
    golden_files = sorted(GOLDEN_DIR.glob("*.txt"))
    assert len(golden_files) == len(libraries)
    for i, library in enumerate(libraries):
        expected = (GOLDEN_DIR / f"{library.title('20260101')}.txt").read_text()
        expected = expected.replace(
            "INFASTQ=fastq.gz\n", f"INFASTQ=/data/{i}.fastq.gz\n"
        )
        assert render_parameter_file(library, "20260101") == expected


def test_overrides_reach_the_file(write_sheet, defaults):
    sheet = write_sheet([make_row(minqual="67", clustermode="virtual")])
    [library] = read_sample_sheet(sheet, defaults, "/data")
    text = render_parameter_file(library, "20260101")
    assert "\nMINQUAL=67\n" in text
    assert "\nCLUSTERMODE=virtual\n" in text


def test_write_parameter_files_and_manifest(tmp_path, write_sheet, defaults):
    rows = [
        make_row(gene="A", fastq="a.fastq.gz"),
        make_row(gene="B", fastq="/abs/b.fastq.gz"),
    ]
    out = tmp_path / "out"
    manifest = write_parameter_files(
        write_sheet(rows), out, defaults, "20260101", "/data"
    )

    with open(manifest) as f:
        entries = list(csv.DictReader(f))
    assert entries == [
        {
            "title": "20260101_A_abundance_na",
            "fastq": "/data/a.fastq.gz",
            "param_file": str(out / "20260101_A_abundance_na.txt"),
        },
        {
            "title": "20260101_B_abundance_na",
            "fastq": "/abs/b.fastq.gz",
            "param_file": str(out / "20260101_B_abundance_na.txt"),
        },
    ]
    assert (
        "INFASTQ=/data/a.fastq.gz\n"
        in (out / "20260101_A_abundance_na.txt").read_text()
    )


def test_bad_run_date(tmp_path, write_sheet, defaults):
    with pytest.raises(SampleSheetError, match="YYYYMMDD"):
        write_parameter_files(
            write_sheet([make_row()]), tmp_path, defaults, "2026-01-01", "/d"
        )


def test_cli_defaults_to_today(tmp_path, write_sheet):
    sheet = write_sheet([make_row()])
    out = tmp_path / "out"
    argv = ["--sample-sheet", str(sheet), "--defaults", json.dumps(DEFAULTS)]
    argv += ["--data-dir", "/data", "--output-dir", str(out)]
    assert main(argv) == 0
    [param_file] = out.glob("*.txt")
    assert re.fullmatch(r"[0-9]{8}_GENE1_abundance_na\.txt", param_file.name)


def test_cli_reports_errors_and_writes_nothing(tmp_path, write_sheet, capsys):
    sheet = write_sheet([make_row(clustermode="bogus")])
    out = tmp_path / "out"
    argv = ["--sample-sheet", str(sheet), "--defaults", json.dumps(DEFAULTS)]
    argv += ["--data-dir", "/data", "--output-dir", str(out)]
    assert main(argv) == 1
    assert "clustermode must be one of" in capsys.readouterr().err
    assert not out.exists()
