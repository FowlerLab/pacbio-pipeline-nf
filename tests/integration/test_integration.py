"""End-to-end `nextflow run -profile local` against a stub pacybara_simplex.sh.

The stub (stubs/pacybara_simplex.sh) parses the parameter file like the
real script and writes a fake clusters_transl_softfilter.csv.gz, so these
tests exercise everything this repo owns -- sample sheet validation,
parameter files, the per-library fan-out, publishing and archiving --
without Pacybara's toolchain. Skipped when nextflow isn't on PATH.
"""

import csv
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
STUBS = Path(__file__).parent / "stubs"

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None, reason="nextflow not on PATH"
)

AMPLICON_PARTS = {
    "fiveprime": "ACGTACGTAC" + "N" * 16 + "GTACGTACGT",
    "target": "ATGGCCATGGCCATGGCC",
    "threeprime": "TTAATTAA",
}


@pytest.fixture
def layout(tmp_path):
    """A miniature Pacbio_processing/ tree."""
    root = tmp_path / "Pacbio_processing"
    dirs = {
        "data": root / "to_process" / "data",
        "output": root / "pacybara_output",
        "completed": root / "completed_fastqs",
        "run": tmp_path / "run",
    }
    for d in dirs.values():
        d.mkdir(parents=True)
    return dirs


def write_fastq(path):
    with gzip.open(path, "wt") as f:
        f.write("@read1\nACGT\n+\nIIII\n")


def write_sheet(path, rows):
    columns = [
        "gene",
        "assay",
        "fragment",
        "fiveprime",
        "target",
        "threeprime",
        "fastq",
    ]
    columns += [c for row in rows for c in row if c not in columns]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {**AMPLICON_PARTS, "assay": "abundance", "fragment": "na", **row}
            )
    return path


def run_pipeline(layout, sheet, *extra):
    env = dict(os.environ)
    # The stub, then this repo's venv (pacybara-make-params).
    env["PATH"] = os.pathsep.join(
        [str(STUBS), str(Path(sys.executable).parent), env["PATH"]]
    )
    cmd = [
        "nextflow",
        "run",
        str(REPO),
        "-profile",
        "local",
        "-params-file",
        str(REPO / "params.yaml"),
        "--sample_sheet",
        str(sheet),
        "--data_dir",
        str(layout["data"]),
        "--output_dir",
        str(layout["output"]),
        "--run_date",
        "20260101",
        *extra,
    ]
    return subprocess.run(
        cmd, cwd=layout["run"], env=env, capture_output=True, text=True
    )


def test_runs_every_library_without_archiving(layout, tmp_path):
    write_fastq(layout["data"] / "gene1.fastq.gz")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    write_fastq(elsewhere / "gene2.fastq.gz")
    sheet = write_sheet(
        tmp_path / "sheet.csv",
        [
            {"gene": "GENE1", "fastq": "gene1.fastq.gz"},
            {
                "gene": "GENE2",
                "fastq": str(elsewhere / "gene2.fastq.gz"),
                "minqual": "67",
            },
        ],
    )

    result = run_pipeline(layout, sheet)
    assert result.returncode == 0, result.stdout + result.stderr

    for gene, fastq, minqual in [
        ("GENE1", layout["data"] / "gene1.fastq.gz", "27"),
        ("GENE2", elsewhere / "gene2.fastq.gz", "67"),
    ]:
        title = f"20260101_{gene}_abundance_na"
        workspace = layout["output"] / title
        prefix = fastq.name.removesuffix(".fastq.gz")
        assert (
            workspace / f"{prefix}_clustering" / "clusters_transl_softfilter.csv.gz"
        ).exists()

        # The published parameter file names the real fastq ...
        published = (workspace / f"{title}.txt").read_text()
        assert f"\nINFASTQ={fastq}\n" in published
        assert f"\nWORKSPACE={title}\n" in published
        assert f"\nMINQUAL={minqual}\n" in published
        # ... while Pacybara ran on the staged copy, in a TITLE-named workspace.
        ran = (workspace / "stub_saw_parameters.txt").read_text()
        assert f"\nINFASTQ={fastq.name}\n" in ran
        assert f"\nWORKSPACE={title}\n" in ran

        assert fastq.exists()  # archive_inputs defaults to false
    assert not any(layout["completed"].iterdir())


def test_archives_only_successful_libraries(layout, tmp_path):
    write_fastq(layout["data"] / "good.fastq.gz")
    write_fastq(layout["data"] / "will_fail.fastq.gz")
    sheet = write_sheet(
        tmp_path / "sheet.csv",
        [
            {"gene": "GOOD", "fastq": "good.fastq.gz"},
            {"gene": "BAD", "fastq": "will_fail.fastq.gz"},
        ],
    )

    result = run_pipeline(
        layout,
        sheet,
        "--archive_inputs",
        "true",
        "--completed_fastq_dir",
        str(layout["completed"]),
    )
    # errorStrategy 'ignore': the failed library doesn't fail the run.
    assert result.returncode == 0, result.stdout + result.stderr

    assert (layout["output"] / "20260101_GOOD_abundance_na").is_dir()
    assert not (layout["output"] / "20260101_BAD_abundance_na").exists()
    assert sorted(p.name for p in layout["completed"].iterdir()) == ["good.fastq.gz"]
    assert sorted(p.name for p in layout["data"].iterdir()) == ["will_fail.fastq.gz"]


def test_invalid_sheet_fails_before_pacybara(layout, tmp_path):
    write_fastq(layout["data"] / "gene1.fastq.gz")
    sheet = write_sheet(
        tmp_path / "sheet.csv",
        [{"gene": "GENE1", "fastq": "gene1.fastq.gz", "clustermode": "sideways"}],
    )

    result = run_pipeline(layout, sheet)
    assert result.returncode != 0
    assert "clustermode must be one of" in result.stdout + result.stderr
    assert not any(layout["output"].iterdir())


def test_missing_fastq_fails_before_pacybara(layout, tmp_path):
    sheet = write_sheet(
        tmp_path / "sheet.csv", [{"gene": "GENE1", "fastq": "absent.fastq.gz"}]
    )

    result = run_pipeline(layout, sheet)
    assert result.returncode != 0
    assert "absent.fastq.gz" in result.stdout + result.stderr
    assert not any(layout["output"].iterdir())
