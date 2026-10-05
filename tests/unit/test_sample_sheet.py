import pytest
from conftest import DATA_DIR, FIVEPRIME, TARGET, make_row

from pacybara_workflow.sample_sheet import SampleSheetError, read_sample_sheet


def test_reads_row_and_applies_defaults(write_sheet, defaults):
    sheet = write_sheet([make_row()])
    [library] = read_sample_sheet(sheet, defaults, data_dir="/data")

    assert library.name == "GENE1_abundance_na"
    assert library.title("20260101") == "20260101_GENE1_abundance_na"
    assert library.fastq == "/data/gene1.fastq.gz"
    assert library.parameters == {
        "barcode": "NNNNNNNNNNNNNNNN",
        "maxqdrops": "4",
        "minbcq": "32",
        "minjaccard": "0.2",
        "minmatches": "1",
        "maxdiff": "1",
        "minqual": "27",
        "clustermode": "uptag",
    }


def test_orf_coordinates_are_one_based_and_inclusive(write_sheet, defaults):
    [library] = read_sample_sheet(write_sheet([make_row()]), defaults, "/data")
    fiveprime_len = len(FIVEPRIME) + 20

    assert library.orfstart == fiveprime_len + 1
    assert library.orfend == fiveprime_len + len(TARGET)
    assert library.amplicon[library.orfstart - 1 : library.orfend] == TARGET


def test_row_overrides_win_and_blank_cells_fall_back(write_sheet, defaults):
    rows = [
        make_row(gene="A", fastq="a.fastq.gz", minqual="67", maxdiff=""),
        make_row(gene="B", fastq="b.fastq.gz", minqual="", maxdiff="2"),
    ]
    a, b = read_sample_sheet(write_sheet(rows), defaults, "/data")

    assert (a.parameters["minqual"], a.parameters["maxdiff"]) == ("67", "1")
    assert (b.parameters["minqual"], b.parameters["maxdiff"]) == ("27", "2")


def test_absolute_fastq_ignores_data_dir(write_sheet, defaults):
    sheet = write_sheet([make_row(fastq="/elsewhere/x/../g.fastq.gz")])
    [library] = read_sample_sheet(sheet, defaults, "/data")
    assert library.fastq == "/elsewhere/g.fastq.gz"


def test_pacybara_prep_sheet_with_bom(tmp_path, defaults):
    # pacybara_prep's own sheet (UTF-8 BOM, no fastq column) plus a fastq
    # column: the BOM must not leak into the first column's name.
    lines = (
        (DATA_DIR / "pacybara_prep_setup.csv").read_text(encoding="utf-8").splitlines()
    )
    assert lines[0].startswith("﻿")
    lines = [lines[0] + ",fastq"] + [
        f"{line},{i}.fastq.gz" for i, line in enumerate(lines[1:])
    ]
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("\n".join(lines), encoding="utf-8")

    libraries = read_sample_sheet(sheet, defaults, "/data")
    assert [lib.name for lib in libraries] == [
        "SOS1_abundance_na",
        "KSR1_abundance_nterm",
        "ARAF_activity_nterm",
    ]


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"gene": ""}, "'gene' is empty"),
        ({"gene": "has space"}, "'gene' may only contain"),
        ({"target": "ATGXYZ"}, "'target' may only contain A, C, G, T or N"),
        ({"fastq": "reads.fastq"}, "must be a .fastq.gz file"),
        ({"maxqdrops": "four"}, "maxqdrops must be a non-negative integer"),
        ({"minjaccard": "1"}, "minjaccard must be written as 0.<digits>"),
        ({"clustermode": "up"}, "clustermode must be one of"),
        ({"barcode": "SWSWSWSW"}, "found 0"),
        ({"clustermode": "downtag"}, "only one barcode"),
    ],
)
def test_invalid_rows(write_sheet, defaults, overrides, message):
    with pytest.raises(SampleSheetError, match=message):
        read_sample_sheet(write_sheet([make_row(**overrides)]), defaults, "/data")


def test_three_barcodes_rejected(write_sheet, defaults):
    row = make_row(threeprime="N" * 32 + "A")
    with pytest.raises(SampleSheetError, match="found 3"):
        read_sample_sheet(write_sheet([row]), defaults, "/data")


def test_relative_fastq_needs_data_dir(write_sheet, defaults):
    with pytest.raises(SampleSheetError, match="no data_dir"):
        read_sample_sheet(write_sheet([make_row()]), defaults)


def test_duplicates_rejected(write_sheet, defaults):
    rows = [make_row(), make_row(fastq="other.fastq.gz"), make_row(gene="G2")]
    with pytest.raises(SampleSheetError) as excinfo:
        read_sample_sheet(write_sheet(rows), defaults, "/data")
    assert excinfo.value.errors == [
        "row 3 (GENE1): duplicates row 2: gene/assay/fragment 'GENE1_abundance_na' "
        "must be unique (it names the output directory)",
        "row 4 (G2): duplicates row 2: fastq '/data/gene1.fastq.gz' is already used",
    ]


def test_missing_and_unknown_columns(write_sheet, defaults):
    row = make_row(extra="x")
    del row["fastq"]
    with pytest.raises(SampleSheetError) as excinfo:
        read_sample_sheet(write_sheet([row]), defaults, "/data")
    assert excinfo.value.errors == [
        "missing required column(s): fastq",
        "unknown column(s): extra",
    ]


def test_incomplete_defaults(write_sheet, defaults):
    del defaults["minqual"]
    with pytest.raises(SampleSheetError, match="pacybara_defaults is missing: minqual"):
        read_sample_sheet(write_sheet([make_row()]), defaults, "/data")


def test_all_errors_reported_together(write_sheet, defaults):
    rows = [make_row(gene="A", fastq="a.fq"), make_row(gene="B", clustermode="x")]
    with pytest.raises(SampleSheetError) as excinfo:
        read_sample_sheet(write_sheet(rows), defaults, "/data")
    assert len(excinfo.value.errors) == 2
