import csv
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"

# params.yaml's pacybara_defaults -- pacybara_prep's hard-coded values.
DEFAULTS = {
    "barcode": "NNNNNNNNNNNNNNNN",
    "maxqdrops": 4,
    "minbcq": 32,
    "minjaccard": 0.2,
    "minmatches": 1,
    "maxdiff": 1,
    "minqual": 27,
    "clustermode": "uptag",
}

FIVEPRIME = "acgt" * 5
BARCODE_REGION = "AC" + "N" * 16 + "GT"
TARGET = "ATGGCC" * 4
THREEPRIME = "TTAA" * 3


@pytest.fixture
def defaults():
    return dict(DEFAULTS)


@pytest.fixture
def write_sheet(tmp_path):
    """Write rows (dicts) to a CSV, with columns in first-seen order."""

    def _write(rows, name="sheet.csv", columns=None):
        if columns is None:
            columns = []
            for row in rows:
                columns += [c for c in row if c not in columns]
        path = tmp_path / name
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        return path

    return _write


def make_row(**overrides):
    row = {
        "gene": "GENE1",
        "assay": "abundance",
        "fragment": "na",
        "fiveprime": FIVEPRIME + BARCODE_REGION,
        "target": TARGET,
        "threeprime": THREEPRIME,
        "fastq": "gene1.fastq.gz",
    }
    row.update(overrides)
    return row
