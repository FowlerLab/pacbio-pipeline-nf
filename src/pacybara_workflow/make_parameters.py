"""`pacybara-make-params`: sample sheet -> one Pacybara parameter file per row.

Writes `<TITLE>.txt` for every library plus `manifest.csv`
(`title,fastq,param_file`), which MAKE_PARAMETER_FILES hands to Nextflow to
fan out one RUN_PACYBARA task per row. Validation failures print every
problem and exit 1 before anything is written.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import re
import sys
from pathlib import Path

from pacybara_workflow.parameter_file import render_parameter_file
from pacybara_workflow.sample_sheet import SampleSheetError, read_sample_sheet

MANIFEST_NAME = "manifest.csv"
_RUN_DATE_RX = re.compile(r"^[0-9]{8}$")


def write_parameter_files(
    sample_sheet: str | os.PathLike,
    output_dir: str | os.PathLike,
    defaults: dict[str, object],
    run_date: str,
    data_dir: str | os.PathLike | None = None,
) -> Path:
    """Validate `sample_sheet` and write its parameter files and manifest.

    Returns the manifest's path. `param_file` entries in the manifest are
    absolute, so the manifest can be read from anywhere.
    """
    if not _RUN_DATE_RX.match(run_date):
        raise SampleSheetError([f"run_date must be YYYYMMDD, got '{run_date}'"])
    libraries = read_sample_sheet(sample_sheet, defaults, data_dir)

    output_dir = Path(output_dir).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / MANIFEST_NAME
    with open(manifest_path, "w", newline="") as manifest:
        writer = csv.writer(manifest)
        writer.writerow(["title", "fastq", "param_file"])
        for library in libraries:
            title = library.title(run_date)
            param_file = output_dir / f"{title}.txt"
            param_file.write_text(render_parameter_file(library, run_date))
            writer.writerow([title, library.fastq, str(param_file)])
            print(f"Generated: {param_file.name}")
    return manifest_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--sample-sheet", required=True, help="CSV, one library per row"
    )
    parser.add_argument(
        "--defaults",
        required=True,
        help="JSON object of Pacybara arguments for blank cells (params.yaml's pacybara_defaults)",
    )
    parser.add_argument(
        "--data-dir", help="base directory for relative `fastq` entries"
    )
    parser.add_argument(
        "--run-date",
        help="YYYYMMDD prefix for every TITLE (default: today, as pacybara_prep does)",
    )
    parser.add_argument("--output-dir", default=".")
    args = parser.parse_args(argv)

    run_date = args.run_date or datetime.date.today().strftime("%Y%m%d")
    try:
        write_parameter_files(
            args.sample_sheet,
            args.output_dir,
            json.loads(args.defaults),
            run_date,
            args.data_dir,
        )
    except SampleSheetError as e:
        print(e, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
