"""Read and validate a Pacybara sample sheet.

The sheet is `pacybara_prep`'s (github.com/dlholmes2117/pacybara_prep) CSV --
`gene,assay,fragment,fiveprime,target,threeprime`, one library per row --
plus a required `fastq` column naming that library's `.fastq.gz`, and
optional per-row overrides for any Pacybara argument in `PARAMETER_KEYS`.

Every problem in the sheet is collected and raised together as one
`SampleSheetError`, so a bad sheet fails the run before any Pacybara job
starts, and with every fix needed listed at once.
"""

from __future__ import annotations

import csv
import dataclasses
import os
import re
from collections.abc import Mapping

REQUIRED_COLUMNS = (
    "gene",
    "assay",
    "fragment",
    "fiveprime",
    "target",
    "threeprime",
    "fastq",
)

# Pacybara's ARGUMENTS-section keys (lowercased) that a row may override and
# params.yaml's `pacybara_defaults` must supply. Order is the parameter
# file's.
PARAMETER_KEYS = (
    "barcode",
    "maxqdrops",
    "minbcq",
    "minjaccard",
    "minmatches",
    "maxdiff",
    "minqual",
    "clustermode",
)

# The same checks pacybara_simplex.sh's validateInteger/validateFloat and its
# CLUSTERMODE case apply -- run here so a typo fails up front, not hours into
# a cluster job.
_INTEGER_RX = re.compile(r"^[0-9]+$")
_FLOAT_RX = re.compile(r"^0\.[0-9]+$")
_INTEGER_KEYS = ("maxqdrops", "minbcq", "minmatches", "maxdiff", "minqual")
CLUSTER_MODES = ("uptag", "downtag", "virtual")

_SEQUENCE_RX = re.compile(r"^[ACGTNacgtn]+$")
# TITLE becomes a directory and file name; keep each component to characters
# that are safe unquoted in Pacybara's bash.
_NAME_RX = re.compile(r"^[A-Za-z0-9._-]+$")


class SampleSheetError(ValueError):
    """A sample sheet failed validation. `errors` lists every problem found."""

    def __init__(self, errors: list[str]):
        self.errors = errors
        super().__init__(
            "Invalid sample sheet:\n" + "\n".join(f"  - {e}" for e in errors)
        )


@dataclasses.dataclass(frozen=True)
class Library:
    """One sample sheet row: a single Pacybara run."""

    gene: str
    assay: str
    fragment: str
    fiveprime: str
    target: str
    threeprime: str
    fastq: str
    parameters: Mapping[str, str]

    @property
    def name(self) -> str:
        """`<gene>_<assay>_<fragment>`: the amplicon's FASTA header."""
        return f"{self.gene}_{self.assay}_{self.fragment}"

    def title(self, run_date: str) -> str:
        """`<run_date>_<gene>_<assay>_<fragment>`, pacybara_prep's TITLE."""
        return f"{run_date}_{self.name}"

    @property
    def amplicon(self) -> str:
        return f"{self.fiveprime}{self.target}{self.threeprime}"

    @property
    def orfstart(self) -> int:
        """1-based start of `target` within the amplicon."""
        return len(self.fiveprime) + 1

    @property
    def orfend(self) -> int:
        """1-based, inclusive end of `target` within the amplicon."""
        return len(self.fiveprime) + len(self.target)


def read_sample_sheet(
    path: str | os.PathLike,
    defaults: Mapping[str, object],
    data_dir: str | os.PathLike | None = None,
) -> list[Library]:
    """Parse and validate `path`, returning one `Library` per row.

    `defaults` supplies every Pacybara argument a row leaves blank (or has no
    column for). A relative `fastq` is resolved against `data_dir`; an
    absolute one is kept as is. Whether the fastq exists is not checked here
    -- the caller may not see the host filesystem (Nextflow checks it when it
    stages the file).
    """
    errors: list[str] = []

    missing_defaults = [k for k in PARAMETER_KEYS if defaults.get(k) in (None, "")]
    if missing_defaults:
        errors.append(f"pacybara_defaults is missing: {', '.join(missing_defaults)}")
    unknown_defaults = sorted(set(defaults) - set(PARAMETER_KEYS))
    if unknown_defaults:
        errors.append(
            f"pacybara_defaults has unknown keys: {', '.join(unknown_defaults)}"
        )

    # utf-8-sig: pacybara_prep's own example sheet starts with a byte-order
    # mark, which would otherwise end up glued to the first column's name.
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        columns = [c.strip() for c in reader.fieldnames or []]
        rows = [
            {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            for row in reader
        ]

    missing_columns = [c for c in REQUIRED_COLUMNS if c not in columns]
    if missing_columns:
        errors.append(f"missing required column(s): {', '.join(missing_columns)}")
    unknown_columns = [c for c in columns if c not in REQUIRED_COLUMNS + PARAMETER_KEYS]
    if unknown_columns:
        errors.append(f"unknown column(s): {', '.join(unknown_columns)}")
    if errors:
        raise SampleSheetError(errors)
    if not rows:
        raise SampleSheetError(["sample sheet has no rows"])

    libraries = []
    seen_names: dict[str, int] = {}
    seen_fastqs: dict[str, int] = {}
    for line_number, row in enumerate(rows, start=2):

        def error(message: str) -> None:
            errors.append(f"row {line_number} ({row.get('gene') or '?'}): {message}")

        for column in REQUIRED_COLUMNS:
            if not row[column]:
                error(f"'{column}' is empty")
        for column in ("gene", "assay", "fragment"):
            if row[column] and not _NAME_RX.match(row[column]):
                error(f"'{column}' may only contain letters, digits, '.', '_' and '-'")
        for column in ("fiveprime", "target", "threeprime"):
            if row[column] and not _SEQUENCE_RX.match(row[column]):
                error(f"'{column}' may only contain A, C, G, T or N")

        parameters = {
            k: row.get(k) or ("" if defaults.get(k) is None else str(defaults[k]))
            for k in PARAMETER_KEYS
        }
        for key in _INTEGER_KEYS:
            if parameters[key] and not _INTEGER_RX.match(parameters[key]):
                error(f"{key} must be a non-negative integer, got '{parameters[key]}'")
        if parameters["minjaccard"] and not _FLOAT_RX.match(parameters["minjaccard"]):
            error(
                f"minjaccard must be written as 0.<digits> (Pacybara's own check), "
                f"got '{parameters['minjaccard']}'"
            )
        if parameters["clustermode"] and parameters["clustermode"] not in CLUSTER_MODES:
            error(
                f"clustermode must be one of {', '.join(CLUSTER_MODES)}, "
                f"got '{parameters['clustermode']}'"
            )

        fastq = row["fastq"]
        if fastq:
            if not fastq.endswith(".fastq.gz"):
                error(f"fastq must be a .fastq.gz file, got '{fastq}'")
            if not os.path.isabs(fastq):
                if data_dir is None:
                    error(
                        f"fastq '{fastq}' is a relative path but no data_dir was given"
                    )
                else:
                    fastq = os.path.join(os.fspath(data_dir), fastq)
            fastq = os.path.normpath(fastq)

        library = Library(
            gene=row["gene"],
            assay=row["assay"],
            fragment=row["fragment"],
            fiveprime=row["fiveprime"],
            target=row["target"],
            threeprime=row["threeprime"],
            fastq=fastq,
            parameters=parameters,
        )

        # pacybara_simplex.sh's findBarcodPos: the (upper-cased) amplicon
        # must contain the barcode exactly once (uptag) or twice
        # (uptag + downtag), matched literally as `grep -o` would.
        if parameters["barcode"] and library.amplicon:
            count = library.amplicon.upper().count(parameters["barcode"].upper())
            if count not in (1, 2):
                error(
                    f"the amplicon must contain the barcode '{parameters['barcode']}' "
                    f"once or twice, found {count}"
                )
            elif parameters["clustermode"] == "downtag" and count != 2:
                error("clustermode is downtag but the amplicon has only one barcode")

        if library.name in seen_names:
            error(
                f"duplicates row {seen_names[library.name]}: gene/assay/fragment "
                f"'{library.name}' must be unique (it names the output directory)"
            )
        else:
            seen_names[library.name] = line_number
        # One library per fastq: archiving moves each fastq once, after its
        # own run succeeds.
        if fastq and fastq in seen_fastqs:
            error(
                f"duplicates row {seen_fastqs[fastq]}: fastq '{fastq}' is already used"
            )
        elif fastq:
            seen_fastqs[fastq] = line_number

        libraries.append(library)

    if errors:
        raise SampleSheetError(errors)
    return libraries
