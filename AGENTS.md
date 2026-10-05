# AGENTS.md — pacbio-pipeline-nf

## Project overview

This is a Nextflow pipeline that automates the Fowler lab's PacBio → Pacybara
mapping step. It takes a sample sheet and validates it, writes one Pacybara
parameter file per library, runs `pacybara_simplex.sh` once per library in a
container, publishes the results, and optionally archives the input fastqs.
`README.md` documents the user-facing contract: the sample sheet columns,
parameters, cluster config and outputs.

Upstream code it depends on:

- **`rothlab/pacybara`** is cloned into the image at
  `Dockerfile`'s `PACYBARA_COMMIT`. It is never vendored or patched.
- **`dlholmes2117/pacybara_prep`** (BSD-3) is the origin of the
  parameter-file template in `src/pacybara_workflow/parameter_file.py`.
  `tests/unit/data/golden/` holds that script's own output for
  `tests/unit/data/pacybara_prep_setup.csv`. The golden test requires our
  output to be byte-identical apart from `INFASTQ`. Don't edit the golden
  files to make a test pass. Regenerate them with the upstream script.

## Repo conventions

These match the sibling repos `fisseq-data-pipeline` and
`fisseq-embeddings-pipeline`:

- **Python owns validation and rendering.** All checks on the sample sheet
  belong in `src/pacybara_workflow/sample_sheet.py`, which reports every
  error at once. Groovy only fans the manifest out with `splitCsv`.
  - The one exception is `file(..., checkIfExists: true)` on each fastq.
    That check runs on the host, which the containerized Python can't see.
- **No runtime Python dependencies.** The package uses stdlib `csv` and
  `argparse`. Pacybara's toolchain lives in the image's conda env, not in
  `pyproject.toml`.
- **Nextflow processes:** one per directory under
  `modules/local/<name>/main.nf`, wired together in `workflows/pacybara.nf`.
  Each process has a `label`, `container "${params.container_image}"` and,
  if it publishes anything, `publishDir ..., mode: 'copy'`.
  - `RUN_PACYBARA` has `errorStrategy 'ignore'`, so one failed library
    doesn't stop the others.
  - Booleans are read with `.toString().toBoolean()`.
- **Config:** defaults belong in `params.yaml`, never in `nextflow.config`
  or a profile. Don't add a `params {}` block to `nextflow.config` (see the
  comment there).
- **No scheduler-specific code.** Cluster settings belong in the user's
  `-c site.config`. `README.md` has an SGE example.

## Pacybara gotchas

- `pacybara_simplex.sh` doesn't copy its parameter file into the workspace,
  although `pacybara.sh` does. `RUN_PACYBARA` copies it explicitly.
- It writes `tmp/` into the current directory. `RUN_PACYBARA` deletes it so
  it isn't confused with output.
- It runs under `set -u` and reads `$CONDA_DEFAULT_ENV`, so it crashes
  if that variable is unset. The image puts the env's `bin/` on PATH instead
  of activating it, so the Dockerfile sets `CONDA_DEFAULT_ENV=base`. The
  "No conda environment detected" message this produces is harmless.
- The barcode is matched literally against the upper-cased amplicon
  (`findBarcodPos`, using `grep -o`). It must occur once (uptag only) or
  twice. `sample_sheet.py` mirrors this check.
- **Upstream bug, not patched:** `pacybara_simplex.sh` compares
  `$CLUSTERMODE == "DOWNTAG"` in upper case. As a result, `clustermode:
  downtag` still feeds the **uptag** barcode file into
  `pacybara_cluster.py`. Only the pre-clustering step uses the downtag
  reads.
- `muscle` is required by Pacybara's dependency check, but no step calls
  it.
- Two of Pacybara's R steps crash on unusually clean libraries. Real
  libraries won't trigger them, but small synthetic test data will
  (`container_smoke.sh` builds its data around both):
  - `pacybara_softfilter.R` fails ("second argument must be a list") if no
    barcode is shared by two clusters.
  - `pacybara_qc.R`'s jackpot plot fails ("need finite 'xlim' values") if
    no clone is wild-type.

## Testing

```bash
uv run pytest tests/unit
uv run pytest tests/integration    # needs nextflow + java; skips otherwise
uv run ruff check --fix . && uv run ruff format .
```

The integration suite runs `nextflow run -profile local` with
`tests/integration/stubs/pacybara_simplex.sh` first on PATH. The stub parses
the parameter file the way the real script does. A fastq whose name contains
`fail` makes the stub fail, which exercises the archive-only-successes
path. Nothing in CI runs real Pacybara.

To smoke-test the real image, build it with `docker build --platform
linux/amd64 -t pacybara-workflow:latest .`. Then
`tests/integration/container_smoke.sh` runs the dependency check and a
synthetic library through real Pacybara.

## Git workflow

Branch off `main` for each unit of work. Run the unit tests and ruff
before every commit, and the integration suite before merging. Merge with
`--no-ff`. Don't rewrite history that's already on `main`.
