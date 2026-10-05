# pacbio-pipeline-nf

A Nextflow pipeline that maps variants to barcodes for a batch of PacBio
libraries with [Pacybara](https://github.com/rothlab/pacybara). You give it
one sample sheet, and it does what used to be done by hand under
`Pacbio_processing/`:

1. Validates the sample sheet and writes one Pacybara parameter file per
   library. This is
   [`pacybara_prep`](https://github.com/dlholmes2117/pacybara_prep)'s
   `make_parameter.py`, with `INFASTQ` filled in for you.
2. Runs Pacybara on every library in parallel, each as its own job
   (`pacybara_simplex.sh` in a container). It stops before the barseqPro
   steps.
3. Publishes each library's results to `<output_dir>/<TITLE>/`, together
   with the parameter file the library was run with.
4. Optionally moves each successfully mapped `.fastq.gz` into
   `completed_fastqs/`.

```text
sample_sheet.csv ──► MAKE_PARAMETER_FILES ──► RUN_PACYBARA (per library) ──► pacybara_output/<TITLE>/
                                                        │
                                                        └──► ARCHIVE_FASTQ (--archive_inputs true) ──► completed_fastqs/
```

## Sample sheet

The columns are `pacybara_prep`'s `pacybara_setup.csv` columns plus `fastq`:

| column | meaning |
|---|---|
| `gene`, `assay`, `fragment` | Name the library. TITLE is `<run_date>_<gene>_<assay>_<fragment>`. |
| `fiveprime`, `target`, `threeprime` | The amplicon, split around the mutated region. ORFSTART and ORFEND are the positions of `target`. |
| `fastq` | The library's `.fastq.gz`, either absolute or relative to `--data_dir`. |
| *(optional)* `barcode`, `maxqdrops`, `minbcq`, `minjaccard`, `minmatches`, `maxdiff`, `minqual`, `clustermode` | Override `params.yaml`'s `pacybara_defaults` for that row. A blank cell uses the default. |

```csv
gene,assay,fragment,fiveprime,target,threeprime,fastq
SOS1,abundance,na,CTAGCG...,ATGCAG...,TGA...,m84000_SOS1_abundance.fastq.gz
KSR1,activity,nterm,CTAGCG...,ATGGAT...,GGC...,m84000_KSR1_activity.fastq.gz
```

The whole sheet is checked before any job starts. Every problem is reported
at once, with its row number. The checks are:

- every required column is present and filled in
- sequences contain only ACGTN
- each `fastq` ends in `.fastq.gz` and exists
- no gene/assay/fragment combination or fastq appears twice
- Pacybara's own checks on its arguments pass
- the barcode appears once or twice in each amplicon

The defaults in `params.yaml` are `pacybara_prep`'s values. Pacybara's README
recommends `minqual` around 67 for CCS reads and around 27 for DeepConsensus
reads.

## Running

You need [Nextflow](https://www.nextflow.io/) (Java 17+) plus Docker or
Apptainer. Point the pipeline at the shared `Pacbio_processing/` tree:

```bash
P=/net/fowler/vol1/shared/Pacbio_processing
nextflow run . -params-file params.yaml -profile apptainer -c site.config \
    --sample_sheet $P/to_process/parameter_files/sample_sheet.csv \
    --data_dir $P/to_process/data \
    --output_dir $P/pacybara_output \
    --archive_inputs true --completed_fastq_dir $P/completed_fastqs
```

- **One failed library doesn't stop the rest.** Nextflow logs it as a warning
  and the run still exits 0. That library's fastq stays in `data/`, so read
  the run summary.
- **`run_date`** defaults to today, matching `pacybara_prep`. A `-resume` on a
  later day keeps the original titles.
- **Parameter files** are generated inside Nextflow's work directory. Each
  library's copy is published into its output folder, so there is nothing
  left over to delete.

### Cluster settings (`site.config`)

This repo contains no scheduler settings, so put your own in `-c site.config`.
Pacybara's README recommends 12 cores and 24 GB per library. Simplex mode
runs the whole library on one node, so give it a generous walltime. Example
for SGE:

```groovy
process {
    executor = 'sge'
    penv = 'serial'
    clusterOptions = '-S /bin/bash'
    withLabel: 'process_single' { cpus = 1;  memory = '2 GB';  time = '1h' }
    withLabel: 'process_high'   { cpus = 12; memory = '24 GB'; time = '72h' }
}
apptainer.cacheDir = '/net/fowler/vol1/shared/apptainer_cache'   // any shared dir
```

`ARCHIVE_FASTQ` always runs on the head node, outside the container.

## Container

`Dockerfile` builds the single image every task runs in:

- Pacybara at a pinned commit (`PACYBARA_COMMIT`)
- Pacybara's conda env from its own `pacybara_env.yml`
- Pacybara's R packages, with the GitHub ones pinned
- this package, in a separate Python 3.13 venv

```bash
docker build --platform linux/amd64 -t pacybara-workflow:latest .
```

On every push to `main`, `.github/workflows/docker.yml` publishes it to
`ghcr.io/fowlerlab/pacbio-pipeline-nf` with a `:latest` tag and a
`:<short-sha>` tag. Set `container_image` in `params.yaml` (or pass
`--container_image`) to the tag you want. On the cluster, Apptainer pulls the
same image as `docker://...`.

## Development

```bash
uv sync --group dev
uv run pytest tests/unit          # sample sheet + parameter files
uv run pytest tests/integration   # real `nextflow run -profile local`, stub Pacybara
uv run pre-commit run --all-files
```

To run only the parameter-file step, outside Nextflow:

```bash
uv run pacybara-make-params --sample-sheet sheet.csv --data-dir /path/to/data \
    --defaults '{"barcode": "NNNNNNNNNNNNNNNN", "maxqdrops": 4, "minbcq": 32, "minjaccard": 0.2, "minmatches": 1, "maxdiff": 1, "minqual": 27, "clustermode": "uptag"}'
```

## License

[MIT](LICENSE.txt). The parameter-file template is adapted from
`pacybara_prep` (BSD-3-Clause, Dan Holmes). Pacybara itself (GPL-3.0) is
cloned into the image at build time and is not vendored here.
