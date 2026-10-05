#!/bin/bash
# Smoke-test the real image: Pacybara's dependency check, then a synthetic
# two-library run of the whole pipeline under the default docker profile
# (real pacybara_simplex.sh, real staging and publishing, real archiving).
# Not part of CI -- the amd64 image is slow to build and to emulate.
#
#   docker build --platform linux/amd64 -t pacybara-workflow:latest .
#   tests/integration/container_smoke.sh [image]
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
IMAGE=${1:-pacybara-workflow:latest}
WORK=$(mktemp -d)
echo "Working in $WORK"

docker run --rm --platform linux/amd64 "$IMAGE" bash -c '
  set -e
  for BIN in muscle bwa bowtie2 samtools seqret Rscript python3 pacybara_simplex.sh pacybara-make-params; do
    command -v $BIN >/dev/null || { echo "missing $BIN"; exit 1; }
  done
  Rscript -e "for (p in c(\"argparser\",\"hash\",\"bitops\",\"pbmcapply\",\"yogitools\",\"yogiseq\",\"hgvsParseR\")) library(p, character.only=TRUE)"
  echo "dependencies OK"'

mkdir -p "$WORK"/{data,output,completed}
python3 - "$WORK" <<'EOF'
import gzip, random, sys
work = sys.argv[1]
rng = random.Random(0)
seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
fiveprime_flank, spacer, threeprime = seq(60), seq(20), seq(30)
orf = "ATG" + "".join(rng.choice(["GCC", "GAA", "CTG", "AAG", "TCC"]) for _ in range(40)) + "TAA"
fiveprime = fiveprime_flank + "N" * 16 + spacer
rows = []
for name in ("ALPHA", "BETA"):
    with gzip.open(f"{work}/data/{name.lower()}.fastq.gz", "wt") as f:
        # Two features every real library has and Pacybara's R steps
        # assume (both crash on data without them):
        #   - a barcode collision: clone 30 reuses clone 0's barcode with 2
        #     reads (else pacybara_softfilter.R's sapply() simplifies to a
        #     vector and its do.call() fails);
        #   - wild-type clones: every 5th but 30 (else pacybara_qc.R's
        #     `cvarTally[-which(names(cvarTally)=="=")]` drops every variant
        #     and its jackpot plot fails).
        barcodes = [seq(16) for _ in range(30)]
        for clone in range(31):
            barcode = barcodes[clone % 30]
            variant = list(orf)
            pos = rng.randrange(3, len(orf) - 3)
            if clone % 5 or clone == 30:
                variant[pos] = rng.choice([b for b in "ACGT" if b != orf[pos]])
            read = fiveprime_flank + barcode + spacer + "".join(variant) + threeprime
            for r in range(6 if clone < 30 else 2):
                f.write(f"@m00000/{clone}{r}/ccs\n{read}\n+\n{'~' * len(read)}\n")
    rows.append(f"{name},abundance,na,{fiveprime},{orf},{threeprime},{name.lower()}.fastq.gz")
with open(f"{work}/sheet.csv", "w") as f:
    f.write("gene,assay,fragment,fiveprime,target,threeprime,fastq\n" + "\n".join(rows) + "\n")
EOF

cd "$WORK"
nextflow run "$REPO" -params-file "$REPO/params.yaml" \
    --container_image "$IMAGE" \
    --sample_sheet "$WORK/sheet.csv" --data_dir "$WORK/data" \
    --output_dir "$WORK/output" --run_date 20260101 \
    --archive_inputs true --completed_fastq_dir "$WORK/completed"

for NAME in ALPHA BETA; do
  OUT="$WORK/output/20260101_${NAME}_abundance_na"
  LOWER=$(echo "$NAME" | tr '[:upper:]' '[:lower:]')
  for F in "${LOWER}_clustering/clusters_transl_softfilter.csv.gz" "20260101_${NAME}_abundance_na.txt"; do
    [[ -s "$OUT/$F" ]] || { echo "FAIL: missing $OUT/$F"; exit 1; }
  done
  echo "$NAME clusters: $(( $(zcat < "$OUT/${LOWER}_clustering/clusters_transl_softfilter.csv.gz" | wc -l) - 1 ))"
  [[ -e "$WORK/completed/${LOWER}.fastq.gz" ]] || { echo "FAIL: ${LOWER}.fastq.gz not archived"; exit 1; }
done
echo "container smoke test OK ($WORK)"
