#!/bin/bash
# Stand-in for Pacybara's pacybara_simplex.sh in the integration suite:
# reads the parameter file the way the real script does (sourcing its
# ARGUMENTS section), checks INFASTQ is readable, and writes a fake
# clusters_transl_softfilter.csv.gz into WORKSPACE. A fastq whose name
# contains "fail" makes it exit 1, standing in for a failed Pacybara run.
set -euo pipefail

CPUS=""
while (( "$#" )); do
  case "$1" in
    -c|--cpus) CPUS=$2; shift 2;;
    *) PARAMETERS=$1; shift;;
  esac
done

ARGS=$(mktemp)
sed -n '/#BEGIN ARGUMENTS/,/#END ARGUMENTS/p' "$PARAMETERS" > "$ARGS"
source "$ARGS"
[[ -r "$INFASTQ" ]] || { echo "cannot read $INFASTQ" >&2; exit 1; }
if [[ "$INFASTQ" == *fail* ]]; then
  echo "stub pacybara failure for $INFASTQ" >&2
  exit 1
fi

mkdir -p tmp "${WORKSPACE}"
OUTPREFIX=$(basename "${INFASTQ%.fastq.gz}")
CLUSTERDIR="${WORKSPACE}/${OUTPREFIX}_clustering"
mkdir -p "$CLUSTERDIR/qc"
printf 'upBarcode,size\nACGT,%s\n' "$CPUS" | gzip -c > "$CLUSTERDIR/clusters_transl_softfilter.csv.gz"
cp "$PARAMETERS" "${WORKSPACE}/stub_saw_parameters.txt"
