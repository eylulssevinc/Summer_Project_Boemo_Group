#!/bin/bash
#
# Smoke test for the --tail-source {none,source1,source2,both} option.
# Dataset-agnostic (DATASET=gem or DATASET=control via --export=ALL,DATASET=...).
#
# MODES defaults to all four (none source1 source2 both) so mode-switching
# logic itself gets validated; pass --export=ALL,DATASET=...,MODES="both" to
# only run the one mode that actually matters once mode-switching has
# already been validated on another dataset (saves ~3/4 of the runtime).
#
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --qos=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:30:00
#SBATCH --job-name=dnascent_tailsource_smoketest
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

echo "=== job start: $(date) ==="
echo "host: $(hostname)"
echo "job ID: ${SLURM_JOB_ID:-unknown}"
echo "dataset: ${DATASET:-gem}"
echo "modes: ${MODES:-none source1 source2 both}"

if [[ -z "${DATASET:-}" ]]; then
    DATASET=gem
fi
if [[ "${DATASET}" != "gem" && "${DATASET}" != "control" ]]; then
    echo "ERROR: DATASET must be 'gem' or 'control', got '${DATASET}'" >&2
    exit 1
fi
if [[ -z "${MODES:-}" ]]; then
    MODES="none source1 source2 both"
fi

REPO=/scratch/esevinc22/Summer_Project_Boemo_Group
DATA=${REPO}/Eylul_gem_data
BIN=${REPO}/bin/DNAscent
REF=${DATA}/ligation_product.fasta
BAM=${DATA}/${DATASET}/basecalled.sorted.bam
INDEX=${DATA}/${DATASET}/dnascent_align/index.dnascent
OUTDIR=${DATA}/smoketest_by_source

THREADS=${SLURM_CPUS_PER_TASK:-16}

for f in "$BIN" "$REF" "$BAM" "$INDEX"; do
    if [[ ! -e "$f" ]]; then
        echo "ERROR: required input missing: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"

for mode in ${MODES}; do
    OUT="${OUTDIR}/${DATASET}.full.${mode}.align"
    echo "=== mode: ${mode} ==="
    "$BIN" align -b "$BAM" -r "$REF" -i "$INDEX" -o "$OUT" -t "$THREADS" --tail-source "$mode"
    echo
done

echo "=== job end: $(date) ==="

echo "=== line counts ==="
for mode in ${MODES}; do
    OUT="${OUTDIR}/${DATASET}.full.${mode}.align"
    echo "-- ${mode} --"
    echo -n "  TAIL_SOURCE1 lines: "; grep -c "^TAIL_SOURCE1" "$OUT" || true
    echo -n "  TAIL_SOURCE2 lines: "; grep -c "^TAIL_SOURCE2" "$OUT" || true
    echo -n "  legacy TAIL lines:  "; grep -c "^TAIL[[:space:]]" "$OUT" || true
done
