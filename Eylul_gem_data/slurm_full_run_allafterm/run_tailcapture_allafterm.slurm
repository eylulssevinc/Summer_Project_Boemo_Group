#!/bin/bash
#
# DNAscent align (with --tail-source both) on the full gem/control datasets,
# using the ALLAFTERM-capable binary. Deliberately writes into
# slurm_full_run_allafterm/ (NOT slurm_full_run/) so the pre-existing
# canonical .align files, figures, and TSVs are left untouched until the new
# ALLAFTERM output has been validated -- see TAIL_SIGNAL_CAPTURE.md's
# validation section for the byte-identical-modulo-ALLAFTERM diff this run's
# output is used for.
#
# Same -l 600 as the original slurm_full_run/run_tailcapture.slurm -- see
# that script's header comment for why -l 600 is load-bearing. Do not change
# it here either.
#
# Submit with (from Eylul_gem_data/slurm_full_run_allafterm/):
#   sbatch --job-name=dnascent_allafterm_gem     --output=logs/gem_%j.out     --error=logs/gem_%j.err     --export=ALL,DATASET=gem     run_tailcapture_allafterm.slurm
#   sbatch --job-name=dnascent_allafterm_control --output=logs/control_%j.out --error=logs/control_%j.err --export=ALL,DATASET=control run_tailcapture_allafterm.slurm
#
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --qos=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=dnascent_allafterm
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

echo "=== job start: $(date) ==="
echo "host: $(hostname)"
echo "job ID: ${SLURM_JOB_ID:-unknown}"
echo "cpus-per-task: ${SLURM_CPUS_PER_TASK:-unset}"
echo "dataset: ${DATASET:-unset}"

if [[ -z "${DATASET:-}" ]]; then
    echo "ERROR: DATASET env var must be set to 'gem' or 'control' (via --export=ALL,DATASET=...)" >&2
    exit 1
fi
if [[ "${DATASET}" != "gem" && "${DATASET}" != "control" ]]; then
    echo "ERROR: DATASET must be 'gem' or 'control', got '${DATASET}'" >&2
    exit 1
fi

REPO=/scratch/esevinc22/Summer_Project_Boemo_Group
DATA=${REPO}/Eylul_gem_data
OUTDIR=${DATA}/slurm_full_run_allafterm/${DATASET}
BIN=${REPO}/bin/DNAscent
REF=${DATA}/ligation_product.fasta
BAM=${DATA}/${DATASET}/basecalled.sorted.bam
INDEX=${DATA}/${DATASET}/dnascent_align/index.dnascent
OUT=${OUTDIR}/${DATASET}.full.align

THREADS=${SLURM_CPUS_PER_TASK:-16}

for f in "$BIN" "$REF" "$BAM" "$INDEX"; do
    if [[ ! -e "$f" ]]; then
        echo "ERROR: required input missing: $f" >&2
        exit 1
    fi
done

echo "binary: $BIN"
echo "reference: $REF"
echo "bam: $BAM"
echo "index: $INDEX"
echo "output: $OUT"
echo "threads: $THREADS"
echo "=== running DNAscent align ==="

"$BIN" align \
    -b "$BAM" \
    -r "$REF" \
    -i "$INDEX" \
    -o "$OUT" \
    -t "$THREADS" \
    -l 600 \
    --tail-source both

echo "=== DNAscent align finished: $(date) ==="
echo "output size: $(du -h "$OUT" 2>/dev/null || echo 'missing')"
echo "TAIL_SOURCE1 line count: $(grep -c '^TAIL_SOURCE1' "$OUT" 2>/dev/null || echo 0)"
echo "TAIL_SOURCE2 line count: $(grep -c '^TAIL_SOURCE2' "$OUT" 2>/dev/null || echo 0)"
echo "NEXT100EVT line count: $(grep -c '^NEXT100EVT' "$OUT" 2>/dev/null || echo 0)"
echo "ALLAFTERM line count: $(grep -c '^ALLAFTERM' "$OUT" 2>/dev/null || echo 0)"
echo "=== job end: $(date) ==="
