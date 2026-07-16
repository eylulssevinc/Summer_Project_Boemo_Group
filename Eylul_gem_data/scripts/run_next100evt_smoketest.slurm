#!/bin/bash
#
# Small, fast smoke test for the new NEXT100EVT pool (added alongside
# --tail-source, anchored at the same boundary Source 2 uses but capped at
# 100 events and kept at event-mean granularity instead of expanded raw
# samples). Uses -m to cap reads processed so this finishes in well under a
# minute, with -l 600 to match the canonical filter.
#
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --qos=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --job-name=dnascent_next100evt_smoketest
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

echo "=== job start: $(date) ==="

REPO=/scratch/esevinc22/Summer_Project_Boemo_Group
DATA=${REPO}/Eylul_gem_data
BIN=${REPO}/bin/DNAscent
REF=${DATA}/ligation_product.fasta
OUTDIR=${DATA}/smoketest_by_source
THREADS=${SLURM_CPUS_PER_TASK:-4}

mkdir -p "$OUTDIR"

for DATASET in gem control; do
    BAM=${DATA}/${DATASET}/basecalled.sorted.bam
    INDEX=${DATA}/${DATASET}/dnascent_align/index.dnascent
    OUT=${OUTDIR}/${DATASET}.next100evt_smoketest.align

    echo "=== dataset: ${DATASET} ==="
    "$BIN" align -b "$BAM" -r "$REF" -i "$INDEX" -o "$OUT" -t "$THREADS" -l 600 -m 500 --tail-source both

    echo "-- ${DATASET} line counts --"
    echo -n "  TAIL_SOURCE1 lines: "; grep -c "^TAIL_SOURCE1" "$OUT" || true
    echo -n "  TAIL_SOURCE2 lines: "; grep -c "^TAIL_SOURCE2" "$OUT" || true
    echo -n "  NEXT100EVT lines:   "; grep -c "^NEXT100EVT" "$OUT" || true
    echo
done

echo "=== job end: $(date) ==="
