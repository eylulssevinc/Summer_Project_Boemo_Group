#!/bin/bash
#
# Same as run_before_after_500_kutem_gpu.slurm but for +/-100 events. Needs
# its own re-parse (a cache built at before-window=500 only retains reads
# that had >=500 pre-anchor events; re-extracting at 100 recovers the many
# additional reads that qualify at this much shorter, easier-to-satisfy
# threshold).
#
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --job-name=dnascent_before_after_100
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

REPO=/scratch/esevinc22/Summer_Project_Boemo_Group
DATA=$REPO/Eylul_gem_data
PYTHON=/opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python
SCRIPT=$DATA/scripts/before_after_match_plot.py
GEM=$DATA/slurm_full_run_allafterm/gem/gem.full.align
CONTROL=$DATA/slurm_full_run_allafterm/control/control.full.align
AFTER_MATRIX=$DATA/slurm_full_run_allafterm/figures_dtw_first100_banded/dtw_warped_matrix.tsv
OUTDIR=$DATA/slurm_full_run_allafterm/figures_before_after_100

for f in "$SCRIPT" "$GEM" "$CONTROL" "$AFTER_MATRIX"; do
    [ -e "$f" ] || { echo "Missing required input: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

"$PYTHON" "$SCRIPT" \
    --gem "$GEM" \
    --control "$CONTROL" \
    --before-window 100 \
    --after-warped-matrix "$AFTER_MATRIX" \
    --before-cache "$OUTDIR/cache_before100_matrix.tsv" \
    --outdir "$OUTDIR"

echo "Output files:"
ls -la "$OUTDIR"
