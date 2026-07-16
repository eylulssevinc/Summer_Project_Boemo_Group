#!/bin/bash
#
# Same banded/regularized DTW pipeline as run_dtw_first100_kutem_gpu.slurm,
# but for the first 500 ALLAFTERM (post-last-Match) events per read. Needs
# its own re-parse of the full .align files (no cache yet at this window)
# and its own barycenter reference / validation, since both scale with
# window size.
#
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --job-name=dnascent_dtw_first500_kutemgpu
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

REPO=/scratch/esevinc22/Summer_Project_Boemo_Group
DATA=$REPO/Eylul_gem_data
PYTHON=/opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python
SCRIPT=$DATA/scripts/dtw_first100_align.py
GEM=$DATA/slurm_full_run_allafterm/gem/gem.full.align
CONTROL=$DATA/slurm_full_run_allafterm/control/control.full.align
OUTDIR=$DATA/slurm_full_run_allafterm/figures_dtw_first500_banded

for f in "$SCRIPT" "$GEM" "$CONTROL"; do
    [ -e "$f" ] || { echo "Missing required input: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

"$PYTHON" "$SCRIPT" \
    --gem "$GEM" \
    --control "$CONTROL" \
    --outdir "$OUTDIR" \
    --window 500 \
    --band-radius 20 \
    --step-penalty 0.3 \
    --barycenter-iters 3 \
    --shuffle-test-n 60 \
    --force-large-window \
    --cache "$OUTDIR/cache_first500_matrix.tsv"

echo "Output files:"
ls -la "$OUTDIR"
