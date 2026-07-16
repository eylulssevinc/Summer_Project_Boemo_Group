#!/bin/bash
#
# Run plot_tail_signal_end_focused.py against the full-dataset gem/control
# tail-capture output. Single-threaded Python (pandas/numpy/matplotlib), but
# parses ~3GB of .align data across two passes per file, so this runs via
# SLURM rather than on the shared login node.
#
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --qos=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --job-name=dnascent_endfocused_plots
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

REPO=/scratch/esevinc22/Summer_Project_Boemo_Group
DATA=$REPO/Eylul_gem_data
PYTHON=/opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python
SCRIPT=$DATA/scripts/plot_tail_signal_end_focused.py
GEM=$DATA/slurm_full_run/gem/gem.full.align
CONTROL=$DATA/slurm_full_run/control/control.full.align
OUTDIR=$DATA/slurm_full_run/figures_end_focused

for f in "$SCRIPT" "$GEM" "$CONTROL"; do
    [ -e "$f" ] || { echo "Missing required input: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

"$PYTHON" "$SCRIPT" \
    --gem "$GEM" \
    --control "$CONTROL" \
    --outdir "$OUTDIR"

echo "Output files:"
ls -la "$OUTDIR"
