#!/bin/bash
#
# Follow-up to Part F: how many ALLAFTERM events are actually "rescued" per
# read, where does per-offset read coverage collapse, and which window is
# therefore trustworthy to compare gem vs control in at the read level.
# Re-parses the same full .align files as run_allafterm_plots_v2.slurm, so
# budget similarly.
#
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --qos=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=dnascent_rescue_window
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

REPO=/scratch/esevinc22/Summer_Project_Boemo_Group
DATA=$REPO/Eylul_gem_data
PYTHON=/opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python
SCRIPT=$DATA/scripts/analyze_rescue_window.py
GEM=$DATA/slurm_full_run_allafterm/gem/gem.full.align
CONTROL=$DATA/slurm_full_run_allafterm/control/control.full.align
OUTDIR=$DATA/slurm_full_run_allafterm/figures_window_analysis

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
