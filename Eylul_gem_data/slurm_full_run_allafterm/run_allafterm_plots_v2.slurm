#!/bin/bash
#
# Rerun plot_allafterm_analysis.py after adding Part F (full, untruncated
# ALLAFTERM per-read boxplots + full HMM-anchored trace). More memory than
# the original run since Part F builds a dense (n_reads x max_len) matrix
# per sample with no cap -- control's max ALLAFTERM length is ~725k events,
# so its matrix alone is ~1.9GB; budget generously.
#
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --qos=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --job-name=dnascent_allafterm_plots_v2
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

REPO=/scratch/esevinc22/Summer_Project_Boemo_Group
DATA=$REPO/Eylul_gem_data
PYTHON=/opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python
SCRIPT=$DATA/scripts/plot_allafterm_analysis.py
GEM=$DATA/slurm_full_run_allafterm/gem/gem.full.align
CONTROL=$DATA/slurm_full_run_allafterm/control/control.full.align
OUTDIR=$DATA/slurm_full_run_allafterm/figures_allafterm

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
