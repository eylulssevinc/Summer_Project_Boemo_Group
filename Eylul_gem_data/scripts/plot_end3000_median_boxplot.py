#!/usr/bin/env python3
"""
Standalone boxplot: per-read median scaled signal over the last 3000 raw
samples before the physical 3' end (column `end_3000_median`), gem vs
control -- the same end-anchored window shown as a band/IQR-over-position
plot in tail_end_anchored_3000.png, but here summarized to one number per
read and compared as a boxplot instead.

Reads the read-level summary TSV that plot_tail_signal_end_focused.py
already produces (../slurm_full_run/figures_end_focused/tail_end_read_level_summary.tsv)
-- no re-parsing of the .align files needed, since end_3000_median is already
computed there (max_lags default includes 3000). See that TSV's column list
and TAIL_SIGNAL_CAPTURE.md sec 7 for what end_3000_median means exactly.

Run with the gmconda_py368 conda env:
    /opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python \\
        plot_end3000_median_boxplot.py \\
        --summary ../slurm_full_run/figures_end_focused/tail_end_read_level_summary.tsv \\
        --outpath ../slurm_full_run/figures_end_focused/tail_end_3000_median_boxplot.png

Same caveats as plot_tail_signal_end_focused.py apply: raw-sample positions
are not base positions, and this is read-level but still a small,
exploratory sample -- don't read a Gem-vs-control conclusion off box overlap
alone.
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", required=True, help="path to tail_end_read_level_summary.tsv")
    ap.add_argument("--outpath", required=True, help="output PNG path")
    args = ap.parse_args()

    df = pd.read_csv(args.summary, sep="\t")
    gem_vals = df.loc[df["sample"] == "gem", "end_3000_median"].dropna().to_numpy()
    control_vals = df.loc[df["sample"] == "control", "end_3000_median"].dropna().to_numpy()

    fig, ax = plt.subplots(figsize=(5, 5))
    data, labels, colors = [], [], []
    if len(gem_vals):
        data.append(gem_vals); labels.append(f"gem\n(n={len(gem_vals)})"); colors.append("tab:red")
    if len(control_vals):
        data.append(control_vals); labels.append(f"control\n(n={len(control_vals)})"); colors.append("tab:gray")

    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=True)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)

    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_ylabel("Per-read median scaled signal, last 3000 raw samples\nbefore physical 3′ end")
    ax.set_title("Read-level median tail signal (last 3000 samples): gem vs control")
    fig.tight_layout()
    fig.savefig(args.outpath, dpi=150)
    plt.close(fig)

    print(f"wrote {args.outpath}")
    print(f"gem:     n={len(gem_vals)}, median of medians={np.median(gem_vals):.4f}" if len(gem_vals) else "gem: no data")
    print(f"control: n={len(control_vals)}, median of medians={np.median(control_vals):.4f}" if len(control_vals) else "control: no data")


if __name__ == "__main__":
    main()
