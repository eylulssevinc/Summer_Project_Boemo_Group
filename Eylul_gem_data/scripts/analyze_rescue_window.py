#!/usr/bin/env python3
"""Follow-up analysis, not a new pipeline: answers "how many ALLAFTERM events
are rescued, and which window is trustworthy enough to plot the gem-vs-control
comparison in?"

Motivation: allafterm_hmm_anchored_full_trace.png (Part F) showed the visible
gem-vs-control difference seen in the first-100/500-event bands (Parts A/B)
apparently vanish once the full ALLAFTERM range is plotted. The coverage
subpanel in that figure already hints why -- the number of unique reads
contributing collapses as event offset grows, and does so at a different rate
for gem (50 reads total) than control (321 reads total). This script makes
that explicit:

  1. Per-read ALLAFTERM length distributions for gem and control (how many
     events are "rescued" per read, and in total), plus explicit outlier
     flagging (longest reads per sample, and what fraction of each sample's
     TOTAL rescued events they alone account for).
  2. The exact event-offset at which each sample's unique-read coverage first
     drops below 100%/95%/90%/80%/50% of its starting count -- this defines,
     from the data itself rather than an arbitrary round number, which window
     is "fully powered" (no attrition), "lightly attrited", or "small-n".
  3. Read-level (one value per read, not pooled events) median/mean
     Mann-Whitney U tests at several windows -- including the coverage-backed
     windows from (2) -- so we can see whether the first-100-events
     difference is a read-level effect, and at what window size (if any) it
     stops being distinguishable from the two groups' overlap.
  4. An outlier-sensitivity check on the FULL (untruncated) per-read
     mean/median comparison: same test, with the most extreme control
     read(s) excluded, to see whether a couple of pathological reads are
     driving (or hiding) any full-range effect.
  5. Two new boxplot figures (windows chosen from step 2's coverage
     boundaries) and one new zoomed coverage-decay figure, so the "which
     events would gem's effect be seen most clearly in" question has a
     concrete, defensible answer to plot from.

Reuses load_blocks / read_level_metrics / plot_read_level_boxplots /
plot_focused_median_boxplot from plot_allafterm_analysis.py rather than
reimplementing the dedup or windowing logic. Writes to a NEW output
directory; never touches the canonical figures_allafterm/ outputs.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_allafterm_analysis import (
    load_blocks,
    read_level_metrics,
    plot_read_level_boxplots,
    plot_focused_median_boxplot,
)

COVERAGE_FRACTIONS = (1.0, 0.95, 0.90, 0.80, 0.50)
STAT_WINDOWS = (100, 500, 1000, 1500, 2000, 3000)
BIG_WINDOW = 10**9  # sentinel: arr[:BIG_WINDOW] == arr, i.e. "untruncated"


def mwu_row(window, metric, gem_vals, control_vals, note=""):
    gem_vals = np.asarray(gem_vals, dtype=float)
    control_vals = np.asarray(control_vals, dtype=float)
    gem_vals = gem_vals[~np.isnan(gem_vals)]
    control_vals = control_vals[~np.isnan(control_vals)]
    row = {
        "window": window, "metric": metric,
        "n_gem": len(gem_vals), "n_control": len(control_vals),
        "U": np.nan, "p_value": np.nan, "auc": np.nan, "note": note,
    }
    if len(gem_vals) and len(control_vals):
        U, p = mannwhitneyu(gem_vals, control_vals, alternative="two-sided")
        row["U"] = float(U)
        row["p_value"] = float(p)
        row["auc"] = float(U / (len(gem_vals) * len(control_vals)))
    return row


def coverage_crossings(lengths, fractions=COVERAGE_FRACTIONS, max_offset=None):
    """lengths: 1D array of per-read ALLAFTERM lengths for one sample.
    Returns dict fraction -> first event-offset (1-based) at which the
    number of reads with length >= offset drops below fraction * n0, plus
    the full per-offset coverage curve up to max_offset (for plotting)."""
    lengths = np.sort(np.asarray(lengths))[::-1]  # descending
    n0 = len(lengths)
    if max_offset is None:
        max_offset = int(lengths.max()) if n0 else 0
    offsets = np.arange(1, max_offset + 1)
    # n_unique_reads(x) = count of reads with length >= x
    # lengths sorted descending; searchsorted on ascending copy is simpler:
    asc = np.sort(lengths)
    n_at_offset = n0 - np.searchsorted(asc, offsets, side="left")
    crossings = {}
    for frac in fractions:
        thresh = frac * n0
        idx = np.argmax(n_at_offset < thresh) if np.any(n_at_offset < thresh) else -1
        crossings[frac] = int(offsets[idx]) if idx >= 0 else None
    return crossings, offsets, n_at_offset, n0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gem", required=True)
    ap.add_argument("--control", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--zoom-max-offset", type=int, default=5000,
                     help="max x for the zoomed coverage-decay figure (default 5000)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("Loading ALLAFTERM blocks (this re-parses the full .align files; "
          "expect this to take a while for control)...")
    main_gem, _, _, _ = load_blocks(args.gem)
    main_control, _, _, _ = load_blocks(args.control)
    gem_allafterm = {rid: v[1] for rid, v in main_gem.items()}
    control_allafterm = {rid: v[1] for rid, v in main_control.items()}

    # ---------------- 1. length distributions + outlier flagging ----------
    gem_lengths = {rid: len(a) for rid, a in gem_allafterm.items()}
    control_lengths = {rid: len(a) for rid, a in control_allafterm.items()}
    print("\n=== ALLAFTERM per-read length summary ===")
    len_rows = []
    for sample, d in (("gem", gem_lengths), ("control", control_lengths)):
        arr = np.array(list(d.values()))
        total = int(arr.sum())
        print(f"{sample}: n_reads={len(arr)} total_events={total} "
              f"min={arr.min()} median={np.median(arr):.0f} mean={arr.mean():.1f} "
              f"max={arr.max()}")
        for rid, L in d.items():
            len_rows.append({"read_id": rid, "sample": sample, "n_allafterm_events": L})
        top5 = sorted(d.items(), key=lambda kv: -kv[1])[:5]
        for rid, L in top5:
            print(f"    top: {rid}  len={L}  ({100*L/total:.1f}% of {sample}'s total rescued events)")
    pd.DataFrame(len_rows).to_csv(os.path.join(args.outdir, "rescue_length_recheck.tsv"), sep="\t", index=False)

    # ---------------- 2. coverage-decay crossings --------------------------
    print("\n=== Coverage-decay crossings (event offset at which unique-read "
          "count first drops below X% of the sample's starting count) ===")
    cov_rows = []
    crossings_by_sample = {}
    for sample, d in (("gem", gem_lengths), ("control", control_lengths)):
        lengths = np.array(list(d.values()))
        crossings, offsets, n_at_offset, n0 = coverage_crossings(lengths, max_offset=args.zoom_max_offset)
        crossings_by_sample[sample] = crossings
        print(f"{sample} (n0={n0}): " + ", ".join(
            f"{int(f*100)}%->offset {crossings[f]}" for f in COVERAGE_FRACTIONS))
        for x, n in zip(offsets, n_at_offset):
            cov_rows.append({"sample": sample, "event_offset_after_last_M": int(x), "n_unique_reads": int(n)})
    pd.DataFrame(cov_rows).to_csv(os.path.join(args.outdir, "coverage_decay_zoomed.tsv"), sep="\t", index=False)

    # ---------------- 3. read-level MWU across coverage-backed windows -----
    print("\n=== Read-level Mann-Whitney U (per-read median / mean), by window ===")
    stat_rows = []
    all_metric_rows = []
    for window in STAT_WINDOWS:
        rows = read_level_metrics(gem_allafterm, "gem", window) + read_level_metrics(control_allafterm, "control", window)
        all_metric_rows += rows
        df = pd.DataFrame(rows)
        for metric in ("median", "mean"):
            g = df.loc[df["sample"] == "gem", metric].to_numpy()
            c = df.loc[df["sample"] == "control", metric].to_numpy()
            row = mwu_row(window, metric, g, c)
            stat_rows.append(row)
            print(f"  window={window:>5} metric={metric:<6} n_gem={row['n_gem']:<4} "
                  f"n_control={row['n_control']:<4} U={row['U']:.1f} p={row['p_value']:.4g} auc={row['auc']:.3f}")

    summary_df = pd.DataFrame(all_metric_rows)
    summary_df.to_csv(os.path.join(args.outdir, "rescue_window_read_level.tsv"), sep="\t", index=False)

    # ---------------- 4. full (untruncated) +/- outlier sensitivity --------
    print("\n=== Full (untruncated) per-read median/mean: all reads vs outliers excluded ===")
    KNOWN_OUTLIERS = {
        "4eeceb33-dda0-4dfc-bf1c-431b804c7623",  # control, 725,085 events
        "e1cd3406-9f89-4862-9055-ca6e242e4e69",  # control, 83,223 events
    }
    full_rows = read_level_metrics(gem_allafterm, "gem", BIG_WINDOW) + read_level_metrics(control_allafterm, "control", BIG_WINDOW)
    full_df = pd.DataFrame(full_rows)
    full_df.to_csv(os.path.join(args.outdir, "rescue_full_read_level.tsv"), sep="\t", index=False)
    for metric in ("median", "mean"):
        g = full_df.loc[full_df["sample"] == "gem", metric].to_numpy()
        c_all = full_df.loc[full_df["sample"] == "control", metric].to_numpy()
        c_excl = full_df.loc[(full_df["sample"] == "control") & (~full_df["read_id"].isin(KNOWN_OUTLIERS)), metric].to_numpy()
        row_all = mwu_row(BIG_WINDOW, metric, g, c_all, note="full, all reads")
        row_excl = mwu_row(BIG_WINDOW, metric, g, c_excl, note="full, 2 known control outliers excluded")
        stat_rows.append(row_all)
        stat_rows.append(row_excl)
        print(f"  metric={metric:<6} ALL:     n_gem={row_all['n_gem']:<4} n_control={row_all['n_control']:<4} "
              f"U={row_all['U']:.1f} p={row_all['p_value']:.4g} auc={row_all['auc']:.3f}")
        print(f"  metric={metric:<6} EXCL-2:  n_gem={row_excl['n_gem']:<4} n_control={row_excl['n_control']:<4} "
              f"U={row_excl['U']:.1f} p={row_excl['p_value']:.4g} auc={row_excl['auc']:.3f}")

    pd.DataFrame(stat_rows).to_csv(os.path.join(args.outdir, "rescue_window_stats.tsv"), sep="\t", index=False)

    # ---------------- 5. figures -------------------------------------------
    print("\n=== Figures ===")
    for window in (1000, 2000):
        plot_read_level_boxplots(summary_df, window, os.path.join(args.outdir, f"rescue_read_level_boxplots_{window}.png"))
        plot_focused_median_boxplot(summary_df, window, os.path.join(args.outdir, f"rescue_focused_median_boxplot_{window}.png"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5.5))
    cov_df = pd.DataFrame(cov_rows)
    colors = {"gem": "tab:red", "control": "tab:gray"}
    for sample in ("gem", "control"):
        sub = cov_df[cov_df["sample"] == sample].sort_values("event_offset_after_last_M")
        ax.plot(sub["event_offset_after_last_M"], sub["n_unique_reads"], color=colors[sample], label=sample, linewidth=1.5)
        n0 = sub["n_unique_reads"].iloc[0]
        for frac, ls in zip((0.5, 0.25), ("--", ":")):
            ax.axhline(frac * n0, color=colors[sample], linewidth=0.5, linestyle=ls, alpha=0.5)
    for frac in COVERAGE_FRACTIONS:
        for sample in ("gem", "control"):
            x = crossings_by_sample[sample][frac]
            if x is not None:
                ax.axvline(x, color=colors[sample], linewidth=0.4, alpha=0.3)
    ax.set_xlabel("event offset after final HMM Match (x=1 is the first rescued event)")
    ax.set_ylabel("number of unique reads still contributing")
    ax.set_title(f"Coverage decay, zoomed to the first {args.zoom_max_offset} rescued events\n"
                 "(thin vertical lines: 100/95/90/80/50% coverage crossings, per sample)")
    ax.legend()
    fig.tight_layout()
    outpath = os.path.join(args.outdir, "coverage_decay_zoomed.png")
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")

    print("\nDone.")


if __name__ == "__main__":
    main()
