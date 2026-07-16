#!/usr/bin/env python3
"""Read-level Mann-Whitney U tests on the banded/regularized DTW-aligned
first-100-event matrix (dtw_warped_matrix.tsv from dtw_first100_align.py),
gem vs control, alongside the same test on the naive (unwarped) values at
the identical offsets for direct comparison.

Two kinds of test are reported, and they are NOT interchangeable:

1. PRE-REGISTERED windows (0..W for W in PREREGISTERED_WINDOWS): chosen
   from the geometry worked out earlier in this project (k=9 last-Match
   9-mer centred at base 611, gemcitabine at base 614 -- already inside
   that 9-mer -- fully clears the pore's sensing window by ~base 618, i.e.
   ~7 bases / ~35 events at ~5 events/base past the anchor) BEFORE looking
   at where the DTW-aligned curves happened to diverge most. These are
   legitimate hypothesis tests.

2. The "peak-divergence" offsets (9, 13, 14, 18 -- picked previously
   because they showed the largest |gem-control| median gap in this same
   DTW-aligned dataset): reported for completeness, but this selection is
   circular (offsets chosen to maximise the apparent effect, then tested on
   the same data) and must NOT be read as confirmatory. Flagged as such in
   both stdout and the output TSV.
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

PREREGISTERED_WINDOWS = (10, 20, 35, 50, 100)
PEAK_DIVERGENCE_OFFSETS = (9, 13, 14, 18)  # 0-based, from prior band_summary analysis


def mwu_row(window_label, method, metric, gem_vals, control_vals, note=""):
    gem_vals = np.asarray(gem_vals, dtype=float)
    control_vals = np.asarray(control_vals, dtype=float)
    gem_vals = gem_vals[~np.isnan(gem_vals)]
    control_vals = control_vals[~np.isnan(control_vals)]
    row = {"window": window_label, "method": method, "metric": metric,
           "n_gem": len(gem_vals), "n_control": len(control_vals),
           "U": np.nan, "p_value": np.nan, "auc": np.nan, "note": note}
    if len(gem_vals) and len(control_vals):
        U, p = mannwhitneyu(gem_vals, control_vals, alternative="two-sided")
        row["U"] = float(U)
        row["p_value"] = float(p)
        row["auc"] = float(U / (len(gem_vals) * len(control_vals)))
    return row


def per_read_summary(df, value_col, offsets):
    sub = df[df["ref_offset"].isin(offsets)]
    g = sub.groupby(["read_id", "sample"])[value_col].agg(["median", "mean"]).reset_index()
    return g


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--warped-matrix", required=True, help="dtw_warped_matrix.tsv from dtw_first100_align.py")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.warped_matrix, sep="\t")
    n_gem = df.loc[df["sample"] == "gem", "read_id"].nunique()
    n_control = df.loc[df["sample"] == "control", "read_id"].nunique()
    print(f"Loaded {args.warped_matrix}: gem={n_gem} unique reads, control={n_control} unique reads, "
          f"{df['ref_offset'].max()+1} offsets (0-based, 0..{df['ref_offset'].max()})")

    rows = []
    print("\n=== Pre-registered windows (chosen from geometry, before seeing this data) ===")
    for window in PREREGISTERED_WINDOWS:
        offsets = range(window)  # 0-based, first `window` events
        for value_col, method in (("dtw_value", "dtw"), ("raw_value", "naive")):
            summ = per_read_summary(df, value_col, offsets)
            for metric in ("median", "mean"):
                g = summ.loc[summ["sample"] == "gem", metric].to_numpy()
                c = summ.loc[summ["sample"] == "control", metric].to_numpy()
                row = mwu_row(f"0-{window}", method, metric, g, c)
                rows.append(row)
                print(f"  window=0-{window:<4} method={method:<6} metric={metric:<6} "
                      f"n_gem={row['n_gem']:<4} n_control={row['n_control']:<4} "
                      f"U={row['U']:.1f} p={row['p_value']:.4g} auc={row['auc']:.3f}")

    n_tests = len(rows)
    bonf = 0.05 / n_tests
    print(f"\n{n_tests} pre-registered tests run; Bonferroni threshold = 0.05/{n_tests} = {bonf:.4g}")
    sig = [r for r in rows if r["p_value"] < bonf]
    if sig:
        print("  surviving Bonferroni correction:")
        for r in sig:
            print(f"    window=0-{r['window'].split('-')[1]} method={r['method']} metric={r['metric']} p={r['p_value']:.4g}")
    else:
        print("  none survive Bonferroni correction.")

    print("\n=== Peak-divergence offsets {9,13,14,18} -- EXPLORATORY / CIRCULAR, NOT CONFIRMATORY ===")
    print("  (these offsets were chosen because they showed the largest gem-vs-control gap in this")
    print("   same DTW-aligned dataset -- testing on them again is selection bias, reported for the")
    print("   record only, not as evidence.)")
    for value_col, method in (("dtw_value", "dtw"), ("raw_value", "naive")):
        summ = per_read_summary(df, value_col, PEAK_DIVERGENCE_OFFSETS)
        for metric in ("median", "mean"):
            g = summ.loc[summ["sample"] == "gem", metric].to_numpy()
            c = summ.loc[summ["sample"] == "control", metric].to_numpy()
            row = mwu_row("peak_offsets_9_13_14_18", method, metric, g, c,
                           note="EXPLORATORY/CIRCULAR -- offsets selected post-hoc from this same data")
            rows.append(row)
            print(f"  method={method:<6} metric={metric:<6} n_gem={row['n_gem']:<4} n_control={row['n_control']:<4} "
                  f"U={row['U']:.1f} p={row['p_value']:.4g} auc={row['auc']:.3f}")

    out = pd.DataFrame(rows)
    outpath = os.path.join(args.outdir, "dtw_read_level_stats.tsv")
    out.to_csv(outpath, sep="\t", index=False)
    print(f"\nwrote {outpath}")


if __name__ == "__main__":
    main()
