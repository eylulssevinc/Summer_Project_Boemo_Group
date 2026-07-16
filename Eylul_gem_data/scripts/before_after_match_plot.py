#!/usr/bin/env python3
"""Extends the naive-vs-DTW ALLAFTERM (post-last-Match) figures with a
"before" panel: the last `--before-window` genuinely reference-ALIGNED
(Match-labelled) events immediately preceding the anchor, for direct visual
comparison against the post-anchor region.

WHY THIS NEEDS NO NEW C++ OUTPUT OR RE-RUN, AND NO SEPARATE INDEX:
the "before" side is already fully determined by output this project's
alignment.cpp has been printing all along -- the ordinary per-position
Match lines (`event_coord \t kmerRef \t scaledEvent \t kmerStrand \t
meanStd.first`) emitted for every window throughout each read, well before
the terminal window's TAIL_SOURCE*/NEXT100EVT/ALLAFTERM lines. Each such
line is one RAW SAMPLE of an event (a Match-labelled event's `.raw[]`
expanded one line per sample, all sharing the same `event_coord`); grouping
consecutive same-`event_coord` lines and averaging their `scaledEvent`
recovers exactly the event's true `.mean`, scaled the identical way as
ALLAFTERM (scaling is affine, so it commutes with averaging). Insertion
lines are already distinguishable (an all-'N' `kmerStrand`, see
alignment.cpp) and are excluded here since they don't correspond to a
confirmed reference base and would break position correspondence across
reads. `event_coord` is monotonic (strictly increasing for forward-strand
reads, strictly decreasing for reverse-strand -- either way, never
revisited) throughout a whole read because alignment.cpp's outer loop only
ever advances `reference_index` forward, so consecutive-line grouping is
safe without needing to track any additional per-event raw index.

Net effect: unlike the post-Match side (which needed DTW because those
events have no reference anchor and reads translocate at different speeds),
the pre-anchor side is genuinely, individually reference-locked across
reads already -- by construction of the HMM alignment itself. No warping is
applied or needed here; "naive" is not an approximation on this side, it is
the ground truth.

CAVEAT DISCOVERED WHILE BUILDING THIS: the reference amplicon here is only
~600-610bp total (the last-Match anchor sits at ~base 611), so asking for
500 pre-anchor events uses nearly the entire aligned region. Not every read
aligns that far: only 28/50 gem reads and the large majority of control
reads have >=500 distinct aligned positions before their own anchor. Reads
below the requested window are dropped from the "before" panel (reported
below), so its n is smaller than the "after" panel's.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def iter_before_after(path):
    """Yield (read_id, before_arr, allafterm_arr) per '>'-delimited block
    that has at least one ALLAFTERM line (matching load_blocks' gate
    elsewhere in this project). before_arr: ALL genuinely Match-aligned
    per-event means (scaled) for the read, in file order (index 0 =
    earliest in the read, last index = immediately before the anchor) --
    callers slice the last N themselves. Insertion lines (all-'N'
    kmerStrand) are excluded."""
    read_id = None
    match_means = []
    cur_coord = None
    cur_vals = []
    allafterm = []

    def _flush():
        nonlocal cur_coord, cur_vals
        if cur_coord is not None and cur_vals:
            match_means.append(sum(cur_vals) / len(cur_vals))
        cur_coord, cur_vals = None, []

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if read_id is not None:
                    _flush()
                    if allafterm:
                        yield read_id, np.asarray(match_means, dtype=float), np.asarray(allafterm, dtype=float)
                parts = line[1:].split()
                read_id = parts[0] if parts else None
                match_means = []
                cur_coord, cur_vals = None, []
                allafterm = []
            elif line.startswith("ALLAFTERM"):
                allafterm.append(float(line.rstrip("\n").split("\t")[2]))
            elif line.startswith(("NEXT100EVT", "TAIL_SOURCE")):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                coord, kmer_strand, scaled = fields[0], fields[3], fields[2]
                if kmer_strand and kmer_strand.count("N") == len(kmer_strand):
                    continue  # insertion placeholder, not a confirmed reference base
                if coord != cur_coord:
                    _flush()
                    cur_coord = coord
                cur_vals.append(float(scaled))
    if read_id is not None:
        _flush()
        if allafterm:
            yield read_id, np.asarray(match_means, dtype=float), np.asarray(allafterm, dtype=float)


def extract_before_matrix(gem_path, control_path, before_window, cache_path):
    if cache_path and os.path.exists(cache_path):
        df = pd.read_csv(cache_path, sep="\t")
        if int(df["offset"].max()) + 1 >= before_window:
            print(f"Loading cached before-matrix from {cache_path}")
            return df
        print(f"  cache at {cache_path} only covers offsets up to {df['offset'].max()} "
              f"(< requested {before_window}) -- re-extracting instead of reusing it.")

    print(f"Re-parsing full .align files to extract the last {before_window} pre-anchor "
          "aligned (Match) events per unique read...")
    t0 = time.time()
    rows = []
    n_total = {"gem": 0, "control": 0}
    n_qualifying = {"gem": 0, "control": 0}
    for sample, path in (("gem", gem_path), ("control", control_path)):
        best = {}  # read_id -> before_arr, longer-ALLAFTERM-wins dedup to match load_blocks elsewhere
        best_allafterm_len = {}
        for read_id, before_arr, allafterm_arr in iter_before_after(path):
            n_total[sample] += 1
            if read_id not in best or allafterm_arr.size > best_allafterm_len[read_id]:
                best[read_id] = before_arr
                best_allafterm_len[read_id] = allafterm_arr.size
        for read_id, before_arr in best.items():
            if before_arr.size < before_window:
                continue
            n_qualifying[sample] += 1
            last_n = before_arr[-before_window:]  # closest-to-anchor last
            for i, val in enumerate(last_n[::-1]):  # offset 0 = closest to anchor
                rows.append({"read_id": read_id, "sample": sample, "offset": i, "value": float(val)})
    df = pd.DataFrame(rows)
    print(f"  gem: {n_total['gem']} qualifying (ALLAFTERM-bearing) blocks -> "
          f"{n_qualifying['gem']} with >= {before_window} pre-anchor aligned events")
    print(f"  control: {n_total['control']} qualifying blocks -> "
          f"{n_qualifying['control']} with >= {before_window} pre-anchor aligned events")
    print(f"  extracted in {time.time()-t0:.1f}s")
    if cache_path:
        df.to_csv(cache_path, sep="\t", index=False)
        print(f"  cached to {cache_path}")
    return df


def median_iqr(mat):
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
    return med, q1, q3


def mean_sd(mat):
    with np.errstate(all="ignore"):
        m = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0, ddof=1)
    return m, m - sd, m + sd


BAND_STATS = {"median": (median_iqr, "median +/- IQR"), "mean": (mean_sd, "mean +/- 1 SD")}


def to_matrix(df, value_col, window, offset_col="offset"):
    piv = df.pivot(index="read_id", columns=offset_col, values=value_col)
    return piv.reindex(columns=range(window)).to_numpy()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gem", required=True)
    ap.add_argument("--control", required=True)
    ap.add_argument("--before-window", type=int, default=500)
    ap.add_argument("--after-warped-matrix", required=True,
                     help="dtw_warped_matrix.tsv from dtw_first100_align.py, at the same window "
                          "as --before-window (so the after-side is already banded-DTW-aligned)")
    ap.add_argument("--before-cache", default=None)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    before_df = extract_before_matrix(args.gem, args.control, args.before_window, args.before_cache)
    after_df = pd.read_csv(args.after_warped_matrix, sep="\t")
    after_window = int(after_df["ref_offset"].max()) + 1
    if after_window != args.before_window:
        print(f"NOTE: after-side window ({after_window}) != before-window ({args.before_window}); "
              "the two sides of the figure will have different x-extents.", file=sys.stderr)

    before_ids_by_sample = {s: set(before_df.loc[before_df["sample"] == s, "read_id"]) for s in ("gem", "control")}
    after_ids_by_sample = {s: set(after_df.loc[after_df["sample"] == s, "read_id"]) for s in ("gem", "control")}
    common_ids = set()
    for sample in ("gem", "control"):
        common = before_ids_by_sample[sample] & after_ids_by_sample[sample]
        common_ids |= common
        print(f"{sample}: {len(before_ids_by_sample[sample])} with >={args.before_window} pre-anchor aligned "
              f"events, {len(after_ids_by_sample[sample])} with a full {after_window}-event after-side, "
              f"{len(common)} with BOTH (this n is what's plotted below)")
    before_df = before_df[before_df["read_id"].isin(common_ids)]
    after_df = after_df[after_df["read_id"].isin(common_ids)]

    for statistic in ("median", "mean"):
        band_fn, band_label = BAND_STATS[statistic]
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
        panel_defs = [
            ("Fully naive (before: aligned positions; after: fixed raw-event offset, no warping)", "raw_value"),
            (f"Before: aligned positions (ground truth) + After: banded DTW-aligned", "dtw_value"),
        ]
        for ax, (title, after_col) in zip(axes, panel_defs):
            for sample, color in (("gem", "tab:red"), ("control", "tab:gray")):
                b_sub = before_df[before_df["sample"] == sample]
                b_mat = to_matrix(b_sub, "value", args.before_window)
                b_center, b_lo, b_hi = band_fn(b_mat)
                x_before = -np.arange(args.before_window, 0, -1)  # ..., -2, -1
                ax.plot(x_before, b_center, color=color, linewidth=1.5)
                ax.fill_between(x_before, b_lo, b_hi, color=color, alpha=0.2)

                a_sub = after_df[after_df["sample"] == sample]
                a_mat = to_matrix(a_sub, after_col, after_window, offset_col="ref_offset")
                a_center, a_lo, a_hi = band_fn(a_mat)
                x_after = np.arange(1, after_window + 1)
                n_reads = b_mat.shape[0]
                ax.plot(x_after, a_center, color=color, linewidth=1.5, label=f"{sample} (n={n_reads})")
                ax.fill_between(x_after, a_lo, a_hi, color=color, alpha=0.2)
            ax.axvline(0, color="black", linewidth=0.8, linestyle="-")
            ax.axhline(0, color="black", linewidth=0.4, linestyle="--")
            ax.set_title(title, fontsize=9)
            ax.set_xlabel("event offset relative to final HMM Match\n"
                          "(negative: aligned reference position; positive: raw post-Match event)")
            ax.legend(fontsize=8)
        axes[0].set_ylabel(f"scaled event signal ({band_label})")
        fig.suptitle(f"Before (last {args.before_window} aligned events) vs after (first {after_window} "
                     f"post-Match events), gem vs control [{statistic}]", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        outpath = os.path.join(args.outdir, f"before_after_match_{args.before_window}_{statistic}.png")
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
        print(f"wrote {outpath}")

    print("\nDone.")


if __name__ == "__main__":
    main()
