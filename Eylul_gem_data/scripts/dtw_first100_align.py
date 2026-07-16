#!/usr/bin/env python3
"""Dynamic-time-warping alignment of the post-last-Match ("ALLAFTERM") event
signal, restricted to a small leading window (default: first 100 events).

Why this exists (PI's "elephant in the room"): events after the last HMM
Match have no reference coordinate at all, and different reads translocate
through the pore at different speeds. So "event offset 7 after the last
Match" in one read and in another read do not necessarily correspond to the
same physical base position -- naively averaging by fixed offset (as
Parts A/B of plot_allafterm_analysis.py do) can blur or wash out a real,
localised signal, or can even manufacture a spurious one. Dynamic time
warping (DTW) addresses this by finding a monotonic correspondence between
two signals based on where their *shapes* actually line up, rather than
assuming a fixed events-per-base rate.

IMPORTANT -- DTW is inherently pairwise (it aligns exactly two signals to
each other; there's no native N-way DTW). To get every read onto one shared,
comparable axis (needed to compute a per-position median/IQR band across
hundreds of reads), something has to be fixed as the common target that
every other read is aligned to. Two ways to pick that target:
  (a) an arbitrary real read (this project's first attempt: the pooled read
      closest, by L1, to the per-offset median) -- simple, but the target
      is exactly as noisy as any other single read.
  (b) a refined consensus ("barycenter"): start from the pooled median
      curve, align every read to it, recompute the target as the median of
      the *resampled* reads, and repeat a few times (DTW Barycenter
      Averaging, Petitjean et al. 2011). This project now uses (b).
Either way, the actual DTW step run against each individual read is the
same pairwise operation -- what changed is only what that fixed target is.

CAUGHT AND FIXED BUG (2026-07-16): the first version of this script used
plain, *unconstrained* DTW (no limit on how far the warping path can stray
from the diagonal). A control test -- DTW-aligning temporally SHUFFLED
(i.e. pure noise, real values but scrambled order) reads to the reference --
showed shuffled noise reached the same or higher post-DTW correlation to the
reference as genuinely-ordered real reads (e.g. one read: 0.55 real vs 0.71
shuffled). That means unconstrained DTW was warping *anything* to resemble
the reference; the near-perfect gem/control overlap it produced was an
artifact, not a real result. Fixed with a Sakoe-Chiba band: the warping path
may not deviate from the diagonal by more than `--band-radius` events. This
also forces the alignment to stay physically plausible (dwell-time
variability is a local speed change, not free reordering across the whole
window) and is much cheaper to compute. `main()` now runs the same
real-vs-shuffled control test automatically every time and prints a
pass/fail verdict, rather than requiring a one-off manual check.

Deliberately restricted to a short leading window (ALLAFTERM only -- no
Source 1/2 -- ~100 events by default): even banded, DTW cost is O(n*band)
per pair of reads, and doing this across the full ALLAFTERM range (up to
~725,000 events for the longest control read) is computationally infeasible.

Pipeline:
  1. Re-parse the full .align files once (cached to a TSV afterwards).
  2. Build a barycenter reference: init at the pooled (gem+control) per-
     offset median curve, then iterate {banded-DTW-align every read to the
     current reference -> resample -> median of resampled = new reference}
     for `--barycenter-iters` rounds.
  3. Run the real-vs-shuffled-noise control test against the FINAL
     reference and band radius; print PASS/FAIL. Figures are annotated with
     the verdict so it travels with the plot, not just the log.
  4. Final banded DTW pass: align every read to the final reference,
     resample onto its grid, band-plot (median +/- IQR) gem vs control,
     side by side with the naive (fixed-offset, no warping) band.
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
from plot_allafterm_analysis import load_blocks

MAX_WINDOW_WITHOUT_FORCE = 300


def extract_first_n_matrix(gem_path, control_path, window, cache_path):
    if cache_path and os.path.exists(cache_path):
        print(f"Loading cached first-{window}-event matrix from {cache_path}")
        df = pd.read_csv(cache_path, sep="\t")
        if int(df["offset"].max()) + 1 >= window and set(df["sample"].unique()) <= {"gem", "control"}:
            return df
        print("  cache does not cover the requested window -- re-extracting.")

    print("Re-parsing full .align files to extract the first "
          f"{window} ALLAFTERM events per unique read (this is the slow part)...")
    t0 = time.time()
    main_gem, _, _, _ = load_blocks(gem_path)
    main_control, _, _, _ = load_blocks(control_path)
    rows = []
    for sample, main in (("gem", main_gem), ("control", main_control)):
        for read_id, (_, allafterm, _has_s1) in main.items():
            w = allafterm[:window]
            if w.size < window:
                continue  # keep the matrix rectangular; report any dropped reads below
            for offset, val in enumerate(w):
                rows.append({"read_id": read_id, "sample": sample, "offset": offset, "value": float(val)})
    df = pd.DataFrame(rows)
    print(f"  extracted {df['read_id'].nunique()} reads x {window} events in {time.time()-t0:.1f}s")
    if cache_path:
        df.to_csv(cache_path, sep="\t", index=False)
        print(f"  cached to {cache_path}")
    return df


def to_wide(df, window):
    """Return {read_id: (sample, np.array length window)}."""
    out = {}
    for read_id, g in df.groupby("read_id"):
        g = g.sort_values("offset")
        if len(g) != window:
            continue
        out[read_id] = (g["sample"].iloc[0], g["value"].to_numpy())
    return out


def dtw_align(query, ref, band_radius=None, step_penalty=0.0):
    """Sakoe-Chiba-banded DTW under absolute-difference cost. `band_radius`
    caps how far the warping path may stray from the diagonal (|i-j| <=
    band_radius); None means unconstrained (kept only for the control test /
    documentation -- main() always passes a finite radius).

    `step_penalty` adds a fixed cost to non-diagonal (horizontal/vertical)
    steps -- i.e. mapping more than one query event onto the same reference
    event or vice versa -- exactly like a gap penalty in sequence alignment.
    Banding alone still leaves the path free to repeatedly hunt, within the
    band, for whichever nearby reference value happens to be closest at
    each step; that many-to-one flexibility is itself enough for pure noise
    to fake a decent match (found empirically: even the best band_radius
    alone left only a marginal real-vs-shuffled margin). The step penalty
    makes that hunting cost something, so a non-diagonal step is only taken
    when it reduces amplitude mismatch by more than the penalty -- pushing
    the path back towards proportional (~1 query event per ~1 reference
    event) unless the data actually supports local stretching.

    Returns (path, total_cost); path is a list of (query_idx, ref_idx)
    pairs, monotonic non-decreasing in both, from (0,0) to (n-1, m-1)."""
    n, m = len(query), len(ref)
    if band_radius is None:
        band_radius = max(n, m)
    D = np.full((n + 1, m + 1), np.inf)
    D[0, 0] = 0.0
    for i in range(1, n + 1):
        j_lo = max(1, i - band_radius)
        j_hi = min(m, i + band_radius)
        prev = D[i - 1]
        cur = D[i]
        qi = query[i - 1]
        for j in range(j_lo, j_hi + 1):
            c = abs(qi - ref[j - 1])
            cur[j] = c + min(prev[j] + step_penalty, cur[j - 1] + step_penalty, prev[j - 1])
    # backtrace -- neighbours outside the band are still +inf, so argmin
    # naturally stays inside it
    path = []
    i, j = n, m
    while i > 0 and j > 0:
        path.append((i - 1, j - 1))
        choices = (D[i - 1, j - 1], D[i - 1, j], D[i, j - 1])
        move = int(np.argmin(choices))
        if move == 0:
            i, j = i - 1, j - 1
        elif move == 1:
            i, j = i - 1, j
        else:
            i, j = i, j - 1
    path.reverse()
    return path, float(D[n, m])


def resample_onto_reference(query, ref_len, path):
    buckets = [[] for _ in range(ref_len)]
    for qi, ri in path:
        buckets[ri].append(query[qi])
    return np.array([np.mean(b) for b in buckets])


def build_barycenter_reference(read_arrays, window, band_radius, n_iters, step_penalty=0.0):
    """Returns the final reference curve, refined from the pooled median by
    `n_iters` rounds of {banded-DTW-align every read to current reference,
    resample, take the median of the resampled reads as the next
    reference}. This is DTW Barycenter Averaging (Petitjean et al. 2011),
    used here instead of one arbitrary real read specifically because a
    single real read's own noise is what let unconstrained DTW produce a
    spurious "match" (see module docstring) -- a median-of-many consensus
    is far harder to chase with spurious warping."""
    ids = list(read_arrays.keys())
    mat = np.stack([read_arrays[rid][1] for rid in ids])
    reference = np.median(mat, axis=0)
    print(f"\nBuilding DTW-barycenter reference (band_radius={band_radius}, step_penalty={step_penalty}, "
          f"{n_iters} iterations):")
    for it in range(n_iters):
        resampled = np.empty_like(mat)
        for k in range(mat.shape[0]):
            path, _ = dtw_align(mat[k], reference, band_radius, step_penalty)
            resampled[k] = resample_onto_reference(mat[k], window, path)
        new_reference = np.median(resampled, axis=0)
        shift = float(np.mean(np.abs(new_reference - reference)))
        print(f"  iter {it + 1}: mean |change| in reference vs previous iteration = {shift:.4f}")
        reference = new_reference
    return reference


def shuffle_control_validation(read_arrays, reference, band_radius, window, seed=0, n_test=20, step_penalty=0.0):
    """The check that caught the original bug: for n_test randomly-chosen
    real reads, DTW-align both the real read and a temporally-shuffled
    (pure noise, same values, scrambled order) copy of it to `reference`,
    and compare post-DTW correlation to the reference. If shuffled noise
    matches about as well as real data, the DTW setup is still pathological
    and results should not be trusted -- regardless of band_radius or
    reference choice, this must be checked, not assumed fixed."""
    rng = np.random.default_rng(seed)
    ids = list(read_arrays.keys())
    test_ids = rng.choice(len(ids), size=min(n_test, len(ids)), replace=False)
    real_corrs, shuf_corrs = [], []
    for k in test_ids:
        arr = read_arrays[ids[k]][1]
        path, _ = dtw_align(arr, reference, band_radius, step_penalty)
        warped = resample_onto_reference(arr, window, path)
        real_corrs.append(np.corrcoef(warped, reference)[0, 1])

        shuffled = arr.copy()
        rng.shuffle(shuffled)
        path_s, _ = dtw_align(shuffled, reference, band_radius, step_penalty)
        warped_s = resample_onto_reference(shuffled, window, path_s)
        shuf_corrs.append(np.corrcoef(warped_s, reference)[0, 1])
    real_corrs, shuf_corrs = np.array(real_corrs), np.array(shuf_corrs)
    margin = float(real_corrs.mean() - shuf_corrs.mean())
    passed = margin > 0.15
    print(f"\nShuffle-control validation (n={len(test_ids)} reads, band_radius={band_radius}, "
          f"step_penalty={step_penalty}):")
    print(f"  real reads:     mean post-DTW corr-to-reference = {real_corrs.mean():.3f}")
    print(f"  shuffled noise: mean post-DTW corr-to-reference = {shuf_corrs.mean():.3f}")
    verdict = "PASS" if passed else "FAIL"
    print(f"  {verdict} (margin={margin:+.3f}; need >0.15 for real reads to clearly beat noise)")
    if not passed:
        print("  WARNING: shuffled noise aligns nearly as well as real data -- the DTW result "
              "below is likely still an artifact. Do not use it as evidence of anything.")
    return passed, margin, real_corrs, shuf_corrs


def median_iqr(mat):
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
    return med, q1, q3


def mean_sd(mat):
    """Mean +/- 1 SD band -- the mean-based analogue of median_iqr's
    median +/- IQR. SD (spread of individual reads' values), not SEM
    (uncertainty of the mean), to keep the same visual meaning as IQR: both
    bands show where most individual reads' values fall, not how precisely
    the central estimate is known."""
    with np.errstate(all="ignore"):
        m = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0, ddof=1)
    return m, m - sd, m + sd


BAND_STATS = {
    "median": (median_iqr, "median +/- IQR"),
    "mean": (mean_sd, "mean +/- 1 SD"),
}


def make_band_figure(warped_df, statistic, window, band_radius, verdict_str, real_corr_mean, shuf_corr_mean,
                      margin, outdir):
    band_fn, band_label = BAND_STATS[statistic]
    band_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, (value_col, title) in zip(axes, (("raw_value", "Naive (fixed-offset, no warping)"),
                                              ("dtw_value", f"Banded DTW-aligned (radius={band_radius}) "
                                                            "onto barycenter reference"))):
        for sample, color in (("gem", "tab:red"), ("control", "tab:gray")):
            sub = warped_df[warped_df["sample"] == sample]
            piv = sub.pivot(index="read_id", columns="ref_offset", values=value_col).to_numpy()
            center, lo, hi = band_fn(piv)
            x = np.arange(1, window + 1)
            ax.plot(x, center, color=color, label=f"{sample} (n={piv.shape[0]})", linewidth=1.5)
            ax.fill_between(x, lo, hi, color=color, alpha=0.2)
            for offset in range(window):
                band_rows.append({"method": "dtw" if value_col == "dtw_value" else "naive", "statistic": statistic,
                                   "sample": sample, "ref_offset": offset + 1,
                                   "center": center[offset], "lo": lo[offset], "hi": hi[offset],
                                   "n": int(piv.shape[0])})
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("event offset after final HMM Match\n(DTW panel: position on barycenter reference's grid)")
        ax.axhline(0, color="black", linewidth=0.4, linestyle="--")
        ax.legend(fontsize=8)
    axes[0].set_ylabel(f"scaled event signal ({band_label})")
    fig.suptitle(f"First {window} post-Match (ALLAFTERM) events: naive vs banded-DTW-aligned, gem vs control "
                 f"[{statistic}]\n{verdict_str} (real corr={real_corr_mean:.2f}, shuffled corr={shuf_corr_mean:.2f}, "
                 f"margin={margin:+.2f})", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    outpath = os.path.join(outdir, f"dtw_vs_naive_first{window}_{statistic}.png")
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"wrote {outpath}")
    return band_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gem", required=True)
    ap.add_argument("--control", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--window", type=int, default=100,
                     help="ALLAFTERM events per read to DTW-align (default 100, per PI's "
                          "quadratic-cost caution; see MAX_WINDOW_WITHOUT_FORCE)")
    ap.add_argument("--band-radius", type=int, default=20,
                     help="Sakoe-Chiba band radius in events (default 20 = 20%% of the default "
                          "100-event window): the warping path cannot map a read-event more than "
                          "this many positions away from the diagonal. Prevents the pathological "
                          "any-to-any matching an unconstrained DTW is prone to on short, noisy "
                          "series (see module docstring). Physically: at ~5 events/base, this "
                          "allows roughly +/-4 bases of local dwell-time drift.")
    ap.add_argument("--barycenter-iters", type=int, default=3,
                     help="refinement rounds for the consensus reference curve (default 3)")
    ap.add_argument("--step-penalty", type=float, default=0.0,
                     help="fixed cost added to non-diagonal (many-to-one) DTW steps, like a gap "
                          "penalty in sequence alignment. Banding alone still leaves the path free "
                          "to hunt locally for whichever nearby reference value happens to match "
                          "by chance, which is itself enough for shuffled noise to fake a good "
                          "score (found empirically). A positive penalty discourages that unless "
                          "it genuinely reduces mismatch. Tune with the shuffle-control margin "
                          "printed below -- see module docstring.")
    ap.add_argument("--shuffle-test-n", type=int, default=20,
                     help="reads used in the automatic real-vs-shuffled-noise validation (default 20)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cache", default=None,
                     help="TSV path to cache the extracted first-window matrix "
                          "(avoids re-parsing the .align files on re-runs)")
    ap.add_argument("--force-large-window", action="store_true",
                     help=f"required if --window > {MAX_WINDOW_WITHOUT_FORCE}")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.window > MAX_WINDOW_WITHOUT_FORCE and not args.force_large_window:
        print(f"Refusing to run with --window={args.window} > {MAX_WINDOW_WITHOUT_FORCE} "
              "without --force-large-window: DTW cost is O(n*band) per read pair, so this "
              f"would cost ~{(args.window/100.0):.0f}x more per-read work than the default "
              "100-event window, on top of however many reads you have. Re-run with "
              "--force-large-window if you've confirmed the expected wall-clock/memory cost "
              "is acceptable.", file=sys.stderr)
        sys.exit(1)

    df = extract_first_n_matrix(args.gem, args.control, args.window, args.cache)
    read_arrays = to_wide(df, args.window)
    n_gem = sum(1 for s, _ in read_arrays.values() if s == "gem")
    n_control = sum(1 for s, _ in read_arrays.values() if s == "control")
    print(f"\nReads with a full {args.window}-event ALLAFTERM window: gem={n_gem} control={n_control}")

    reference = build_barycenter_reference(read_arrays, args.window, args.band_radius, args.barycenter_iters,
                                            step_penalty=args.step_penalty)

    passed, margin, real_corrs, shuf_corrs = shuffle_control_validation(
        read_arrays, reference, args.band_radius, args.window, seed=args.seed, n_test=args.shuffle_test_n,
        step_penalty=args.step_penalty)

    print(f"\nRunning final DTW pass: {len(read_arrays)} reads x window={args.window}, "
          f"band_radius={args.band_radius}...")
    t0 = time.time()
    warped_rows = []
    cost_rows = []
    for rid, (sample, arr) in read_arrays.items():
        path, cost = dtw_align(arr, reference, args.band_radius, args.step_penalty)
        warped = resample_onto_reference(arr, args.window, path)
        for offset in range(args.window):
            warped_rows.append({"read_id": rid, "sample": sample, "ref_offset": offset,
                                 "raw_value": float(arr[offset]), "dtw_value": float(warped[offset])})
        cost_rows.append({"read_id": rid, "sample": sample, "dtw_cost": cost, "dtw_cost_per_event": cost / args.window})
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s ({elapsed/len(read_arrays)*1000:.2f} ms/read)")

    warped_df = pd.DataFrame(warped_rows)
    cost_df = pd.DataFrame(cost_rows)
    warped_df.to_csv(os.path.join(args.outdir, "dtw_warped_matrix.tsv"), sep="\t", index=False)
    cost_df.to_csv(os.path.join(args.outdir, "dtw_cost_per_read.tsv"), sep="\t", index=False)
    pd.DataFrame({"ref_offset": np.arange(args.window), "value": reference}).to_csv(
        os.path.join(args.outdir, "dtw_barycenter_reference.tsv"), sep="\t", index=False)

    print("\nDTW alignment cost (final pass; total path cost / window) by sample:")
    print(cost_df.groupby("sample")["dtw_cost_per_event"].describe())

    verdict_str = "VALIDATION PASSED" if passed else "VALIDATION FAILED -- see log, do not trust this panel"
    all_band_rows = []
    for statistic in ("median", "mean"):
        all_band_rows += make_band_figure(warped_df, statistic, args.window, args.band_radius, verdict_str,
                                           real_corrs.mean(), shuf_corrs.mean(), margin, args.outdir)
    pd.DataFrame(all_band_rows).to_csv(os.path.join(args.outdir, "dtw_vs_naive_band_summary.tsv"), sep="\t", index=False)

    # cost distribution figure -- flags whether a handful of reads warp far more than
    # the rest (which would mean their DTW resampling is on shakier ground than most)
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    for sample, color in (("gem", "tab:red"), ("control", "tab:gray")):
        vals = cost_df.loc[cost_df["sample"] == sample, "dtw_cost_per_event"].to_numpy()
        ax2.hist(vals, bins=30, color=color, alpha=0.5, label=f"{sample} (n={len(vals)})", density=True)
    ax2.set_xlabel("DTW cost per event (total path cost / window)")
    ax2.set_ylabel("density")
    ax2.set_title(f"DTW alignment cost to barycenter reference (band_radius={args.band_radius})")
    ax2.legend()
    fig2.tight_layout()
    outpath2 = os.path.join(args.outdir, "dtw_cost_distribution.png")
    fig2.savefig(outpath2, dpi=150)
    plt.close(fig2)
    print(f"wrote {outpath2}")

    print("\nDone.")


if __name__ == "__main__":
    main()
