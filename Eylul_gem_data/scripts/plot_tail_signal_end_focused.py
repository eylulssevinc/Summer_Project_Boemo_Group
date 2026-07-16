#!/usr/bin/env python3
"""
End-focused visualization of the recovered 3' tail signal produced by DNAscent's
captureTailSignal feature (see ../../TAIL_SIGNAL_CAPTURE.md and
../../PROJECT_OVERVIEW.md for what this data is and how it was generated).

Why this script exists, and how it differs from plot_tail_signal.py
---------------------------------------------------------------------
plot_tail_signal.py's headline comparison plot anchors mostly on the cutoff
point (x=0 = first recovered TAIL sample). That's a natural anchor for "what
did the aligner give up on first", but it's the WRONG anchor for this project's
actual biology: gemcitabine is expected at the second-to-last base of the 615 bp
ligation-product reference, i.e. right at the true physical 3' end of the
molecule -- not necessarily right after wherever the aligner happened to stop.
Tail length varies enormously between reads (roughly 5.6k-59.5k raw samples in
the gem data), so "100 samples after the cutoff" means something completely
different physically from read to read, while "100 samples before the physical
end" is the same physical distance from the one point that actually matters
biologically, for every read. This script anchors everything from that
physical 3' end instead, and adds read-level (not just sample-pooled) summaries
so a handful of reads with unusually large tails can't quietly dominate a plot.

Run with the gmconda_py368 conda env (matplotlib/numpy/pandas):

    /opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python \\
        plot_tail_signal_end_focused.py \\
        --gem     ../slurm_full_run/gem/gem.full.align \\
        --control ../slurm_full_run/control/control.full.align \\
        --outdir  ../slurm_full_run/figures_end_focused

IMPORTANT SCIENTIFIC CAUTIONS (see also the printed terminal summary):
  - End-anchored raw-sample positions are NOT base positions. x=0 is the final
    raw SAMPLE of the read, not "the final base", and x=1 is not "the
    second-to-last base". A single base/k-mer dwell in the pore can produce
    anywhere from a handful to hundreds of raw samples, and that dwell time is
    itself variable and not something this script corrects for.
  - Because gemcitabine is expected right at the molecule's physical 3' end,
    anchoring from that end is more biologically relevant than anchoring from
    the aligner's cutoff point -- but "more relevant" is not "precise": we still
    do not know exactly which raw samples correspond to the modified base.
  - Do NOT read a Gem-vs-control difference off visual band overlap/non-overlap
    alone. Sample sizes here are modest (50 gem / 347 control qualifying reads),
    bands are pooled across reads of very different tail length, and no
    correction for multiple comparisons is applied anywhere in this script. Use
    the read-level summary TSVs for anything beyond a first, exploratory look,
    and treat any apparent difference as a hypothesis to test properly (e.g.
    with an explicit read-level statistical test, or the changepoint analysis
    recommended in TAIL_SIGNAL_CAPTURE.md), not a conclusion.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # no display on a compute/login node
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def iter_read_tails(path):
    """Yield (read_id, tail_values) for every read in a DNAscent .align file
    that has at least one TAIL line.

    read_id is parsed from the '>' header line (first whitespace-delimited
    token, with the leading '>' stripped) -- this is the biological unit for
    everything in this script: TAIL lines are raw *samples* within one read's
    recovered tail, not independent observations, so every downstream
    computation groups by read_id before doing anything statistical.

    tail_values is returned in the file's own chronological order (index 0 =
    the first recovered sample right after the point the normal, reference-
    anchored aligner stopped assigning this read to the genome). Use
    end_align_tail() to re-anchor it from the physical 3' end instead.
    """
    read_id, tail = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if read_id is not None and tail:
                    yield read_id, np.asarray(tail, dtype=float)
                parts = line[1:].split()
                read_id = parts[0] if parts else None
                tail = []
            elif line.startswith("TAIL"):
                tail.append(float(line.rstrip("\n").split("\t")[2]))
    if read_id is not None and tail:
        yield read_id, np.asarray(tail, dtype=float)


def load_tail_dict(path):
    """Return an ordered {read_id: tail_array} dict for every qualifying read
    in path (path may be None, in which case an empty dict is returned).

    IMPORTANT: the same physical read (readID) can produce more than one
    qualifying alignment record/block in a DNAscent .align file -- confirmed
    in the control data, where 25 readIDs each had 2-3 separate qualifying
    blocks (secondary/supplementary alignments of the same molecule against
    this short 615 bp reference, e.g. one block '0 606 fwd' and another
    '0 615 fwd' for the same readID -- only the latter satisfies the
    captureTailSignal gate, but a molecule can occasionally satisfy it twice
    via two distinct records). Counting blocks instead of unique reads is
    exactly the "treating TAIL lines/blocks as independent observations"
    pitfall this project explicitly wants to avoid -- the biological unit is
    the READ, so this function deduplicates by readID, keeping the record
    with the LONGEST recovered tail (the most complete/informative capture of
    that molecule's post-cutoff signal) wherever a collision occurs, and
    prints how many collisions it resolved so this is never silent.
    """
    tails_by_read = {}
    n_blocks = 0
    n_collisions = 0
    if path is None or not os.path.exists(path):
        return tails_by_read
    for read_id, tail in iter_read_tails(path):
        n_blocks += 1
        if read_id in tails_by_read:
            n_collisions += 1
            if len(tail) <= len(tails_by_read[read_id]):
                continue
        tails_by_read[read_id] = tail
    if n_collisions:
        print(f"  [{os.path.basename(path)}] {n_blocks} qualifying alignment blocks -> "
              f"{len(tails_by_read)} unique reads ({n_collisions} duplicate-readID block(s) "
              f"collapsed, longest tail per read kept)")
    return tails_by_read


def extract_regular_end_signal(path, read_ids_of_interest, coord_lo, coord_hi):
    """Single extra pass through the .align file collecting, for each read_id
    in read_ids_of_interest (expected to be the already-known qualifying tail
    reads -- there's no point scanning for reads we won't report on), every
    *regular* (non-TAIL) line's scaled-signal value whose first field (the
    DNAscent reference coordinate) falls in [coord_lo, coord_hi] inclusive.

    Coordinate convention note: this is the literal 0-based coordinate as it
    appears in column 1 of the .align file (e.g. a read with refStart=0 has
    its first centered-9mer coordinate at 4, since k//2=4 -- see
    PROJECT_OVERVIEW.md section 3.7). This is one less than the 1-based
    genomic positions used in the older gem_forward_600_end_median_iqr.ipynb
    analysis (which used genomic positions 600-611 = align coordinates
    599-610). The default range here (600-611) is chosen only to land in the
    same tail-end stretch of the reference immediately upstream of where TAIL
    signal takes over; exact boundary alignment with the older notebook does
    not matter for this read-level, order-of-magnitude comparison.

    Returns {read_id: [scaled_signal, ...]} -- the list may be empty for a
    read if no regular line happened to fall in the requested coordinate
    window (possible if the read has an indel there, though unlikely for
    reads that are already known to reach the true reference end).
    """
    wanted = set(read_ids_of_interest)
    out = {rid: [] for rid in wanted}
    if path is None or not os.path.exists(path) or not wanted:
        return out

    current_id = None
    collecting = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].split()
                current_id = parts[0] if parts else None
                collecting = current_id in wanted
                continue
            if not collecting or line.startswith("TAIL"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                coord = int(fields[0])
            except ValueError:
                continue
            if coord_lo <= coord <= coord_hi:
                out[current_id].append(float(fields[2]))
    return out


# ---------------------------------------------------------------------------
# Core end-anchoring logic
# ---------------------------------------------------------------------------

def end_align_tail(tail):
    """Reverse a chronological tail array so index 0 = the final raw sample of
    the read (the physical 3' terminus), index 1 = one sample before that,
    etc. This is the ONE anchor point that is genuinely comparable across
    reads regardless of how much signal got trimmed before it -- tail length
    varies ~10x across reads, so a fixed cutoff-relative position does not
    mean the same physical thing from read to read, but a fixed
    end-relative position always does."""
    return np.asarray(tail)[::-1]


def _end_matrix(tails, max_lag):
    """Stack end-aligned tails into an (n_reads x max_lag) matrix, NaN where
    a read's tail doesn't reach that far from the end."""
    n = len(tails)
    mat = np.full((n, max_lag), np.nan)
    for i, t in enumerate(tails):
        et = end_align_tail(t)
        L = min(len(et), max_lag)
        if L == 0:
            continue
        mat[i, :L] = et[:L]
    return mat


def compute_end_median_iqr(tails, max_lag, min_n=5):
    """Per-position median/Q1/Q3 across reads, end-anchored, truncated at the
    last position where at least min_n reads still contribute.

    Returns (x, median, q1, q3, n), each an array of the same (possibly zero)
    length. Truncating rather than plotting a band built from 1-2 reads at an
    extreme lag avoids presenting a misleading "confidence band" out where the
    sample size is too small to mean anything.
    """
    if not tails:
        return (np.array([]),) * 5
    mat = _end_matrix(tails, max_lag)
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
    n = np.sum(~np.isnan(mat), axis=0)
    valid = n >= min_n
    if not np.any(valid):
        return (np.array([]),) * 5
    last = int(np.max(np.where(valid)))
    x = np.arange(last + 1)
    return x, med[: last + 1], q1[: last + 1], q3[: last + 1], n[: last + 1]


def reads_reaching(tails, lag):
    """How many reads have a tail at least `lag` raw samples long -- i.e. how
    many reads still contribute at the far edge of an end-anchored window of
    that size, before any min_n truncation is applied."""
    return sum(1 for t in tails if len(t) >= lag)


# ---------------------------------------------------------------------------
# Read-level summaries (the biological unit is the READ, not the TAIL line)
# ---------------------------------------------------------------------------

def compute_read_level_summary(tails_by_read, sample, max_lags=(100, 500, 1000, 3000)):
    """One row per read. For each window size in max_lags, summarize the last
    `lag` raw samples of that read's tail (tail[-lag:] -- if the tail is
    shorter than `lag`, this naturally (via numpy negative-slice semantics)
    just uses the whole tail; it is never padded or extrapolated). Also
    reports, for the 500-sample window specifically, what fraction of samples
    sit above +1 / below -1 in scaled-signal units -- a crude, threshold-based
    way to ask "does this read's near-end signal spend unusual amounts of time
    far from the pore-model baseline" without committing to any particular
    changepoint model.
    """
    rows = []
    for read_id, tail in tails_by_read.items():
        row = {"sample": sample, "read_id": read_id, "tail_length": len(tail)}
        for lag in max_lags:
            window = tail[-lag:]
            row[f"end_{lag}_median"] = float(np.median(window)) if window.size else np.nan
            row[f"end_{lag}_mean"] = float(np.mean(window)) if window.size else np.nan
            row[f"end_{lag}_std"] = float(np.std(window, ddof=1)) if window.size > 1 else np.nan

        window_500 = tail[-500:]
        if window_500.size:
            row["fraction_end_500_above_1"] = float(np.mean(window_500 > 1.0))
            row["fraction_end_500_below_minus1"] = float(np.mean(window_500 < -1.0))
        else:
            row["fraction_end_500_above_1"] = np.nan
            row["fraction_end_500_below_minus1"] = np.nan
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_end_anchored_bands(gem_tails, control_tails, max_lag, min_n, outpath):
    """Gem vs control median +/- IQR, end-anchored, for one max_lag value.
    One of these is produced per requested max_lag (100/500/1000/3000) so the
    near-the-end region isn't visually compressed by also trying to show the
    full multi-thousand-sample range in the same panel."""
    fig, ax = plt.subplots(figsize=(8, 5))
    any_drawn = False

    for tails, label, color in ((gem_tails, "gem", "tab:red"), (control_tails, "control", "tab:gray")):
        if not tails:
            continue
        x, med, q1, q3, n = compute_end_median_iqr(tails, max_lag, min_n)
        if x.size == 0:
            continue
        any_drawn = True
        ax.plot(x, med, color=color, linewidth=1.4, label=f"{label} (median, n={int(n[0])} at x=0)")
        ax.fill_between(x, q1, q3, color=color, alpha=0.25, label=f"{label} IQR")

    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Distance from physical 3′ end (raw samples; 0 = final sample)")
    ax.set_ylabel("Scaled signal")
    ax.set_title("Gem vs control tail signal near physical 3′ end")
    if any_drawn:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[end_anchored max_lag={max_lag}] wrote {outpath}")


def plot_individual_end_examples(read_ids, tails, outpath, sample_label, n_examples=12, window=500, seed=0):
    """Small multiples of individual raw tail traces, end-anchored: the LAST
    `window` samples before the physical 3' end (tail[::-1][:window]), for a
    representative, fixed-seed random sample of reads. Every trace here is
    real recovered signal, not a summary statistic -- this is the "what does
    the very end of one read actually look like" view."""
    if not tails:
        print(f"[individual_end_examples:{sample_label}] no reads to plot, skipping {outpath}")
        return

    rng = np.random.RandomState(seed)
    idx = np.arange(len(tails))
    if len(idx) > n_examples:
        idx = rng.choice(idx, size=n_examples, replace=False)
        idx.sort()

    ncols = 3
    nrows = int(np.ceil(len(idx) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 2.2 * nrows),
                              sharex=True, sharey=True, squeeze=False)

    for ax, i in zip(axes.flat, idx):
        t_end = end_align_tail(tails[i])[:window]
        ax.plot(np.arange(len(t_end)), t_end, linewidth=0.6, color="tab:blue")
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.set_title(f"{read_ids[i]}\n(tail length={len(tails[i])})", fontsize=7)

    for ax in axes.flat[len(idx):]:
        ax.axis("off")

    fig.suptitle(f"Individual recovered tail signal, last {window} samples from physical 3′ end — {sample_label}")
    # fig.supxlabel/supylabel need matplotlib>=3.4; this env has 3.3.4, so use
    # shared fig.text() labels instead for compatibility (same workaround as
    # plot_tail_signal.py).
    fig.text(0.5, 0.02, "Distance from physical 3′ end (raw samples; 0 = final sample)", ha="center", fontsize=9)
    fig.text(0.02, 0.5, "Scaled signal", va="center", rotation="vertical", fontsize=9)
    fig.tight_layout(rect=[0.03, 0.04, 1, 0.96])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[individual_end_examples:{sample_label}] wrote {outpath} ({len(idx)} reads shown)")


def _boxplot_panel(ax, gem_vals, control_vals, title):
    data, labels, colors = [], [], []
    if len(gem_vals):
        data.append(gem_vals); labels.append("gem"); colors.append("tab:red")
    if len(control_vals):
        data.append(control_vals); labels.append("control"); colors.append("tab:gray")
    if not data:
        ax.set_title(f"{title}\n(no data)", fontsize=8)
        return
    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=True)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)
    ax.set_title(title, fontsize=9)


def plot_read_level_boxplots(summary_df, outpath):
    """Read-level (not sample-pooled) gem-vs-control comparison for the
    metrics computed by compute_read_level_summary. Each read contributes
    exactly one point to each panel, so a handful of reads with unusually long
    or extreme tails can't silently dominate the comparison the way a
    per-sample-pooled band plot could."""
    metrics = [
        "end_100_median", "end_500_median", "end_1000_median",
        "end_500_std", "fraction_end_500_above_1", "fraction_end_500_below_minus1",
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    gem_df = summary_df[summary_df["sample"] == "gem"]
    control_df = summary_df[summary_df["sample"] == "control"]

    for ax, metric in zip(axes.flat, metrics):
        gem_vals = gem_df[metric].dropna().to_numpy() if metric in gem_df else np.array([])
        control_vals = control_df[metric].dropna().to_numpy() if metric in control_df else np.array([])
        _boxplot_panel(ax, gem_vals, control_vals, metric)

    fig.suptitle("Read-level tail-signal summary near the physical 3′ end: gem vs control")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[read_level_boxplots] wrote {outpath}")


def plot_regular_vs_tail_boxplot(compare_df, outpath):
    """Read-level comparison of the last regular, reference-anchored signal
    (coordinates just upstream of the tail cutoff) against the recovered
    tail's own near-end signal, and the delta between them -- gem vs control."""
    metrics = [
        "regular_end_median_scaled_signal",
        "tail_end_500_median_scaled_signal",
        "delta_tail_end_vs_regular_end",
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    gem_df = compare_df[compare_df["sample"] == "gem"]
    control_df = compare_df[compare_df["sample"] == "control"]

    for ax, metric in zip(axes, metrics):
        gem_vals = gem_df[metric].dropna().to_numpy() if metric in gem_df else np.array([])
        control_vals = control_df[metric].dropna().to_numpy() if metric in control_df else np.array([])
        _boxplot_panel(ax, gem_vals, control_vals, metric)

    fig.suptitle("Regular (reference-anchored) vs. recovered tail-end signal: gem vs control")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[regular_vs_tail_boxplot] wrote {outpath}")


# ---------------------------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------------------------

def print_summary(label, tails_by_read, max_lags):
    tails = list(tails_by_read.values())
    n_reads = len(tails)
    total_samples = sum(len(t) for t in tails)
    print(f"\n{label}: {n_reads} qualifying tail reads, {total_samples} total TAIL samples")
    if not tails:
        return
    lengths = [len(t) for t in tails]
    print(f"  tail length: min={min(lengths)} median={int(np.median(lengths))} "
          f"mean={np.mean(lengths):.1f} max={max(lengths)}")
    for lag in max_lags:
        print(f"  reads reaching >= {lag} samples from the end: {reads_reaching(tails, lag)}/{n_reads}")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gem", required=True, help="path to gem .align file (from captureTailSignal run)")
    ap.add_argument("--control", default=None, help="path to control .align file (optional)")
    ap.add_argument("--outdir", required=True, help="directory to write figures/TSVs into (created if missing)")
    ap.add_argument("--max-lags", default="100,500,1000,3000",
                    help="comma-separated end-anchored window sizes to analyze/plot (default: 100,500,1000,3000)")
    ap.add_argument("--min-n", type=int, default=5,
                    help="stop drawing the end-anchored median/IQR band once fewer than this many reads remain at a position (default 5)")
    ap.add_argument("--n-examples", type=int, default=12, help="number of individual reads to show per individual-examples figure")
    ap.add_argument("--example-window", type=int, default=500, help="samples from the physical end to show per example read")
    ap.add_argument("--seed", type=int, default=0, help="random seed for selecting example reads (reproducibility)")
    ap.add_argument("--regular-coord-lo", type=int, default=600, help="lower bound (inclusive) of .align reference coordinates treated as the 'regular end' region for part 5 (default 600)")
    ap.add_argument("--regular-coord-hi", type=int, default=611, help="upper bound (inclusive) of .align reference coordinates treated as the 'regular end' region for part 5 (default 611)")
    ap.add_argument("--skip-regular-comparison", action="store_true",
                    help="skip part 5 (regular-vs-tail comparison), which requires an extra pass over each .align file")
    args = ap.parse_args()

    max_lags = [int(x) for x in args.max_lags.split(",") if x.strip()]
    os.makedirs(args.outdir, exist_ok=True)

    print("Loading gem tail reads...")
    gem_by_read = load_tail_dict(args.gem)
    print("Loading control tail reads...")
    control_by_read = load_tail_dict(args.control)

    print_summary("gem", gem_by_read, max_lags)
    print_summary("control", control_by_read, max_lags)

    if not gem_by_read:
        print("ERROR: no qualifying gem reads found in", args.gem, file=sys.stderr)
        sys.exit(1)

    gem_ids, gem_tails = list(gem_by_read.keys()), list(gem_by_read.values())
    control_ids, control_tails = list(control_by_read.keys()), list(control_by_read.values())

    # --- 1. end-anchored median/IQR bands, one figure per max_lag -----------
    for lag in max_lags:
        plot_end_anchored_bands(
            gem_tails, control_tails, max_lag=lag, min_n=args.min_n,
            outpath=os.path.join(args.outdir, f"tail_end_anchored_{lag}.png"),
        )

    # --- 2. individual end-anchored example traces, gem and control separately
    plot_individual_end_examples(
        gem_ids, gem_tails, os.path.join(args.outdir, "gem_individual_end_examples.png"),
        sample_label="gem", n_examples=args.n_examples, window=args.example_window, seed=args.seed,
    )
    plot_individual_end_examples(
        control_ids, control_tails, os.path.join(args.outdir, "control_individual_end_examples.png"),
        sample_label="control", n_examples=args.n_examples, window=args.example_window, seed=args.seed,
    )

    # --- 3. read-level summary TSV ------------------------------------------
    rows = compute_read_level_summary(gem_by_read, "gem", max_lags=max_lags) + \
           compute_read_level_summary(control_by_read, "control", max_lags=max_lags)
    summary_df = pd.DataFrame(rows)
    summary_path = os.path.join(args.outdir, "tail_end_read_level_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"[read_level_summary] wrote {summary_path} ({len(summary_df)} reads)")

    # --- 4. read-level boxplots ----------------------------------------------
    plot_read_level_boxplots(summary_df, os.path.join(args.outdir, "tail_end_read_level_boxplots.png"))

    # --- 5. regular (reference-anchored) vs tail-end comparison, optional ---
    if not args.skip_regular_comparison:
        print("Extracting regular-region signal for gem reads (extra pass over the .align file)...")
        gem_regular = extract_regular_end_signal(args.gem, gem_ids, args.regular_coord_lo, args.regular_coord_hi)
        print("Extracting regular-region signal for control reads (extra pass over the .align file)...")
        control_regular = extract_regular_end_signal(args.control, control_ids, args.regular_coord_lo, args.regular_coord_hi)

        compare_rows = []
        for sample, ids_regular, tails_by_read in (
            ("gem", gem_regular, gem_by_read),
            ("control", control_regular, control_by_read),
        ):
            for read_id, regular_vals in ids_regular.items():
                tail = tails_by_read[read_id]
                regular_median = float(np.median(regular_vals)) if regular_vals else np.nan
                tail_end_500 = tail[-500:]
                tail_median = float(np.median(tail_end_500)) if tail_end_500.size else np.nan
                delta = tail_median - regular_median if (regular_vals and tail_end_500.size) else np.nan
                compare_rows.append({
                    "sample": sample,
                    "read_id": read_id,
                    "n_regular_samples": len(regular_vals),
                    "regular_end_median_scaled_signal": regular_median,
                    "tail_end_500_median_scaled_signal": tail_median,
                    "delta_tail_end_vs_regular_end": delta,
                })

        compare_df = pd.DataFrame(compare_rows)
        compare_path = os.path.join(args.outdir, "tail_end_vs_regular_end_summary.tsv")
        compare_df.to_csv(compare_path, sep="\t", index=False)
        print(f"[regular_vs_tail_summary] wrote {compare_path} ({len(compare_df)} reads)")

        plot_regular_vs_tail_boxplot(compare_df, os.path.join(args.outdir, "tail_end_vs_regular_end_boxplot.png"))
    else:
        print("Skipping part 5 (regular-vs-tail comparison) due to --skip-regular-comparison.")

    print(f"\nDone. Reminder: these are exploratory, read-level-summarized views of a modest sample "
          f"({len(gem_tails)} gem / {len(control_tails)} control qualifying reads, deduplicated to "
          f"unique physical reads -- see any duplicate-readID messages above). See the module "
          f"docstring for scientific caveats on interpreting raw-sample positions and on not "
          f"over-reading band overlap/non-overlap.")


if __name__ == "__main__":
    main()
