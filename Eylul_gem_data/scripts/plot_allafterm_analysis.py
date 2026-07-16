#!/usr/bin/env python3
"""
Analysis of DNAscent's NEXT100EVT and ALLAFTERM tail-capture pools (see
TAIL_SIGNAL_CAPTURE.md sec 9-10 and PROJECT_OVERVIEW.md sec 4.6-4.7 for the
full writeup of what these are and how they're produced). This script does
NOT overwrite plot_tail_signal_by_source.py or plot_tail_signal_end_focused.py
-- it is a new, independent script for the new output.

WHAT NEXT100EVT AND ALLAFTERM ARE, AND HOW THEY RELATE TO EVERYTHING ELSE
---------------------------------------------------------------------------
Both are anchored at `lastMatchRawIdx` -- the raw r.events index of the last
event the reference-anchored HMM actually labelled Match, in the read's final
processed window (src/alignment.cpp; see PROJECT_OVERVIEW.md sec 4.6 for the
full derivation, including why this differs from Source 2's own boundary
whenever a read has real Source 1 signal). NEXT100EVT is exactly the first
up-to-100 valid entries of ALLAFTERM (same anchor, same 0<mean<250 guard, same
shift/scale normalization, same iteration order) -- ALLAFTERM simply removes
the +100 cap and continues to the end of r.events. This equality is asserted
by validate_allafterm.py and is NOT re-derived here; this script parses the
two prefixes independently from the .align file regardless, and uses whichever
one is the natural fit for a given figure (NEXT100EVT for the first-100 view,
ALLAFTERM for anything longer).

Four important distinctions, restated explicitly because they are easy to
blur when just looking at plots:
  1. Source 1 / Source 2 (TAIL_SOURCE1 / TAIL_SOURCE2 lines) are printed at
     EXPANDED RAW-ADC-SAMPLE granularity -- one line per raw sample, several
     samples per detected event. NEXT100EVT and ALLAFTERM are printed at
     EVENT-MEAN granularity -- one line per detected event (its `.mean`),
     never expanded. A "length" in Source 1/2 space (raw samples) and a
     "length" in NEXT100EVT/ALLAFTERM space (events) are NOT the same unit
     and must never be compared directly.
  2. ALLAFTERM can and typically does include BOTH Source 1 territory (if the
     read's terminal window assigned any trailing events as Insertion) AND
     Source 2 territory (events the rough aligner trimmed before Viterbi ever
     ran) -- it does not respect that mechanism-level split at all; it is
     defined purely by position in r.events relative to lastMatchRawIdx.
  3. Neither NEXT100EVT nor ALLAFTERM proves that a recovered event
     corresponds to one of the reference's final bases (e.g. the gemcitabine
     position). This is recovered, past-the-modelled-region raw signal --
     it may include pore adapter signal, motor protein signal, ordinary dwell
     variability, or genuinely non-reference sequence (e.g. extra ligated
     material for reads that would have failed the 3' soft-clip gate --
     though such reads never reach this pool at all, see the gate). Treat
     every plot here as "what is the recovered signal like", not "where is
     the modified base".
  4. Events already labelled Match by the HMM are NOT part of either pool --
     they remain in the ordinary, reference-anchored `.align` output exactly
     as before this feature existed. NEXT100EVT/ALLAFTERM start strictly
     AFTER the last such Match event.

IMPORTANT SCIENTIFIC CAUTIONS
---------------------------------------------------------------------------
  - "Detected event number after final HMM Match" (Part A) and "distance from
    the final detected event" (Part B) are event-index positions, not base
    positions and not raw-sample positions -- see distinction 1 above.
  - Do not read a gem-vs-control difference off band overlap/non-overlap
    alone -- see the closing read-level discussion this script's caller
    should provide; every figure here also reports the number of
    *unique reads* contributing at each position, specifically so a viewer
    can tell when a band's tail is being carried by a handful of reads.
  - Sample sizes are the same modest ones as the rest of this project (order
    of tens to a few hundred qualifying reads) -- treat everything here as
    exploratory.

DEDUPLICATION RULE (read carefully -- differs from plot_tail_signal_by_source.py)
---------------------------------------------------------------------------
A physical read can produce more than one qualifying alignment block in a
.align file (documented throughout this project). For Parts A-D, this script
deduplicates by read_id, and on a collision **keeps the block with the longer
ALLAFTERM array** (not the longer Source 1, not a per-channel merge) -- this
is deliberately simpler than plot_tail_signal_by_source.py's per-channel rule,
because NEXT100EVT is *always* a prefix of its OWN block's ALLAFTERM (they are
computed together, from the same lastMatchRawIdx, in the same eventalign()
call), so there is no complementary-channel scenario to reconcile here the way
there was for Source1 vs Source2. One consequence: a read's *winning* block
for ALLAFTERM-length purposes is not guaranteed to be the block that happens
to carry Source 1 signal (Source 1 can be small even when its block's overall
ALLAFTERM/Source 2 pool is smaller than a sibling block's). Part E's example
selection specifically works around this -- see its own docstring.

Run with the gmconda_py368 conda env (matplotlib/numpy/pandas/scipy):

    /opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python \\
        plot_allafterm_analysis.py \\
        --gem     ../slurm_full_run_allafterm/gem/gem.full.align \\
        --control ../slurm_full_run_allafterm/control/control.full.align \\
        --outdir  ../slurm_full_run_allafterm/figures_allafterm
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib
matplotlib.use("Agg")  # no display on a compute/login node
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def iter_all_blocks(path):
    """Yield (read_id, next100_arr, allafterm_arr, has_source1) for EVERY
    '>'-delimited alignment block in path that has at least one NEXT100EVT or
    ALLAFTERM line. next100/allafterm are parsed independently of one another
    (separate accumulators below) -- neither is derived from the other here,
    even though validate_allafterm.py has already confirmed the expected
    prefix relationship holds. has_source1 is a lightweight boolean (TAIL_SOURCE1
    line present in this block at all, values not retained) used only by Part
    E's example selection.
    """
    read_id = None
    next100, allafterm = [], []
    has_source1 = False

    def _has_data():
        return bool(next100) or bool(allafterm)

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if read_id is not None and _has_data():
                    yield read_id, np.asarray(next100, dtype=float), np.asarray(allafterm, dtype=float), has_source1
                parts = line[1:].split()
                read_id = parts[0] if parts else None
                next100, allafterm = [], []
                has_source1 = False
            elif line.startswith("NEXT100EVT"):
                next100.append(float(line.rstrip("\n").split("\t")[2]))
            elif line.startswith("ALLAFTERM"):
                allafterm.append(float(line.rstrip("\n").split("\t")[2]))
            elif line.startswith("TAIL_SOURCE1"):
                has_source1 = True
    if read_id is not None and _has_data():
        yield read_id, np.asarray(next100, dtype=float), np.asarray(allafterm, dtype=float), has_source1


def load_blocks(path):
    """Return (main_dict, all_blocks_dict, n_blocks, n_collisions).

    main_dict: {read_id: (next100_arr, allafterm_arr, has_source1)} -- one
    entry per unique read, the block with the LONGER ALLAFTERM array kept on
    collision (see module docstring for why this rule, and why it differs
    from plot_tail_signal_by_source.py's per-channel rule).

    all_blocks_dict: {read_id: [(next100_arr, allafterm_arr, has_source1), ...]}
    -- every block, undeduped, kept only so Part E can find a Source-1-bearing
    block for a read even when it isn't the ALLAFTERM-length winner.

    Prints block/unique-read/collision counts (never silent about dedup, per
    this project's established convention).
    """
    main = {}
    all_blocks = {}
    n_blocks = 0
    n_collisions = 0
    if path is None or not os.path.exists(path):
        return main, all_blocks, n_blocks, n_collisions

    for read_id, next100, allafterm, has_s1 in iter_all_blocks(path):
        n_blocks += 1
        all_blocks.setdefault(read_id, []).append((next100, allafterm, has_s1))
        if read_id in main:
            n_collisions += 1
            prev_next100, prev_allafterm, prev_s1 = main[read_id]
            if allafterm.size > prev_allafterm.size:
                main[read_id] = (next100, allafterm, has_s1)
            # else keep the existing (longer-or-equal) entry
        else:
            main[read_id] = (next100, allafterm, has_s1)

    n_unique = len(main)
    print(f"  [{os.path.basename(path)}] {n_blocks} qualifying alignment blocks -> "
          f"{n_unique} unique reads" + (f" ({n_collisions} duplicate-readID block(s) "
          f"collapsed, longer ALLAFTERM kept)" if n_collisions else ""))
    return main, all_blocks, n_blocks, n_collisions


def source1_aware_allafterm(read_id, main, all_blocks):
    """For Part E only: return the ALLAFTERM array to use for this read's
    example trace. If any of this read's blocks has Source 1 signal, use the
    longest-ALLAFTERM block AMONG THOSE (so the trace shown actually reflects
    the Source-1-influenced block), rather than blindly using `main`'s
    overall winner, which the module docstring explains is not guaranteed to
    be the same block. Falls back to `main`'s entry otherwise."""
    candidates = all_blocks.get(read_id, [])
    s1_candidates = [(next100, allafterm) for next100, allafterm, has_s1 in candidates if has_s1]
    if s1_candidates:
        best = max(s1_candidates, key=lambda pair: pair[1].size)
        return best[1]
    return main[read_id][1]


# ---------------------------------------------------------------------------
# Matrix / statistics helpers
# ---------------------------------------------------------------------------

def stack_matrix(arrays, max_len):
    """(n_reads x max_len) matrix, NaN-padded past each array's own length."""
    n = len(arrays)
    mat = np.full((n, max_len), np.nan)
    for i, a in enumerate(arrays):
        L = min(len(a), max_len)
        if L:
            mat[i, :L] = a[:L]
    return mat


def median_iqr(mat):
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
    n = np.sum(~np.isnan(mat), axis=0)
    return med, q1, q3, n


# ---------------------------------------------------------------------------
# Part A / Part B: two-panel (median+IQR, then n-contributing-reads) band plot
# ---------------------------------------------------------------------------

def plot_band_two_panel(gem_arrays, control_arrays, max_len, min_n, x_offset,
                         outpath_png, outpath_tsv, title, xlabel):
    """x_offset=1 for Part A (HMM-Match-anchored, 1-based: x=1 is the first
    event after Match). x_offset=0 for Part B (physical-end-anchored,
    0-based: x=0 is the final detected event). Truncates each sample's curve
    once fewer than min_n unique reads remain at that x, rather than drawing
    a band from a couple of surviving reads at an extreme position. Lower
    panel plots n-contributing-unique-reads vs x directly, so a viewer never
    has to infer coverage from the band's disappearance alone."""
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    tsv_rows = []
    any_drawn = False

    for arrays, label, color in ((gem_arrays, "gem", "tab:red"), (control_arrays, "control", "tab:gray")):
        if not arrays:
            continue
        mat = stack_matrix(arrays, max_len)
        med, q1, q3, n = median_iqr(mat)
        valid = n >= min_n
        if not np.any(valid):
            continue
        last = int(np.max(np.where(valid)))
        x = np.arange(last + 1) + x_offset
        any_drawn = True
        n_at_start = int(n[0])
        ax_top.plot(x, med[: last + 1], color=color, linewidth=1.4,
                    label=f"{label} (median, n={n_at_start} unique reads at x={x_offset})")
        ax_top.fill_between(x, q1[: last + 1], q3[: last + 1], color=color, alpha=0.25, label=f"{label} IQR")
        ax_bot.plot(x, n[: last + 1], color=color, linewidth=1.2)
        for xi, m, qq1, qq3, ni in zip(x, med[: last + 1], q1[: last + 1], q3[: last + 1], n[: last + 1]):
            tsv_rows.append({"x": int(xi), "sample": label, "n_contributing_unique_reads": int(ni),
                              "median": m, "q1": qq1, "q3": qq3})

    ax_top.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax_top.set_ylabel("Scaled event mean")
    ax_top.set_title(title, fontsize=10)
    if any_drawn:
        ax_top.legend(fontsize=8)
    else:
        ax_top.set_title(title + "\n(no data)", fontsize=10)
    ax_bot.set_xlabel(xlabel)
    ax_bot.set_ylabel("N contributing\nunique reads")
    fig.tight_layout()
    fig.savefig(outpath_png, dpi=150)
    plt.close(fig)

    pd.DataFrame(tsv_rows).to_csv(outpath_tsv, sep="\t", index=False)
    print(f"  wrote {outpath_png}")
    print(f"  wrote {outpath_tsv} ({len(tsv_rows)} rows)")


# ---------------------------------------------------------------------------
# Part F: FULL (untruncated) ALLAFTERM -- every event after the final HMM
# Match, for every deduplicated read, no 100/500/3000-event cutoff anywhere.
# All calculations below use ONLY the `allafterm` arrays already loaded from
# ALLAFTERM lines (via load_blocks -> main dict) -- not NEXT100EVT, not
# Source 1/2 raw samples, and not any physical-end-reversed subset.
# ---------------------------------------------------------------------------

def plot_full_read_boxplots(gem_allafterm, control_allafterm, outpath_png, outpath_tsv):
    """Exactly one (mean, median) pair per unique read, computed over that
    read's ENTIRE ALLAFTERM array (no truncation at all -- a read with
    725,085 events contributes one mean and one median, computed from all
    725,085, same as a read with 1,042)."""
    rows = []
    for sample, d in (("gem", gem_allafterm), ("control", control_allafterm)):
        for read_id, arr in d.items():
            rows.append({
                "read_id": read_id, "sample": sample, "n_allafterm_events": int(arr.size),
                "mean": float(np.mean(arr)) if arr.size else np.nan,
                "median": float(np.median(arr)) if arr.size else np.nan,
            })
    df = pd.DataFrame(rows)
    df.to_csv(outpath_tsv, sep="\t", index=False)
    print(f"  wrote {outpath_tsv} ({len(df)} rows, one per unique read)")

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    for ax, metric in zip(axes, ("mean", "median")):
        gem_vals = df.loc[df["sample"] == "gem", metric].dropna().to_numpy()
        control_vals = df.loc[df["sample"] == "control", metric].dropna().to_numpy()
        _boxplot_panel(ax, gem_vals, control_vals, f"Full-ALLAFTERM per-read {metric}")
        ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_ylabel("Scaled event mean")
    fig.suptitle("Per-read mean/median over the COMPLETE post-Match ALLAFTERM sequence\n"
                 "(no 100/500-event cutoff -- every event after the final HMM Match, per read)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(outpath_png, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath_png}")
    print(f"  gem:     n={df[df['sample']=='gem'].shape[0]} unique reads")
    print(f"  control: n={df[df['sample']=='control'].shape[0]} unique reads")


def median_mean_iqr(mat):
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
        mean = np.nanmean(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
    n = np.sum(~np.isnan(mat), axis=0)
    return med, mean, q1, q3, n


def threshold_crossing_x(n_at_x1, n, fractions=(0.5, 0.25, 0.10, 0.05)):
    """First 0-based column index i where n[i] < fraction * n_at_x1, for each
    fraction. n_at_x1 is this sample's total unique-read count (== n[0],
    since every read has at least 1 ALLAFTERM event by construction of the
    gate). Returns {fraction: index_or_None}."""
    out = {}
    for f in fractions:
        below = np.where(n < f * n_at_x1)[0]
        out[f] = int(below[0]) if below.size else None
    return out


def _draw_threshold_lines(ax, crossings, x_offset, color, sample_label, y_text):
    """Vertical lines (with compact text labels) at each coverage-threshold
    crossing, in the given sample's color. y_text controls vertical text
    placement so gem/control labels don't overlap when both are drawn on the
    same axis."""
    linestyles = {0.5: ":", 0.25: "--", 0.10: "-.", 0.05: (0, (1, 1))}
    for frac, idx in crossings.items():
        if idx is None:
            continue
        x = idx + x_offset
        ax.axvline(x, color=color, linestyle=linestyles[frac], linewidth=0.9, alpha=0.6)
        ax.text(x, y_text, f"{sample_label} {int(frac*100)}%", color=color, fontsize=6,
                rotation=90, va="top", ha="right", alpha=0.8)


def decimate_for_plot(x, *arrays, max_points=20000):
    """Stride-based decimation for RENDERING ONLY -- never applied to the
    TSV, which always retains full, one-row-per-event resolution.
    matplotlib's Agg backend cannot reliably draw a line/filled polygon with
    hundreds of thousands of vertices on this environment (raises
    `OverflowError: In draw_markers: Exceeded cell block limit` past roughly
    a few hundred thousand points, hit directly when first plotting control's
    full ~725k-event trace). Keeps the exact first and last point so the
    complete x-range's endpoints are always represented, even though not
    every single intermediate point gets its own vertex -- this is a display
    resolution limit, not a truncation of the plotted range."""
    n = len(x)
    if n <= max_points:
        return (x,) + arrays
    stride = int(np.ceil(n / max_points))
    idx = np.arange(0, n, stride)
    if idx[-1] != n - 1:
        idx = np.append(idx, n - 1)
    return (x[idx],) + tuple(a[idx] for a in arrays)


def plot_full_hmm_anchored_trace(gem_allafterm, control_allafterm, outdir):
    """Gem vs control median+IQR from x=1 (first event after the final HMM
    Match) through the FULL available ALLAFTERM range -- i.e. through
    whichever read (in either sample) has the most events. No 100/500-event
    cutoff. Lower panel: n contributing unique reads vs x, with vertical
    markers at each sample's own 50/25/10/5%-of-starting-count crossing.
    Produces both a linear-x and a log-x version, both showing the complete
    range (log-x is NOT a truncation -- every position is still plotted,
    just with a different x-axis scale for readability given a ~700x span
    between the shortest and longest reads).
    """
    gem_arrays = list(gem_allafterm.values())
    control_arrays = list(control_allafterm.values())
    max_len = max((a.size for a in gem_arrays + control_arrays), default=0)
    if max_len == 0:
        print("  [full_hmm_anchored_trace] no data, skipping")
        return

    print(f"  building full-length matrices (max_len={max_len} events, "
          f"gem n={len(gem_arrays)}, control n={len(control_arrays)}) ...")

    per_sample = {}
    tsv_rows = []
    for arrays, label, color in ((gem_arrays, "gem", "tab:red"), (control_arrays, "control", "tab:gray")):
        if not arrays:
            continue
        mat = stack_matrix(arrays, max_len)
        med, mean, q1, q3, n = median_mean_iqr(mat)
        last_with_data = int(np.max(np.where(n > 0))) if np.any(n > 0) else -1
        per_sample[label] = (med, mean, q1, q3, n, color, last_with_data)
        crossings = threshold_crossing_x(n[0], n)
        for i in range(last_with_data + 1):
            tsv_rows.append({
                "sample": label, "event_offset_after_last_M": i + 1,
                "n_unique_reads": int(n[i]), "median": med[i], "mean": mean[i],
                "q1": q1[i], "q3": q3[i],
            })
        crossing_report = {f"{int(f*100)}%": (idx + 1 if idx is not None else "never")
                            for f, idx in crossings.items()}
        print(f"    {label}: starting n={int(n[0])}, longest read={last_with_data + 1} events, "
              f"coverage drops below {crossing_report}")

    tsv_path = os.path.join(outdir, "allafterm_hmm_anchored_full_trace.tsv")
    pd.DataFrame(tsv_rows).to_csv(tsv_path, sep="\t", index=False)
    print(f"  wrote {tsv_path} ({len(tsv_rows)} rows)")

    for logx, suffix in ((False, ""), (True, "_logx")):
        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, figsize=(9, 7.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
        )
        for label, (med, mean, q1, q3, n, color, last_with_data) in per_sample.items():
            x_full = np.arange(last_with_data + 1) + 1  # x=1 = first event after Match
            x, med_p, q1_p, q3_p, n_p = decimate_for_plot(
                x_full, med[: last_with_data + 1], q1[: last_with_data + 1],
                q3[: last_with_data + 1], n[: last_with_data + 1],
            )
            ax_top.plot(x, med_p, color=color, linewidth=1.0,
                        label=f"{label} (median, n={int(n[0])} unique reads at x=1)")
            ax_top.fill_between(x, q1_p, q3_p, color=color, alpha=0.2, label=f"{label} IQR")
            ax_bot.plot(x, n_p, color=color, linewidth=1.0)
            crossings = threshold_crossing_x(n[0], n)
            _draw_threshold_lines(ax_top, crossings, 1, color, label, ax_top.get_ylim()[1])
            _draw_threshold_lines(ax_bot, crossings, 1, color, label, n[0])

        ax_top.axhline(0, color="black", linewidth=0.5, linestyle="--")
        ax_top.set_ylabel("Scaled event mean")
        title_suffix = " (log-x)" if logx else ""
        ax_top.set_title(
            f"Gem vs control: COMPLETE post-Match ALLAFTERM trace{title_suffix}\n"
            f"(x=1 = first event after Match, full range, no cutoff; vertical\n"
            f"lines = 50/25/10/5% of that sample's starting unique-read count)",
            fontsize=9,
        )
        ax_top.legend(fontsize=7, loc="upper right")
        ax_bot.set_xlabel("Detected event number after final HMM Match (1 = first event after Match)")
        ax_bot.set_ylabel("N contributing\nunique reads")
        if logx:
            ax_top.set_xscale("log")
            ax_bot.set_xscale("log")
        fig.tight_layout()
        outpath = os.path.join(outdir, f"allafterm_hmm_anchored_full_trace{suffix}.png")
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
        print(f"  wrote {outpath}")


# ---------------------------------------------------------------------------
# Part C: read-level windowed metrics + boxplots
# ---------------------------------------------------------------------------

def read_level_metrics(allafterm_by_read, sample, window):
    """One row per unique read, for the first `window` events after the
    final HMM Match (allafterm[:window] -- if the read has fewer events than
    `window`, this uses however many it has, same as numpy slice semantics
    elsewhere in this project; never padded or extrapolated)."""
    rows = []
    for read_id, arr in allafterm_by_read.items():
        w = arr[:window]
        if w.size == 0:
            continue
        row = {"read_id": read_id, "sample": sample, "window": window, "n_events_in_window": int(w.size)}
        row["median"] = float(np.median(w))
        row["mean"] = float(np.mean(w))
        row["std"] = float(np.std(w, ddof=1)) if w.size > 1 else np.nan
        if w.size >= 2:
            x = np.arange(1, w.size + 1, dtype=float)
            slope, _intercept, _r, _p, _se = linregress(x, w)
            row["slope_vs_event_number"] = float(slope)
        else:
            row["slope_vs_event_number"] = np.nan
        row["fraction_above_plus1"] = float(np.mean(w > 1.0))
        row["fraction_below_minus1"] = float(np.mean(w < -1.0))
        rows.append(row)
    return rows


def _boxplot_panel(ax, gem_vals, control_vals, title):
    data, labels, colors = [], [], []
    if len(gem_vals):
        data.append(gem_vals); labels.append(f"gem\n(n={len(gem_vals)})"); colors.append("tab:red")
    if len(control_vals):
        data.append(control_vals); labels.append(f"control\n(n={len(control_vals)})"); colors.append("tab:gray")
    if not data:
        ax.set_title(f"{title}\n(no data)", fontsize=8)
        return
    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=True)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)
    ax.set_title(title, fontsize=9)


def plot_read_level_boxplots(summary_df, window, outpath):
    metrics = ["median", "mean", "std", "slope_vs_event_number", "fraction_above_plus1", "fraction_below_minus1"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    df_w = summary_df[summary_df["window"] == window]
    gem_df = df_w[df_w["sample"] == "gem"]
    control_df = df_w[df_w["sample"] == "control"]
    for ax, metric in zip(axes.flat, metrics):
        gem_vals = gem_df[metric].dropna().to_numpy()
        control_vals = control_df[metric].dropna().to_numpy()
        _boxplot_panel(ax, gem_vals, control_vals, metric)
    fig.suptitle(f"Read-level post-Match signal summary, first {window} events: gem vs control")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")


def plot_focused_median_boxplot(summary_df, window, outpath):
    """Single-panel boxplot of just the per-read median over the first
    `window` post-Match events -- the ALLAFTERM-analogue of
    tail_end_3000_median_boxplot.png."""
    df_w = summary_df[summary_df["window"] == window]
    gem_vals = df_w.loc[df_w["sample"] == "gem", "median"].dropna().to_numpy()
    control_vals = df_w.loc[df_w["sample"] == "control", "median"].dropna().to_numpy()

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
    ax.set_ylabel(f"Per-read median scaled event mean,\nfirst {window} events after final HMM Match")
    ax.set_title(f"Read-level median post-Match signal (first {window} events): gem vs control")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")
    if len(gem_vals):
        print(f"    gem:     n={len(gem_vals)}, median of medians={np.median(gem_vals):.4f}")
    if len(control_vals):
        print(f"    control: n={len(control_vals)}, median of medians={np.median(control_vals):.4f}")


# ---------------------------------------------------------------------------
# Part D: remaining-tail length distribution
# ---------------------------------------------------------------------------

def plot_remaining_length_distributions(gem_lengths, control_lengths, outdir):
    # Histogram
    fig, ax = plt.subplots(figsize=(7, 5))
    bins = 40
    if len(gem_lengths):
        ax.hist(gem_lengths, bins=bins, color="tab:red", alpha=0.5, label=f"gem (n={len(gem_lengths)})")
    if len(control_lengths):
        ax.hist(control_lengths, bins=bins, color="tab:gray", alpha=0.5, label=f"control (n={len(control_lengths)})")
    ax.set_xlabel("ALLAFTERM events per unique read")
    ax.set_ylabel("Count")
    ax.set_title("Remaining post-Match event count per read: gem vs control")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(outdir, "allafterm_remaining_length_hist.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print(f"  wrote {p}")

    # Log-scaled histogram (x-axis log, since tails can be very long)
    fig, ax = plt.subplots(figsize=(7, 5))
    all_pos = np.concatenate([a[a > 0] for a in (gem_lengths, control_lengths) if len(a)]) if (len(gem_lengths) or len(control_lengths)) else np.array([1])
    lo = max(1, all_pos.min()) if all_pos.size else 1
    hi = all_pos.max() if all_pos.size else 10
    logbins = np.logspace(np.log10(lo), np.log10(hi), 40)
    if len(gem_lengths):
        ax.hist(gem_lengths, bins=logbins, color="tab:red", alpha=0.5, label=f"gem (n={len(gem_lengths)})")
    if len(control_lengths):
        ax.hist(control_lengths, bins=logbins, color="tab:gray", alpha=0.5, label=f"control (n={len(control_lengths)})")
    ax.set_xscale("log")
    ax.set_xlabel("ALLAFTERM events per unique read (log scale)")
    ax.set_ylabel("Count")
    ax.set_title("Remaining post-Match event count per read (log-scaled x): gem vs control")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(outdir, "allafterm_remaining_length_hist_log.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print(f"  wrote {p}")

    # ECDF
    fig, ax = plt.subplots(figsize=(7, 5))
    for lengths, label, color in ((gem_lengths, "gem", "tab:red"), (control_lengths, "control", "tab:gray")):
        if not len(lengths):
            continue
        xs = np.sort(lengths)
        ys = np.arange(1, len(xs) + 1) / len(xs)
        ax.step(xs, ys, where="post", color=color, label=f"{label} (n={len(lengths)})")
    ax.set_xscale("log")
    ax.set_xlabel("ALLAFTERM events per unique read (log scale)")
    ax.set_ylabel("Empirical CDF")
    ax.set_title("Remaining post-Match event count per read (ECDF): gem vs control")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(outdir, "allafterm_remaining_length_ecdf.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print(f"  wrote {p}")

    # Box/violin
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    data, labels, colors = [], [], []
    if len(gem_lengths):
        data.append(gem_lengths); labels.append(f"gem\n(n={len(gem_lengths)})"); colors.append("tab:red")
    if len(control_lengths):
        data.append(control_lengths); labels.append(f"control\n(n={len(control_lengths)})"); colors.append("tab:gray")
    if data:
        bp = axes[0].boxplot(data, labels=labels, patch_artist=True, showfliers=True)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color); patch.set_alpha(0.4)
        axes[0].set_yscale("log")
        axes[0].set_ylabel("ALLAFTERM events per read (log scale)")
        axes[0].set_title("Boxplot")
        vp = axes[1].violinplot(data, showmedians=True)
        for body, color in zip(vp["bodies"], colors):
            body.set_facecolor(color); body.set_alpha(0.4)
        axes[1].set_xticks(range(1, len(labels) + 1))
        axes[1].set_xticklabels(labels)
        axes[1].set_yscale("log")
        axes[1].set_title("Violin")
    fig.suptitle("Remaining post-Match event count per read: gem vs control")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(outdir, "allafterm_remaining_length_boxplot.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print(f"  wrote {p}")


# ---------------------------------------------------------------------------
# Part E: individual traces
# ---------------------------------------------------------------------------

def select_examples(read_ids, all_blocks, n_examples, seed):
    """Fixed-seed selection: every read with a Source-1-bearing block is
    force-included (documented as "if practical" -- practical here, since
    both datasets have very few such reads), topped up to n_examples with a
    reproducible random sample of the rest."""
    rng = np.random.RandomState(seed)
    forced = [rid for rid in read_ids if any(has_s1 for _, _, has_s1 in all_blocks.get(rid, []))]
    remaining = [rid for rid in read_ids if rid not in forced]
    n_more = max(0, n_examples - len(forced))
    if len(remaining) > n_more:
        chosen_more = list(rng.choice(remaining, size=n_more, replace=False))
    else:
        chosen_more = remaining
    selected = forced + chosen_more
    return selected


def plot_examples(read_ids, main, all_blocks, outpath, sample_label, window=500):
    """First `window` ALLAFTERM events for each selected example read, one
    subplot per read. x=0 is marked as the final HMM Match (true by
    construction -- ALLAFTERM's own index 0 is defined as the first event
    after that Match). Source 1 presence is annotated in the subplot title
    when known, but the Source1-to-event-count boundary itself is NOT marked
    on the trace: Source 1/Source 2 are printed at expanded raw-ADC-sample
    granularity while ALLAFTERM is event-mean granularity, and the per-event
    raw-sample count needed to convert between them is not recoverable from
    this output -- see module docstring distinction 1. Stating this rather
    than estimating a boundary that can't be reliably reconstructed."""
    if not read_ids:
        print(f"  [examples:{sample_label}] no reads to plot, skipping {outpath}")
        return

    ncols = 3
    nrows = int(np.ceil(len(read_ids) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 2.4 * nrows), squeeze=False)

    for ax, rid in zip(axes.flat, read_ids):
        arr = source1_aware_allafterm(rid, main, all_blocks) if rid in main else np.array([])
        has_s1 = any(has_s1 for _, _, has_s1 in all_blocks.get(rid, []))
        t = arr[:window]
        ax.plot(np.arange(len(t)), t, linewidth=0.6, color="tab:blue")
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.axvline(0, color="black", linewidth=0.8, linestyle=":")
        s1_tag = "Source1: yes" if has_s1 else "Source1: no"
        ax.set_title(f"{rid}\n(n={len(arr)} total, {s1_tag})", fontsize=6.5)

    for ax in axes.flat[len(read_ids):]:
        ax.axis("off")

    fig.suptitle(
        f"{sample_label}: individual post-Match traces (first {window} ALLAFTERM events)\n"
        f"x=0 = final HMM Match (dotted line). Source1/Source2 event-boundary NOT marked "
        f"-- not reliably reconstructable at event granularity from raw-sample-level Source "
        f"output (see module docstring).", fontsize=8
    )
    fig.text(0.5, 0.01, "Detected event number after final HMM Match", ha="center", fontsize=9)
    fig.text(0.01, 0.5, "Scaled event mean", va="center", rotation="vertical", fontsize=9)
    fig.tight_layout(rect=[0.03, 0.05, 1, 0.88])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath} ({len(read_ids)} reads shown)")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gem", required=True, help="path to gem .align file (from an ALLAFTERM-capable --tail-source both run)")
    ap.add_argument("--control", default=None, help="path to control .align file (optional)")
    ap.add_argument("--outdir", required=True, help="directory to write figures/TSVs into (created if missing)")
    ap.add_argument("--min-n", type=int, default=5, help="stop drawing a band once fewer than this many unique reads remain at a position (default 5)")
    ap.add_argument("--n-examples", type=int, default=12, help="individual reads shown per Part E figure, per sample (default 12, >=10 required)")
    ap.add_argument("--example-window", type=int, default=500, help="ALLAFTERM events shown per Part E example trace (default 500)")
    ap.add_argument("--seed", type=int, default=0, help="random seed for Part E example selection (reproducibility)")
    args = ap.parse_args()

    if args.n_examples < 10:
        print("WARNING: --n-examples < 10; the analysis request calls for at least 10 example reads per sample.", file=sys.stderr)

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading gem read blocks...")
    gem_main, gem_all, gem_n_blocks, gem_n_coll = load_blocks(args.gem)
    print("Loading control read blocks...")
    control_main, control_all, control_n_blocks, control_n_coll = load_blocks(args.control)

    if not gem_main:
        print("ERROR: no qualifying gem reads (no NEXT100EVT/ALLAFTERM lines) found in", args.gem, file=sys.stderr)
        sys.exit(1)

    gem_next100 = {rid: v[0] for rid, v in gem_main.items()}
    gem_allafterm = {rid: v[1] for rid, v in gem_main.items()}
    control_next100 = {rid: v[0] for rid, v in control_main.items()}
    control_allafterm = {rid: v[1] for rid, v in control_main.items()}

    print(f"\ngem:     {len(gem_main)} unique reads, {gem_n_blocks} alignment blocks, {gem_n_coll} collision(s)")
    print(f"control: {len(control_main)} unique reads, {control_n_blocks} alignment blocks, {control_n_coll} collision(s)")

    # ==================== Part A: HMM-Match-anchored bands =================
    print("\n=== Part A: HMM-Match-anchored median +/- IQR ===")
    plot_band_two_panel(
        list(gem_next100.values()), list(control_next100.values()),
        max_len=100, min_n=args.min_n, x_offset=1,
        outpath_png=os.path.join(args.outdir, "allafterm_hmm_anchored_first100.png"),
        outpath_tsv=os.path.join(args.outdir, "allafterm_hmm_anchored_first100.tsv"),
        title="Gem vs control: first 100 events after final HMM Match (NEXT100EVT)",
        xlabel="Detected event number after final HMM Match (1 = first event after Match)",
    )
    plot_band_two_panel(
        list(gem_allafterm.values()), list(control_allafterm.values()),
        max_len=500, min_n=args.min_n, x_offset=1,
        outpath_png=os.path.join(args.outdir, "allafterm_hmm_anchored_first500.png"),
        outpath_tsv=os.path.join(args.outdir, "allafterm_hmm_anchored_first500.tsv"),
        title="Gem vs control: first 500 events after final HMM Match (ALLAFTERM)",
        xlabel="Detected event number after final HMM Match (1 = first event after Match)",
    )

    # ==================== Part B: physical-end-anchored bands ==============
    print("\n=== Part B: physical-end-anchored median +/- IQR ===")
    gem_reversed = [a[::-1] for a in gem_allafterm.values()]
    control_reversed = [a[::-1] for a in control_allafterm.values()]
    for lag in (100, 500, 3000):
        plot_band_two_panel(
            gem_reversed, control_reversed,
            max_len=lag, min_n=args.min_n, x_offset=0,
            outpath_png=os.path.join(args.outdir, f"allafterm_physical_end_last{lag}.png"),
            outpath_tsv=os.path.join(args.outdir, f"allafterm_physical_end_last{lag}.tsv"),
            title=f"Gem vs control: last {lag} detected events before physical 3′ end (ALLAFTERM, reversed)",
            xlabel="Detected events from physical 3′ end (0 = final detected event)",
        )

    # ==================== Part F: FULL (untruncated) ALLAFTERM ==============
    print("\n=== Part F: full-ALLAFTERM per-read boxplots and full HMM-anchored trace (no cutoff) ===")
    plot_full_read_boxplots(
        gem_allafterm, control_allafterm,
        outpath_png=os.path.join(args.outdir, "allafterm_full_read_mean_median_boxplots.png"),
        outpath_tsv=os.path.join(args.outdir, "allafterm_full_read_mean_median.tsv"),
    )
    plot_full_hmm_anchored_trace(gem_allafterm, control_allafterm, args.outdir)

    # ==================== Part C: read-level windowed metrics ==============
    print("\n=== Part C: read-level boxplots (first 25/50/100 post-Match events) ===")
    rows = []
    for window in (25, 50, 100):
        rows += read_level_metrics(gem_allafterm, "gem", window)
        rows += read_level_metrics(control_allafterm, "control", window)
    summary_df = pd.DataFrame(rows)
    summary_path = os.path.join(args.outdir, "allafterm_read_level_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"  wrote {summary_path} ({len(summary_df)} rows)")

    for window in (25, 50, 100):
        plot_read_level_boxplots(summary_df, window, os.path.join(args.outdir, f"allafterm_read_level_boxplots_{window}.png"))
    plot_focused_median_boxplot(summary_df, 100, os.path.join(args.outdir, "allafterm_first100_median_boxplot.png"))

    # ==================== Part D: remaining-tail length =====================
    print("\n=== Part D: remaining post-Match event count per read ===")
    gem_lengths = np.array([len(a) for a in gem_allafterm.values()])
    control_lengths = np.array([len(a) for a in control_allafterm.values()])
    plot_remaining_length_distributions(gem_lengths, control_lengths, args.outdir)

    length_rows = (
        [{"read_id": rid, "sample": "gem", "n_allafterm_events": len(a)} for rid, a in gem_allafterm.items()]
        + [{"read_id": rid, "sample": "control", "n_allafterm_events": len(a)} for rid, a in control_allafterm.items()]
    )
    length_path = os.path.join(args.outdir, "allafterm_remaining_length.tsv")
    pd.DataFrame(length_rows).to_csv(length_path, sep="\t", index=False)
    print(f"  wrote {length_path} ({len(length_rows)} rows)")

    # ==================== Part E: individual traces =========================
    print("\n=== Part E: individual traces (Source1/Source2 boundary NOT estimated -- see caveats) ===")
    gem_examples = select_examples(list(gem_main.keys()), gem_all, args.n_examples, args.seed)
    control_examples = select_examples(list(control_main.keys()), control_all, args.n_examples, args.seed)
    print(f"  gem: {len(gem_examples)} example reads selected "
          f"({sum(1 for r in gem_examples if any(h for _, _, h in gem_all.get(r, [])))} with Source1)")
    print(f"  control: {len(control_examples)} example reads selected "
          f"({sum(1 for r in control_examples if any(h for _, _, h in control_all.get(r, [])))} with Source1)")
    plot_examples(gem_examples, gem_main, gem_all,
                  os.path.join(args.outdir, "allafterm_examples_gem.png"), "gem", window=args.example_window)
    plot_examples(control_examples, control_main, control_all,
                  os.path.join(args.outdir, "allafterm_examples_control.png"), "control", window=args.example_window)

    print(f"\nDone. Reminder: {len(gem_main)} gem / {len(control_main)} control unique qualifying reads "
          f"(deduplicated -- see collision counts above). This is exploratory, event-count-anchored, "
          f"read-level-summarized analysis -- see module docstring for what NEXT100EVT/ALLAFTERM can "
          f"and cannot tell you, and do not read a gem-vs-control conclusion off band/box overlap alone.")


if __name__ == "__main__":
    main()
