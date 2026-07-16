#!/usr/bin/env python3
"""
Plot recovered 3' tail signal from DNAscent's captureTailSignal feature,
split by the mechanism that produced it (see TAIL_SIGNAL_CAPTURE.md for the
full writeup of the --tail-source option this script depends on).

Source 1 = trailing insertion-labelled events from the terminal, reference-
anchored Viterbi HMM window (printed as TAIL_SOURCE1 lines). Directly tied to
the final reference-anchored window, so may be the cleaner signal to inspect.

Source 2 = events trimmed by the rough, basecall-level aligner before
eventalign ever saw them (printed as TAIL_SOURCE2 lines). Real signal, but
less specifically reference-anchored, and typically much larger/more
variable in length than Source 1.

This script does NOT itself run DNAscent -- it reads .align output already
generated with `DNAscent align --tail-source {source1,source2,both}`. Old
.align files that still use the combined, unlabelled `TAIL\\t<idx>\\t<val>`
format (from before the --tail-source option existed) are also readable here
and treated as a single "legacy" pseudo-source, since Source 1 and Source 2
can no longer be told apart in that format -- that ambiguity is exactly the
reason --tail-source was added.

IMPORTANT SCIENTIFIC CAUTIONS (see also the printed terminal summary):
  - x-axis positions are raw signal sample indices, not base coordinates.
  - Do NOT treat TAIL_SOURCE1/TAIL_SOURCE2 lines as independent observations
    -- the biological/statistical unit is the READ. Every plot and summary
    here is built per-read first (see load_read_blocks's deduplication) and
    only pooled across reads for the band plots.
  - Source 1 and Source 2 may have very different length distributions (in
    the current gem/control data, Source 1 is bounded by one small HMM
    window while Source 2 can run to tens of thousands of raw samples), so
    comparing their absolute tail lengths directly is not meaningful; compare
    within a source, across gem vs control, instead.
  - A "combined" view (Source 1 + Source 2 concatenated, chronological) is
    provided for continuity with the older, unsplit TAIL output, but mixes
    two different mechanisms -- prefer the per-source views for anything
    beyond a sanity check that recoverable signal exists at all.

Run with the gmconda_py368 conda env (matplotlib/numpy/pandas):

    /opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python \\
        plot_tail_signal_by_source.py \\
        --gem     ../smoketest_by_source/gem.full.both.align \\
        --control ../smoketest_by_source/control.full.both.align \\
        --outdir  ../smoketest_by_source/figures_by_source
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

def iter_read_blocks(path):
    """Yield (read_id, source1_values, source2_values, legacy_values) for
    every '>'-delimited read block in a .align file that has at least one
    tail line of any kind.

    Source-specific prefixes (TAIL_SOURCE1/TAIL_SOURCE2) are checked before
    the bare "TAIL" prefix, since both source lines also start with "TAIL" --
    only a line that is neither of the specific sources but still starts with
    "TAIL" is legacy combined output from before --tail-source existed.
    """
    read_id = None
    s1, s2, legacy = [], [], []

    def _has_data():
        return bool(s1) or bool(s2) or bool(legacy)

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if read_id is not None and _has_data():
                    yield read_id, s1, s2, legacy
                parts = line[1:].split()
                read_id = parts[0] if parts else None
                s1, s2, legacy = [], [], []
            elif line.startswith("TAIL_SOURCE1"):
                s1.append(float(line.rstrip("\n").split("\t")[2]))
            elif line.startswith("TAIL_SOURCE2"):
                s2.append(float(line.rstrip("\n").split("\t")[2]))
            elif line.startswith("TAIL"):
                legacy.append(float(line.rstrip("\n").split("\t")[2]))
    if read_id is not None and _has_data():
        yield read_id, s1, s2, legacy


def load_read_blocks(path):
    """Return {read_id: {'source1': arr, 'source2': arr, 'legacy': arr}} for
    every qualifying read in path, deduplicated by read_id.

    Same duplicate-readID situation documented for plot_tail_signal_end_focused.py
    applies here: a physical read can produce more than one qualifying
    alignment block against this short reference (secondary/supplementary
    alignments). On collision, each channel (source1/source2/legacy) is kept
    independently from whichever colliding block has the longer array for
    THAT channel -- picking one whole block by total signal is not safe here,
    because colliding blocks can be complementary rather than redundant: in
    control's dataset, one read had two blocks, {source1=990, source2=405202}
    and {source1=0, source2=409401} -- the second block's total is only 0.8%
    larger, but picking it whole (as an earlier version of this function did)
    silently discarded the single largest Source 1 signal in the entire
    dataset. Keeping the longer array per channel avoids that failure mode.
    Prints how many collisions it resolved.
    """
    blocks = {}
    n_blocks = 0
    n_collisions = 0
    if path is None or not os.path.exists(path):
        return blocks

    for read_id, s1, s2, legacy in iter_read_blocks(path):
        n_blocks += 1
        s1a = np.asarray(s1, dtype=float)
        s2a = np.asarray(s2, dtype=float)
        lega = np.asarray(legacy, dtype=float)
        if read_id in blocks:
            n_collisions += 1
            prev = blocks[read_id]
            blocks[read_id] = {
                "source1": s1a if s1a.size > prev["source1"].size else prev["source1"],
                "source2": s2a if s2a.size > prev["source2"].size else prev["source2"],
                "legacy": lega if lega.size > prev["legacy"].size else prev["legacy"],
            }
        else:
            blocks[read_id] = {"source1": s1a, "source2": s2a, "legacy": lega}

    if n_collisions:
        print(f"  [{os.path.basename(path)}] {n_blocks} qualifying alignment blocks -> "
              f"{len(blocks)} unique reads ({n_collisions} duplicate-readID block(s) "
              f"collapsed, longest array per channel kept independently across "
              f"colliding blocks -- see this function's docstring)")
    return blocks


def extract_source(blocks, key):
    """{read_id: array} for reads that have non-empty signal under the given
    channel key ('source1', 'source2', or 'legacy')."""
    return {rid: b[key] for rid, b in blocks.items() if b[key].size > 0}


def extract_combined(blocks):
    """{read_id: array} combining source1+source2 (chronological, matching
    the pre-split combined TAIL ordering) where either is present, falling
    back to legacy TAIL data for reads/files that only have that. Also
    returns whether any legacy data was used, so callers can label the plot
    accordingly."""
    combined = {}
    used_legacy = False
    for rid, b in blocks.items():
        if b["source1"].size or b["source2"].size:
            combined[rid] = np.concatenate([b["source1"], b["source2"]])
        elif b["legacy"].size:
            combined[rid] = b["legacy"]
            used_legacy = True
    return combined, used_legacy


# ---------------------------------------------------------------------------
# Building comparable matrices across reads of very different tail length
# ---------------------------------------------------------------------------

def position_matrix(tails, max_len, anchor):
    """Stack reads into an (n_reads x max_len) matrix, NaN-padded, so
    per-position stats can be computed across reads of very different tail
    length. anchor="cutoff": column 0 = this source's first recovered
    sample. anchor="end": column 0 = this source's last recovered sample
    (tail[::-1])."""
    n = len(tails)
    mat = np.full((n, max_len), np.nan)
    for i, t in enumerate(tails):
        L = min(len(t), max_len)
        if L == 0:
            continue
        if anchor == "cutoff":
            mat[i, :L] = t[:L]
        elif anchor == "end":
            mat[i, :L] = t[::-1][:L]
        else:
            raise ValueError("anchor must be 'cutoff' or 'end'")
    return mat


def median_iqr(mat):
    """Per-column median/Q1/Q3/n, ignoring NaN."""
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
    n = np.sum(~np.isnan(mat), axis=0)
    return med, q1, q3, n


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_band(gem_tails, control_tails, anchor, outpath, title, max_lag, min_n):
    """Gem vs control median +/- IQR band for one source and one anchor
    (cutoff or end). Truncated once fewer than min_n reads remain at a
    position, rather than drawing a misleading band from a couple of
    surviving reads at an extreme lag."""
    fig, ax = plt.subplots(figsize=(8, 5))
    any_drawn = False

    for tails, label, color in ((gem_tails, "gem", "tab:red"), (control_tails, "control", "tab:gray")):
        if not tails:
            continue
        mat = position_matrix(tails, max_len=max_lag, anchor=anchor)
        med, q1, q3, n = median_iqr(mat)
        valid = n >= min_n
        if not np.any(valid):
            continue
        last = int(np.max(np.where(valid)))
        x = np.arange(last + 1)
        any_drawn = True
        ax.plot(x, med[: last + 1], color=color, linewidth=1.4, label=f"{label} (median, n={int(n[0])} at x=0)")
        ax.fill_between(x, q1[: last + 1], q3[: last + 1], color=color, alpha=0.25, label=f"{label} IQR")

    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
    xlabel = ("Position from cutoff (raw samples; 0 = first recovered sample)" if anchor == "cutoff"
              else "Distance from physical 3′ end (raw samples; 0 = final sample)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Scaled signal")
    ax.set_title(title, fontsize=10)
    if any_drawn:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[band:{anchor}] wrote {outpath}")


def plot_examples(read_ids, tails, outpath, title, n_examples=12, window=500, seed=0):
    """Small multiples of individual raw traces, first `window` samples of
    this source's own signal (index 0 = this source's own cutoff -- for
    Source 1 that's essentially the reference-anchored cutoff itself; for
    Source 2 it's wherever the rough aligner gave up, which may be well
    downstream of Source 1's samples). Fixed-seed random sample of reads."""
    if not tails:
        print(f"[examples] no reads to plot, skipping {outpath}")
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
        t = tails[i][:window]
        ax.plot(np.arange(len(t)), t, linewidth=0.6, color="tab:blue")
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.set_title(f"{read_ids[i]}\n(length={len(tails[i])})", fontsize=7)

    for ax in axes.flat[len(idx):]:
        ax.axis("off")

    fig.suptitle(title, fontsize=10)
    # fig.supxlabel/supylabel need matplotlib>=3.4; this env has 3.3.4, so use
    # shared fig.text() labels instead for compatibility (same workaround as
    # the other plotting scripts in this directory).
    fig.text(0.5, 0.02, "Position from this source's own cutoff (raw samples)", ha="center", fontsize=9)
    fig.text(0.02, 0.5, "Scaled signal", va="center", rotation="vertical", fontsize=9)
    fig.tight_layout(rect=[0.03, 0.04, 1, 0.95])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[examples] wrote {outpath} ({len(idx)} reads shown)")


# ---------------------------------------------------------------------------
# Read-level summary
# ---------------------------------------------------------------------------

def read_level_rows(tails_by_read, sample, source):
    rows = []
    for read_id, tail in tails_by_read.items():
        rows.append({
            "sample": sample,
            "source": source,
            "read_id": read_id,
            "tail_length": len(tail),
            "median_scaled_signal": float(np.median(tail)) if tail.size else np.nan,
            "mean_scaled_signal": float(np.mean(tail)) if tail.size else np.nan,
            "std_scaled_signal": float(np.std(tail, ddof=1)) if tail.size > 1 else np.nan,
        })
    return rows


def print_summary(label, tails_by_read):
    tails = list(tails_by_read.values())
    n_reads = len(tails)
    total_samples = sum(len(t) for t in tails)
    print(f"  {label}: {n_reads} reads, {total_samples} total samples", end="")
    if tails:
        lengths = [len(t) for t in tails]
        print(f" (length min={min(lengths)} median={int(np.median(lengths))} "
              f"mean={np.mean(lengths):.1f} max={max(lengths)})")
    else:
        print()


# ---------------------------------------------------------------------------

def process_source(source_key, source_title, gem_tails, control_tails, args, summary_rows):
    print(f"\n=== {source_title} ===")
    print_summary("gem", gem_tails)
    print_summary("control", control_tails)

    if not gem_tails and not control_tails:
        print(f"  no data for this source -- skipping its plots")
        return

    plot_band(
        list(gem_tails.values()), list(control_tails.values()), anchor="cutoff",
        outpath=os.path.join(args.outdir, f"tail_by_source_{source_key}_cutoff_band.png"),
        title=f"{source_title}: gem vs control, anchored at cutoff",
        max_lag=args.max_lag, min_n=args.min_n,
    )
    plot_band(
        list(gem_tails.values()), list(control_tails.values()), anchor="end",
        outpath=os.path.join(args.outdir, f"tail_by_source_{source_key}_end_band.png"),
        title=f"{source_title}: gem vs control, anchored at physical 3′ end",
        max_lag=args.max_lag, min_n=args.min_n,
    )
    # Individual traces for BOTH gem and control (not gem-only): a source can be
    # sparse or entirely absent in one arm and only present in the other (e.g.
    # Source 1 turned out to be gem-empty/control-only in this project's actual
    # data), so gem-only examples would silently show nothing for exactly the
    # reads worth looking at.
    for tails, label in ((gem_tails, "gem"), (control_tails, "control")):
        plot_examples(
            list(tails.keys()), list(tails.values()),
            outpath=os.path.join(args.outdir, f"tail_by_source_{source_key}_{label}_examples.png"),
            title=f"{source_title}: individual {label} read traces",
            n_examples=args.n_examples, window=args.example_window, seed=args.seed,
        )

    summary_rows.extend(read_level_rows(gem_tails, "gem", source_key))
    summary_rows.extend(read_level_rows(control_tails, "control", source_key))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gem", required=True, help="path to gem .align file (from a --tail-source run)")
    ap.add_argument("--control", default=None, help="path to control .align file (optional)")
    ap.add_argument("--outdir", required=True, help="directory to write figures/TSV into (created if missing)")
    ap.add_argument("--max-lag", type=int, default=3000, help="max position (either anchor) for the band plots (default 3000)")
    ap.add_argument("--min-n", type=int, default=5, help="stop drawing a band once fewer than this many reads remain at a position (default 5)")
    ap.add_argument("--n-examples", type=int, default=12, help="number of individual reads to show per examples figure")
    ap.add_argument("--example-window", type=int, default=500, help="samples from this source's own cutoff to show per example read")
    ap.add_argument("--seed", type=int, default=0, help="random seed for selecting example reads (reproducibility)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading gem read blocks...")
    gem_blocks = load_read_blocks(args.gem)
    print("Loading control read blocks...")
    control_blocks = load_read_blocks(args.control)

    if not gem_blocks:
        print("ERROR: no qualifying gem reads (no TAIL_SOURCE1/TAIL_SOURCE2/legacy TAIL lines) found in", args.gem, file=sys.stderr)
        sys.exit(1)

    summary_rows = []

    process_source(
        "source1", "Source 1 (terminal-window trailing insertions)",
        extract_source(gem_blocks, "source1"), extract_source(control_blocks, "source1"),
        args, summary_rows,
    )
    process_source(
        "source2", "Source 2 (rough-aligner-trimmed events)",
        extract_source(gem_blocks, "source2"), extract_source(control_blocks, "source2"),
        args, summary_rows,
    )

    gem_combined, gem_used_legacy = extract_combined(gem_blocks)
    control_combined, control_used_legacy = extract_combined(control_blocks)
    used_legacy = gem_used_legacy or control_used_legacy
    combined_title = ("Combined (legacy, unsplit TAIL lines)" if used_legacy
                       else "Combined (Source 1 + Source 2, chronological)")
    process_source("combined", combined_title, gem_combined, control_combined, args, summary_rows)

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(args.outdir, "tail_by_source_read_level_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\n[read_level_summary] wrote {summary_path} ({len(summary_df)} rows)")

    print(f"\nDone. Reminder: Source 1 and Source 2 are different recovery mechanisms with "
          f"different length distributions (see module docstring) -- don't compare their raw "
          f"tail lengths directly, and treat the 'combined' view as a continuity check against "
          f"older unsplit output, not as the primary analysis.")


if __name__ == "__main__":
    main()
