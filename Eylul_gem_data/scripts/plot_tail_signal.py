#!/usr/bin/env python3
"""
Plot the recovered 3' tail signal produced by DNAscent's captureTailSignal
feature (see ../../TAIL_SIGNAL_CAPTURE.md for the full writeup of what this
data is and how it was generated).

Run with the gmconda_py368 conda env — it has matplotlib/numpy and matches
the Python version the earlier *_median_iqr.ipynb notebooks in this repo used.
The bare system `python3` does NOT have matplotlib/numpy installed.

    /opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python \
        plot_tail_signal.py \
        --gem     ../slurm_full_run/gem/gem.full.align \
        --control ../slurm_full_run/control/control.full.align \
        --outdir  ../slurm_full_run/figures

`--control` is optional — if omitted (or the file has no qualifying reads
yet), the script just skips the gem-vs-control comparison and still produces
the gem-only example plot.

Outputs (into --outdir):
  tail_signal_examples.png        small multiples of individual raw tail
                                   traces for a handful of gem reads
  tail_signal_mean_band.png       gem (and control, if given) median +/- IQR
                                   band, one panel anchored from the cutoff
                                   point and one anchored from the true 3'
                                   end of the read (see notes below on why
                                   both views exist)
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")  # no display on a compute/login node
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def iter_read_tails(path):
    """Yield (header, tail_values) for every read in a DNAscent .align file
    that has at least one TAIL line.

    TAIL lines look like: "TAIL\t<idx>\t<scaled_value>" — idx is a 0-based
    sequential position within this read's recovered tail, NOT a reference
    coordinate (there is no valid reference coordinate for this signal; that
    is the whole point of the feature). tail_values[0] is therefore the first
    raw sample immediately after the point where the normal, reference-
    anchored aligner stopped assigning this read to the genome.

    A read with zero TAIL lines either didn't qualify for capture (see
    TAIL_SIGNAL_CAPTURE.md section 4 for the exact gate: forward-mapped,
    alignment reaches the true reference end, no 3' soft clip) or, in a rare
    edge case, qualified but had nothing left to recover. Either way it's
    uninformative here, so such reads are simply not yielded.
    """
    header, tail = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None and tail:
                    yield header, tail
                header, tail = line.strip(), []
            elif line.startswith("TAIL"):
                tail.append(float(line.rstrip("\n").split("\t")[2]))
    if header is not None and tail:
        yield header, tail


def load_tails(path):
    """Return (headers, tails) as parallel lists. tails[i] is a 1D float
    numpy array, in chronological order (position 0 = right after cutoff)."""
    headers, tails = [], []
    if path is None or not os.path.exists(path):
        return headers, tails
    for header, tail in iter_read_tails(path):
        headers.append(header)
        tails.append(np.asarray(tail, dtype=float))
    return headers, tails


# ---------------------------------------------------------------------------
# Building comparable matrices across reads of very different tail length
# ---------------------------------------------------------------------------

def position_matrix(tails, max_len, anchor):
    """Stack reads into a (n_reads x max_len) matrix so per-position stats
    can be computed across reads of very different tail length.

    anchor="cutoff": column 0 = first recovered sample (right after the
        normal aligner stopped). Reads shorter than max_len simply leave
        the remaining columns as NaN (they don't get padded/extrapolated).
    anchor="end": column 0 = the very last raw sample of the read (the true
        physical 3' terminus). This is the one point that's genuinely
        comparable across reads regardless of how much signal got trimmed
        before it — see TAIL_SIGNAL_CAPTURE.md's caveat that tail length
        itself varies hugely (5.6k-59.5k raw samples in the gem data), so
        "position from cutoff" mixes together very different physical
        distances-from-the-end depending on the read.

    Either way, NaN means "this read doesn't reach this position" — it is
    never filled in, so per-column stats below are computed only over reads
    that actually have data there (see median_iqr).
    """
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
    """Per-column median/IQR/n, ignoring NaN. n is returned so callers can
    stop drawing (or fade) the band once too few reads remain to be
    meaningful, rather than plotting a "confidence band" built from one or
    two reads out at the far tail."""
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
    n = np.sum(~np.isnan(mat), axis=0)
    return med, q1, q3, n


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_examples(headers, tails, outpath, n_examples=12, window=500, seed=0):
    """Small multiples of individual raw tail traces, first `window` samples
    from the cutoff point, for a representative (fixed-seed random) sample of
    reads. Individual traces, not summarized — this is the "what does one
    read actually look like" view."""
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
        read_id = headers[i].split()[0].lstrip(">")
        ax.set_title(f"{read_id}\n(tail length={len(tails[i])})", fontsize=7)

    for ax in axes.flat[len(idx):]:
        ax.axis("off")

    fig.suptitle(f"Individual recovered tail signal, first {window} samples from cutoff")
    # fig.supxlabel/supylabel need matplotlib>=3.4; this env has 3.3.4, so use
    # shared fig.text() labels instead for compatibility.
    fig.text(0.5, 0.02, "Position from cutoff (raw samples)", ha="center", fontsize=9)
    fig.text(0.02, 0.5, "Scaled signal", va="center", rotation="vertical", fontsize=9)
    fig.tight_layout(rect=[0.03, 0.04, 1, 0.96])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[examples] wrote {outpath} ({len(idx)} reads shown)")


def plot_mean_band(gem_tails, control_tails, outpath, max_lag=3000, min_n=5):
    """Gem (and control, if available) median +/- IQR band, two panels:
    left anchored from the cutoff, right anchored from the true 3' end.

    Bands are truncated once the number of reads contributing at that
    position drops below `min_n` (default 5) — with sample sizes this small
    (order of tens of reads), a "band" built from 1-2 reads at an extreme
    lag is not a meaningful confidence interval and would be misleading to
    show without qualification.
    """
    have_control = bool(control_tails)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    for ax, anchor, title in zip(
        axes,
        ("cutoff", "end"),
        ("Anchored at cutoff (pos 0 = first recovered sample)",
         "Anchored at true 3' end (pos 0 = last raw sample of read)"),
    ):
        for tails, label, color in (
            (gem_tails, "gem", "tab:red"),
            (control_tails if have_control else [], "control", "tab:gray"),
        ):
            if not tails:
                continue
            mat = position_matrix(tails, max_len=max_lag, anchor=anchor)
            med, q1, q3, n = median_iqr(mat)
            valid = n >= min_n
            if not np.any(valid):
                continue
            last = np.max(np.where(valid))
            x = np.arange(last + 1)
            ax.plot(x, med[: last + 1], color=color, label=f"{label} (median, n up to {int(n[0])})")
            ax.fill_between(x, q1[: last + 1], q3[: last + 1], color=color, alpha=0.25, label=f"{label} IQR")

        ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Position (raw samples)")
        ax.set_ylabel("Scaled signal")
        ax.legend(fontsize=7)

    fig.suptitle(
        f"Gem vs control recovered tail signal (median ± IQR, truncated below n={min_n} reads)"
        if have_control else
        f"Gem recovered tail signal (median ± IQR, truncated below n={min_n} reads)"
    )
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[mean_band] wrote {outpath}")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gem", required=True, help="path to gem .align file (from captureTailSignal run)")
    ap.add_argument("--control", default=None, help="path to control .align file (optional)")
    ap.add_argument("--outdir", required=True, help="directory to write figures into (created if missing)")
    ap.add_argument("--n-examples", type=int, default=12, help="number of individual gem reads to show in the small-multiples plot")
    ap.add_argument("--example-window", type=int, default=500, help="samples from cutoff to show per example read")
    ap.add_argument("--max-lag", type=int, default=3000, help="max position (either anchor) considered for the median/IQR band plot")
    ap.add_argument("--min-n", type=int, default=5, help="stop drawing the median/IQR band once fewer than this many reads remain at a position")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    gem_headers, gem_tails = load_tails(args.gem)
    print(f"gem: {len(gem_tails)} qualifying reads with recovered tail signal")
    if gem_tails:
        lengths = [len(t) for t in gem_tails]
        print(f"  tail length: min={min(lengths)} median={int(np.median(lengths))} max={max(lengths)}")

    control_headers, control_tails = load_tails(args.control)
    if args.control is not None:
        print(f"control: {len(control_tails)} qualifying reads with recovered tail signal")
        if not control_tails:
            print("  (no qualifying control reads found -- gem-vs-control comparison will be gem-only)")

    if not gem_tails:
        print("ERROR: no qualifying gem reads found in", args.gem, file=sys.stderr)
        sys.exit(1)

    plot_examples(
        gem_headers, gem_tails,
        outpath=os.path.join(args.outdir, "tail_signal_examples.png"),
        n_examples=args.n_examples, window=args.example_window,
    )

    plot_mean_band(
        gem_tails, control_tails,
        outpath=os.path.join(args.outdir, "tail_signal_mean_band.png"),
        max_lag=args.max_lag, min_n=args.min_n,
    )


if __name__ == "__main__":
    main()
