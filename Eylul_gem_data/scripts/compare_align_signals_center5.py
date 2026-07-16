#!/usr/bin/env python3
"""Compare DNAscent align signal residuals using center-coordinate labels.

DNAscent align reports signal against a 9-mer pore model, and the expected
signal column is the expected unmodified signal for that full 9-mer context.
This script keeps the residual calculation unchanged:

    residual = observed_scaled_signal - expected_unmodified_9mer_signal

The only difference from compare_align_signals.py is labeling/grouping: each
9-mer is reported at its central base coordinate (coord + 4) with the middle
5-mer context (kmer[2:7]).

A center-coordinate signal should still be interpreted as nanopore signal from
the surrounding sequence context, not as proof that a single base is modified.
"""

import argparse
import csv
import math
from collections import defaultdict


class RunningStats:
    def __init__(self):
        self.n = 0
        self.mean = 0.0
        self.m2 = 0.0

    def add(self, x):
        self.n += 1
        delta = x - self.mean
        self.mean += delta / self.n
        self.m2 += delta * (x - self.mean)

    @property
    def var(self):
        if self.n < 2:
            return float("nan")
        return self.m2 / (self.n - 1)


def read_align(path, min_events_per_read_coord=1):
    stats = defaultdict(RunningStats)
    original_9mer_examples = {}
    current_read = None
    per_read_values = defaultdict(list)

    def flush_read():
        for (center_coord, center_5mer), vals in per_read_values.items():
            if len(vals) >= min_events_per_read_coord:
                stats[(center_coord, center_5mer)].add(sum(vals) / len(vals))

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith(">"):
                if current_read is not None:
                    flush_read()
                current_read = line
                per_read_values = defaultdict(list)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            coord = int(fields[0])
            kmer = fields[1]
            if len(kmer) < 9:
                continue

            center_coord = coord + 4
            center_5mer = kmer[2:7]
            observed = float(fields[2])
            expected = float(fields[4])
            residual = observed - expected
            key = (center_coord, center_5mer)

            per_read_values[key].append(residual)
            original_9mer_examples.setdefault(key, kmer)

    if current_read is not None:
        flush_read()

    return stats, original_9mer_examples


def normal_p_two_sided(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def benjamini_hochberg(rows):
    ordered = sorted(
        [row for row in rows if not math.isnan(row["p_approx"])],
        key=lambda row: row["p_approx"],
    )
    m = len(ordered)
    prev = 1.0
    for rank, row in reversed(list(enumerate(ordered, start=1))):
        q = min(prev, row["p_approx"] * m / rank)
        row["q_approx_bh"] = q
        prev = q
    for row in rows:
        row.setdefault("q_approx_bh", float("nan"))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare DNAscent align scaled-signal residuals between gem and "
            "control, reporting the 9-mer center coordinate and middle 5-mer."
        )
    )
    parser.add_argument("--gem", required=True)
    parser.add_argument("--control", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-n", type=int, default=20)
    args = parser.parse_args()

    gem, gem_9mers = read_align(args.gem)
    control, control_9mers = read_align(args.control)

    rows = []
    for key in sorted(set(gem) & set(control)):
        center_coord, center_5mer = key
        g = gem[key]
        c = control[key]
        if g.n < args.min_n or c.n < args.min_n:
            continue

        g_var = g.var
        c_var = c.var
        delta = g.mean - c.mean
        se = math.sqrt(g_var / g.n + c_var / c.n) if g.n > 1 and c.n > 1 else float("nan")
        z = delta / se if se and not math.isnan(se) else float("nan")
        p = normal_p_two_sided(z) if not math.isnan(z) else float("nan")
        pooled = math.sqrt(((g.n - 1) * g_var + (c.n - 1) * c_var) / (g.n + c.n - 2))
        cohen_d = delta / pooled if pooled else float("nan")

        rows.append(
            {
                "center_coord": center_coord,
                "center_5mer": center_5mer,
                "original_9mer_example": gem_9mers.get(key, control_9mers.get(key, "")),
                "gem_n_reads": g.n,
                "control_n_reads": c.n,
                "gem_mean_residual": g.mean,
                "control_mean_residual": c.mean,
                "gem_minus_control": delta,
                "cohen_d": cohen_d,
                "z_approx": z,
                "p_approx": p,
            }
        )

    benjamini_hochberg(rows)

    fields = [
        "center_coord",
        "center_5mer",
        "original_9mer_example",
        "gem_n_reads",
        "control_n_reads",
        "gem_mean_residual",
        "control_mean_residual",
        "gem_minus_control",
        "cohen_d",
        "z_approx",
        "p_approx",
        "q_approx_bh",
    ]
    with open(args.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} tested center positions to {args.out}")


if __name__ == "__main__":
    main()
