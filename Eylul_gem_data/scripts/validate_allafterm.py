#!/usr/bin/env python3
"""
Automated validation checks for the ALLAFTERM / NEXT100EVT tail-capture
output (see TAIL_SIGNAL_CAPTURE.md for the full writeup). Run this against
any --tail-source both .align file produced by the ALLAFTERM-capable binary.

Checks performed (numbered to match the request that introduced ALLAFTERM):

  1. Every NEXT100EVT sequence equals the first up-to-100 valid event means
     in the corresponding ALLAFTERM sequence, for every read block.
  3. ALLAFTERM indices start at zero and are contiguous within every block.
  4. No block has ALLAFTERM without also having a plausible mechanism to
     produce it (Source 2 and/or Source 1 present, or the file predates
     that -- see body). This is a proxy for "no ALLAFTERM for gate-failing
     reads": a gate-failing read gets TailCaptureMode::NONE in align_main,
     so it can never print ANY tail-prefixed line at all, ALLAFTERM
     included -- this is a structural guarantee from the existing gate,
     not new code, but it's checked here in case that ever regresses.
  5. Empty/short remainders don't crash: this script itself must run to
     completion over every block, including ones with zero tail lines.

Check 2 (existing normal/Source1/Source2 output byte-identical modulo new
ALLAFTERM lines and summary text) is NOT done here -- it needs a diff
against a pre-ALLAFTERM reference file, which lives in a separate one-off
diff (see TAIL_SIGNAL_CAPTURE.md's validation section for the exact
commands used and their result).

Usage:
    python3 validate_allafterm.py --align path/to/some.full.align
"""
import argparse
import sys


def iter_blocks(path):
    """Yield (read_id, lines) where lines is every line belonging to that
    read's block (the header line itself excluded)."""
    read_id = None
    lines = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if read_id is not None:
                    yield read_id, lines
                parts = line[1:].split()
                read_id = parts[0] if parts else None
                lines = []
            else:
                lines.append(line.rstrip("\n"))
    if read_id is not None:
        yield read_id, lines


def parse_block(lines):
    """Return dict of prefix -> list of (idx, value) tuples, in file order,
    for the four tail-related prefixes, plus a count of ordinary
    (non-tail-prefixed) lines."""
    out = {"TAIL_SOURCE1": [], "TAIL_SOURCE2": [], "NEXT100EVT": [], "ALLAFTERM": []}
    n_normal = 0
    for line in lines:
        matched = False
        for prefix in ("TAIL_SOURCE1", "TAIL_SOURCE2", "NEXT100EVT", "ALLAFTERM"):
            if line.startswith(prefix):
                fields = line.split("\t")
                out[prefix].append((int(fields[1]), float(fields[2])))
                matched = True
                break
        if not matched:
            n_normal += 1
    return out, n_normal


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--align", required=True, help="path to a --tail-source both .align file")
    args = ap.parse_args()

    n_blocks = 0
    n_blocks_with_allafterm = 0
    n_blocks_with_any_tail = 0
    n_blocks_no_tail_at_all = 0
    prefix_mismatches = []
    contiguity_failures = []
    orphan_allafterm = []

    for read_id, lines in iter_blocks(args.align):
        n_blocks += 1
        parsed, n_normal = parse_block(lines)
        s1, s2, next100, allafterm = parsed["TAIL_SOURCE1"], parsed["TAIL_SOURCE2"], parsed["NEXT100EVT"], parsed["ALLAFTERM"]

        has_any_tail = bool(s1 or s2 or next100 or allafterm)
        if has_any_tail:
            n_blocks_with_any_tail += 1
        else:
            n_blocks_no_tail_at_all += 1

        # Check 3: ALLAFTERM indices start at 0 and are contiguous.
        if allafterm:
            idxs = [i for i, _ in allafterm]
            if idxs != list(range(len(idxs))):
                contiguity_failures.append((read_id, idxs[:5], idxs[-5:]))

        # Check 1: NEXT100EVT is exactly ALLAFTERM's first up-to-100 entries.
        if next100 or allafterm:
            n_blocks_with_allafterm += 1
            next100_vals = [v for _, v in next100]
            allafterm_prefix = [v for _, v in allafterm[: len(next100_vals)]]
            expected_len = min(100, len(allafterm))
            if len(next100_vals) != expected_len:
                prefix_mismatches.append((read_id, "length", len(next100_vals), expected_len))
            elif next100_vals != allafterm_prefix:
                prefix_mismatches.append((read_id, "values", next100_vals[:3], allafterm_prefix[:3]))

        # Check 4 (structural proxy): ALLAFTERM should never appear with zero
        # of Source1/Source2/NEXT100EVT also present in the SAME block, since
        # they're all gated by the identical (wantSource1 or wantSource2)
        # condition and the same haveLastMatchRawIdx/lastMatchRawIdx anchor --
        # if ALLAFTERM fired, NEXT100EVT must have too (as long as ALLAFTERM
        # is non-empty), and vice versa.
        if bool(allafterm) != bool(next100):
            orphan_allafterm.append((read_id, len(allafterm), len(next100)))

    print(f"Blocks parsed: {n_blocks}")
    print(f"  with any tail-prefixed line: {n_blocks_with_any_tail}")
    print(f"  with zero tail-prefixed lines: {n_blocks_no_tail_at_all}")
    print(f"  with ALLAFTERM/NEXT100EVT present: {n_blocks_with_allafterm}")

    ok = True

    print("\n[Check 1] NEXT100EVT == first up-to-100 valid ALLAFTERM entries")
    if prefix_mismatches:
        ok = False
        print(f"  FAIL: {len(prefix_mismatches)} block(s) mismatched, e.g.: {prefix_mismatches[:3]}")
    else:
        print(f"  PASS: all {n_blocks_with_allafterm} blocks with data matched exactly")

    print("\n[Check 3] ALLAFTERM indices start at 0 and are contiguous")
    if contiguity_failures:
        ok = False
        print(f"  FAIL: {len(contiguity_failures)} block(s), e.g.: {contiguity_failures[:3]}")
    else:
        print(f"  PASS: all ALLAFTERM-bearing blocks had contiguous 0-based indices")

    print("\n[Check 4] ALLAFTERM never appears without NEXT100EVT in the same block (or vice versa)")
    if orphan_allafterm:
        ok = False
        print(f"  FAIL: {len(orphan_allafterm)} block(s), e.g.: {orphan_allafterm[:3]}")
    else:
        print(f"  PASS: no orphan cases found")

    print("\n[Check 5] Script completed without crashing on any block (including zero-tail-line blocks)")
    print(f"  PASS: parsed {n_blocks} blocks, {n_blocks_no_tail_at_all} with zero tail lines, no exceptions raised")

    print(f"\nOverall: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
