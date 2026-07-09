# Gemcitabine tail-signal capture: what changed, why, and how to use it

## 1. Background and motivation

DNAscent's `align` executable normally stops assigning signal to a read as soon as
the reference runs out of a full k-mer window to call against. Concretely, this
codebase uses `k = 9` (`src/config.h:45`) and an HMM chunking window of
`windowLength_align = 50` bp (`src/config.h:46`). For BrdU/EdU analysis this is a
non-issue — those analogues are never expected right at a read's 3' terminus. For
the gemcitabine (Gem) experiment, the opposite is true: Gem is expected to be the
**second-to-last base** of the ligation-product reference
(`Eylul_gem_data/ligation_product.fasta`, 615 bp, ends `...C-A-C-T-C-G`, i.e. the
`C` immediately before the final `G` is the Gem position). That position sits
inside the last few bases of the reference, where no full centred 9-mer exists —
exactly the region the original code silently discards.

The physical DNA molecule still passes through the pore past that point, so the
raw signal exists; DNAscent just never looked at it. This effort adds the ability
to **recover and print that discarded signal** so it can be inspected for anything
that looks like Gem incorporation.

Everything below documents the code as it stands in this repo
(`/scratch/esevinc22/Summer_Project_Boemo_Group`), diffed against the untouched
upstream copy at `/scratch/esevinc22/DNAscent/DNAscent`.

## 2. File-by-file changes

`eventalign()` is a function inside
`src/alignment.cpp`, and that's where nearly all of the new logic lives.

### `src/alignment.cpp`

**New function `hasSoftClipAtReferenceThreePrime` (line 519)**
Original behavior: didn't exist. BAM CIGAR strings are always stored in
reference/genome coordinate order, regardless of the read's mapped strand. So the
reference's 3' end corresponds to the *last* CIGAR operation for a forward-mapped
read, and the *first* CIGAR operation for a reverse-mapped read. This function
checks that operation and returns `true` if it's a soft clip (`BAM_CSOFT_CLIP`).

Why: the synthetic construct here is a ligation product (500 bp + 100 bp), and
some molecules may have extra material ligated on past the designed 3' end. If a
read has a 3' soft clip, whatever signal follows the aligned region belongs to
that unknown extra sequence, not to "a few bases past position 615 of our
construct." This function is how the gate (below) tells the two cases apart.

**`eventalign()` signature change (line 563)**
Original: `eventalign(DNAscent::read &r, unsigned int totalWindowLength)`.
Now: takes a third parameter, `bool captureTailSignal`, controlling whether the
new recovery logic (below) runs at all. `alignment.h` was updated to match.

**New local buffer `tailSignal` (line 577)**
A `std::vector<double>` that collects scaled raw signal samples that have no
reference coordinate. Two independent things feed into it:

**Source 1 — trailing insertions in the terminal HMM window (lines 783–792).**
In the existing per-window Viterbi decoding, insertion-labelled (`I`) events past
the window's last match (`evIdx >= lastM_ev`) were previously always dropped,
because normally the *next* window picks them up. For the terminal window there
is no next window. `reachedFinalKmer` (line 712) detects that this window is the
last one the outer loop will ever run for this read. When it's both the final
window and `captureTailSignal` is true, those trailing-insertion events are
pushed into `tailSignal` instead of being discarded.

**Source 2 — events past the end of the rough event alignment (lines 802–825).**
This is the bigger source in practice (see §4's caveats). The adaptive-banded aligner
(`event_handling.cpp`) picks a best-scoring final event index and silently trims
anything after it — it does *not* force-match trailing events if that scores
worse than just stopping. Those trimmed events never enter `r.eventAlignment` at
all, so nothing above this point in the function would ever see them. This block
walks `r.events` from `r.eventAlignment.back().first + 1` to the end of the read
and appends every raw sample (after the same `0 < mean < 250` sanity guard used
elsewhere) to `tailSignal`.

**Output (lines 827–833)**
`tailSignal` is emitted in chronological order (source 1 first, then source 2) as
lines `TAIL\t<idx>\t<scaled_value>` appended to `r.humanReadable_eventalignOut` —
the same string that becomes the `.align` output for this read. `idx` is just a
0-based sequential position within the tail, **not** a reference coordinate.

### `src/alignment.h`
Signature updated to match: `eventalign(DNAscent::read &, unsigned int, bool);`
plus the declaration for `hasSoftClipAtReferenceThreePrime(bam1_t *)`.

### `src/alignment.cpp`, `align_main()` (lines 967–977)
This is the gate. Before calling `eventalign()`:
```cpp
bool captureTailSignal = (r.strand == "fwd")
                          and (r.refEnd == (int) reference.at(r.referenceMappedTo).size())
                          and not hasSoftClipAtReferenceThreePrime(r.record);
eventalign(r, Pore_Substrate_Config.windowLength_align, captureTailSignal);
```
All three conditions must hold:
1. **Forward-mapped.** The project only cares about the 3' end; for a
   reverse-mapped read the physical 3' end is at the reference's 5' side, which
   this code doesn't handle (explicitly out of scope — see the original prompt:
   "on 5' end is fine because our goal is to look at Gem at the 3' end").
2. **Alignment reaches the true end of the reference contig** (`refEnd ==
   reference length`). If the read's alignment stops short of position 615,
   there's no "past the end" signal to recover for this read.
3. **No 3' soft clip.** This is the ligation-artifact guard described above.

### `src/detect.cpp` (line ~893) and `src/trainCNN.cpp` (two call sites, ~line 328
and ~line 334)
No logic changes — both files' `eventalign(...)` calls simply gained a literal
`false` as the third argument, since `detect` and `trainCNN` have no use for tail
signal and must keep their existing behavior identical.

## 3. Build/environment issue and how to reproduce a working build

Neither this repo nor the untouched upstream copy had ever linked successfully
here — `bin/` only contained an empty `.keep` file. The compile step works fine;
the final link fails:
```
/usr/bin/ld: cannot find -llzma
collect2: error: ld returned 1 exit status
```
This machine has the *runtime* library (`/usr/lib64/liblzma.so.5`, found via
`ldconfig -p`) but no `-dev` package, so there's no unversioned `liblzma.so`
symlink for the linker to find. This is an environment gap, not a code bug — the
`DNAscent.def` Singularity recipe installs `liblzma-dev` via `apt`, which isn't
available on this host outside a container.

**Fix used during this session** (fragile — depends on a specific conda package
still being present in the shared pkg cache):
```bash
make LDFLAGS="-L/opt/ohpc/pub/groups/kuiscid/conda3/latest/pkgs/xz-5.2.10-h5eee18b_1/lib -ldl -llzma -lbz2 -lm -lz"
```

**More durable fix (recommended going forward):** a personal symlink was created
at `$HOME/lib/liblzma.so -> /usr/lib64/liblzma.so.5` (this doesn't touch any
shared/system path, just `$HOME`). Rebuild with:
```bash
make LDFLAGS="-L$HOME/lib -ldl -llzma -lbz2 -lm -lz"
```
This has been verified to resolve correctly (`gcc ... -L$HOME/lib -llzma` links
clean) and doesn't depend on any particular conda package staying installed.
If the underlying HPC gap is ever fixed properly (i.e. `liblzma-dev` installed
system-wide), plain `make` with no `LDFLAGS` override should also start working.

After linking, `bin/DNAscent --version` should print `Version: 4.2.1` (a
TensorFlow warning about `libcudart.so.11.0` is expected/harmless on a CPU-only
run — it's just probing for a GPU that isn't there).

## 4. Output interpretation

### File locations

| Run | Path | Status |
|---|---|---|
| Interactive validation (gem) | `Eylul_gem_data/tailtest/gem.full.align` | Complete — 3488/3488 reads, 50 qualifying, used to validate the feature before the SLURM run (superseded by the SLURM output below, kept only as a cross-check — the two are byte-for-byte consistent, same 50 qualifying reads and 1,094,466 TAIL lines) |
| Interactive validation (control) | `Eylul_gem_data/tailtest/control.full.align` | **Incomplete** — killed early due to login-node contention, only ~3 qualifying reads. Superseded by the SLURM output below. |
| **Full SLURM run (gem)** | `Eylul_gem_data/slurm_full_run/gem/gem.full.align` | **Complete.** Job 1326383, 3m51s. 3488/3488 reads processed (95 failed rough event alignment), **50 qualifying reads**, 1,094,466 TAIL lines, 1.2 GB. |
| **Full SLURM run (control)** | `Eylul_gem_data/slurm_full_run/control/control.full.align` | **Complete.** Job 1326384, 12m43s. 4778/4778 reads processed (197 failed), **347 qualifying reads**, 12,586,861 TAIL lines, 1.8 GB. |
| SLURM stdout/stderr | `Eylul_gem_data/slurm_full_run/logs/{gem,control}_<jobid>.{out,err}` | Both clean — no errors, no crashes, no recurrence of the one-off segfault from §5. |
| Figures | `Eylul_gem_data/slurm_full_run/figures/*.png` | Generated by `plot_tail_signal.py` — see §6. |

The `tailtest/` directory is disposable session scratch — safe to delete now
that the SLURM outputs above are the confirmed, complete source of truth.

**Notable asymmetry**: control's qualifying rate is much higher than gem's
(347/4778 ≈ 7.3% vs 50/3488 ≈ 1.4%). Both datasets are gated by the same
3' soft-clip / full-length criteria, so this ~5x difference suggests the gem
prep specifically has more ligation artifacts (or a different degree of them)
than the control prep — worth asking your PI whether that's expected from the
library prep protocol, or whether it hints at something gem-specific (e.g. if
gemcitabine incorporation itself perturbs the ligation/adapter step). This is
an observation, not a conclusion — flagging it because it's the kind of thing
that's easy to miss if you only look at aggregate qualifying counts.

### Line format

Normal `.align` lines (unchanged): `<ref_coord>\t<kmer>\t<scaled_signal>\t<model_kmer>\t<flag>`.

New tail lines: `TAIL\t<idx>\t<scaled_signal>` — only 3 fields, first field is the
literal string `TAIL`, not a coordinate. `idx` is 0-based and sequential *within
that read's tail*, starting right after the last position the normal aligner
assigned. `<scaled_signal>` is in the same shift/scale-normalized units as the
`<scaled_signal>` column in normal lines (not raw pA) — i.e. directly comparable
across reads and to the rest of the file, but not physical current.

Existing coordinate-based parsers are already safe by construction if they check
field count before doing `int(fields[0])` (e.g.
`Eylul_gem_data/scripts/compare_align_signals.py` skips any line with fewer than
5 fields) — but don't rely on that being accidental; explicitly skip `TAIL` lines
in any new code that expects a numeric first field.

### Extracting tail signal per read (concrete example)

```python
def iter_read_tails(path):
    """Yield (read_header, [tail_values]) for every read in a .align file
    that has at least one TAIL line."""
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

for header, tail in iter_read_tails("Eylul_gem_data/slurm_full_run/gem/gem.full.align"):
    print(header, len(tail), tail[:5])
```
This is the same logic used by `Eylul_gem_data/scripts/plot_tail_signal.py`
(§6 / task 3) — see that file for the full, documented version.

### What "qualifying read" means, and how to filter for it

A read "qualifies" for tail capture if, at the point `eventalign()` was called,
`captureTailSignal` was `true` for it (§2, gate in `align_main`): forward-mapped,
alignment reaches the true reference end, no 3' soft clip. The `.align` file
itself doesn't record this flag explicitly — the practical, reliable proxy is:
**a read qualified and had recoverable signal if and only if it has at least one
`TAIL` line before the next `>` header** (exactly what `iter_read_tails` above
filters for via `if tail:`).

One edge case worth knowing: it's theoretically possible for a read to satisfy
all three gate conditions and still emit zero `TAIL` lines, if the rough event
alignment happened to already reach the very last raw event with nothing left to
recover. This has not been observed in the gem data (all 50 qualifying reads had
thousands of samples), but means "zero TAIL lines" is not proof a read *failed*
the gate — just that no extra signal was found either way. It doesn't matter in
practice for filtering (both cases are equally uninformative), but don't
over-interpret a `TAIL`-line count of exactly 0 as "definitely disqualified."

### Caveats and limitations

- **Qualifying-read yield is low: ~1.5% of gem reads (50/3425 that reached rough
  event alignment).** The dominant reason, confirmed by inspecting CIGAR strings
  directly, is that most reads carry a substantial 3' soft clip (one example: 107
  bp soft-clipped even though the aligned portion spans the full 615 bp
  reference) — i.e. most physical molecules in this prep have extra ligated
  material past the designed construct. The gate is working as designed; the low
  yield reflects the data, not a bug. See §5 for a suggested follow-up.
- **Tail length is large and highly variable**: median 17,081 raw samples, mean
  21,889, range 5,673–59,543 (gem, n=50). This is on the order of seconds of
  signal — much more than "a few extra bases' worth." The large majority of this
  comes from Source 2 (events the adaptive-banded aligner trimmed before ever
  attempting to align them), not Source 1 (terminal-window trailing insertions,
  which is bounded by one small HMM window and should only ever contribute a
  small fraction of a typical tail). Don't assume the whole tail is
  "post-terminus" signal in the same sense — the far end of a long tail is
  further from the reference-anchored region than the near end.
- **One segfault, not reproduced.** The very first test invocation of the freshly
  built binary crashed (`-t 4 -m 5`). Bisection (forcing `captureTailSignal =
  false`, then computing-but-discarding the real value, then restoring the real
  code path) never reproduced it, and the identical command succeeded 3/3 times
  afterward, plus in every subsequent larger-scale run (up to and including the
  full SLURM jobs). No correlation with the new code was found. Best guess: a
  transient artifact from first-execution disk-cache warming on a heavily
  contended shared node. Flagging rather than hiding it — if it recurs (check
  the SLURM `.err` logs for `Segmentation fault`), the next step would be
  installing `gdb` (not currently available) or enabling core dumps for a real
  backtrace.
- **Values are scale-normalized, not raw pA.** Fine for relative/statistical
  comparison, but if you ever want physical current values you'd need to reverse
  the per-read shift/scale (`r.scalings`), which isn't stored per-line in the
  output.
- The interactive `tailtest/control.full.align` is incomplete (killed early due
  to login-node contention, ~3 qualifying reads only) — this has since been
  superseded by the complete `slurm_full_run/control/` output (347 qualifying
  reads); don't use the `tailtest/` copy for anything beyond a sanity check.
- **Sample sizes are still modest for a real statistical claim**: 50 qualifying
  gem reads, 347 qualifying control reads. Enough for a first-look comparison
  (§6), not enough to claim a validated difference either way.

## 5. Recommended next steps

1. **Changepoint / segmentation analysis on the tail itself.** With tails this
   long and variable, a natural next question is *where within the tail* (if
   anywhere) does the signal look different from a stable baseline — right at
   the cutoff, or only near the true physical end of the read. A simple
   moving-window mean/variance shift detector (or reusing the vendored `scrappie`
   event-detection segmenter already in this repo, `src/scrappie/`) would be a
   reasonable first pass before anything more sophisticated.
2. **Investigate the 3' soft-clip rate.** ~98% of reads are currently excluded by
   the ligation-artifact gate. Worth checking whether this is inherent to the
   ligation chemistry, specific to this prep batch, or partially an artifact of
   dorado's adapter/primer trimming settings (the BAM's `@RG` line shows
   `tm:adapter,primer`, i.e. dorado already tried to trim some of this). If a
   meaningful fraction of the soft-clipped tail turns out to be a fixed-length
   adapter rather than genuine extra ligation, it might be possible to safely
   recover more reads.
3. **Anchor comparisons from the true 3' end, not from the cutoff.** Tail length
   varies enormously per read, but the physical 3' terminus (the very last raw
   sample of the read) is a fixed, comparable reference point across reads. The
   plotting script in §6 does this ("position from the end") in addition to
   the naive "position from the cutoff" view.
4. Once gem vs. control profiles are characterized well enough to see (or rule
   out) a consistent shift, revisit the original question of whether a simple
   classifier could detect Gem incorporation from this tail signal.
5. **Follow up on the qualifying-rate asymmetry** noted in §4 (control 7.3% vs
   gem 1.4%) — if that turns out to be gem-specific rather than a prep-batch
   effect, it's potentially interesting on its own.

## 6. Plotting script

`Eylul_gem_data/scripts/plot_tail_signal.py` reads the `.align` output(s) above
and produces two figures. It is self-contained and documented inline (module
docstring + a comment on every non-obvious function) — the summary here is
just enough to run it; read the script itself for the details.

**Environment**: the system `python3` (3.6.8) has neither matplotlib nor
numpy. Use the `gmconda_py368` conda env instead (matplotlib 3.3.4, numpy
1.19.5) — this is the same Python version the earlier `*_median_iqr.ipynb`
notebooks in this repo were built with:
```bash
/opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python \
    Eylul_gem_data/scripts/plot_tail_signal.py \
    --gem     Eylul_gem_data/slurm_full_run/gem/gem.full.align \
    --control Eylul_gem_data/slurm_full_run/control/control.full.align \
    --outdir  Eylul_gem_data/slurm_full_run/figures
```
`--control` is optional. Useful flags: `--n-examples` (default 12, how many
individual reads to show), `--example-window` (default 500 samples),
`--max-lag` (default 3000, how far out the median/IQR band plot goes),
`--min-n` (default 5, the band stops being drawn once fewer reads than this
remain at a given position — see below).

**Outputs**:
- `tail_signal_examples.png` — small multiples of individual raw tail traces
  (first `--example-window` samples from the cutoff) for a representative,
  fixed-seed random sample of gem reads. This is the "what does one read
  actually look like" view — every trace here is real recovered signal, not a
  summary statistic, and they show clear multi-level structure (consistent
  with genuine ionic-current state transitions, not noise).
- `tail_signal_mean_band.png` — two panels, gem vs control median ± IQR:
  left panel anchored at the cutoff (position 0 = first recovered sample),
  right panel anchored at the true 3' end of the read (position 0 = the very
  last raw sample). Both are shown because they answer different questions —
  "cutoff" shows what immediately follows the point the normal aligner gave
  up, while "end" is the one axis that's physically comparable across reads
  regardless of how much got trimmed before it (tail length varies ~10x
  across reads, so a fixed cutoff-relative position means something different
  read to read). The band is deliberately cut off wherever fewer than
  `--min-n` reads still contribute at that position, rather than drawing a
  misleading "confidence interval" from a couple of surviving reads at an
  extreme lag.

On the current full-dataset run (50 gem / 347 control qualifying reads), the
two medians track each other closely with substantial IQR overlap across most
of the range — no dramatic, obvious separation yet, though there are local
patches (e.g. gem dips further below baseline than control in the
cutoff-anchored view, roughly 100-500 samples in) that may be worth a closer,
statistically-grounded look via the changepoint analysis suggested above
rather than reading anything into by eye at this stage.

To rerun on updated data (e.g. after a bigger sequencing run), just re-point
`--gem`/`--control` at the new `.align` files — nothing else needs to change.
