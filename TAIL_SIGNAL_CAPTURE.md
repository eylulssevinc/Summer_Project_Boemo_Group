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

**Superseded 2026-07-10 — see §8.** Everything below describes the original
feature as a single on/off `bool captureTailSignal` that captured and printed
Source 1 and Source 2 together as combined `TAIL` lines. That's no longer how
the code works: `captureTailSignal` was replaced by a `TailCaptureMode` enum
and a `--tail-source {none,source1,source2,both}` CLI option, and Source 1/
Source 2 are now collected in separate buffers and printed as separate
`TAIL_SOURCE1`/`TAIL_SOURCE2` lines. This section is kept as-is for historical
context on *why* each piece of logic exists (the gate conditions, the two
recovery mechanisms, the reasoning behind each), which is all still accurate —
just read `bool captureTailSignal` below as "tail capture is not disabled"
and mentally replace the described single-buffer/single-line-format output
with the split version in §8.

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
| **Full SLURM run (gem)** | `Eylul_gem_data/slurm_full_run/gem/gem.full.align` | **Complete, regenerated 2026-07-10** (job 1337782, 2m50s, `-l 600 --tail-source both`) — see note below. 3488/3488 records processed (95 failed), **50 qualifying reads**, 1,094,466 Source 2 samples (0 Source 1), 1.2 GB. Byte-identical in substance to the original 2026-07-08 run (job 1326383); now source-labelled. |
| **Full SLURM run (control)** | `Eylul_gem_data/slurm_full_run/control/control.full.align` | **Complete, regenerated 2026-07-10** (job 1337783, 10m57s, `-l 600 --tail-source both`) — see note below. 4778/4778 records processed (197 failed), **347 qualifying alignment blocks / 321 unique reads** (see correction note below), 12,585,799 Source 2 + 1,062 Source 1 samples (12,586,861 combined — matches the original run's total exactly), 1.9 GB. |
| Pre-`--tail-source` backups | `Eylul_gem_data/slurm_full_run/pre_tailsource_backup/{gem,control}.full.align.combined_format_backup` | The original 2026-07-08 combined-`TAIL`-format files, kept for reference. Not needed for anything — the regenerated files above are an exact reproduction, just source-labelled — but kept rather than deleted. |
| SLURM stdout/stderr | `Eylul_gem_data/slurm_full_run/logs/{gem,control}_<jobid>.{out,err}` | Both clean — no errors, no crashes, no recurrence of the one-off segfault from §5. |
| Figures | `Eylul_gem_data/slurm_full_run/figures/*.png` | Generated by `plot_tail_signal.py` — see §6. |
| Figures (end-focused) | `Eylul_gem_data/slurm_full_run/figures_end_focused/*.png`, `*.tsv` | Generated by `plot_tail_signal_end_focused.py` — see §7. Regenerated 2026-07-10 from the files above; numbers unchanged (50 gem / 321 control, same sample counts) since the regeneration was an exact reproduction. |

**Important — do not confuse this with the §8 smoke test's 59/339 figures.**
The `Eylul_gem_data/smoketest_by_source/` files from §8 were generated with
the *default* `-l 100` (minimum aligned length), not this pipeline's `-l 600`.
That's a materially looser filter that admits thousands of short
supplementary-alignment fragments this pipeline is designed to exclude — see
§8's investigation. The smoke test's 59 gem / 339 control are **not** a
correction to 50/321; they're a different, non-comparable read population.
50 gem / 321 control (`-l 600`) remain the authoritative figures for this
project's actual biological question, confirmed by directly regenerating them
with the current codebase and getting an exact match.

The `tailtest/` directory is disposable session scratch — safe to delete now
that the SLURM outputs above are the confirmed, complete source of truth.

**Correction (found 2026-07-09, while building `plot_tail_signal_end_focused.py`):
"347 qualifying control reads" is actually a count of qualifying *alignment
blocks*, not unique physical reads.** The same physical molecule can produce
more than one BAM alignment record against this short 615 bp reference
(secondary/supplementary alignments — e.g. one block reaching only position
606 and another block for the *same* readID reaching the true end at 615,
both forward-mapped). 25 control readIDs each had 2-3 separate qualifying
blocks, accounting for exactly the 347 - 321 = 26 gap. **The correct unique-read
count for control is 321** (gem is unaffected: 50 blocks = 50 unique readIDs,
zero collisions). This matters because the read, not the alignment block or
TAIL line, is the meaningful biological/statistical unit here — see
`Eylul_gem_data/scripts/plot_tail_signal_end_focused.py`'s `load_tail_dict`,
which deduplicates by readID (keeping the longest recovered tail per read) and
prints a warning whenever this happens. Any earlier numbers below quoting
"347 qualifying control reads" refer to the block count; treat 321 as the
authoritative per-read figure going forward.

**Notable asymmetry**: control's per-read qualifying rate is still much higher
than gem's (321/4778 ≈ 6.7% vs 50/3488 ≈ 1.4%). Both datasets are gated by the
same 3' soft-clip / full-length criteria, so this ~5x difference suggests the
gem prep specifically has more ligation artifacts (or a different degree of
them) than the control prep — worth asking your PI whether that's expected
from the library prep protocol, or whether it hints at something gem-specific
(e.g. if gemcitabine incorporation itself perturbs the ligation/adapter step).
This is an observation, not a conclusion — flagging it because it's the kind
of thing that's easy to miss if you only look at aggregate qualifying counts.

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
  superseded by the complete `slurm_full_run/control/` output (321 unique
  qualifying reads, see the correction note above); don't use the `tailtest/`
  copy for anything beyond a sanity check.
- **Sample sizes are still modest for a real statistical claim**: 50 qualifying
  gem reads, 321 qualifying control reads (unique-read counts — see the
  correction note above). Enough for a first-look comparison (§6), not enough
  to claim a validated difference either way.

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

On the current full-dataset run (50 gem / 347 control qualifying *alignment
blocks* — see the correction note in §4: this script pools by block, not by
deduplicated unique read, so control's true unique-read count is 321; the
difference doesn't change the qualitative picture below, but see §7 for the
version of this analysis that dedupes properly), the two medians track each
other closely with substantial IQR overlap across most of the range — no
dramatic, obvious separation yet, though there are local patches (e.g. gem
dips further below baseline than control in the cutoff-anchored view, roughly
100-500 samples in) that may be worth a closer, statistically-grounded look
via the changepoint analysis suggested above rather than reading anything into
by eye at this stage.

To rerun on updated data (e.g. after a bigger sequencing run), just re-point
`--gem`/`--control` at the new `.align` files — nothing else needs to change.

## 7. End-focused plotting script (read-level, deduplicated)

`Eylul_gem_data/scripts/plot_tail_signal_end_focused.py` is a second, later
analysis script that supersedes §6's cutoff-anchored comparison for anything
biological: since gemcitabine is expected right at the physical 3' end of the
molecule (not necessarily near wherever the aligner happened to stop), every
plot here anchors from the true 3' end, and — unlike `plot_tail_signal.py` —
explicitly deduplicates by readID (keeping the longest recovered tail per read
when a readID has multiple qualifying alignment blocks, see §4's correction
note), so the read is genuinely the unit of analysis throughout, not the
alignment block or the individual TAIL line.

```bash
/opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python \
    Eylul_gem_data/scripts/plot_tail_signal_end_focused.py \
    --gem     Eylul_gem_data/slurm_full_run/gem/gem.full.align \
    --control Eylul_gem_data/slurm_full_run/control/control.full.align \
    --outdir  Eylul_gem_data/slurm_full_run/figures_end_focused
```

Run via SLURM (`Eylul_gem_data/scripts/run_end_focused_plots.slurm`, `cosbi`
partition/account/QOS — a second CPU partition available to this project
alongside `kutem`, see `reference-dnascent-infra` in memory) rather than
directly on the login node: it parses the full ~3 GB of `.align` data across
two passes per file. In practice this finishes in under a minute once a node
is available.

**Outputs** (in `Eylul_gem_data/slurm_full_run/figures_end_focused/`):
- `tail_end_anchored_{100,500,1000,3000}.png` — gem vs control median ± IQR,
  end-anchored, at four different zoom levels near the physical 3' end.
- `gem_individual_end_examples.png`, `control_individual_end_examples.png` —
  individual raw traces, last 500 samples before the physical end, for a
  fixed-seed random sample of reads from each group.
- `tail_end_read_level_summary.tsv` — one row per unique read (50 gem + 321
  control = 371 rows) with per-read median/mean/std over the last
  100/500/1000/3000 samples, plus the fraction of the last 500 samples
  above +1 / below -1 in scaled-signal units.
- `tail_end_read_level_boxplots.png` — read-level (not sample-pooled) gem vs
  control boxplots for six of those summary metrics.
- `tail_end_vs_regular_end_summary.tsv`, `tail_end_vs_regular_end_boxplot.png`
  — per read, the median scaled signal from the last regular,
  reference-anchored `.align` lines (coordinates 600-611, immediately
  upstream of the tail cutoff) vs. the median of the last 500 tail samples,
  and their delta, gem vs control.

Same scientific caveats as §6 apply, plus one specific to this script: an
end-anchored raw-sample position is not a base position (a single base/k-mer
dwell can span many raw samples, and dwell time itself varies), so "x=1" is
not "the second-to-last base" even though gemcitabine is expected there — see
the script's module docstring for the full caution, and don't read a
gem-vs-control difference off band overlap/non-overlap alone.

## 8. Source-selectable tail capture (`--tail-source`, added 2026-07-10)

Everything in §2-§7 above treated Source 1 and Source 2 as one combined pool,
because that's how the original feature worked: both fed a single `tailSignal`
buffer and were printed together as undifferentiated `TAIL` lines. That made
it impossible to ask whether the two mechanisms behave differently — which
matters here because they're not equivalent:

- **Source 1** (terminal-window trailing insertions) comes directly out of
  the final reference-anchored Viterbi HMM window (§2/§3's `builtinViterbi`
  fine-alignment stage). It's the more directly reference-tied signal.
- **Source 2** (rough-aligner-trimmed events) never reached that stage at all
  — it's whatever the rough, basecall-level aligner (`event_handling.cpp`,
  `adaptive_banded_simple_event_align`) gave up on before eventalign ever ran.
  Real signal, but one step further removed from the reference. It's also
  typically much larger and more variable in length than Source 1 (Source 1
  is bounded by one small HMM window; Source 2 can run to tens of thousands
  of raw samples — see the gem/control length distributions in §4).

Since the current dataset's Source 1 signal may turn out to be too sparse to
say much about on its own, but a larger dataset is expected soon, the code now
supports selecting which source(s) get printed, so this doesn't need
revisiting later.

### What changed

**`src/alignment.h`** — added:
```cpp
enum class TailCaptureMode { NONE, SOURCE1, SOURCE2, BOTH };
```
and changed `eventalign()`'s declaration from
`void eventalign(DNAscent::read &, unsigned int, bool)` to
`void eventalign(DNAscent::read &, unsigned int, TailCaptureMode)`.

**`src/alignment.cpp`**:
- New CLI option **`--tail-source <none|source1|source2|both>`**, parsed by a
  new `parseTailSourceArg()` and stored on `Arguments.tailSourceMode`.
  **Default is `none`** — running `DNAscent align` without this flag now
  captures no tail signal at all and behaves like stock DNAscent, unlike the
  old code where tail capture was unconditionally on for any read passing the
  gate. This was a deliberate choice: normal `align` behavior should not
  change unless a user explicitly opts in.
- The read-level gate in `align_main()` is unchanged in substance (still
  forward-mapped + reaches the true reference end + no 3' soft clip), just
  renamed to `readPassesTailGate` for clarity, and combined with the
  requested mode:
  ```cpp
  bool readPassesTailGate = (r.strand == "fwd")
                            and (r.refEnd == (int) reference.at(r.referenceMappedTo).size())
                            and not hasSoftClipAtReferenceThreePrime(r.record);
  TailCaptureMode effectiveTailMode = readPassesTailGate ? args.tailSourceMode : TailCaptureMode::NONE;
  eventalign(r, Pore_Substrate_Config.windowLength_align, effectiveTailMode);
  ```
  So a read that fails the gate gets `NONE` regardless of `--tail-source`,
  and a read that passes it gets whatever the user asked for.
- Inside `eventalign()`, the single `tailSignal` buffer was split into
  `tailSignalSource1` and `tailSignalSource2`. Two local booleans,
  `wantSource1`/`wantSource2` (derived from `tailMode`), gate both collection
  and printing per source — e.g. under `--tail-source source1`, Source 2's
  collection loop (the one that can walk tens of thousands of trailing raw
  events) never runs at all, rather than collecting-then-discarding.
- **`src/detect.cpp`** and **`src/trainCNN.cpp`**: their three `eventalign(...)`
  call sites now pass `TailCaptureMode::NONE` instead of the old literal
  `false` — behavior unchanged, just updated to match the new signature.

### Output format

```
TAIL_SOURCE1\t<idx>\t<scaled_value>
TAIL_SOURCE2\t<idx>\t<scaled_value>
```
`idx` is 0-based and sequential *within that source, for that read* (i.e.
Source 1 and Source 2 each restart at 0). Under `--tail-source both`, Source 1
lines are printed before Source 2 lines within a read's block — the same
chronological order the old combined `TAIL` output used. Under `none`,
`source1`, or `source2`, only the relevant line type (or none) appears.
Existing parsers that skip any line not starting with a numeric first field
are unaffected — `TAIL_SOURCE1`/`TAIL_SOURCE2` are just as non-numeric as the
old `TAIL` prefix was. Old `.align` files generated before this change still
use the unlabelled `TAIL` format; nothing here rewrites them, and
`plot_tail_signal_by_source.py` (below) reads that format too, as a "legacy"
pseudo-source (Source 1 and Source 2 can't be told apart in it).

### Terminal summary

At the end of every `align` run (both the normal end-of-BAM exit and the
early exit under `-m`/`--maxReads`), a summary now prints:
```
Tail-source capture summary (--tail-source <mode>):
  reads passing tail gate (forward-mapped, reaches true ref end, no 3' soft clip): <n>
  reads with Source 1 signal printed this run: <n>, total Source 1 samples: <n>
  reads with Source 2 signal printed this run: <n>, total Source 2 samples: <n>
  Source 1 tail length per read: <n> reads, min=<> median=<> mean=<> max=<>
  Source 2 tail length per read: <n> reads, min=<> median=<> mean=<> max=<>
```
"reads passing tail gate" reflects the gate only, independent of `--tail-source`.
The Source 1/Source 2 counts, however, reflect *what this run actually printed*
— under `--tail-source source1`, the Source 2 numbers will always read zero
even for reads that do have real Source 2 signal, because it was never
collected that run. This is intentional (collection is skipped entirely for
an unrequested source, for efficiency — see above), but means these numbers
are not a substitute for a `--tail-source both` run if you want to know what
exists in principle for a read.

### Smoke test results

**Note (added 2026-07-10, after this section was first written): the numbers
below used `-l 100` (this tool's default), NOT this project's `-l 600`.**
That difference alone explains why these figures (59 gem / 339 control) don't
match §4's authoritative 50 gem / 321 control — it is not a correction to
those numbers, just a different, looser read population. §4 has since been
regenerated directly with `-l 600 --tail-source both` and reproduces 50/321
exactly, with **0 Source 1 reads for gem and 6 for control** (1,062 samples,
max 990 in one read) — the authoritative Source 1 figures for this project,
slightly different from this section's 0/7 below because the `-l 600`
population differs slightly from the looser `-l 100` one. Everything below is
still useful as a *mode-switching correctness*
test (that part doesn't depend on `-l`) and as the investigation that found
the Source 1 dedup bug (fixed, see below) — just don't read the 59/339 counts
as datapoints about the actual gem/control comparison.

Run via `Eylul_gem_data/scripts/run_tailsource_smoketest.slurm` (`cosbi`
partition): all four modes against the full gem BAM (job 1336986, 38m25s),
then `--tail-source both` only against the full control BAM (job 1337344,
39m52s — the four-mode sweep was already validated on gem, so only the mode
that matters was re-run for control).

**What "qualifying records" means, precisely.** `DNAscent align`'s BAM buffer
accepts *any* record — primary, secondary, or supplementary — that
independently satisfies `-q/--quality` (mapping quality ≥ 20 by default) and
`-l/--length` (aligned length ≥ 100bp) with nonzero query length; this is
pre-existing behavior, unrelated to this session's changes. A small
standalone flagstat tool (linked against this repo's vendored htslib) was
built to check the composition directly:

```
              total    unmapped  secondary  supplementary  primary-mapped
gem BAM:      35,849    4,854      56         14,675          16,264
control BAM:  39,303    5,720      91         18,771          14,721
```

Both BAMs have roughly as many supplementary records as primary-mapped ones —
this dataset has a high rate of split/chimeric alignment (dorado representing
one physical read as multiple BAM records, e.g. when it can't cleanly align
the whole read as one contiguous span). The quality/length filter admits a
large fraction of these independently of whether they're primary — hence
countRecords reporting ~18,608 (gem) / ~23,967 (control) qualifying
*records*, not qualifying *reads*. All four gem modes processed the exact
same 18,608/18,608 records with an identical 2,578 failed (rough-alignment
QC) count, confirming `--tail-source` has zero effect on which records get
processed or which fail — nothing is silently skipped differently between
modes.

**Results (unique reads, after deduplication by read_id — see the dedup note
below):**

```
                    reads passing tail gate   Source 1 reads (samples)     Source 2 reads (samples)
gem   (--tail-source both):   59 blocks = 59 unique reads (0 collisions)   0 (0)                        59 (1,172,252; len min=2331 median=16243 max=59543)
control (--tail-source both): 365 blocks -> 339 unique reads (26 collisions)  7 (1,084; len min=5 median=22 max=990)   339 (11,796,270; len min=4294 median=19966 max=3618689)
```

`grep -c` counts on the raw (block-level, pre-dedup) `.align` files matched
the terminal summary exactly in every case (gem: `TAIL_SOURCE2`=1,172,252,
`TAIL_SOURCE1`=0; control: `TAIL_SOURCE2`=13,031,431 raw /
11,796,270 after dedup, `TAIL_SOURCE1`=1,084 raw, unaffected by dedup since
none of the 26 collisions involved two blocks both carrying Source 1 data).

**Is Source 1 = 0 for gem a bug?** No — verified two ways. (1) Non-trailing
`I`-state (insertion) lines — printed as ordinary `.align` lines with a
`NNNNNNNNN` placeholder k-mer field, unaffected by `--tail-source` — occur
**3,282,037 times** in gem's output alone, so the HMM's insertion mechanism
plainly fires constantly throughout the alignment; it is specifically the
*trailing*-insertion-in-the-*terminal*-window pattern (Source 1's trigger)
that's rare. (2) The HMM's fixed transition probabilities
(`src/config.h`'s `HMM_TransitionProbs_DNA_R10`) set `internalM2I = 0.001`
against a "stay" (dwell) probability that's typically close to 1 — i.e. the
model is tuned to explain extra events at a fixed position via repeated
"stay" almost always, and only ever prefers "insert" when stay's emission
score gets bad enough. Whether that happens in a read's specific terminal
window is data-dependent per read, not something the code can force. Control
independently having 7 real Source 1 reads (up to 990 samples) confirms the
mechanism is fully functional; gem's 59 gate-passing reads simply didn't
happen to produce that pattern in this dataset. Larger datasets may or may
not change this — it isn't guaranteed to scale with read count the way
Source 2 does.

**Dedup bug found and fixed while investigating this.** `plot_tail_signal_by_source.py`'s
`load_read_blocks` originally deduplicated colliding blocks by picking the
one with the larger *total* signal (source1+source2+legacy summed) and
discarding the other whole. For control read `e1cd3406-9f89-4862-9055-ca6e242e4e69`,
one colliding block had `{source1=990, source2=405202}` and the other had
`{source1=0, source2=409401}` — only 0.8% more total signal — so the old rule
silently discarded the single largest Source 1 signal in the entire dataset.
Fixed to keep the longer array **per channel independently** across colliding
blocks, rather than picking one whole block (see the function's current
docstring). All 7 control Source 1 reads (not 6) are now retained correctly.

Manually inspected one qualifying gem read's block in the `both` output
(`ea67543f-9e40-4e6c-bcb9-ae15ba02d967`): normal reference-anchored lines are
untouched up through coordinate 610, immediately followed by
`TAIL_SOURCE2\t0\t...` through `TAIL_SOURCE2\t14654\t...` (14,655 lines,
0-based and sequential), no `TAIL_SOURCE1` lines, no interleaving/corruption.

Note: this smoke test's full-dataset output files are large (gem ~4.5 GB × 4
modes, control ~6 GB for the one mode run, ~35 GB total) and live in
`Eylul_gem_data/smoketest_by_source/` — disposable once the numbers above are
trusted, not wired into any other script's default paths.

### Plotting script

`Eylul_gem_data/scripts/plot_tail_signal_by_source.py` is a new script (does
not overwrite `plot_tail_signal.py` or `plot_tail_signal_end_focused.py`) that
reads `--tail-source both` (or `source1`/`source2`-only, or legacy combined)
`.align` output and, for each of Source 1, Source 2, and a combined view,
produces: a cutoff-anchored median±IQR band, an end-anchored median±IQR band,
individual example traces for both gem and control, and rows in a shared
read-level summary TSV (`tail_by_source_read_level_summary.tsv`, columns:
sample, source, read_id, tail_length, median/mean/std scaled signal). Like
`plot_tail_signal_end_focused.py`, it deduplicates by read_id — but unlike
that script, on a collision it keeps the longer array **per channel**
(source1/source2/legacy) independently rather than picking one whole
colliding block, since two blocks for the same read can be complementary
rather than redundant (see the smoke-test dedup bugfix above for the concrete
case that motivated this). Run it the same way as the other plotting scripts,
with the `gmconda_py368` conda env:
```bash
/opt/ohpc/pub/compiler/conda3/latest/envs/gmconda_py368/bin/python \
    Eylul_gem_data/scripts/plot_tail_signal_by_source.py \
    --gem     <path to a --tail-source both .align file> \
    --control <path to a --tail-source both .align file> \
    --outdir  <output dir>
```
Combined-view plots mix two different recovery mechanisms and should be
treated as a continuity check against the older, unsplit output rather than
the primary analysis — prefer the per-source views.

## 9. `NEXT100EVT`: a fixed-width, mechanism-agnostic tail window (added 2026-07-15)

### Motivation

The PI proposed a simpler, implementation-agnostic recovery method, independent
of the Source 1 / Source 2 split above: take the read's full per-read event
vector (`r.events`, downsampled raw signal, one entry per `event{mean, raw[]}`
— *not* to be confused with individual raw ADC samples), find the index of the
last event the alignment pipeline (adaptive-banded rough alignment, then the
reference-anchored Viterbi HMM) ever aligned to anything, and grab the next 100
events after that index directly from the vector — regardless of *why* the
pipeline stopped there.

**Correction (same day, after the PI reviewed a first version of this
feature)**: my first implementation anchored this at
`lastAlignedEvent = r.eventAlignment.back().first` — the *rough aligner's own*
last-touched event, i.e. literally Source 2's boundary — reasoning that since
Viterbi can never see an event the rough aligner didn't already hand it, this
was "the last event the whole pipeline aligned." The PI's actual, more precise
ask (with exact line-number references into upstream DNAscent, confirmed by
fetching `reads.h`/`alignment.cpp` at the cited commit) is **not** that boundary
— it's the raw index of the last event the HMM specifically **Match**-labelled
against a reference k-mer, which is `lastM_ev` (§3.7/§4.4's existing
last-Match-finding loop) translated back into a `r.events` index. These two
anchors coincide *only* when a read has zero Source 1 signal (true for all 50
gem reads); for the 6 control reads that do have Source 1 signal, the old
anchor was strictly *later* than the PI's — it skipped straight past whatever
got labelled a trailing insertion, when it should have started with it. Fixed
by adding a new parallel-tracking array and reusing the existing PI-provided
`eventalign()` per-window computation. See `PROJECT_OVERVIEW.md` §4.6 for the
full mechanism explanation (`eventIndeces`, `lastMatchRawIdx`) — this section
now only records the operational/testing side.

### What changed

`eventalign()` in `src/alignment.cpp` builds a new parallel array,
`eventIndeces` (same length as `eventSnippet`/`eventSnippet_means`, pushed to
at the same point, holding each entry's raw index into `r.events`), and uses
`eventIndeces[lastM_ev]` from the **last** window processed to set
`lastMatchRawIdx` (a variable declared before the outer per-window loop,
unconditionally overwritten every iteration, so it holds the final window's
value once the loop ends). It then collects a third, independent pool,
`tailSignalNext100Events`, whenever any `--tail-source` mode other than `none`
is requested (reuses the existing `wantSource1 or wantSource2` gate — no new
CLI flag):

```cpp
size_t stop = std::min( r.events.size(), lastMatchRawIdx + 1 + 100 );
for ( size_t evi = lastMatchRawIdx + 1; evi < stop; evi++ ){
    double event_mean = r.events[evi].mean;
    if ( not (0. < event_mean and event_mean < 250.) ) continue;
    tailSignalNext100Events.push_back( (event_mean - r.scalings.shift) / r.scalings.scale );
}
```

One value per **event** (its `.mean`, not its expanded `.raw` samples), capped
at the next 100 events after `lastMatchRawIdx`, fewer only if the read's own
event vector runs out first. Printed as `NEXT100EVT\t<idx>\t<val>` lines,
immediately after Source 1/Source 2 for that read (`lastMatchRawIdx` can be
*earlier* than Source 2's own boundary for the rare reads with real Source 1
signal — see verification below — even though the printed line order is always
Source 1 → Source 2 → NEXT100EVT regardless).

**Deliberately does not start with `TAIL`**: `plot_tail_signal.py` and
`plot_tail_signal_end_focused.py` both key off a bare `line.startswith("TAIL")`
to mean "legacy combined tail signal," and would otherwise silently fold this
new pool into that bucket, double-counting samples already present in Source 1/2.
Because of this, regenerating the canonical `.align` files with the rebuilt
binary required **no changes and no reruns** of `figures_end_focused/` — those
scripts simply don't see `NEXT100EVT` lines at all.

The terminal summary (`printTailCaptureSummary`) gained a matching third block
(`countTailLines` now returns `n3` alongside `n1`/`n2`).

### Smoke test

`Eylul_gem_data/scripts/run_next100evt_smoketest.slurm` (`-l 600 -m 500
--tail-source both`, both datasets): 9 gem / 59 control reads passed the tail
gate; **every one** got exactly 100 `NEXT100EVT` samples (min=median=mean=max=100
in both datasets) — every gate-passing read had at least 100 events of
recoverable signal available past `lastMatchRawIdx`, consistent with Source 2
alone always running into the thousands of samples per read. Unaffected by the
anchor correction below (total sample counts don't reveal *which* raw events
were selected — only content-level inspection does).

**Anchor-shift verification** (the check that actually matters, since file
line order is always Source 1 → Source 2 → NEXT100EVT regardless of which
underlying events NEXT100EVT draws from — so *position in the file* doesn't
prove *which* anchor was used). Added a temporary debug print comparing
`lastMatchRawIdx` (the corrected anchor) against `r.eventAlignment.back().first`
(Source 2's anchor, i.e. what the first, since-corrected version of this
feature used) for control read `e1cd3406-9f89-4862-9055-ca6e242e4e69` (the
same duplicate-block readID from the §8 dedup bugfix, still present as two
separate blocks):

| block | Source 1 samples | `lastMatchRawIdx` | Source 2's anchor | difference |
|---|---|---|---|---|
| 1 | 0 | 50434 | 50434 | 0 |
| 2 | 990 | 51082 | 51275 | 193 events |

Block 1 (no Source 1 signal): both anchors identical, as expected — nothing
for them to disagree about. Block 2 (990 Source 1 *samples*): the two anchors
differ by 193 raw *events* (990 samples / 193 events ≈ 5.1 samples/event, a
plausible dwell time) — confirming `lastMatchRawIdx` sits *before* Source 2's
boundary by exactly the span of that block's Source 1 signal, exactly as
predicted. Debug print removed after verification; binary rebuilt clean.

### Canonical regeneration

Full-dataset canonical files (`slurm_full_run/{gem,control}/*.full.align`)
regenerated **twice** with the rebuilt binary (same `-l 600 --tail-source
both` as before — see §8's canonical-vs-smoke-test caveat, which still applies
unchanged): once right after the first (since-corrected) version of this
feature, and again after the anchor fix above — the numbers below are from the
final, corrected regeneration. Pre-regeneration copies backed up to
`slurm_full_run/pre_next100evt_backup/` (taken before the *first*
regeneration, i.e. the true pre-`NEXT100EVT` baseline). Confirmed **exact
reproduction** of every pre-existing number, in both regenerations (50 gem /
347 control qualifying reads, 0/6 Source 1 reads, identical Source 1/2 sample
counts down to the last digit) — this only *adds* `NEXT100EVT` lines; every
previously-existing line in these files is
byte-identical to before.

New in this regeneration: **every one** of the 50 gem and 347 control
qualifying reads got exactly 100 `NEXT100EVT` samples (5,000 and 34,700 total
respectively) — no read ran out of events before reaching the 100-event cap,
consistent with Source 2 alone never dropping below ~5,000 raw samples per
read in either dataset.

## 10. `ALLAFTERM`: the uncapped pool, plus a full read-level analysis (added 2026-07-15)

### What changed and why

`NEXT100EVT` (§9) is a fixed, 100-event window — useful for a uniform,
comparable-across-reads view, but it discards everything past event 100, and
every qualifying read in this dataset turned out to have well over 100 events
available (`NEXT100EVT` never got truncated by running out of data — see §9).
`ALLAFTERM` removes the cap: it is **every valid event in `r.events` from
`lastMatchRawIdx + 1` to `r.events.size()`**, i.e. `NEXT100EVT`'s full,
uncapped superset. `NEXT100EVT` itself is **unchanged** — same anchor, same
collection code, same output lines as before.

`src/alignment.cpp` changes:
- **Collection** (lines 996–1016), immediately after `NEXT100EVT`'s own
  collection block, reusing the identical `lastMatchRawIdx`/`haveLastMatchRawIdx`
  computed once per read (§9) — no new anchor logic:
  ```cpp
  std::vector<double> tailSignalAllAfterM;
  if ( (wantSource1 or wantSource2) and haveLastMatchRawIdx and lastMatchRawIdx < r.events.size() ){
      for ( size_t evi = lastMatchRawIdx + 1; evi < r.events.size(); evi++ ){
          double event_mean = r.events[evi].mean;
          if ( not (0. < event_mean and event_mean < 250.) ) continue;
          double scaledEvent = (event_mean - r.scalings.shift) / r.scalings.scale;
          tailSignalAllAfterM.push_back(scaledEvent);
      }
  }
  ```
  The extra `lastMatchRawIdx < r.events.size()` guard (beyond the existing
  `haveLastMatchRawIdx`) is defensive: it makes the invariant "never compute
  `lastMatchRawIdx + 1` from an out-of-range value, never index an empty
  `r.events`" explicit, rather than relying on the for-loop's own bound
  (`evi < r.events.size()`) to make an already-out-of-range start a silent
  zero-iteration no-op. Gated on the exact same `(wantSource1 or wantSource2)`
  condition as `NEXT100EVT` — no new CLI flag, still off by default under
  `--tail-source none`.
- **Output** (lines 1043–1051), printed immediately after `NEXT100EVT` for the
  same read: `ALLAFTERM\t<idx>\t<val>` lines. Deliberately does **not** start
  with `TAIL`, for the identical reason `NEXT100EVT` doesn't (§9) — existing
  scripts keying off a bare `"TAIL"` prefix would otherwise silently fold this
  pool in as more "legacy" signal.
- **Terminal summary**: `countTailLines` (line 199) gained a fourth counter
  `n4`; `printTailCaptureSummary` (lines 250–277) gained block-level count,
  **unique-read-ID count** (`std::set<std::string> uniqueReadIDsWithAllAfterM`,
  line 1115, accumulated inside the same `#pragma omp critical` section as
  every other counter — thread safety unchanged), total event count, and
  min/median/mean/max events-per-alignment-block.
- **Unchanged, verified**: ordinary reference-anchored output, Source 1/2
  collection, `NEXT100EVT` itself, `r.QCpassed` / read QC, the forward/full-
  reference/no-3′-soft-clip gate (§8), and `src/detect.cpp`/`src/trainCNN.cpp`
  (both still pass `TailCaptureMode::NONE` literally, so `wantSource1`/
  `wantSource2` are always false there and this new code never executes for
  either tool) — see validation below for how this was actually checked, not
  just asserted.

### Output semantics — read this before using the data

- **`NEXT100EVT` is exactly the first up-to-100 valid entries of `ALLAFTERM`**,
  unless the read has fewer than 100 remaining valid events (in this dataset,
  none did — every qualifying read had ≥100). Verified exactly, per read
  block, not just claimed — see validation §1 below.
- **`ALLAFTERM` can and typically does span both Source 1 and Source 2
  territory.** It is defined purely by position relative to `lastMatchRawIdx`,
  with no reference to which mechanism (terminal-window trailing insertion vs.
  rough-aligner trim) produced a given event — unlike Source 1/Source 2, which
  are split by mechanism.
- **Granularity differs by pool.** Source 1/Source 2 are printed at *expanded
  raw-ADC-sample* granularity (several lines per detected event, from each
  event's `.raw[]`). `NEXT100EVT` and `ALLAFTERM` are printed at *event-mean*
  granularity (one line per detected event, its `.mean`, never expanded). A
  "length" in one space is not a length in the other; never compare them
  directly.
- **Neither `NEXT100EVT` nor `ALLAFTERM` proves biological origin.** This is
  recovered signal past the last position the reference-anchored HMM modelled
  — it may include pore adapter signal, motor-protein signal, ordinary dwell
  variability, or (for reads that pass the 3′ gate, which already excludes
  soft-clipped/ligation-artifact reads) genuine non-reference sequence. A
  plot showing recovered signal is not itself evidence that recovered signal
  corresponds to the gemcitabine position specifically.
- **Events already labelled Match remain excluded from both pools** and stay
  exactly where they always were: ordinary, reference-anchored `.align`
  output, untouched by any of this.

### Validation

Ran `Eylul_gem_data/scripts/validate_allafterm.py` (new script, checks 1/3/4/5
below) plus two one-off checks (2/6/7) not built into that script:

1. **`NEXT100EVT` == first up-to-100 valid `ALLAFTERM` entries, every block.**
   PASS on the smoke test (9 gem / 59 control blocks with data) and the full
   canonical run (50 gem / 347 control blocks with data) — zero mismatches in
   either case.
2. **Existing normal/Source 1/Source 2 output byte-identical to pre-change
   output, modulo new `ALLAFTERM` lines.** A naive whole-file `diff` on the
   full-dataset output was NOT usable for this — it never finished, because
   **the order read blocks are written in is not deterministic across
   separate runs** (an existing property of `align_main`'s OpenMP-parallel
   loop: whichever thread's `#pragma omp critical` section runs first writes
   first, and thread scheduling isn't guaranteed identical run-to-run — this
   is pre-existing, unrelated to `ALLAFTERM`, and was confirmed directly:
   both runs' files contain the exact same *set* of read headers, just in a
   different *order*). Verified instead with a small order-independent
   script (per-read-ID, sorted multiset of SHA-256 block hashes, ALLAFTERM
   lines stripped from the new file before hashing — handles the
   documented duplicate-block-per-read-ID case too, since a read_id's several
   blocks are compared as a set, not by position): **zero content mismatches**
   across all 3,002 gem / 3,764 control unique read_ids common to both the
   pre-`ALLAFTERM` canonical files (`slurm_full_run/{gem,control}/*.full.align`)
   and the new ones (`slurm_full_run_allafterm/{gem,control}/*.full.align`).
3. **`ALLAFTERM` indices start at 0 and are contiguous per block.** PASS,
   both smoke test and full run.
4. **No `ALLAFTERM` for gate-failing reads.** Structurally guaranteed by
   `effectiveTailMode = readPassesTailGate ? args.tailSourceMode :
   TailCaptureMode::NONE` (§8, unchanged) — a gate-failing read gets `NONE`,
   which makes `wantSource1`/`wantSource2` both false, which is `ALLAFTERM`'s
   own collection gate. Checked as a proxy (every block with `ALLAFTERM` also
   has `NEXT100EVT` and vice versa, since they share one gate) rather than
   re-deriving the gate itself: PASS, zero orphan cases in either dataset.
5. **Empty/short remainders don't crash.** 3,343/3,393 gem blocks and
   4,234/4,581 control blocks in the full canonical files have *zero*
   tail-prefixed lines at all (reads that never reached `eventalign` with a
   passing gate, or extremely short remainders) — `validate_allafterm.py`
   parsed all of them without exception.
6. **Duplicate read IDs reported and handled deterministically by the
   plotting script.** `plot_allafterm_analysis.py`'s dedup rule (keep the
   block with the longer `ALLAFTERM` array) is applied and reported — see
   its run log below for exact block/unique-read/collision counts.
7. **Plotted sample sizes refer to unique read IDs, not alignment blocks.**
   Every legend/title in the new script reports `n` from the deduplicated
   per-read dictionaries, and the terminal-summary unique-read-ID count
   (independently computed in C++, via `uniqueReadIDsWithAllAfterM`) matches
   the Python-side dedup count exactly (321 control unique reads, both sides)
   — an incidental but reassuring cross-check that both independently-written
   dedup mechanisms agree.

### Commands used

```bash
# Build (cosbi partition, see Eylul_gem_data/scripts/... smoke test scripts for the pattern)
make -j8 LDFLAGS="-L$HOME/lib -ldl -llzma -lbz2 -lm -lz"

# Smoke test (500-read cap, both datasets)
Eylul_gem_data/scripts/run_next100evt_smoketest.slurm-style invocation, -l 600 -m 500 --tail-source both

# Full canonical run -- NEW directory, existing slurm_full_run/ untouched
Eylul_gem_data/slurm_full_run_allafterm/run_tailcapture_allafterm.slurm
  sbatch --export=ALL,DATASET=gem     run_tailcapture_allafterm.slurm   # job 1359705
  sbatch --export=ALL,DATASET=control run_tailcapture_allafterm.slurm   # job 1359706

# Validation
python3 Eylul_gem_data/scripts/validate_allafterm.py --align <file>

# Plotting
Eylul_gem_data/slurm_full_run_allafterm/run_allafterm_plots.slurm       # job 1359779
```

### Counts (full canonical run, `-l 600 --tail-source both`, identical filter to §8-9)

| | gem | control |
|---|---|---|
| reads passing tail gate | 50 | 347 (321 unique) |
| Source 1 reads (samples) | 0 (0) | 6 (1,062) |
| Source 2 reads (samples) | 50 (1,094,466) | 347 (12,585,799) |
| NEXT100EVT reads (samples) | 50 (5,000) | 347 (34,700) |
| ALLAFTERM blocks (unique reads) | 50 (50) | 347 (321) |
| ALLAFTERM total events | 212,484 | 2,476,953 |
| ALLAFTERM events/block: min/median/mean/max | 1,134 / 3,346 / 4,249.7 / 11,499 | 1,042 / 3,950 / 7,138.2 / 725,085 |

Every pre-existing number (reads passing gate, Source 1/2, NEXT100EVT) is an
**exact match** to §8-9's figures, confirming the ALLAFTERM addition changed
nothing about existing behavior.

### Figures and TSVs

See `PROJECT_OVERVIEW.md` §4.7 for the full walkthrough of what each figure
shows and the scientific-limitations discussion; paths below, all under
`Eylul_gem_data/slurm_full_run_allafterm/figures_allafterm/`:

- Part A (HMM-Match-anchored): `allafterm_hmm_anchored_first100.{png,tsv}`,
  `allafterm_hmm_anchored_first500.{png,tsv}`
- Part B (physical-end-anchored): `allafterm_physical_end_last{100,500,3000}.{png,tsv}`
- Part C (read-level): `allafterm_read_level_summary.tsv`,
  `allafterm_read_level_boxplots_{25,50,100}.png`,
  `allafterm_first100_median_boxplot.png`
- Part D (remaining-length distribution): `allafterm_remaining_length_hist.png`,
  `allafterm_remaining_length_hist_log.png`, `allafterm_remaining_length_ecdf.png`,
  `allafterm_remaining_length_boxplot.png`, `allafterm_remaining_length.tsv`
- Part E (individual traces): `allafterm_examples_gem.png`,
  `allafterm_examples_control.png`
- **Part F (added 2026-07-15, full/untruncated), see below**:
  `allafterm_full_read_mean_median_boxplots.png` + `allafterm_full_read_mean_median.tsv`,
  `allafterm_hmm_anchored_full_trace.png` + `_logx.png` + `.tsv`

Disposable smoke-test output: `Eylul_gem_data/smoketest_allafterm/` — not
wired into any other script's default paths, same convention as §8-9's
smoke-test directories.

### Part F: full, untruncated ALLAFTERM (added 2026-07-15)

Every earlier figure (Parts A-E) caps at 100/500/3000 events for readability
and comparability across reads. Part F removes every cap: **every event in
every deduplicated read's `ALLAFTERM` array is used, in full** — not
`NEXT100EVT`, not Source 1/2 raw samples, not any physical-end-truncated
subset. Confirmed directly from the run: `plot_full_read_boxplots` and
`plot_full_hmm_anchored_trace` both read from the `gem_allafterm`/
`control_allafterm` dicts (the same ones Parts A-D already use, populated
solely from `ALLAFTERM` lines via `load_blocks`), and both are called with no
`window`/`max_len` argument — the mean/median in
`allafterm_full_read_mean_median.tsv` are computed over `arr` directly (no
slicing at all), and `allafterm_hmm_anchored_full_trace`'s matrix width is
`max(len(a) for a in gem_arrays + control_arrays)` — i.e. exactly as long as
the single longest read in either sample (725,085 events, a control read).

**Counts**: 50 gem / 321 control unique reads in both new outputs (identical
to every other Part in this analysis — same dedup, same `main` dicts).

**`allafterm_full_read_mean_median_boxplots.png`**: two panels (per-read mean,
per-read median), one point per unique read regardless of tail length — a
read with 725,085 events and a read with 1,042 each contribute exactly one
number to each panel. `allafterm_full_read_mean_median.tsv` (371 rows:
read_id, sample, n_allafterm_events, mean, median) has the underlying values.

**`allafterm_hmm_anchored_full_trace.png`** (+ `_logx.png`): gem vs control
median±IQR from x=1 (first event after the final HMM Match) through
whichever read goes furthest in either sample, with a coverage (n
contributing unique reads) subpanel and vertical markers (dotted/dashed/
dashdot/finely-dotted, colored per sample) at each sample's own 50%/25%/10%/5%-
of-starting-count crossing — printed in the run log too (gem crosses 50% at
event 3,347, 25% at 5,973, 10% at 8,527, 5% at 9,086, out of 50 starting
reads/11,499-event longest read; control crosses 50% at 3,854, 25% at 5,519,
10% at 8,105, 5% at 10,433, out of 321 starting reads/725,085-event longest
read). The linear-x version is dominated by the ~14,500:1 span between the
shortest and longest contributing reads (gem's entire curve is compressed
into the leftmost sliver) — exactly the "unreadable" case anticipated, so the
log-x version is the one to actually read; both retain the complete x=1..725,085
range (log-x is a different x-axis *scale*, not a truncation — every position
is still plotted in both).

**One real rendering limitation hit and fixed**: `fig.savefig` initially
crashed with `OverflowError: In draw_markers: Exceeded cell block limit` —
matplotlib's Agg backend on this environment cannot draw a line/filled-polygon
with control's full ~725,085 vertices. This is a *display* limitation, not a
data problem — `allafterm_hmm_anchored_full_trace.tsv` (736,584 rows) already
has full, one-row-per-event resolution and was written successfully before
the crash. Fixed with `decimate_for_plot()`: stride-based decimation *for
rendering only*, capped at 20,000 points per line, always keeping the exact
first and last point so the plotted range's endpoints are exact — the TSV is
never decimated. This does not affect any other Part's plots (their
100/500/3000-point windows were always well under this limit).

**Companion TSV columns** (`allafterm_hmm_anchored_full_trace.tsv`): `sample`,
`event_offset_after_last_M` (1-based), `n_unique_reads`, `median`, `mean`,
`q1`, `q3` — one row per (sample, offset) pair where at least one read
contributes, i.e. up to that sample's own longest read (not padded with
placeholder rows beyond it).

### Scientific findings: is the post-Match signal actually different, gem vs control?

Pooled-event median±IQR bands (Parts A/B, all zoom levels) show substantial
overlap throughout, consistent with every earlier analysis in this project —
**but IQR-band overlap/non-overlap is not a significance test**, so a direct
read-level check was run instead: for each unique read, compute the summary
statistic (median/mean/std/slope-vs-event-number/fraction above +1/fraction
below −1) over its own first 25, 50, and 100 post-Match events, then compare
gem's 50 values against control's 321 values per metric with a two-sided
Mann-Whitney U test (nonparametric, doesn't assume normal per-read summary
statistics; one test per read, not per pooled event — the read is the unit).

At the 100-event window: read-level **mean** (p≈0.012), **slope vs. event
number** (p≈0.0015), and **fraction of events below −1** (p≈0.0046) all show a
nominally significant difference, with control trending more negative than
gem as the window widens (25→50→100 events); the **median** itself is only
borderline (p≈0.062, common-language/AUC≈0.58 — a randomly picked gem read
exceeds a randomly picked control read about 58% of the time, a modest, not
large, effect). **Robustness check**: dropping the 3 most extreme reads (by
|median|) from each group at the 100-event window *tightens* the median's
p-value (0.062 → 0.0087) rather than weakening it — the apparent difference is
not being carried by a handful of outlier reads; if anything, a few atypical
reads in each group were adding noise that partially masked it.

**Caveats that must accompany this, not follow it as an afterthought:**
- ~15 metric/window combinations were tested with no multiple-comparison
  correction — at nominal p<0.05 alone, a false positive or two among 15
  correlated tests is expected by chance. The *consistency of direction*
  across window sizes and correlated metrics (mean, slope, and
  fraction-below-−1 all pointing the same way) is reassuring but not a
  substitute for a pre-registered or corrected test.
- Sample sizes remain modest and asymmetric (50 gem vs. 321 control).
- This region is entirely **past** the reference-anchored, HMM-modelled part
  of the read (§4.4/§4.7's "does not prove biological origin" caveat applies
  in full) — a genuine, reproducible gem-vs-control difference *here* still
  cannot, by this analysis alone, be attributed to gemcitabine specifically
  rather than to some other systematic difference between the two library
  preps (the same open question already flagged for the qualifying-rate
  asymmetry, `TAIL_SIGNAL_CAPTURE.md` §4/§8).

**Conclusion**: there is a modest, read-level, outlier-robust signal — not
just a pooled-event or band-overlap artifact — but it is not large, not
correction-adjusted, and not yet distinguishable from a prep-batch effect.
Treat as a specific, well-defined hypothesis for follow-up (e.g. a
pre-registered test on independent/held-out reads, or matched-prep controls),
not as evidence of gemcitabine detection.
