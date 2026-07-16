# Gemcitabine detection with DNAscent: project overview and eventalign deep-dive

This document is a from-scratch narrative of the whole effort so far: what question
we're trying to answer, the first analysis that was run (and where it hits a hard
wall), a ground-up explanation of how DNAscent's `eventalign` pipeline actually
works internally (event detection, adaptive-banded alignment, the per-window
Viterbi HMM, all the way to the lines written to the `.align` file), and finally
the tail-signal-capture extension that was built to get past that wall. It is meant
to be readable on its own, without having read `TAIL_SIGNAL_CAPTURE.md` first —
that file is still the authoritative, more operational reference for build steps,
exact output file locations, and the plotting script; this document points to it
rather than repeating it wherever the two overlap.

---

## 1. The question, and the experimental setup

Gemcitabine (Gem, brand name Gemzar; chemical name 2′,2′-difluorodeoxycytidine,
dFdC) is a cytidine-analogue chemotherapy drug. Like BrdU/EdU (thymidine
analogues DNAscent was originally built to detect), it gets incorporated into DNA
in place of the natural base during replication — in Gem's case, in place of
cytidine. The scientific question driving this whole project is: **can nanopore
sequencing signal be used to detect gemcitabine incorporation**, the same way
DNAscent already detects BrdU/EdU?

To test this in a controlled way, a synthetic construct was made rather than
using genomic DNA directly:

- **Reference**: `Eylul_gem_data/ligation_product.fasta` — a 615 bp "ligation
  product" (500 bp + 100 bp pieces ligated together), sequence ending
  `...GACGTAGCAGCTGCGATAGTAGTCACTCG`.
- **Gemcitabine is expected at the second-to-last base**, i.e. 1-based position
  614 (the `C` immediately before the final `G`) — consistent with dFdC being a
  cytidine analogue.
- **Two samples were sequenced**: a **gem** sample (expected to carry the
  modified base at position 614) and a **control** sample (unmodified, unlabelled
  DNA), both nanopore-sequenced and basecalled/aligned the same way.

Both samples were run through DNAscent's `align` executable to turn raw nanopore
signal into a per-reference-position table of observed signal (see §3 for exactly
how that table is built). The question then becomes: **does the signal at or near
position 614 look different between the gem and control samples?**

## 2. First analysis: reference-anchored median/IQR comparison, and where it stops

Before any code was modified, the existing (stock) `DNAscent align` output was
analysed directly with a set of Jupyter notebooks in `Eylul_Code/` /
`Eylul_gem_data/`:

- `gem_forward_600_end_median_iqr.ipynb` (identical copies exist in both
  directories — confirmed byte-for-byte identical) — **this is the analysis the
  referenced figure comes from**:
  `Eylul_Code/gem_vs_control.forward.600_to_end.median_iqr.png`
  (regenerated copy also at
  `Eylul_gem_data/outputs/gem_vs_control.forward.600_to_end.median_iqr.png`).
- Two sibling notebooks trying alternative coordinate-labelling conventions:
  `gem_forward_600_end_directcoord_median_iqr.ipynb` (label by raw `align` file
  coordinate, no k-mer-centering) and
  `gem_forward_600_end_centered9_median_iqr.ipynb` (center every plotted genomic
  position on its own 9-mer, rather than only the last window). Both explore the
  same underlying data and hit the same wall described below; they aren't
  discussed further here.
- Input data: `Eylul_gem_data/{gem,control}/dnascent_align/{gem,control}.l600.test.align`
  — these are **standard, unmodified `DNAscent align` output** (produced with
  `-l 600`, i.e. minimum read length 600 bp, well before any of the tail-capture
  code existed).

### 2.1 What the notebook actually computes

DNAscent's normal `.align` output gives one line per raw signal sample that could
be assigned a reference coordinate, in the format (see §3.6 for the full
derivation):

```
<ref_coord 0-based>   <kmer on read strand>   <scaled observed signal>   <kmer on reference strand>   <expected model mean>   [BrdU_call  EdU_call]
```

The notebook:

1. Reads the 615 bp reference and, for every 1-based **9-mer-center** position
   from 600 up to the last position that still has a full, in-bounds centered
   9-mer, builds a mapping table (`build_center_position_mapping`). Because
   `align` coordinates behave as the 0-based center coordinate of a 9-mer window,
   position `pos` (1-based) maps to `align_coord_0based = pos - 1`, with the
   window spanning `pos-4 .. pos+4` (`CENTER_OFFSET = k//2 = 4`).
2. For a 615 bp reference, the **last position with a full centered 9-mer is
   611** (window `607-615`, the last 9 bases of the reference). Positions
   **612-615 have no valid centered 9-mer at all** — there's no reference beyond
   615 to complete the window. Position 614 (the actual Gem site) is one of
   these — **it is structurally impossible for this method to report a value
   at the Gem site itself.**
3. For each `.align` file, walks every **forward-strand only** read
   (`strand == 'fwd'`, since only the forward mapping direction is used in this
   project — see §4.1 for why), and for every line whose coordinate falls in the
   mapping table and whose observed k-mer matches the expected reference k-mer
   (mismatches are counted and asserted to be zero — a consistency check that the
   `.align` file and the reference agree), computes a **residual**:
   `residual = observed_scaled_signal - expected_model_mean`. Multiple raw-sample
   rows at the same position within one read are averaged first
   (`per_read_values`), so each read contributes at most one value per position.
4. Aggregates residuals **across reads** at each position into a **median** and
   **interquartile range (IQR: 25th-75th percentile)** — median/IQR rather than
   mean/std specifically because they're robust to the occasional outlier read
   (mean/std are also computed and saved to the TSV, just not the headline plot).
5. Plots gem (red) vs control (blue) median ± IQR across positions 600-611, with
   the windows that *contain* the Gem site (606-614 and 607-615, the two terminal
   windows) highlighted in yellow for visual reference even though the Gem base
   itself isn't a plotted x-position.

### 2.2 What the result shows

Looking at `gem_vs_control.forward.600_to_end.median_iqr.png`: **the gem and
control medians track each other closely across the whole 600-611 range, with
substantial IQR overlap at every position** — there is no visually obvious,
consistent separation between the two samples in this reference-anchored,
fully-analyzable region. The read counts backing this comparison are large
(from the underlying TSV,
`Eylul_gem_data/outputs/gem_vs_control.forward.600_to_end.median_iqr.tsv`):
**755 gem / 816 control reads contribute at position 600**, dropping steadily to
**294 gem / 628 control reads by position 611** as fewer reads' alignments
extend that far without being trimmed.

That drop-off is itself an important, easy-to-miss clue: **even under this much
looser criterion (a read just needs *some* aligned coverage at position 611, not
a full end-to-end alignment), control reads outnumber gem reads by roughly 2:1
at the last analyzable position (628 vs 294).** This is the same asymmetry that
reappears, much more starkly, in the tail-capture qualifying-read rates in §5
(control 7.3% vs gem 1.4%) — worth keeping in mind as likely the same underlying
phenomenon (probably a difference in how much of each prep's molecules carry
extra ligated material past the reference), not two unrelated observations.

### 2.3 Why this analysis can't see the Gem site — and what that motivated

The core limitation is structural, not statistical: **DNAscent's k=9 pore model
needs a full, centered 9-mer of reference sequence to make any call at all**, and
a 615 bp reference simply doesn't have one centered on position 614 (the window
would need to span 610-618, four bases past the end of the molecule). So no
matter how much data you collect, this method — and the standard, unmodified
`DNAscent align` output it's built on — **cannot ever report a value at the Gem
site itself.** The physical DNA molecule does keep moving through the pore past
that point, generating real raw signal, but the standard aligner discards it
because it has no reference k-mer to hang it on (see §3.6 for exactly where and
why this discarding happens in the code).

This is the direct motivation for the tail-signal-capture extension (§4-5):
**recover that discarded signal instead of throwing it away**, so it can be
inspected directly rather than inferred indirectly from neighbouring, unmodified
k-mers that never actually contain the Gem base in their window.

---

## 3. How DNAscent's `align` / `eventalign` pipeline works, from raw signal to output

This section explains the full pipeline as it exists in this repository
(`/scratch/esevinc22/Summer_Project_Boemo_Group`), independent of the
tail-capture changes (which are layered on top and described separately in §4).
Line numbers refer to `src/alignment.cpp` and `src/event_handling.cpp` unless
stated otherwise.

### 3.1 Top-level flow (`align_main`, `src/alignment.cpp:839`)

1. Parse CLI args (`-b` BAM, `-r` reference FASTA, `-i` DNAscent index mapping
   readID → pod5/fast5 file, `-o` output, plus `-t` threads / `-q` min mapping
   quality / `-l` min read length / `-m` max reads to process).
2. Load the index and the reference FASTA into memory (`import_reference_pfasta`).
3. Open the BAM file, read records one at a time, and **buffer** them
   (`maxBufferSize = threads`, or `4*threads` for `>4` threads) — records that
   don't pass `-q`/`-l` filtering, or that have zero query length, are discarded
   immediately without ever reaching the buffer.
4. Once a buffer fills (or the BAM is exhausted), process the whole buffer **in
   parallel** with `#pragma omp parallel for schedule(dynamic)`, one `DNAscent::read`
   object per BAM record, across `args.threads` OpenMP threads. Each read goes
   through: construction → signal fetch → event normalisation/rough alignment →
   `eventalign()` (the fine per-window HMM pass) → (if it passed QC) appended to
   the output file inside a `#pragma omp critical` block (so writes from
   different threads never interleave).
5. Reads fail (counted in `failed`) at exactly two points: if the rough
   event alignment QC rejects the read (`r.eventAlignment.size() == 0`, see
   §3.4's QC criteria), or if `eventalign()` itself sets `r.QCpassed = false`
   (only currently possible via the insertion-rate check in §3.6).

### 3.2 Constructing a `DNAscent::read` (`src/reads.h:211`, the `read` constructor)

Each BAM record becomes one `DNAscent::read` object, which is essentially the
per-read "workspace" that every later stage reads from and writes into. On
construction:

- Extracts the readID from the BAM query name, and (if this is a **split read**
  — Dorado can split one physical read into several BAM records, tagged with a
  `pi` parent-ID and `sp` start-coordinate aux tag) resolves which parent pod5/fast5
  file and signal offset to actually fetch raw signal from.
- **`parseCigar`** (`src/htsInterface.cpp:59`) walks the BAM CIGAR string and
  builds two coordinate maps used constantly downstream: `refToQuery` (reference
  0-based index → basecall/query 0-based index) and `queryToRef` (the inverse).
  Match/`=`/`X` operations advance both coordinates together; deletions/`N`
  advance only the reference (and are additionally flagged in `refToDel`);
  insertions/soft-clips advance only the query; hard clips advance neither.
  Crucially, **CIGAR is always stored in reference-coordinate order regardless of
  mapping strand** — for a reverse-mapped read, the walk over CIGAR operations
  happens back-to-front (`for (int i = n_cigar-1; i >= 0; i--)`) so that
  `refToQuery`/`queryToRef` still come out in a consistent, increasing-reference
  orientation. This same reference-order-regardless-of-strand fact is exactly
  what `hasSoftClipAtReferenceThreePrime` (§4.2) relies on.
- Looks up the reference contig this read mapped to and slices out
  `referenceSeqMappedTo` = the exact reference substring this read aligns
  against (`reference[refStart:refEnd]`).
- Fetches the basecalled sequence from the BAM record (`getQuerySequence`).
- **If the read mapped to the reverse strand**, both `basecall` and
  `referenceSeqMappedTo` are reverse-complemented in place, and `isReverse`/
  `strand="rev"` are set. This means that from this point on, **everything in the
  rest of the pipeline operates in "sequencing 5'→3'" orientation**, not
  reference-forward orientation — reference coordinates are only reconstructed
  right before being written to output (§3.6), by walking `reference_coord`
  forwards for `fwd` reads or backwards for `rev` reads.

### 3.3 Signal extraction and event detection (`normaliseEvents`, `src/event_handling.cpp:545`)

1. **Raw signal fetch**: `pod5_getSignal`/`fast5_getSignal` (dispatched on file
   extension) populate `r.raw` — the complete raw ADC/pA trace for this read,
   full length, nothing trimmed yet.
2. **Event detection**: `detect_events(...)` (from the vendored `scrappie`
   library, `src/scrappie/event_detection.*`) segments the raw trace into a
   table of candidate "events" — an *event* here means a maximal run of raw
   samples believed to represent one step of DNA through the pore (i.e. one
   translocation step, roughly but not exactly one base), detected via a
   changepoint/t-test-style statistic on the raw signal, independent of any
   sequence information. This is a purely signal-level segmentation — it knows
   nothing about the reference or the basecall yet.
3. The scrappie event table is converted into DNAscent's own `r.events` vector
   (each entry an `event{ double mean; std::vector<double> raw; }`) by merging
   consecutive scrappie sub-events that share the same reported mean (a
   scrappie implementation detail) and recording, for each final event, both its
   mean level and the actual slice of raw samples that produced it — this raw
   slice is exactly what's needed later to emit one output line per raw sample
   rather than one per event (§3.6), and is also what tail-capture (§4.4)
   collects when there are no more events but there are still raw samples.
4. **Pore-model k-mer ranks are precomputed** for both the basecalled query
   sequence and the reference sequence the read mapped to
   (`kmer_ranks_query`, `kmer_ranks_ref`), via `kmer2index` — a simple base-4
   positional encoding (A=0,T=1,G=2,C=3, most-significant base first) that maps
   any 9-mer string to an integer index into the flat `pore_model` lookup table
   (`Pore_Substrate_Config.pore_model`, loaded once at startup from
   `pore_models/r10.4.1_400bps.nucleotide.9mer.model` — one `(mean, std)` pair
   per possible 9-mer).

### 3.4 Scaling estimate, pass 1: quantile regression (`estimateScaling_quantiles`)

Raw nanopore signal isn't in the same numeric units as the pore model out of the
box — every read (and even every pore/channel) has its own affine
shift/scale, mostly from pore-to-pore and run-to-run current variation. Before
any alignment can compare observed signal to expected model levels, this needs
correcting.

The first-pass estimate doesn't require any alignment at all: it sorts the
observed event means and the reference k-mers' expected model means
**independently**, splits each into 10 quantile-median buckets
(`quantileMedians`), and fits a simple linear regression
(`linear_regression`) between the two sets of decile medians. This gives a rough
`(shift, scale)` such that `(observed - shift) / scale` should land roughly in
pore-model units — good enough to bootstrap the next, alignment-based stage, but
not accurate per-event since it never actually pairs up *which* event
corresponds to *which* k-mer.

### 3.5 Rough alignment: adaptive-banded event alignment (`adaptive_banded_simple_event_align`, `src/event_handling.cpp:149`)

This is the first of **two distinct alignment stages** in the whole pipeline —
easy to conflate, so it's worth being explicit about the division of labour:
this stage works in **query/basecall coordinates** (not reference coordinates)
and produces a rough, whole-read event↔k-mer correspondence; the second stage
(§3.6, `eventalign()`/`builtinViterbi`) re-does a much more careful alignment
**per ~50 bp reference window**, using this rough alignment only to know
roughly which events belong roughly where.

This function is adapted from nanopolish's classic adaptive-banded DTW-style
aligner:

- **Why banded, and why adaptive:** a full dynamic-programming alignment matrix
  between every event and every k-mer would be `O(n_events × n_kmers)` — too
  slow for whole reads. Instead, only a diagonal **band** of fixed `bandwidth`
  (100, from `AdaptiveBanded_Params_DNA_R10`, `config.h:41`) around the expected
  alignment path is computed at any time. The band's position is **adaptive**:
  at each step (`band_idx`), it decides whether the next band should shift "down"
  (consume another event) or "right" (consume another k-mer) by comparing the
  best scores at the lower-left vs. upper-right corners of the previous band
  (Suzuki's rule) — this lets the band track the true alignment path even though
  the overall event-rate-per-base isn't perfectly uniform.
- **Three move types per cell**, each a distinct HMM-like transition:
  - **Diagonal ("step"/D)**: consume one event **and** one k-mer — a normal
    1-event-1-base match. Scored as `previous-diagonal-cell + lp_step + emission`.
  - **Up ("stay"/U)**: consume one event but stay on the same k-mer — multiple
    events per base (translocation was slower than one event per base at this
    point). Scored as `previous-up-cell + lp_stay + emission`.
  - **Left ("skip"/L)**: consume one k-mer with **no** event — heavily
    penalized (`lp_skip = log(1e-30)`, i.e. essentially forbidden except when
    truly necessary) since it means a base produced no distinguishable event at
    all.
  - `p_stay` (and hence `lp_stay`) is derived once per read from
    `events_per_kmer = n_events / n_kmers`, so the stay/step balance adapts to
    how fast or slow this particular read translocated on average.
- **Emission probability** (`logProbabilityMatch`): the log of a Gaussian PDF —
  `(observed_event_mean - shift)/scale` compared against that k-mer's
  `(mean, std)` from the pore model. This is the only place actual signal-level
  information enters the DP; transitions are otherwise fixed per-read constants.
- **Initialization and termination are deliberately asymmetric**: band 0/1 are
  seeded assuming the alignment starts at k-mer 0 (`start_cell_offset`,
  `first_trim_offset`), but the **end** of the alignment is found by scanning
  every event against the *last* k-mer and picking whichever gives the best
  score **after adding a `lp_trim` penalty for every trailing event not
  included** (`+ (n_events - event_idx) * lp_trim`, `lp_trim = log(0.01)`). This
  is the single most important design detail for the tail-capture work: **the
  aligner is explicitly allowed to decide that it's better to stop early and
  "trim" the remaining events than to force-align them to the last k-mer** if
  those events don't look like that k-mer at all. Those trimmed events are never
  added to `r.eventAlignment` and are otherwise invisible to the rest of the
  pipeline — this is exactly the population that tail-capture Source 2 (§4.4)
  recovers.
- **Backtrace**: starting from the chosen best end cell, walk `trace[][]`
  backwards (`FROM_D`/`FROM_U`/`FROM_L`), pushing `(event_idx, kmer_idx)` pairs
  into `r.eventAlignment` at every D or U step (L/skip steps consume a k-mer but
  emit nothing), until reaching kmer_idx < 0. The result is reversed into
  chronological order. Along the way, consecutive same-k-mer events are also
  averaged together into `cleanedSignals`/`cleanedRanks` (used only for QC and
  for scaling refinement below, not for the final output).
- **QC gate** (`r.alignmentQCs.recordQCs`, checked right after): the read is
  **rejected outright** (`r.eventAlignment.clear()`, which `align_main` detects
  as `r.eventAlignment.size() == 0` and counts as `failed`) if:
  - `avg_log_emission < min_average_log_emission` (-6.0) — the alignment's
    average per-event fit to the model was too poor, or
  - `!spanned` — the alignment didn't reach from k-mer 0 to the last k-mer
    (`n_kmers - 1`) at all, or
  - `max_gap > max_gap_threshold` (5) — too many consecutive skip/L moves in a
    row at some point, or
  - fewer than 100 cleaned signal/rank pairs survived.

  This is a **whole-read, coordinate-free** QC check — importantly, it's
  unrelated to (and stricter in a different way than) the *reference-end*
  reachability that the tail-capture gate checks later in `align_main` (§4.3).

### 3.6 Scaling estimate, pass 2: Theil-Sen refinement (`estimateScaling_theilSen`)

Once a real event↔k-mer correspondence exists (`cleanedSignals`/`cleanedRanks`
from the rough alignment above), the shift/scale estimate is refined using a
**Theil-Sen estimator** — the median of all pairwise slopes between sampled
(signal, model-mean) points, which is far more outlier-robust than an ordinary
least-squares fit. This produces the `r.scalings` values (`shift`, `scale`) used
for every subsequent normalisation in the pipeline, plus `eventsPerBase`
(`et.n / (basecall.size() - k)`, the empirical average number of raw events per
reference base for this specific read — used to set the Match-state self-loop
probability in the fine HMM below). If this step fails to find a valid slope
(`slope_median == 0`), scalings are set to `-1` as a sentinel and the read is
failed (`r.eventAlignment.clear()`).

### 3.7 The fine per-window realignment: `eventalign()` (`src/alignment.cpp:563`)

**Why a second alignment stage exists at all:** the rough alignment above
(§3.5) works entirely in **query/basecall coordinates**, and the basecall can
disagree with the reference at any given position (sequencing/basecalling
errors, or genuine analogue-induced basecalling artifacts). DNAscent's whole
purpose, though, is to report per-**reference**-position statistics that are
comparable *across many different reads* — which requires re-expressing
everything in reference coordinates, and requires actually knowing, precisely,
whether each event should be labelled a Match, Insertion, or Deletion against
the *reference* k-mer at that position (not just "roughly aligned to it"). That
precise, reference-anchored labelling is what this function computes, one
manageable ~50 bp window at a time, using a full (small) Viterbi HMM rather than
the coarser banded heuristic.

**Outer loop — walking the reference in windows** (the `while` loop starting at
line 580):

- `windowLength_align = 50` bp (`config.h:46`) is the nominal window size, capped
  by however many bases are actually left (`basesToEnd`).
- **Breakpoint search** (lines 594-625): if there's substantially more than one
  window's worth of reference left (`basesToEnd > 1.5 * windowLength`), the code
  looks for a **natural boundary** to end the window at, rather than always
  cutting at a fixed length: it scans candidate positions in the next 1.5
  windows for a spot where the expected pore-model signal level changes sharply
  (`> 0.75` model-mean units) on **both** sides relative to the neighbouring
  k-mer. Ending a window right at such a natural "cliff" in the expected signal,
  rather than at an arbitrary length cutoff, is a heuristic to avoid crossing an
  ambiguous stretch of similar-looking k-mers exactly at a window boundary.
- **`referenceDefined`** (line 535): skips windows containing anything other
  than `A/T/G/C` (e.g. reference `N`s) — the pore model has no entry for
  ambiguous bases, so such a window is simply advanced past rather than aligned.
- **Gathering the events for this window** (lines 636-662): walks
  `r.eventAlignment` starting from `readHead` (carried over from the previous
  window, so this is a single forward sweep across the whole read, not
  re-scanned from the start each time), and collects every event whose rough
  `(event, kmer)` pair falls within this window's `refToQuery` range into
  `eventSnippet`/`eventSnippet_means`. Events are additionally sanity-guarded
  (`0 < mean < 250`, an implausible-signal-level filter used identically at
  every point raw signal is consumed in this file). If fewer than 2 events
  survive (typically because the rough alignment placed a deletion across this
  entire window), the window is skipped (`reference_index += windowLength;
  continue`) rather than aligned.
- **`indelScore`** (lines 664-668): `querySpan - referenceSpan`, a simple
  diagnostic of how much net insertion/deletion the rough alignment implies over
  this window — recorded per-event in the output (via `r.addSignal`'s `quality`
  argument) but not used to skip or reject the window itself.

**The fine alignment itself: `builtinViterbi`** (line 193, called at line 682)

This is a **from-scratch, from-scratch-per-window profile HMM Viterbi decoder**
(not reusing any external library), whose job is: given the (short, ~50 bp)
reference sequence for this window and the events roughly known to belong to
it, find the single most likely labelling of every event as
Match/Insertion/Deletion against a specific reference k-mer position, along with
a full backtrace so that labelling can actually be recovered (not just its
score).

- **State space**: for a window with `n_states = window_length - k + 1`
  reference k-mer positions, there are `3 * n_states` HMM states — one
  Match (M), Insertion (I), and Deletion (D) state **per reference k-mer
  position** in the window. This is the standard "profile HMM" topology used for
  nanopore signal alignment (the same structure nanopolish/DNAscent's ancestor
  tools use): M states emit an event whose level should match that k-mer's pore
  model; I states emit an event that doesn't correspond to advancing the
  reference at all (extra signal); D states are silent (advance the reference
  by one k-mer without consuming any event, i.e. "this k-mer produced no
  distinguishable signal", or more commonly here, just "the window boundary
  logic skipped past it").
- **Transition probabilities** — a mix of fixed constants and one per-read
  quantity:
  - Fixed ("external"/"internal", from `HMM_TransitionProbs_DNA_R10`,
    `config.h:42`): `externalD2D=0.3`, `externalD2M=0.7`, `externalI2M=0.999`,
    `externalM2D=0.0025`, `internalM2I=0.001`, `internalI2I=0.001`. ("external"
    = a transition that also advances to the next k-mer position; "internal" =
    a transition that stays on the same k-mer position.)
  - **Per-read**: `internalM12M1 = eln(1 - 1/eventsPerBase)` — the Match state's
    self-loop probability (stay on the same k-mer for another event), directly
    parameterized by this read's own average events-per-base from §3.6. A read
    that translocated slowly (high eventsPerBase) gets a higher self-loop
    probability, i.e. the model expects to see more consecutive events per base
    for that read.
  - All arithmetic is done in **log-space** throughout
    (`eln`/`lnSum`/`lnProd`/`lnGreaterThan`, `src/probability.cpp`) purely to
    avoid numerical underflow over long products of small probabilities — `eln`
    maps `0 → NaN` (log of zero, used as "impossible"/uninitialized) and is
    otherwise just `log(x)`.
- **Emission probabilities**: for a Match state at k-mer index `i`,
  `matchProb = eln(normalPDF(model_mean_i, model_std_i, (observed - shift)/scale))`
  — i.e. exactly the same kind of Gaussian log-likelihood used in the rough
  alignment (§3.5), just now scored against a single, exactly-known reference
  k-mer rather than a banded region. Insertion states have **no** emission
  penalty (`insProb = 0.0`, i.e. probability 1) — there's no specific k-mer an
  insertion "should" look like, so only the transition probability penalizes
  using an I state at all.
- **DP recursion** (the `for t in observations` / `for i in n_states` loops,
  lines 257-433): standard profile-HMM forward recursion computing
  `I_curr[i]`, `M_curr[i]`, `D_curr[i]` as the log-max over every valid
  predecessor transition (I from I-stay or M-to-I; M from I-to-M, M-to-M
  external, M-to-M internal-stay, or D-to-M; D from M-to-D or D-to-D — D states
  are filled in a second inner loop since they're silent and must be resolved
  after all M/I states at the same time step are known). **Backtrace arrays**
  (`backtraceS`: which state, `backtraceT`: which event/time index) are stored
  at *every* cell, not just kept as running maxima — this is what allows a full
  labelled path to be reconstructed afterward, not just a final score.
- **Termination and traceback**: the best of ending on D / M(+transition-to-end)
  / I(+transition-to-end) at the last k-mer position is chosen, then the
  backtrace is walked from there until reaching the sentinel `-1` (the start
  state), producing a list of `"<kmer_pos>_M"` / `"<kmer_pos>_I"` /
  `"<kmer_pos>_D"` labels which is then reversed into chronological
  (`stateLabels`, used back in `eventalign()`).

**Back in `eventalign()`, using the Viterbi output** (lines 684-803):

- A first pass over `stateLabels` finds `lastM_ev`/`lastM_ref`: the event index
  and k-mer position of the **last Match-labelled state** in this window's path.
  This can fall short of the nominal window end if Viterbi ends the window on a
  D or a run of I's — in which case the *next* iteration of the outer loop
  simply resumes from there (`reference_index += lastM_ref + 1`, `readHead +=
  lastM_ev + 1`), effectively re-processing the small leftover as part of a new
  window. `reachedFinalKmer` (line 712, part of the tail-capture work, §4.3)
  checks whether this Match position is at or past the very end of the whole
  read's reference span — i.e. whether the outer loop is about to terminate for
  good.
- A **second pass** over the same `stateLabels` actually emits output:
  - **Match (M)** states: compute the reference coordinate this event/k-mer
    corresponds to (`reference_coord + pos` for forward reads, `reference_coord
    - pos - 1` with the k-mer reverse-complemented for reverse reads — recall
    `referenceSeqMappedTo` was already put in sequencing-orientation back in
    §3.2, so this is where reference-genome orientation gets reconstructed for
    output), then write **one output line per raw sample** belonging to that
    event (`eventSnippet[evIdx].raw`, not just one line per event) — this is why
    a single 9-mer position typically has *many* lines in the `.align` file, all
    sharing the same reference coordinate and k-mer but each holding a
    different individual raw sample's scaled value. Each line is:
    `<ref_coord>\t<kmer_on_read_strand>\t<scaled_signal>\t<kmer_on_ref_strand>\t<model_mean>`,
    with two extra tab-separated fields (existing BrdU/EdU calls at that
    coordinate) appended if `r.refCoordToCalls` already has an entry there
    (relevant for `detect`'s reuse of this code path, not for plain `align`).
    `r.addSignal(...)` is also called here to populate `r.refCoordToAP`, the
    structure `detect.cpp`'s neural-network feature tensor is built from later —
    irrelevant for plain `align` output but why this bookkeeping exists at all.
  - **Insertion (I)** states: ordinarily written as a line with no real
    reference coordinate association beyond the position where the insertion
    occurred, kmer field replaced with `NNNNNNNNN` and model-mean field `0` —
    **except** insertions at or past `lastM_ev` (i.e. *trailing* insertions with
    no following Match in this window) are normally silently dropped, on the
    assumption that the *next* window will pick up whatever they represent. For
    the **truly final** window of a read, there is no next window — this is
    exactly the gap the tail-capture Source 1 addresses (§4.4).
  - **Deletion (D)** states emit nothing (silent, no event) but do increment
    `runningDeletions`, which combined with `runningInsertions` (incremented on
    every kept I) drives the whole-read `insRate` QC check at the very top of
    the outer loop — if net insertion rate ever exceeds 20%, the read is
    abandoned mid-way (`r.QCpassed = false; return;`), the only way `eventalign`
    itself fails a read.

### 3.8 Summary: the two alignments, side by side

| | Rough alignment (§3.5) | Fine alignment (§3.7) |
|---|---|---|
| Function | `adaptive_banded_simple_event_align` | `builtinViterbi` (called from `eventalign`) |
| Scope | whole read, once | one ~50 bp reference window at a time |
| Coordinate space | query/basecall | reference |
| Algorithm | banded DP, adaptive band placement | exact Viterbi, small enough to be exhaustive |
| States | 3 per k-mer (step/stay/skip moves) | 3 per k-mer (M/I/D, full profile HMM) |
| What it produces | `r.eventAlignment`: rough (event, query-kmer) pairs | the actual `.align` output lines, per reference coordinate |
| What it explicitly allows | **trimming** trailing events that don't fit | dropping trailing insertions *unless* this is the read's last window |
| Relevance to tail-capture | trimmed events = **Source 2** (§4.4) | dropped trailing insertions in the last window = **Source 1** (§4.4) |

---

## 4. The tail-signal-capture extension

Section 2.3 established the wall: the standard pipeline in §3 structurally
cannot report anything at reference position 614, because no full centered
9-mer exists there. But §3.5 and §3.7 both showed that real, raw signal *does*
exist past that point in many reads — it's just discarded rather than emitted.
This section documents, in detail, the code added to recover and print it
instead. (This is the same feature summarized more operationally in
`TAIL_SIGNAL_CAPTURE.md`; here the goal is to explain *why* each piece of code
does what it does, building on the full pipeline description in §3.)

### 4.1 What "past the end" means, concretely, and the strand restriction

The project only cares about the read's 3' end, since that's where Gem is
expected (position 614 of 615). For a **forward-mapped** read, "3' end of the
molecule" and "3' end / high-coordinate end of the reference" are the same
direction, so "signal past the last reference k-mer" directly means "signal
from the true 3' end of this molecule." For a **reverse-mapped** read, the
physical 3' end of the molecule corresponds to the reference's **5'** (low
coordinate) side — handling that case would require an entirely separate,
mirrored set of checks and was explicitly out of scope from the start (this
project only ever looks at forward-mapped reads for this reason). This is why
every gate described below starts with `r.strand == "fwd"`.

### 4.2 `hasSoftClipAtReferenceThreePrime` (`src/alignment.cpp:519`, new function)

As established in §3.2, BAM CIGAR strings are stored in reference-coordinate
order **regardless of mapping strand** — so the CIGAR operation corresponding to
the reference's 3' end is the **last** operation for a forward-mapped read, and
the **first** operation for a reverse-mapped read (`bam_is_rev(record) ?
cigar[0] : cigar[n_cigar-1]`). This function checks that one operation and
returns `true` if it's a soft-clip (`BAM_CSOFT_CLIP`).

**Why this matters**: the reference is a synthetic **ligation product** (500 bp
+ 100 bp pieces joined together) — real physical DNA molecules in this prep can
plausibly have **extra material ligated on past the designed 615 bp construct**.
If a read has a 3' soft-clip, the aligner is explicitly saying "there was more
sequence here that I couldn't align to the reference at all" — meaning whatever
raw signal follows the aligned region most likely corresponds to that unknown
extra ligated sequence, **not** to "a few more bases of our known construct
past position 615." Recovering tail signal from such a read would be
scientifically meaningless (or actively misleading) since there's no reference
sequence to say what it *should* look like. This function is precisely how the
gate below (§4.3) tells the two cases apart. This was confirmed in practice by
directly inspecting a real CIGAR string during development: one read had a
**107 bp 3' soft-clip** despite its aligned portion already spanning the *full*
615 bp reference — exactly the ligation-artifact scenario this check is meant
to exclude.

### 4.3 The gate (`align_main`, `src/alignment.cpp:1158-1167`)

```cpp
bool readPassesTailGate = (r.strand == "fwd")
                          and (r.refEnd == (int) reference.at(r.referenceMappedTo).size())
                          and not hasSoftClipAtReferenceThreePrime(r.record);
TailCaptureMode effectiveTailMode = readPassesTailGate ? args.tailSourceMode : TailCaptureMode::NONE;
eventalign(r, Pore_Substrate_Config.windowLength_align, effectiveTailMode);
```

All three conditions must hold simultaneously for `readPassesTailGate`:

1. **`r.strand == "fwd"`** — per §4.1, reverse-mapped reads are out of scope.
2. **`r.refEnd == reference length`** — the read's alignment must reach the
   *true* end of the 615 bp reference contig. If a read's alignment stops short
   (say, at position 550), there's no "signal past the reference end" to
   recover for it at all — it simply didn't get that far.
3. **`not hasSoftClipAtReferenceThreePrime(r.record)`** — the ligation-artifact
   guard from §4.2.

(This gate was originally a plain `bool captureTailSignal` passed straight to
`eventalign`; as of 2026-07-10 it's `readPassesTailGate`, a hard read-level
requirement combined with a **user-selectable** `--tail-source
{none,source1,source2,both}` CLI option — a gate-passing read still gets
`TailCaptureMode::NONE` if the user didn't ask for tail capture at all, and the
default is `none` for safety. `TailCaptureMode` is declared in
`src/alignment.h`; `eventalign`'s signature changed accordingly
(`eventalign(DNAscent::read &, unsigned int, TailCaptureMode)`), and the two
*other* call sites — `src/detect.cpp` and `src/trainCNN.cpp` (two call sites) —
simply pass `TailCaptureMode::NONE` literally: those tools have no use for tail
signal, and passing `NONE` guarantees bit-for-bit-identical behaviour to before
this feature existed.)

### 4.4 Where the recovered signal actually comes from: three independent pools

As of 2026-07-10 (`--tail-source` split) and 2026-07-15 (`NEXT100EVT`, §4.6),
what started as one combined buffer is now three separate pools, each gated by
`wantSource1 = (tailMode == SOURCE1 or tailMode == BOTH)` /
`wantSource2 = (tailMode == SOURCE2 or tailMode == BOTH)` — a read that failed
the gate in §4.3 gets `TailCaptureMode::NONE` and none of this runs at all.

**Source 1 — trailing insertions in the read's truly final HMM window**
(`src/alignment.cpp:919-931`, inside the existing per-window loop from §3.7,
pushed to `tailSignalSource1`). As described in §3.7, insertion-labelled events
past a window's last Match (`evIdx >= lastM_ev`) are normally always dropped,
because the *next* window is expected to pick up whatever they represent. For
the read's **last** window, there is no next window — `reachedFinalKmer` (line
848) is precisely the flag that detects this ("is the outer `while` loop about
to terminate for good after this window"). When both `reachedFinalKmer` and
`wantSource1` are true, these otherwise-doomed trailing-insertion events are
pushed into `tailSignalSource1` instead of being silently discarded. Note
carefully (documented in the source comment right above `reachedFinalKmer`):
reaching what *looks* like the window boundary isn't sufficient on its own to
guarantee no more windows will run — if Viterbi D-labels some trailing
reference positions, `lastM_ref` can fall short and the outer loop genuinely
does run once more on the leftover, so `reachedFinalKmer` has to check the
*actual* Match position reached, not just whether this window's nominal length
reached the end. **Empirically rare**: 0 of 50 gem reads, 6 of 347 control
reads (max 990 samples) have any Source 1 signal at all — see
`TAIL_SIGNAL_CAPTURE.md` §8.

**Source 2 — events the rough aligner trimmed before `eventalign()` ever ran**
(`src/alignment.cpp:941-961`, *after* the entire outer per-window loop
finishes, i.e. a whole-read-level step, not per-window, pushed to
`tailSignalSource2`). As detailed in §3.5, the adaptive-banded rough alignment
(`adaptive_banded_simple_event_align`) explicitly allows itself to "trim"
trailing events rather than force-match them to the last k-mer, if doing so
scores better — exactly the behaviour you'd expect if a modified base right at
the molecule's end produces signal that doesn't look like any unmodified
k-mer's expected level. Those trimmed events never enter `r.eventAlignment` at
all, so **nothing** in §3.7's window-processing loop ever sees them (that loop
only ever walks forward through `r.eventAlignment` via `readHead`). This block
walks `r.events` directly, starting right after `r.eventAlignment.back().first`
(the last event index the rough aligner *did* use) through to the true end of
`r.events`, applying the same `0 < mean < 250` signal-sanity guard used
everywhere else in this file, and appends every raw sample from each such event
to `tailSignalSource2`. **In practice this is the dominant source of tail
signal** — see the caveats in `TAIL_SIGNAL_CAPTURE.md` §8 for the actual
measured tail-length statistics (median 17,081 raw samples per qualifying gem
read), which are far larger than what Source 1 alone (bounded by one small HMM
window) could ever produce.

**NEXT100EVT — a third, mechanism-agnostic pool, see §4.6.**

### 4.5 Output format

```cpp
if (wantSource1){
    for ( size_t i = 0; i < tailSignalSource1.size(); i++ )
        r.humanReadable_eventalignOut += "TAIL_SOURCE1\t" + std::to_string(i) + "\t" + std::to_string(tailSignalSource1[i]) + "\n";
}
if (wantSource2){
    for ( size_t i = 0; i < tailSignalSource2.size(); i++ )
        r.humanReadable_eventalignOut += "TAIL_SOURCE2\t" + std::to_string(i) + "\t" + std::to_string(tailSignalSource2[i]) + "\n";
}
```

Each recovered raw sample becomes a line `TAIL_SOURCE1\t<idx>\t<scaled_value>`
or `TAIL_SOURCE2\t<idx>\t<scaled_value>`: three tab-separated fields (vs. 5 or 7
for a normal line), where the first field is never a numeric reference
coordinate (there isn't one — that's the entire point), `idx` restarts at 0
*per source, per read* (not a single running index across both), and
`<scaled_value>` is in the exact same shift/scale-normalized units as the
`<scaled_signal>` column of a normal line (i.e. directly comparable to the rest
of the file and across reads, but not raw pA current). Source 1 is always
written before Source 2 for a given read, matching the old combined buffer's
chronological construction. (Pre-2026-07-10 `.align` files use a single
unlabelled `TAIL` line instead — still parseable, treated as "legacy" by the
newer plotting script.)

See `TAIL_SIGNAL_CAPTURE.md` §8 for: the full `--tail-source` CLI option, the
file-location table, a concrete Python parsing example, the precise
operational definition of a "qualifying read", and the full list of caveats.

### 4.6 `NEXT100EVT`: the PI's proposed anchor, independent of the Source 1/2 split (added 2026-07-15)

§4.4 above is fundamentally an implementation-driven split: Source 1 and Source
2 exist as two separate pools because they're produced by two different parts
of the *existing* code (the terminal Viterbi window vs. the rough aligner's own
trim point), not because there's anything biologically meaningful about that
distinction. The PI proposed a simpler, implementation-agnostic alternative,
working directly from `r.events` (`src/reads.h:190` — the read's full
downsampled-signal vector, **always stored 5'→3'**, one entry per detected
event with a `.mean` and a `.raw[]`) and `r.eventAlignment`
(`src/reads.h:194` — the rough alignment's `vector<pair<event index, query
index>>`, used e.g. at `src/alignment.cpp:779` to look up `r.events[...].mean`
for a given aligned pair): find the index, in the read's full event vector,
of the **last event the HMM actually matched to a reference k-mer**, then take
the next 100 events after that index directly from `r.events` — regardless of
which of the two mechanisms above is responsible for whatever comes next.

**The problem this solves**: `eventSnippet`/`eventSnippet_means` (§3.7,
`src/alignment.cpp:758-759`) are rebuilt fresh for every ~50 bp window and
*only* ever hold event **means**, not their index into `r.events` — so once
`lastM_ev` (the *local*, within-this-window* index of the last Match state,
§3.7) is computed, there is no way to translate it back into a position in the
full per-read event vector. The fix (also proposed by the PI, "hacky, but the
data structures weren't set up for this either way") is a third array, built
in lockstep with `eventSnippet`/`eventSnippet_means`:

```cpp
// src/alignment.cpp, inside the per-window loop, right after eventSnippet is declared:
std::vector< size_t > eventIndeces;
// ...inside the same `if (0. < event_mean and event_mean < 250.)` block that
// pushes to eventSnippet_means/eventSnippet:
eventIndeces.push_back((r.eventAlignment)[j].first);
```

Since `eventIndeces` is pushed to at exactly the same point as
`eventSnippet_means`/`eventSnippet`, it's guaranteed to stay the same length and
in step with them — `eventIndeces[i]` is "where `eventSnippet_means[i]` lives
in `r.events`." Once the existing lastM_ev-finding loop (§3.7,
`src/alignment.cpp:817-831`) has run for this window, `eventIndeces[lastM_ev]`
is exactly the number the PI asked for: the raw `r.events` index of the last
event genuinely Match-labelled against the reference, in **this** window. That
value is captured into `lastMatchRawIdx`, a variable declared once *before*
the outer per-window `while` loop (`src/alignment.cpp:698`) and unconditionally
overwritten every iteration — so once the loop finishes, whatever the **last**
window computed is what's left in `lastMatchRawIdx`. This is a direct
implementation of "the number you want is `eventIndeces[lastM_ev]` for the last
time you do this in the read": no extra bookkeeping about the outer loop's
termination condition is needed, since a variable declared outside a loop and
overwritten inside it naturally ends up holding the last iteration's value.

**Why this is *not* simply "Source 2 again"**: it's tempting to assume this
anchor is the same as Source 2's (`r.eventAlignment.back().first`, §4.4) —
after all, Viterbi can never see an event the rough aligner didn't already hand
it, so both anchors are upper-bounded by the same rough-alignment boundary.
But they can genuinely differ: if the terminal window's HMM labelled some
trailing events as **Insertion** (i.e. Source 1 fired for this read) rather
than extending the Match, those events sit strictly *between*
`lastMatchRawIdx` and `r.eventAlignment.back().first` — Viterbi consumed them
(they're in `eventSnippet`, hence in `eventIndeces`), just didn't label them
`M`. So for the 0 gem / 6 control reads that have real Source 1 signal
(§4.4), `lastMatchRawIdx` is *earlier* than Source 2's boundary by exactly the
length of that read's Source 1 pool — meaning a "next 100 events from
`lastMatchRawIdx`" window folds any Source 1 signal into its first entries,
rather than skipping straight past it into Source 2 territory the way an
anchor at `r.eventAlignment.back().first` would.

**Collection and output** (`src/alignment.cpp:963-984`, after the outer
`while` loop, alongside Source 1/2):

```cpp
std::vector<double> tailSignalNext100Events;
if ( (wantSource1 or wantSource2) and haveLastMatchRawIdx ){
    size_t stop = std::min( r.events.size(), lastMatchRawIdx + 1 + 100 );
    for ( size_t evi = lastMatchRawIdx + 1; evi < stop; evi++ ){
        double event_mean = r.events[evi].mean;
        if ( not (0. < event_mean and event_mean < 250.) ) continue;
        tailSignalNext100Events.push_back( (event_mean - r.scalings.shift) / r.scalings.scale );
    }
}
```

Reused the existing `--tail-source != none` gate rather than adding new CLI
surface — one value is pushed **per event** (`.mean`), not per raw ADC sample
the way Source 1/2 do (they expand each event's `.raw[]`), since the PI
specifically asked for "the next 100 *events*." Printed as
`NEXT100EVT\t<idx>\t<val>` lines, after Source 1/Source 2 for that read.
Deliberately **not** prefixed `TAIL` — `plot_tail_signal.py` and
`plot_tail_signal_end_focused.py` both key off a bare `line.startswith("TAIL")`
to mean "legacy combined signal," and would otherwise silently fold this pool
into that bucket, double-counting samples already present in Source 1/2.

**Verified on the full canonical dataset** (`-l 600 --tail-source both`,
`slurm_full_run/{gem,control}/*.full.align`): every one of the 50 gem and 347
control qualifying reads got exactly 100 `NEXT100EVT` samples — no read ran
out of events before the cap.

**Important note on how this was actually verified** (correcting an earlier,
weaker check): printed line *position* in the `.align` file is always
`TAIL_SOURCE1` → `TAIL_SOURCE2` → `NEXT100EVT` regardless of which raw events
`NEXT100EVT` draws from (Source 1/2 print their full buffers unconditionally
when requested; `NEXT100EVT` is a separately-collected pool printed after
them) — so "does `NEXT100EVT` appear right after `TAIL_SOURCE1`/`TAIL_SOURCE2`
in the file" does **not** actually prove *which* raw events it contains. The
real proof is a direct comparison of `lastMatchRawIdx` (this anchor) against
`r.eventAlignment.back().first` (Source 2's own boundary), which was done with
a temporary debug print (see `TAIL_SIGNAL_CAPTURE.md` §9's table): for control
read `e1cd3406-9f89-4862-9055-ca6e242e4e69`'s block with 0 Source 1 signal,
both anchors were identical (50434); for its other block, with 990 Source 1
*samples*, `lastMatchRawIdx` (51082) was 193 raw *events* earlier than Source
2's boundary (51275) — confirming the anchor genuinely differs, by exactly the
span of that block's Source 1 signal, exactly as this section predicts. Full
details, exact sample counts, and the smoke-test history (including an earlier,
since-corrected version of this feature that anchored at Source 2's boundary
instead) are in `TAIL_SIGNAL_CAPTURE.md` §9.

### 4.7 `ALLAFTERM`: the uncapped pool, and a full read-level analysis (added 2026-07-15)

`NEXT100EVT` (§4.6) is fixed at 100 events — useful for a uniform,
comparable-across-reads view, but every qualifying read in this dataset had
well over 100 events available past `lastMatchRawIdx`, so 100 events discards
most of the recoverable signal. `ALLAFTERM` is the same pool, uncapped:
**every valid event in `r.events` from `lastMatchRawIdx + 1` to
`r.events.size()`**, i.e. `NEXT100EVT`'s full superset, using the identical
anchor, guard, and scaling. `NEXT100EVT` itself did **not** change.

```cpp
// src/alignment.cpp:996-1016, immediately after NEXT100EVT's own collection,
// reusing the same lastMatchRawIdx/haveLastMatchRawIdx computed once per read:
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

Printed as `ALLAFTERM\t<idx>\t<val>` lines (`src/alignment.cpp:1043-1051`,
immediately after `NEXT100EVT` for the same read) — deliberately not
`TAIL`-prefixed, same reasoning as `NEXT100EVT` (§4.6). Terminal summary
(`printTailCaptureSummary`) gained block count, a genuinely **unique-read-ID**
count (a `std::set<std::string>`, accumulated inside the same
`#pragma omp critical` section as every other counter — thread safety
unchanged), total event count, and min/median/mean/max events per alignment
block. Gated on the identical `(wantSource1 or wantSource2)` condition as
everything else in this feature — no new CLI surface, still off by default.

**Output semantics — four distinctions that matter when interpreting the new
plots (§4.7's plotting script, below):**
1. **`NEXT100EVT` is exactly `ALLAFTERM`'s first up-to-100 valid entries** —
   proven per-block (not just asserted) by `validate_allafterm.py`, PASS on
   both the smoke test and the full canonical run, zero mismatches.
2. **Granularity differs by pool.** Source 1/Source 2 print *expanded raw-ADC
   samples* (several lines per event). `NEXT100EVT`/`ALLAFTERM` print one line
   per *event mean*, never expanded. A length in one space is not a length in
   the other.
3. **`ALLAFTERM` can span both Source 1 and Source 2 territory** — it has no
   concept of which mechanism produced a given event, unlike the §4.4 split.
4. **Neither pool proves the recovered signal is the modelled construct's
   final bases specifically** — it may be adapter, motor-protein, ordinary
   dwell variability, or (for reads already past the 3′ soft-clip gate)
   genuine non-reference sequence. Match-labelled events remain excluded from
   both pools and untouched in ordinary output.

**Validation** (full detail, exact numbers, and the one genuine surprise
encountered — pre-existing OpenMP block-order nondeterminism unrelated to this
change — are in `TAIL_SIGNAL_CAPTURE.md` §10): all 7 requested checks passed,
including a byte-identical comparison of every pre-existing line against the
pre-`ALLAFTERM` canonical files (order-independent, per-read-ID multiset
hash comparison, since raw `diff` doesn't handle across-run block reordering).

**New plotting script**: `Eylul_gem_data/scripts/plot_allafterm_analysis.py`
(does not overwrite the two existing plotting scripts). Parses `NEXT100EVT`
and `ALLAFTERM` independently; deduplicates by read_id keeping the block with
the **longer `ALLAFTERM` array** (documented as deliberately simpler than
§4.4's per-channel Source 1/2 rule, since `NEXT100EVT` is always that same
block's own prefix — no complementary-channel scenario to reconcile). Part E's
example-read selection specifically works around one consequence of that rule:
a read's ALLAFTERM-length-winning block is not guaranteed to be the block that
happens to carry Source 1 signal, so example traces for Source-1-bearing reads
use that read's Source-1-bearing block specifically, not blindly the winner.

Six groups of output, all in
`Eylul_gem_data/slurm_full_run_allafterm/figures_allafterm/`:
- **A — HMM-Match-anchored median±IQR**: first 100 events (from `NEXT100EVT`)
  and first 500 (from `ALLAFTERM`), each with a lower panel of unique-read
  coverage vs. x and a companion TSV (x, sample, n, median, Q1, Q3).
- **B — physical-end-anchored median±IQR**: last 100/500/3000 detected events
  (reversed `ALLAFTERM`, x=0 = final detected event), same panel/TSV structure.
- **C — read-level boxplots**: median/mean/std/slope-vs-event-number/fraction
  above +1/fraction below −1, over the first 25/50/100 post-Match events, plus
  a focused single boxplot of the first-100 median (the `ALLAFTERM` analogue
  of `tail_end_3000_median_boxplot.png`).
- **D — remaining-length distribution**: histogram, log-scaled histogram,
  ECDF, box/violin of `ALLAFTERM` event count per unique read.
- **E — individual traces**: first 500 `ALLAFTERM` events for ≥10 gem + ≥10
  control reads (fixed seed 0), every Source-1-bearing read force-included.
  x=0 marked as the final HMM Match (true by construction). The Source
  1/Source 2 boundary is **deliberately not marked** — Source 1/2 are
  raw-sample-granularity and `ALLAFTERM` is event-mean-granularity, and the
  per-event raw-sample count needed to convert between them is not
  recoverable from this output at all, so estimating that boundary on the
  plot would be fabricating precision the data doesn't support. States this
  explicitly rather than guessing.
- **F — full, untruncated `ALLAFTERM`** (added 2026-07-15, follow-up request):
  no 100/500/3000-event cap anywhere. Per-read mean and median over the
  **entire** `ALLAFTERM` array (`allafterm_full_read_mean_median_boxplots.png`,
  one value per unique read regardless of tail length — a 725,085-event read
  and a 1,042-event read each contribute exactly one point), and a
  median±IQR trace from x=1 through whichever read goes furthest in either
  sample (`allafterm_hmm_anchored_full_trace.png` + a log-x version, since
  the ~700x span between shortest and longest reads makes the linear version
  unreadable in the informative sub-10,000-event region — both retain the
  complete range, log-x is a different axis scale, not a truncation), with a
  coverage subpanel and explicit 50/25/10/5%-of-starting-count crossing
  markers so a viewer can't mistake a thinly-covered tail for a well-covered
  one. Hit and worked around a genuine `matplotlib`/Agg rendering limit
  (`OverflowError: In draw_markers: Exceeded cell block limit` on a
  ~725,085-vertex polygon) with display-only stride decimation — the
  companion TSV (`allafterm_hmm_anchored_full_trace.tsv`, 736,584 rows) is
  never decimated and keeps full per-event resolution. Full detail in
  `TAIL_SIGNAL_CAPTURE.md` §10's Part F addendum.

Counts (full canonical run): 50 gem / 321 control unique reads (347 blocks,
26 collisions — cross-checked: C++'s independently-computed unique-read-ID
count and the plotting script's independent dedup both landed on 321 for
control). `ALLAFTERM` totals 212,484 events (gem) / 2,476,953 events
(control); events per alignment block range 1,134–11,499 (gem) /
1,042–725,085 (control) — see `TAIL_SIGNAL_CAPTURE.md` §10 for the full table.

**Scientific limitations, restated for this specific analysis**: sample sizes
remain modest and asymmetric (50 gem vs. 321 control); event-index positions
are not base positions; sequencing/prep differences between the gem and
control BAMs (already flagged in §4.4/`TAIL_SIGNAL_CAPTURE.md` as an
unexplained qualifying-rate asymmetry) could in principle produce systematic
signal differences unrelated to gemcitabine chemistry, and this analysis
cannot distinguish that from a genuine gemcitabine-specific effect. See the
closing discussion below (and the conversation this document accompanies) for
whether the observed post-Match differences are visually apparent and whether
they hold up under an actual read-level statistical check, rather than being
read off IQR-band overlap.

### 4.8 Time-warped post-Match analysis: which window to trust, dynamic time warping, and a ground-truth before/after comparison (added 2026-07-16)

Everything below works from Source 1/2 unchanged and, following the PI's
instruction, deliberately sets Source 1/2 aside — it is entirely about the
`ALLAFTERM` ("events after the last M") pool from §4.7. No C++ code changed
in this section; everything here is new Python analysis of the already-
generated `.align` output.

**4.8.1 Verifying the implementation against the PI's exact description.**
The PI separately re-described the mechanism in more precise terms, citing
specific upstream lines (commit `e830fe7b`, `reads.h:190,194`,
`alignment.cpp:629,634`). Fetched and compared directly:
- The reference↔read map ("genome alignment") is `r.refToQuery`/`queryToRef`
  (`src/reads.h:192`, populated by `parseCigar` from the BAM record) and the
  read↔events map ("rough event alignment") is `r.eventAlignment`
  (`src/event_handling.cpp:149-422`) — confirmed to match the PI's
  description of "a map between reference and read given by genome
  alignment" and "a map between read and events given by the rough event
  alignment" exactly.
- The scaling formula the PI pointed to (upstream `alignment.cpp#L718`,
  `scaledEvent = (raw - r.scalings.shift) / r.scalings.scale`) is identical,
  verbatim, to what this codebase already applies to every rescued pool
  (Source 1, Source 2, `NEXT100EVT`, `ALLAFTERM`) — same per-read
  `PoreParameters` (`src/reads.h:188`), same formula, no divergence.
- One deliberate divergence from a specific PI *suggestion* (not a
  requirement): the PI floated "don't label these as insertions, label them
  as something else, add it to the human-readable string." What's actually
  implemented instead routes trailing insertions to a **separate buffer**
  (`tailSignalSource1`, still internally Viterbi-labelled `"I"` — that part
  wasn't changed) rather than inventing a new in-line state label, and prints
  all four rescued pools as their own tagged blocks (`TAIL_SOURCE1`,
  `TAIL_SOURCE2`, `NEXT100EVT`, `ALLAFTERM`) after the normal per-position
  output, rather than embedding them into the reference-coordinate-anchored
  string (which has no valid coordinate to anchor to past the last Match in
  the first place). Functionally equivalent for the stated goal ("make it
  clear these come after the last aligned event"); noted since the PI said
  the exact implementation doesn't matter as long as it's correct.

**4.8.2 How many events are actually "rescued," and which window is trustworthy.**
`Eylul_gem_data/scripts/analyze_rescue_window.py` (SLURM:
`slurm_full_run_allafterm/run_rescue_window_analysis.slurm`; output in
`figures_window_analysis/`) answers three practical questions that came up
when deciding what to plot next:
- **Total/per-read `ALLAFTERM` length** (post-dedup, longest-`ALLAFTERM`-wins,
  matching §4.7's rule): gem 50 reads / 212,484 events total (min 1,134,
  median 3,328, mean 4,250, max 11,499); control 321 reads / 2,233,201 events
  total (min 1,042, median 3,853, mean 6,957, max **725,085**). Note this
  2,233,201 is *post-dedup* — §4.7's 2,476,953 figure is the pre-dedup,
  all-347-blocks total; the two aren't the same denominator.
- **One read dominates control's tail**: the single 725,085-event control
  read alone accounts for 32.5% of *all* of control's pooled `ALLAFTERM`
  events combined — 8.7× longer than control's next-longest read (83,223).
  Gem has nothing comparable (max 11,499). This is why control's *mean*
  tail length (6,957) looks so different from its *median* (3,853) while
  gem's mean/median are close (4,250 / 3,328) — the difference in "events
  rescued" between conditions is mostly this one outlier plus control simply
  having 6.4× more qualifying reads, not a shifted population of ordinary
  reads.
- **Coverage decay** (event-offset at which each sample's unique-read count
  first drops below X% of its own starting N): both samples hold effectively
  100% coverage only out to ~1,000-1,100 events (gem: 100% to offset 1,135,
  50% by offset 3,347; control: 100% to offset 1,043, 50% by offset 3,854).
  Past a few thousand events the pooled median at large offsets is
  increasingly just whichever 1-2 reads happen to still be that long — by
  offset 80,000+ in control that's literally the one 725,085-event read.
- **Read-level Mann-Whitney U by window** (median/mean per read, gem vs.
  control, no pooled-event fallacy): window 100 → p≈0.06 (median) / 0.012
  (mean); window **500 → p≈0.0037 (median) / 0.0036 (mean), the strongest
  result, AUC≈0.63**; window 1,000 onward → null (p 0.25-0.91); full
  untruncated comparison → null (p≈0.29) whether or not the two largest
  control outlier reads are excluded (excluding them doesn't materially
  change the null result, so the far tail's non-signal isn't an outlier
  artifact — it's genuinely not there).
- **Geometric argument for which window is biologically sensible**: the
  center of the last aligned 9-mer (k=9) is reference position 611, so that
  9-mer spans bases 607-615 — meaning gemcitabine, at base 614, is **already
  inside the sensing footprint of the very last successfully-aligned event**,
  not somewhere downstream of it. It should fully clear the pore's 9-mer
  window only once the window's leading edge reaches ~base 618 (614+4), i.e.
  ~7 bases past the anchor — at the PI's own ~5-events/base estimate, ~35
  events, generously doubled or tripled for lesion-induced dwell-time slowing
  still well under 100-150 events. **Conclusion: the first ~100-500 events
  make much more sense than 1,000+ as the biologically relevant window** —
  the ~1,000-event full-coverage boundary is about statistical power (every
  read still contributes), not about where gemcitabine's signal should
  actually be.

**4.8.3 Dynamic time warping (DTW) of the first-N post-Match events.**
Motivation (the PI's "elephant in the room"): post-Match events have no
reference coordinate, and reads translocate at different speeds, so "event
offset 7 after the last Match" in two different reads need not be the same
physical base — naive fixed-offset averaging (§4.7 Part A/F) can blur or
manufacture signal. DTW finds a monotonic correspondence between two
signals' *shapes* rather than assuming a fixed events/base rate. DTW is
inherently **pairwise** — there is no native N-way DTW — so aligning many
reads onto one shared, comparable axis requires designating something fixed
as the common target every read is aligned to.

*First attempt and the bug it had.* Picked one real read (from the pooled
gem+control set, closest by L1 distance to the pooled per-offset median) as
that target, and ran plain unconstrained DTW (no limit on how far the
warping path may stray from the diagonal). Caught via a control test before
trusting the result: DTW-aligning **temporally shuffled** (pure noise, same
values, scrambled order) copies of real reads to the same reference achieved
**equal-or-higher** post-DTW correlation-to-reference than the genuinely-
ordered reads (e.g. one read: 0.55 real vs. 0.71 shuffled). Unconstrained DTW
was warping *anything* to resemble the reference; the resulting near-perfect
gem/control overlap was an artifact, not a finding.

*Fix* (`Eylul_gem_data/scripts/dtw_first100_align.py`, generalized past its
filename to take `--window` — validated at both 100 and 500):
- **Sakoe-Chiba band** (`--band-radius`, default 20 events): the warping path
  cannot map a read-event more than this many positions from the diagonal.
- **Step penalty** (`--step-penalty`): a fixed cost on non-diagonal (many-to-
  one) moves, exactly like a gap penalty in sequence alignment. Banding
  *alone* was not enough — swept radii 3-50 and even the best (20) only
  barely passed the validation below (margin +0.15, right at threshold),
  because a band still leaves the path free to hunt locally for whatever
  nearby value happens to match noise by chance.
- **DTW-barycenter reference** instead of one arbitrary read: start from the
  pooled per-offset median, then iterate {banded-DTW-align every read to the
  current reference → resample → median-of-resampled = next reference}
  (`--barycenter-iters`, default 3) — a consensus is far harder to chase with
  spurious warping than one specific read's own noise.
- **The same real-vs-shuffled-noise control is now built into the script
  itself** (`shuffle_control_validation`), run automatically every time and
  printed as an explicit PASS/FAIL, annotated directly onto the figure's
  title (not just the log) so the verdict travels with the plot.
- **Validated parameters**: `band_radius=20, step_penalty=0.3,
  barycenter_iters=3` → real reads reach mean post-DTW correlation-to-
  reference 0.79 (window 100) / 0.77 (window 500) vs. 0.58 for shuffled
  noise in both — margins +0.21 / +0.19, comfortably passing.
- Cost is `O(window × band_radius)` per read pair; the script refuses
  `--window` above 300 without an explicit `--force-large-window`
  acknowledging the cost tradeoff (this is why the pool is only ever used at
  small windows — the full uncapped `ALLAFTERM` range, up to 725,085 events
  for one control read, is computationally infeasible for DTW regardless of
  banding).

*Visual result*: once genuinely pathological warping is ruled out, gem and
control's DTW-aligned median curves largely **track each other** through the
whole window (sharp shared peaks recur at the same warped positions in both
samples, extending across the full 500-event window tested, not just the
first 100). This is expected, not a null result — both samples read the same
reference locus, so the ordinary (non-lesion) sequence context downstream of
the anchor should produce near-identical expected signal; DTW recovers that
shared, sequence-driven shape that fixed-offset averaging had smeared into
noise. The original (buggy, unconstrained) run is kept at
`figures_dtw_first100/` purely as a record of the bug — **do not use it**;
the validated outputs are in `figures_dtw_first100_banded/` and
`figures_dtw_first500_banded/` (`dtw_vs_naive_first{100,500}_{mean,median}.png`
— both a `mean ± 1 SD` and a `median ± IQR` variant were built, since the two
turned out visually near-identical, i.e. this signal is not a few-outlier-
reads artifact; `dtw_warped_matrix.tsv`, `dtw_barycenter_reference.tsv`,
`dtw_cost_per_read.tsv`, `dtw_vs_naive_band_summary.tsv`, a cached extraction
TSV to avoid re-parsing the full `.align` files on every re-run).

**4.8.4 Read-level statistics on the DTW-aligned data**
(`Eylul_gem_data/scripts/dtw_read_level_stats.py`, output
`figures_dtw_first100_banded/dtw_read_level_stats.tsv`). Mann-Whitney U,
gem vs. control, on **per-read** median/mean over pre-registered windows
(10, 20, 35, 50, 100 events — chosen from the §4.8.2 geometric estimate),
comparing the DTW-aligned values against the naive (unwarped) values at the
identical offsets, 20 tests total, Bonferroni threshold 0.0025:
- **Only one result survives correction: DTW-aligned per-read median over
  the first 20 events, p=0.00053, AUC=0.65.** At that same window the naive
  (unwarped) version shows nothing (p=0.52) — a real instance of DTW
  correction recovering a signal that fixed-offset averaging missed.
- At window 0-35 (the geometric "full lesion clearance" estimate) the
  pattern **inverts**: DTW shows nothing (p≈0.9-0.99) while naive shows a
  weak, uncorrected-only signal (p≈0.03, does not survive Bonferroni) that
  disappears under DTW. The two methods disagree at both windows, in
  opposite directions — a sign that any real effect is narrow and fragile
  (roughly events 9-20 after the anchor), not a broad signature across the
  whole predicted footprint.
- Window 0-100 naive-mean result (p=0.012) exactly reproduces §4.8.2's
  independently-computed window-100 result — a useful cross-check that the
  two scripts' per-read aggregation logic agree.
- A separate, explicitly-flagged **exploratory/circular** check on the exact
  offsets {9, 13, 14, 18} (chosen *because* they showed the largest
  gem-control gap in this same dataset) reaches even lower p-values
  (down to 0.001) — reported for the record only; selection bias makes this
  non-confirmatory, and the pre-registered 0-20 window above isn't fully
  independent of that same earlier observation either, so treat 0.00053 as
  a strong lead rather than a clean confirmation.

**4.8.5 A ground-truth before/after comparison**
(`Eylul_gem_data/scripts/before_after_match_plot.py`; SLURM:
`run_before_after_{100,500}_kutem_gpu.slurm`; output in
`figures_before_after_{100,500}/`). Motivated by a natural question: does
gem vs. control track as closely on the *known-aligned* side of the anchor
as DTW claims to recover on the unaligned side? Built without any new C++
output or a persistent raw-index tracker — **the needed data was already
being printed**: the ordinary per-position Match lines emitted throughout
every window of the read (`event_coord\tkmerRef\tscaledEvent\tkmerStrand\t
meanStd.first`, one line per **raw ADC sample**, several consecutive lines
sharing one `event_coord`). Grouping consecutive same-`event_coord` lines
and averaging their `scaledEvent` recovers the true per-event mean at
identical scaling to `ALLAFTERM` (scaling is affine, so it commutes with
averaging); insertion lines are already distinguishable (an all-`N`
`kmerStrand`) and excluded, since they don't correspond to a confirmed
reference base. Safe without tracking any extra index because
`event_coord` only ever moves in one direction through a whole read
(alignment.cpp's outer loop never revisits an earlier `reference_index`) —
consecutive lines sharing a coordinate are guaranteed to be one physical
event, never two coincidentally-revisited ones. Unlike the post-Match side,
**this pre-anchor side needs no DTW at all** — it's genuinely reference-
locked across reads already, by construction of the HMM alignment.

Constraint discovered while building this: the reference amplicon is only
~600-610bp total (the last-Match anchor sits at ~base 611), so requesting
500 pre-anchor events uses nearly the *entire* aligned region — only 28/50
gem reads and 235/321 control reads have that many distinct aligned
positions (window 100 has full coverage: 50/50 gem, 321/321 control).

Result, both window sizes: gem and control track each other **tightly**
before the anchor (real, position-locked ground truth — a positive control
confirming the shared-signal phenomenon isn't a scaling artifact). That
tight co-variation **collapses immediately** in the naive (unwarped)
after-panel, but is **substantially recovered** in the banded-DTW after-
panel — visually a far smoother transition across the anchor than the naive
version's sharp break. The one place gem and control still visibly separate
in the DTW panel is the same narrow bump around offset +10 to +20 already
flagged statistically in §4.8.4 — consistent, but (per that section's own
caveat) not to be read as a significance test by eye.

**Scientific bottom line for §4.8**: there is a real, narrow, DTW-dependent
signal around events 9-20 after the last Match (read-level MWU p=0.00053,
survives Bonferroni across the pre-registered window sweep) that naive
fixed-offset analysis misses, geometrically consistent with where
gemcitabine (base 614) should actually manifest, and visually consistent
across the coverage-decay analysis, the DTW band plots, and the
before/after ground-truth comparison. It is not confirmed at the wider
window the geometry also seemed to predict (0-35), the window-20 boundary
isn't fully independent of an earlier data-driven observation, and — as in
§4.7 — a real, reproducible difference here still needs to be distinguished
from a systematic gem-vs-control library-prep effect before attributing it
to gemcitabine specifically. Treat as the most promising lead so far, worth
independent validation (e.g. a held-out split, or a different barycenter
reference) before presenting it as a finding.

---

## 5. Current state and where to look for what

| Question | Where to look |
|---|---|
| Does gem vs control differ anywhere in the normal, reference-anchored, fully-analyzable region (positions 600-611)? | §2 above; no clear separation seen; `Eylul_Code/gem_vs_control.forward.600_to_end.median_iqr.png` |
| How does the standard `align` pipeline work end-to-end? | §3 above |
| What was changed to recover signal past position 611 (up to and including 614), and why? | §4 above |
| Exact build steps, full-dataset run status/job IDs, output file paths, output line format, parsing examples, qualifying-read definition, and all known caveats/limitations for the tail-capture output | `TAIL_SIGNAL_CAPTURE.md` (§§3-4, 8-9) |
| What is `--tail-source`/`TAIL_SOURCE1`/`TAIL_SOURCE2`, and why a third `NEXT100EVT` pool exists (the PI's proposed `eventIndeces`/`lastMatchRawIdx` anchor) | §4.4-4.6 above; `TAIL_SIGNAL_CAPTURE.md` §§8-9 |
| How to plot/re-plot the tail-capture results | Cutoff-anchored: `TAIL_SIGNAL_CAPTURE.md` §6, `Eylul_gem_data/scripts/plot_tail_signal.py`. End-anchored, read-deduplicated (the more biologically relevant view, since Gem is expected right at the physical 3' end): `TAIL_SIGNAL_CAPTURE.md` §7, `Eylul_gem_data/scripts/plot_tail_signal_end_focused.py`. Event-count-anchored, uncapped (`ALLAFTERM`/`NEXT100EVT`), with read-level statistics: `TAIL_SIGNAL_CAPTURE.md` §10, `Eylul_gem_data/scripts/plot_allafterm_analysis.py` |
| What is `ALLAFTERM`, how does it relate to `NEXT100EVT`, and is the post-Match signal actually different gem vs control at the read level | §4.7 above; `TAIL_SIGNAL_CAPTURE.md` §10 |
| How many events are actually "rescued," which window is trustworthy, why does the first-100-event difference disappear over the full range, and where should gemcitabine's signal actually show up geometrically | §4.8.2 above; `Eylul_gem_data/scripts/analyze_rescue_window.py`; `figures_window_analysis/` |
| What is dynamic time warping being used for here, why did the first attempt need to be thrown out, and what are the validated parameters | §4.8.3 above; `Eylul_gem_data/scripts/dtw_first100_align.py`; `figures_dtw_first{100,500}_banded/` (NOT `figures_dtw_first100/`, which is the invalidated first attempt) |
| Is there an actual statistically-defensible gem-vs-control difference anywhere in the post-Match region, and how strong is it | §4.8.4 above; `Eylul_gem_data/scripts/dtw_read_level_stats.py`; `figures_dtw_first100_banded/dtw_read_level_stats.tsv` |
| Does gem vs control track as closely on the known-aligned (pre-anchor) side as DTW claims to recover on the post-Match side | §4.8.5 above; `Eylul_gem_data/scripts/before_after_match_plot.py`; `figures_before_after_{100,500}/` |
| Recommended next analytical steps | `TAIL_SIGNAL_CAPTURE.md` §5 (changepoint analysis on the tail, investigating the ~98% 3' soft-clip rate, revisiting the qualifying-rate asymmetry noted independently in both §2.2 above and there); §4.8's own bottom line calls for independent validation (held-out split or a different barycenter reference) of the events-9-20 signal before treating it as more than a lead |

**Headline status as of the full SLURM runs**: 50/3488 gem reads (~1.4%) and
321/4778 control reads (~6.7%, unique-read count — the raw qualifying-*block*
count was 347, but 25 control readIDs each had 2-3 separate qualifying
alignment blocks against this short reference; see the correction note in
`TAIL_SIGNAL_CAPTURE.md` §4) qualify for tail capture under the strict gate in
§4.3. The gem-vs-control comparison of the recovered tail signal itself (median
± IQR, anchored both from the cutoff point and from the true 3' end) shows the
two samples tracking closely with substantial overlap — similar to the §2
result in the analyzable region, no dramatic separation has emerged yet, though
some local divergence exists that's better suited to a proper changepoint
analysis than to reading the plot by eye (see `TAIL_SIGNAL_CAPTURE.md` §5-7 for
specifics and the actual figures).

**Update (2026-07-15, `ALLAFTERM` read-level analysis, §4.7)**: pooled-event
median±IQR bands still show substantial overlap, consistent with the above —
but an explicit **read-level** statistical check (Mann-Whitney U on each
read's own summary statistic, not eyeballing band overlap) over the first 100
post-Match events finds a modest, direction-consistent difference: control
reads trend more negative than gem reads as the window grows from 25→50→100
events (read-level mean p≈0.012, slope-vs-event-number p≈0.0015, fraction of
events below −1 p≈0.0046; the median itself is only borderline at p≈0.062).
This pattern **survives trimming the 3 most extreme reads per group** (p gets
tighter, not weaker, after trimming) — so it is not obviously an artifact of
a few outlier reads — but effect size is modest (common-language/AUC≈0.58,
i.e. a randomly picked gem read's value exceeds a randomly picked control
read's about 58% of the time, not a clean separation), no correction for the
~15 metric/window combinations tested was applied, and — critically — this
region is entirely past the reference-anchored, HMM-modelled part of the read,
so a real, reproducible difference here would still need to be distinguished
from a systematic gem-vs-control library-prep effect (the same open question
already flagged for the qualifying-rate asymmetry) before attributing it to
gemcitabine specifically. Treat as a hypothesis worth a proper follow-up, not
a conclusion.

**Update (2026-07-16, time-warped analysis and ground-truth comparison,
§4.8)**: the 2026-07-15 read-level signal above was measured with naive,
fixed-offset averaging over the first 100 events — that same window turns
out to be one where naive and DTW-corrected analysis actively disagree
(§4.8.4), so it should now be read in light of the following, more targeted
result rather than standing alone. Caught and fixed a real bug first: an
initial unconstrained-DTW attempt was pathological (shuffled noise matched
the reference as well as real reads did) and was rebuilt with a banded,
regularized, barycenter-referenced version whose validation is now
automatic and printed with every run. With that fixed implementation, a
read-level Mann-Whitney U on the DTW-aligned per-read median over the first
20 post-Match events is the only result (of 20 pre-registered tests, gem vs
control, several window sizes) that survives Bonferroni correction
(p=0.00053, AUC=0.65) — and it lines up with three independent lines of
evidence: the geometric prediction that gemcitabine (base 614) sits just
inside the last aligned 9-mer's own sensing window rather than far
downstream; the read-level windowed-MWU sweep in §4.8.2, which peaks at
window 500 and goes null by window 1,000; and a new ground-truth
before/after comparison (§4.8.5) showing gem and control track tightly on
the known-aligned side of the anchor, lose that tracking under naive
post-Match averaging, and substantially recover it once DTW-aligned — with
the one visible remaining divergence sitting at the same events-9-20
location. None of this is yet an independently-validated finding (the
window-20 choice isn't fully free of having already seen where the
divergence clustered), and the same library-prep-confound caveat above
still applies — but it is the most specific, best-corroborated lead this
project has produced on the post-Match region so far.
