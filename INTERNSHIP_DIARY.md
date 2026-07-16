# Internship Diary

A day-by-day log of what I actually did on the DNAscent/gemcitabine project,
kept in my own words so I can look back later and remember not just *what*
I built but *why*, and what it found. Technical reference detail (exact code,
file paths, commands) lives in `PROJECT_OVERVIEW.md` and
`TAIL_SIGNAL_CAPTURE.md` — this file is the narrative version, for the
internship report.

Entries are added at the bottom as new days happen.

---

## 16/07/2026

Today was entirely about the region of a nanopore read that comes *after*
the last base DNAscent's HMM could confidently align to the reference — the
"post-Match" or `ALLAFTERM` signal we'd built in an earlier session, and
where gemcitabine incorporation is actually expected to show up. The day
moved through four connected pieces of work: checking our implementation
against my PI's exact description, figuring out *which* window of events is
actually worth looking at, building (and then fixing a serious bug in) a
dynamic-time-warping analysis to correct for the fact that different reads
move through the pore at different speeds, and finally building a
ground-truth comparison to sanity-check all of it. I finished by writing
everything up in `PROJECT_OVERVIEW.md`.

### 1. Checking the implementation against my PI's description

After my meeting with my PI, he gave a much more precise description of how
the alignment mechanism should work, with exact line citations into the
upstream DNAscent GitHub repository. Rather than assume my earlier
implementation matched what he meant, I had it fetch those exact lines and
compare directly:

- The way the code maps reference positions to read positions (from the BAM
  alignment) and read positions to raw signal events (from DNAscent's own
  rough aligner) matched his description exactly.
- The formula used to scale/normalize the raw signal — necessary because
  every physical nanopore is slightly different, so you have to put all
  reads on the same scale before comparing them — was also character-for-
  character identical to what he pointed at upstream, and confirmed to be
  applied consistently to every part of the "rescued" tail signal, not just
  the normal aligned part of the read.
- One thing didn't match his suggested *implementation approach*, though not
  in a way that mattered for correctness: he suggested labelling the
  rescued events as a new category directly in the human-readable alignment
  string. What we'd actually done was print them as separate, clearly
  tagged blocks after the normal output instead. Functionally equivalent,
  just a different way of achieving "make it obvious these events come
  after the last aligned base" — worth flagging honestly since he said the
  exact method didn't matter as long as it was correct.

### 2. How many events are actually "rescued," and which window should I trust?

Before building anything new, I needed a real, numbers-based answer to a
question I'd been going back and forth on: when I plot the first 100 events
after the last aligned base, I see a visible difference between gem and
control — but when I plot *all* of the rescued events (some reads have over
700,000 of them!), that difference disappears. Is the full-range plot
telling me the effect isn't real, or is something else going on?

I had a proper analysis built to check this, and it answered the question
clearly:

- **How much data is there really?** Gem reads have a median of about 3,300
  rescued events each; control reads about 3,850 — not wildly different.
  But control has one single read with over 725,000 rescued events, which
  by itself makes up nearly a third of *all* of control's rescued signal
  combined. That one read is clearly some kind of outlier/artifact, not
  representative of anything biological.
- **Coverage drops off fast.** Both samples have basically every read still
  contributing out to about 1,000 events past the anchor point, but after
  that the number of reads still that long shrinks quickly — by 3,000-4,000
  events in, more than half the reads have already run out. So any
  "median" computed way out at, say, event 50,000 isn't really an average
  of many reads anymore — it's just whatever the one or two longest reads
  happen to look like.
- **A read-level statistical test (not just eyeballing the plots) confirmed
  the same story numerically**: comparing gem vs control read-by-read, the
  difference is weak at the first 100 events, strongest and actually
  statistically significant around the first 500 events, and then
  completely disappears by 1,000 events and stays gone all the way out to
  the full uncapped range — even after removing that one giant outlier
  read. So the full-range plot isn't wrong; the signal really does wash out
  once you include everything past a few hundred events.
- **The geometry backs this up.** I worked out, from how DNAscent reports
  positions (the middle of a 9-base window centred on each aligned
  position), that the last position we can confidently align is base 611,
  and gemcitabine sits at base 614 — only 3 bases later. Since the
  "window" DNAscent looks at is 9 bases wide, gemcitabine is actually
  *already inside* the sensing window of the very last aligned base, and
  should fully pass through that window within about 7 more bases, which
  at the expected pore speed is roughly 35 raw events. That means the first
  100-500 events is exactly where I'd expect gemcitabine to show up
  biologically — going out to 1,000+ events is mostly just adding
  irrelevant, unrelated sequence and diluting any real effect, which
  matches exactly what the statistics showed.

### 3. Dynamic time warping — and a bug I'm glad I caught before trusting the results

My PI's key concern was what he called "the elephant in the room": since
these post-Match events have no reference position at all, and different
reads move through the pore at different speeds, "the 7th event after the
anchor" in one read might correspond to a completely different physical
base than "the 7th event" in another read. Averaging across reads without
correcting for this could either blur out a real signal or, worse,
manufacture a fake one. His suggested fix was dynamic time warping (DTW): a
technique that stretches and compresses one signal to line it up with
another based on where their shapes actually match, rather than assuming
every read moves at the same speed.

I implemented this — picked one representative read as a reference and
aligned every other read to it — and got a result that looked almost too
good: gem and control's aligned signals overlapped almost perfectly. That
immediately made me suspicious, and I'm glad it did. I ran a control test:
I took real reads, scrambled their values into pure random noise, and ran
the same alignment process on the scrambled version. **The scrambled noise
matched the reference just as well as the real data did.** That meant the
DTW algorithm, as I'd first written it, was flexible enough to warp
*anything* — including garbage — into looking like a match. The
beautiful-looking result was an artifact of the algorithm, not a real
finding.

I fixed this properly rather than just tweaking parameters until the
scrambled-noise test happened to pass once:

- Added a constraint (a "band") limiting how far in time the alignment is
  allowed to shift a point, so it can't cheat by matching things that are
  far apart.
- Added a penalty for unnecessary stretching, similar to a gap penalty in
  DNA sequence alignment, so the algorithm only warps when it actually
  improves the match, not just because it can.
- Replaced the single reference read with a "consensus" signal built by
  averaging many reads together over a few rounds — much harder for random
  noise to accidentally resemble than one specific noisy read.
- Built the scrambled-noise control test permanently into the script itself,
  so it runs and prints a pass/fail verdict automatically every time, not
  just as a one-off manual check I might forget to redo later.

With the fix in place, the same test now clearly shows real reads matching
much better than scrambled noise, with a comfortable safety margin. I ran
this for both the first 100 and first 500 events, and made both a version
using each read's median value and one using the mean, since I originally
wasn't sure which one I'd used for the first version of the plot (it was
median) and wanted both available side by side.

### 4. Testing statistically whether the DTW-corrected signal is real

With a validated, trustworthy alignment, I ran a proper statistical test
(read-by-read, not just looking at whether the bands on a plot overlap) at
several different window sizes decided in advance from the geometry above,
so I couldn't be accused of just picking whichever window looked best after
the fact. Correcting for the fact that I ran multiple tests, exactly one
result survived: the first 20 events after the anchor, using the
DTW-corrected median, showed a real, statistically significant difference
between gem and control. Interestingly, at that same 20-event window, the
*uncorrected* version showed nothing at all — a genuine case of the time-
correction recovering a real signal that naive averaging was hiding.
Slightly further out (about 35 events, matching my geometric estimate for
where gemcitabine's influence should fully clear), the pattern flipped: DTW
showed nothing while the uncorrected version showed a weak, borderline
signal. I also ran a version restricted to the exact points that looked
most different in an earlier plot, out of curiosity, but flagged clearly
that testing on points chosen because they looked different is circular
and shouldn't be treated as real evidence — only the windows I picked in
advance count as a fair test.

### 5. A ground-truth sanity check: before vs. after the anchor

As a final check, I built one more comparison: what does the signal look
like for the last several hundred events *before* the anchor — i.e. the
part of the read that IS confidently aligned to the reference, where we
know for certain that two reads at the same position really are looking at
the same physical base? I realized this didn't need any new code changes at
all — the information was already being written out by the program, just
not in a form anyone had extracted before. I wrote a new parser to pull it
out directly from the existing output files.

The result was reassuring: before the anchor, where we have genuine ground
truth, gem and control track each other very tightly, confirming this
shared-signal pattern is real and not some artifact of how the values are
scaled. Right after the anchor, without any time correction, that tight
tracking collapses immediately. With the DTW correction applied, most of
that tight tracking is recovered smoothly across the boundary — and the one
place gem and control still visibly pull apart is that same narrow region
around events 9-20 after the anchor that the statistics had already flagged.
I made this comparison at both the ±500-event and ±100-event scale (the
500-event version has fewer usable reads than I expected, because it turns
out the reference sequence itself is only about 600 bases long, so asking
for 500 aligned bases *before* the anchor uses up almost the entire
reference — a genuine, useful thing to have discovered along the way).

### Where things stand

There's a real, narrow, and now statistically-defensible signal sitting
right around events 9-20 after the last aligned base — right where the
geometry says gemcitabine should be — that only shows up once you correct
for reads moving through the pore at different speeds, and that's now
backed up by four independent angles (the coverage/statistics analysis, the
validated DTW plots, the formal statistical test, and the before/after
ground-truth comparison). It isn't a fully independent, publication-ready
result yet — the window I tested wasn't chosen completely blind to earlier
observations, and I still can't rule out that some of what I'm seeing is a
systematic difference between how the gem and control samples were
prepared rather than gemcitabine chemistry itself. But it's by far the most
specific and well-supported lead this project has produced so far, and a
clear next step for follow-up validation.

Everything from today is written up in detail in `PROJECT_OVERVIEW.md`
section 4.8.

### 6. Making the work portable to a different server

Near the end of the day I realized I'd want to continue this same Claude
session later from a different server, still on the same account, so I
asked what the best way to carry everything over would be — not just the
conversation, but the actual scripts and the reasoning behind them. This
turned into a useful discovery: `PROJECT_OVERVIEW.md` and
`TAIL_SIGNAL_CAPTURE.md` had already been getting committed and pushed to
GitHub as I went (good — that meant all of today's write-up was already
safe), but the entire `Eylul_gem_data/` directory — where every script I
wrote today actually lives — was `.gitignore`'d, so none of the code itself
had ever been tracked. Fixed `.gitignore` to keep the big data (`.align`
files, caches, generated figures) excluded while explicitly allowing the
scripts and SLURM job files to be tracked, and staged (not yet committed —
leaving that step to me) every script from today plus a few earlier ones
that had never been committed either. Also hit a small scare in the
process: this very diary file had somehow disappeared from disk between
being written and being staged, with no trace in git history — had to
recreate it from what was still in the conversation. Worth remembering for
future sessions: commit new files reasonably soon after creating them,
rather than assuming they'll still be there later in the same session.
