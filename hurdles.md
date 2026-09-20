# Hurdles Log — Kushige 2013 vs. Ito 2009 rhythmic-gene comparison

> **This file is the deliverable.** Mike wants the analysis hurdles characterized,
> because whatever breaks on this data will break on the lab's own across-organism
> metagenomic data. Keep this open in a tab while you work and append to the
> ENCOUNTERED section the moment something bites — don't reconstruct it at the end.
> The friction IS the data; the grubby dataset-specific snags you'd otherwise smooth
> over are the most valuable entries.
>
> For each encountered hurdle, three things: **what bit**, **what you did**, and
> **does it transfer** to env-vs-reference and/or env-vs-env (MAG-to-MAG). That third
> line is what turns a lab notebook into what Mike asked for.
>
> Scope reminder: this run is **Version A** (compare orthologs among the RHYTHMIC
> genes; amplitude in ranks; phase as a dawn/dusk 2x2). Version B (all-orthologs,
> unbiased) is phase two — the full datasets ARE available, so it's a real option,
> but not on the few-days timeline.

---

## PART 1 — PREDICTED HURDLES (from the project brief)

These were anticipated before opening the files. Update each `Status:` line as you
actually hit it — confirm, revise, or note it didn't materialize. Predictions being
wrong is itself worth recording.

### H1 — Rhythmicity-threshold selection bias  [DEMONSTRATED]
Conditioning a correlation on "rhythmic in dataset X" is selection on a noisy variable
-> range restriction + regression to the mean -> correlation biased down, plus a
spurious "Anabaena is less rhythmic" story. Version A IS the biased version by
construction; demonstrate the trap by running it and noting that an all-orthologs
analysis (Version B) removes it.
Transfer: recurs identically the first time the lab compares rhythmic genes across its
own conditions. HIGHEST transfer value.
**Status:** _confirmed, sharply — see 2026-07-27 entries below. Restricting to either
paper's own headline rhythmic set collapses the name-matched ortholog overlap from an
ad hoc 16 genes down to 1-2, stable across the entire 50th-90th amplitude-percentile
range. Selection on rhythmicity doesn't just bias a correlation here, it nearly empties
the comparison set outright. Version B (all-RBH-orthologs) mostly fixes this -- 565
usable pairs, no two-sided rhythmicity gate -- but see the "Version B built" entry for a
residual ONE-SIDED version of H1 that survived anyway, driven by data publication, not
methodology: Ito2009 only reports a peak time for genes it calls cycling at all._

### H2 — Orthology mapping  [DEMONSTRATED, easy instance]
5,368 vs. ~2,700 genes, paralog families in Anabaena. Default RBH (clean, 1:1, discards
paralogs). FALLBACK to protect timeline: name-match on gene symbols both papers print
(kaiABC, apcE, pecB, groEL...) — for Version A you only need orthologs among the
rhythmic sets, so name-matching is a legitimate primary method here.
Split transfer: env-vs-reference = well rehearsed (ambiguity one-sided). env-vs-env
(MAG-to-MAG) = only partially rehearsed and strictly harder — no ground truth either
side, and MAG incompleteness is non-random (GC/coverage skewed) so "gene absent"
correlates with gene class correlates with rhythmicity class = a confound shaped like
biology, invisible to RBH. Do NOT read a clean result here as evidence the method
survives two MAGs.
**Status:** _name-match used FIRST (no proteome FASTAs in the repo yet; RBH looked like
its own sub-project). See 2026-07-27 "Version A" entry — name-matching's recall turned
out to be a real limiter, not just a documented caveat: only ~25-50% of hits on either
side carry a gene symbol at all. Went back and actually built RBH afterward (see
"Version B built" entry) -- turned out NOT to be the multi-day swamp the brief warned
about: both reference genomes still carry their classic locus tags as GFF `old_locus_tag`
(96.5% / 90.4% coverage), so the ID-reconciliation step that usually makes RBH painful
was clean. 1,846 genome-wide RBH pairs vs. name-matching's 1-2. Lesson: the brief's own
fallback advice ("use RBH only if it's already familiar and quick") undersold how quick
it actually was once tried -- worth trying RBH first, not last, when both genomes are
this well-curated._

### H3 — CV as a mean-dependent, cross-platform amplitude measure  [DEMONSTRATED]
"Amplitude" = CV = SD/mean, inflates for low-expression genes; two different arrays so
absolute CVs not comparable. Mitigation: compare in RANKS (Spearman), not raw CV.
Transfer: recurs, likely worse, on coverage-dependent metagenomic data.
**Status:** _confirmed and quantified — see 2026-07-27 entry. The shared numeric anchor
10^-0.6 (~0.251) that both papers use in their threshold definitions sits at the 90th
percentile of Ito2009's own amplitude distribution but only the 19th percentile of
Kushige2013's (Kushige's raw CVs run ~2.7x higher across the board: median 0.341 vs
0.124). The "same number" is not remotely the same stringency on the two platforms._

### H4 — Attenuation of correlation by measurement noise  [DESCRIBED — not executed]
Both amplitudes noisy -> any r attenuated. In principle: estimate reliability from
replicate–replicate correlation (both n=2) and disattenuate. NOT done on this timeline
(needs per-replicate values that may not be in the SIs; it's a mini-project). Named as
the correction a rigorous Version B would need.
Transfer: recurs whenever replication is thin — i.e. essentially always in env sampling.
**Status:** _**UPGRADED 2026-07-29 from "blocked by data" to EXECUTABLE.** The earlier
status said this was foreclosed because Ito publishes no per-replicate values. That was
true of Ito and false of the problem. GSE52486 (Markson 2013's WT arm) states in its own
GEO summary that its purpose was "to obtain a replicate of the wild-type LL circadian
timecourse published in Vijayan et al" — same protocol, same platform (GPL9534), 7
timepoints against Vijayan's 16. So the 7942 side now has a genuine biological replicate
pair, and Kushige is already n=2 on the Anabaena side. Both sides have what
disattenuation needs. `src/07_fetch_vijayan_markson.py` computes the Vijayan-vs-Markson
gene-wise reliability and prints sqrt(rho) as the attenuation factor. Lesson worth
keeping: "blocked by data" was really "blocked by THIS dataset", and nobody re-checked
until the dataset was up for replacement anyway._

### H5 — Phase is circular  [DEMONSTRATED]
Peak time wraps at ~24 h. Do NOT Pearson peak times. Version A uses a dawn/dusk 2x2
contingency table (Fisher + odds ratio). Fuller circular stats (Rayleigh, mean
resultant length, circular correlation coefficient) are in Schmelling 2017 if wanted.
Transfer: recurs identically.
**Status:** _pending_

### H6 — FDR / rhythmicity-call calibration  [DESCRIBED — not executed]
Both papers' rhythmicity = cosine-correlation P on a short, autocorrelated series ->
effective df below timepoint count -> P anticonservative; amplitude thresholds
stringent; gene counts (78, 97) are loose-test x tight-filter products. Use the papers'
published calls as-is. Re-calling with JTK_CYCLE / RAIN is parked.
Transfer: named as the caller-choice decision the lab faces for its own data.
**Status:** _using published calls as-is_

### H7 — Data provenance / recoverability  [DEMONSTRATED — it already bit]
Manifested during scoping: SIs unverifiable from an automated context; Kushige's raw
arrays in a dead database (KEGG Expression DB, ex0001892-ex001947). NOTE: both full SIs
were successfully downloaded manually — so the gate cleared, but the lesson stands.
Lesson for the lab's pipeline: archive every intermediate per-gene table yourself,
durable and machine-readable — journal SIs and project databases are not reliable
five-year stores.
Transfer: recurs as "did our pipeline preserve the intermediates we'll need later."
**Status:** _confirmed — full tables obtained by hand; log the source URLs + dates in
data/raw/SOURCES.txt_

---

## PART 2 — ENCOUNTERED HURDLES (append live)

The valuable part. Every time something is annoying, ambiguous, or forces a judgment
call, add an entry. Template below — copy it, don't overthink it, one line each.

<!--
### YYYY-MM-DD — <short title>
**What bit:** <the specific snag, with the actual column names / values / IDs>
**What I did:** <the fix or workaround>
**Transfers?:** <env-vs-ref? env-vs-env? recurs identically / worse / not at all — why>
-->

### 2026-07-27 — Properly-scoped Version A overlap is too small to support phase statistics, at any threshold

**What bit:** An earlier ad hoc pass (matching gene names across the FULL 2515/5336-gene
tables using each paper's loosest own "cycling vs. AR" call) produced a workable-looking
16-gene overlap and a flashy result (narM near-antiphase between species). That set does
NOT correspond to either paper's actual headline rhythmic definition:
- Kushige2013's "78 cycling genes" (Table S2, in the SI PDF -- NOT in the Data Set S1
  spreadsheet, which has no amplitude/P columns at all, only raw time series) = tier
  a+b under N+ (P<0.01-0.05, amplitude>10^-0.6 to 10^-0.7). The looser tier c (+114 more
  genes) is explicitly excluded from the paper's own "78" figure.
- Ito2009's "97 higher-amplitude genes" (per the SI PDF's Table S2 footnote, not tabulated
  by name anywhere in the SI) = P<0.001 AND amplitude>10^-0.6 in the wild-type. Filtering
  Dataset S1 directly on those numeric columns reproduces 95, not exactly 97 -- close
  enough to confirm the threshold, off by a small rounding/boundary-inclusivity gap not
  worth chasing further.

Re-matching gene names between these two PROPERLY-scoped sets: only **2 pairs** survive
(`gap1`, `kaiB`). narM is in neither set -- it's Kushige tier c (excluded from the 78) and
its Ito values (P=0.017, amplitude=0.117) are nowhere near the 97-gene bar. The 16-gene
"finding" doesn't survive contact with the papers' own definitions.

Suspecting the amplitude thresholds themselves weren't comparable across platforms
(H3), checked directly: the shared numeric anchor 10^-0.6 sits at the 90th percentile of
Ito's full amplitude distribution but only the 19th percentile of Kushige's (Kushige's
raw CVs run ~2.7x higher across the board -- different array platform, different noise
floor). Re-derived both rhythmic sets by AMPLITUDE PERCENTILE RANK within each platform's
own distribution instead of a shared raw number (P-value criteria left as each paper's
own). Result: WORSE, not better -- 0 overlapping names at the (matched) 90th percentile.
Swept the whole percentile range as a sensitivity check instead of cherry-picking one
number:

| percentile | Ito n (named) | Kushige n (named) | name-matched overlap |
|---|---|---|---|
| 90th | 96 (44) | 105 (28) | 0 |
| 85th | 102 (45) | 159 (40) | 1 (gap1) |
| 80th | 107 (47) | 197 (52) | 1 (gap1) |
| 75th | 111 (50) | 243 (63) | 1 (gap1) |
| 70th | 114 (52) | 305 (77) | 1 (gap1) |
| 60th | 117 (53) | 412 (99) | 2 (gap1, rpaC) |
| 50th | 122 (56) | 508 (114) | 2 (gap1, rpaC) |

Overlap sits at n=1-2 across the ENTIRE 50th-90th range. Not a threshold-calibration
problem -- moving the cutoff doesn't rescue it. n=1-2 cannot support a circular
phase-offset histogram, a Rayleigh non-uniformity test, or an amplitude-vs-|Δφ|
correlation; none of those were built.

**What I did:** Built the actual `data/` pipeline the brief specifies (`src/01_parse_kushige.py`,
`02_parse_ito.py` -- parsing the papers' own curated tables from the SI PDFs via
`pdftotext -layout`, not recomputing from raw arrays where avoidable; `03_orthologs_namematch.py`,
`04_percentile_rhythmic_sets.py`), with every ortholog row explicitly flagged
`ortholog_method="name-match", rbh_confirmed=False` rather than letting that distinction
get silently lost. Stopped short of forcing a percentile choice that "works" -- reported
the sensitivity table instead and let the small-n result stand as the finding.
One process note: `pdftotext`-extracted text from a PowerPoint-derived PDF contained
stray `\r` characters from overstruck text that Python's universal-newline `readlines()`
silently split on (but `wc`/`grep` don't), desynchronizing line-indexed parsing by ~30
lines with no error -- fixed by reading bytes and splitting on `\n` only. Also: writing
scratch files to plain `/tmp` in this environment was unstable (content differed between
back-to-back reads of the same path) -- moved to the session scratchpad dir instead.

**Open question, not resolved:** whether the small overlap reflects real biology (Kushige's
highest-amplitude genes are heavily nitrogen-fixation/heterocyst-related --
nifB/nifS/nifX/nifV1/ctaC/hupC -- and phycobilisome genes; Ito's are core clock machinery,
transporters, and phototaxis genes -- almost disjoint functional programs) or is a
name-matching recall artifact (only ~25-50% of hits on either side carry a gene symbol at
any percentile; a real ortholog with an unnamed partner on either side is invisible to
this method). Cannot be resolved without RBH.

**Transfers?:** Both directions, and this is probably the single highest-value entry in
this register for that reason. env-vs-reference: identical mechanism -- if the lab's own
rhythmicity threshold and the reference paper's threshold aren't calibrated to the same
percentile of each platform's/method's own noise distribution, any "shared rhythmic set"
comparison is comparing differently-selected populations before any biology enters, and a
small resulting n may look like a clean negative result when it's actually an artifact of
the intersection procedure. env-vs-env (MAG-to-MAG): recurs and compounds -- two
non-reference genomes means the H3 percentile-mismatch problem AND the H2 name-matching
recall problem both apply on both sides simultaneously, with no anchor to calibrate
against. Recommend for the lab's pipeline: (1) never compare two platforms' rhythmicity
calls on shared raw thresholds without checking where those thresholds land in each
platform's own full distribution first; (2) treat a small name-matched overlap as
ambiguous (biology vs. method-recall) by default, not as a negative result, unless RBH or
equivalent backs it up.


---

### 2026-07-27 — Version B built: RBH orthologs, no rhythmicity gate (mostly)

**What bit:** After the Version A writeup above, decided to actually build Version B
rather than just recommend it. Fetched both organisms' own reference-genome proteomes
from NCBI (RefSeq GCF_000012525.1 for 7942, GCF_000009705.1 for Anabaena/Nostoc PCC
7120 -- the exact genomes each paper's own locus tags come from). Confirmed the ID
chain is clean before committing: GFF `old_locus_tag` -> shared gene `ID` -> CDS
`protein_id` -> FASTA header, on both genomes, no ambiguous (duplicated)
`old_locus_tag` mappings. Ran reciprocal `blastp` (evalue 1e-5, best hit by bitscore,
not `-max_target_seqs` which isn't a strict best-hit guarantee) -> **1,846 RBH ortholog
pairs**, no rhythmicity filter anywhere in that step.

Joining those pairs against full per-gene tables (not just the rhythmic subsets) hit an
asymmetry that only showed up at execution time: Ito2009's spreadsheet has a literal
`"-"` placeholder instead of a peak time for every gene it doesn't call cycling (1,715
of 2,515 genes) -- there's no raw time series in Ito's public data to fall back on and
compute an alternative phase estimate the way Kushige's data allows (Kushige2013 DOES
publish per-timepoint raw values for all 5,336 genes, letting a peak-time proxy be
computed for every one of them regardless of significance). Dropping rows with an
undefined peak time therefore reimposes a rhythmicity filter -- but only on the 7942
side, and it's a DATA-AVAILABILITY constraint from what Ito chose to publish, not a
methodology choice made here. n: 1,846 RBH pairs -> 1,787 with expression data on both
sides -> **565 with a defined peak time on both sides** (7942's own gene population
is 800/2515 = 31.8% cycling by its own call, and 565/1787 = 31.6% of the ortholog
subset survives the same filter -- consistent, so this isn't distorting the *rate*,
just re-narrowing the set to "cycling-in-7942" specifically).

Then ran the actual requested analysis on those 565 pairs: circular Δφ (wrapped to
±12h), 4h-binned histogram (`figures/06_phase_offset_distribution.png`), a from-scratch
Rayleigh test (no astropy/pycircstat available -- implemented directly, Zar's
higher-order p-value correction), and a Spearman check of |Δφ| against each pair's
weaker-side amplitude PERCENTILE RANK (not raw CV -- H3 again: raw CV isn't comparable
across platforms, so ranked each gene within its own platform's full distribution
first, took the worse rank per pair).

**Results:**
- Rayleigh: p=1.4e-7 (n=565) rejects uniformity, but weakly -- mean resultant length
  R-bar=0.167 out of 1. Mean direction -4.6h.
- Shape: closest to **(a) unimodal with a long tail**, not (b) bimodal dawn/dusk-swap.
  Histogram peaks at [-4,0)h (24.1%); antipodal mass (|Δφ|>=8h, 31.3%) sits BELOW the
  uniform baseline (33.3%) -- the opposite of what a real dawn<->dusk cluster-swap
  would produce.
- Amplitude: Spearman rho=+0.014, p=0.74 between |Δφ| and weaker-side amplitude
  percentile -- **no relationship**. Bottom-quartile-amplitude pairs have essentially
  the same mean |Δφ| (5.61h) as everyone else (5.81h). At n=565 the "big offsets are
  just low-amplitude noise" hypothesis (motivated by the earlier ad hoc narM look) does
  NOT hold up.
- narM: Δφ=-12.0h, landing exactly on the antipode (the one point on the wrap where the
  sign is inherently ambiguous -- noted explicitly in the histogram's axis label).
  Amplitude confidence solid on both sides (44th/92nd percentile) -- not a weak-signal
  gene, so this reads as a genuine large offset, not measurement noise. But n=1: can't
  call one gene "representative of a class" either way.
- kaiB: Δφ=-8.8h, 98th-percentile amplitude in 7942 (as expected -- it's the clock
  gene) but only 27th-percentile in Anabaena (borderline low). Independently confirms
  the brief's own §4 stated prior ("kai: top-amplitude in 7942, flat in Anabaena")
  without having gone looking for it.

**What I did:** Extended `VersionB_AllOrthologs.ipynb` through genome fetch -> locus/protein
mapping -> RBH -> full-table join -> Δφ histogram -> Rayleigh test -> amplitude overlay
-> narM/gap1/kaiB lookup, all executed end-to-end (`data/processed/rbh_orthologs.csv`,
`data/processed/all_orthologs_rbh.csv`, `figures/06_phase_offset_distribution.png`).
Genome provenance logged in `data/raw/genomes/SOURCES.txt`. Also hit and fixed a
separate git issue along the way: `origin/master` had diverged (someone created, then
deleted, an empty stub file at the notebook's exact path directly on GitHub) --
resolved with a standard merge (net-zero on the remote side, no real conflict), and
hardened `git_workflow.py`'s `save_and_push` to surface the actual git error on push
failure instead of a bare traceback.

**Deliverable summary (paste-ready):**

> The circular phase-offset distribution across 565 confirmed RBH orthologs (defined
> peak time on both sides) is closest to **unimodal with a long tail** (Rayleigh
> p=1.4×10⁻⁷ rejects uniformity, but R-bar=0.167 is weak; mean direction −4.6h), not
> the bimodal dawn/dusk cluster-swap hypothesis — antipodal mass sits below, not above,
> the uniform baseline. |Δφ| does **not** track amplitude confidence (Spearman
> ρ=+0.014, p=0.74; bottom-quartile-amplitude pairs show the same mean |Δφ| as
> everyone else), so large offsets in general are not simply low-amplitude measurement
> noise. narM itself sits at the extreme antipodal edge (Δφ=−12.0h) with solid
> amplitude confidence on both sides (44th/92nd percentile) — not obviously
> noise-driven, but a lone outlier at n=1 rather than representative of a broader
> class, while kaiB (Δφ=−8.8h, 98th/27th percentile) independently reproduces the
> brief's own predicted kai-gene divergence.

**Transfers?:** The RBH-was-actually-easy result transfers directly to env-vs-reference
(same "well-curated anchor genome" situation the brief already expected). The
one-sided-filter-from-data-availability finding transfers to BOTH cases and is worth
generalizing: whenever two data sources are joined, check whether either side silently
drops non-significant/non-called records rather than reporting a null value for them --
that silently reimposes exactly the kind of selection filter this whole exercise was
built to avoid, and it will not announce itself with an error.

### 2026-07-29 — H10: selecting a dataset on availability metadata instead of design

**What bit:** Asked to find a replacement for the Ito 2009 side, I ranked candidates on
"does it have archived per-timepoint data" and recommended **Piechura, Amarnath & O'Shea
2017 (GSE104204)** as the strongest option — in a slide deck, to be sent to two PIs.
Reading the actual GEO record afterwards: GSE104204's design is *simulated natural light*,
"Low Light" and "Clear Day" arms sampled at 0.5, 2, 4, 6, 8, 9, 10 and 12 h **since dawn**,
plus high-light/shade pulse perturbations and RpaA/RpaB/RNAP ChIP-seq. That is a 12-hour
diel window under entraining light — **not a free-running LL time course, and less than
one full cycle.** A 24 h cosine cannot be fitted to it and a peak time cannot be estimated
from it. The recommendation was wrong, and wrong in a way that would have survived review
by anyone who also only read the abstract and the data-availability statement.

The failure mode is specific and worth naming separately from H7/H8: those are about
whether data *exists* and whether a join silently drops it. This is about a dataset being
fully archived, well documented, correctly described — and structurally unable to answer
the question, for reasons that appear nowhere in the metadata fields you would filter on.
"RNA-seq, circadian, S. elongatus, raw reads in GEO" is true of GSE104204 and useless.

**What I did:** Pulled the GEO records for every candidate and tabulated *design* —
entrainment, free-run vs. entrained, sampling interval, duration, cycles covered,
replicates — before comparing anything else. That immediately reordered the list and
surfaced GSE18902 and GSE52486, neither of which had been top-ranked on the availability
criterion. Rewrote the deck slide.

**Transfers?:** Both cases, and it is the cheapest hurdle in the register to avoid. Any
reference dataset the lab pulls for an env-vs-reference comparison has to be screened on
experimental design first — for rhythms specifically: was it free-running or entrained,
and does it span at least one full cycle at a sampling interval finer than the feature
you are trying to resolve. A dataset can be perfectly archived and still not contain the
axis you need. Recommend a written design checklist applied *before* any download, since
the download-and-inspect step is where sunk cost starts accumulating.

---

### 2026-07-29 — Two silent bugs in the new cosinor, caught only by synthetic ground truth

**What bit:** Wrote `src/07_fetch_vijayan_markson.py` to fit one identical cosinor to both
organisms. Two defects, neither of which raised an error:

1. **Peak time offset by a constant.** Wrote the acrophase as `arctan2(-g, b)/w` instead of
   `arctan2(g, b)/w`. For `y = M + A·cos(w(t−φ))` the design matrix gives
   `b = A·cos(wφ)`, `g = A·sin(wφ)`, so `φ = atan2(g, b)/w`. The sign error shifts every
   peak time by a fixed amount that depends on the true phase — recovered peak was 14.0 h
   for a true 10.0 h, and 4.0 h for a true 20.0 h. **Invisible within one dataset** (the
   distribution of peak times just rotates) and it corrupts every cross-dataset Δφ, which
   is the entire deliverable. Genes with a true peak at 0 h came back correct, which is
   exactly enough to make a spot-check pass.
2. **ZeroDivisionError on a perfect fit.** `r2` rounds to exactly 1.0 on noiseless input,
   so the F-statistic divided by zero. Guarded on `ss_res > 0`, which does not catch it —
   `ss_res` is tiny-but-nonzero while `1 - r2` is exactly zero.

**What I did:** Tested the fit against synthetic series with known peak times and
amplitudes rather than eyeballing real output. All seven test phases now recover to
<1e-14 h. Also ran a 2,000-draw pure-noise null: false-positive rate 4.2% at nominal 5%,
0.90% at nominal 1% — calibrated on *independent* noise, which does not clear H6's
autocorrelation concern (the p-value is still anticonservative on a real 4 h-sampled
series; kept it for ranking, not as a calibrated FDR). Added a bin-edge snap so a peak at
0 h cannot come back as 23.999 and land in the wrong 4 h bin.

**Transfers?:** Directly and to everything. The general lesson is that a phase or
periodicity estimator must be validated against simulated data with a known answer before
it touches real data, because the characteristic failure is a *plausible* number, not a
crash or an obvious outlier. A constant phase offset applied to both sides of a comparison
would have cancelled and hidden; applied to one side, as here, it becomes the result. For
the metagenomics pipeline: any Δφ, period, or amplitude estimator gets a synthetic-recovery
test committed alongside it, including the degenerate cases (perfect fit, all-NaN, fewer
points than parameters, peak exactly on a bin boundary).

---

### 2026-07-29 — Hardcoded sample→timepoint maps are an unchecked assumption

**What bit:** `src/07_fetch_vijayan_markson.py` maps GEO sample accessions to timepoints
with two hardcoded dicts. The Markson one was written from sample titles read off the GEO
record. The Vijayan one was **inferred** — from the accession range on the GPL9534 platform
page plus the paper's stated design — without ever opening GSE18902. It happened to be
right (GSM468463 = T=24 h through GSM468478 = T=84 h, with 468479–485 the novobiocin arm),
but "happened to be right" is not a property you can rely on twice.

The failure mode if it had been wrong is the dangerous kind: a permuted sample→time
assignment produces a completely plausible cosinor fit, a plausible amplitude, and a
plausible peak time, with no error and nothing visibly odd. Same silent shape as the
acrophase sign bug logged above.

**What I did:** Added `verify_sample_map()`, which parses `!Sample_title` out of the series
matrix and asserts the hardcoded map against what GEO actually ships, aborting on any
mismatch. Tested against the verbatim real titles from both series, plus four negative
controls: map shifted by one timepoint, two samples transposed, a novobiocin sample pulled
in by mistake, and an accession absent from the series. All four abort.

Writing the check surfaced a bug in the check: the hour regex used `\b` after the `h`,
which never fires in `WT_LL_36h_replicate` because underscore is a word character — so it
rejected all seven of Markson's genuinely-correct samples. Fixed with a negative lookahead,
which also keeps `T = 5 minutes` from matching as an hour.

**Transfers?:** Both cases, and it generalises past GEO. Any time sample identity, ordering,
or condition labelling is asserted in code rather than read from the data, it needs a
runtime check against the source — because the wrong answer looks exactly like the right
one. For the metagenomics pipeline this is sample-sheet validation: assert that the
sample→timepoint, sample→site, and sample→treatment maps agree with the sequencing
metadata before anything downstream runs, and fail loudly rather than proceeding. Note also
that the validator itself needed testing; a check you have not tried to make fail is not
yet a check.

---

### 2026-07-30 — H11: a merged spreadsheet header silently mislabelled a whole column

**What bit:** `VersionB_AllOrthologs.ipynb` §3 builds `corr_p_anabaena` from
`kushige_raw["Unnamed: 7"]`, treating it as Kushige's cosine-correlation p-value for
rhythmicity. It is not. It is **the t-test p-value for Ri, the heterocyst-enrichment
ratio** — an entirely unrelated quantity.

The cause is a merged header cell. Row 5 of `Kushige2013_so2.xlsx` labels column 6
`Ri value (t-test: p-value)`, and that single label spans two columns: col 6 holds the
ratio, col 7 holds the p in parentheses (`0.8638...`, `(0.301854...)`). Pandas names the
unlabelled second column `Unnamed: 7`, which reads like an anonymous overflow column
rather than half of a two-cell header. Nothing errors; the values are plausible
p-values in [0, 1]; and the column then travels under a name asserting it means
something else.

Four checks that identify it, none of which require the SI PDF:
- Correlation with Ri: rho = -0.237. Correlation with N+ amplitude: rho = +0.046. A
  rhythmicity p should track amplitude and not the enrichment ratio; this does the
  opposite.
- Only **44%** of tier-a genes fall below 0.05, though tier a is *defined* as p < 0.01.
- **19%** of arrhythmic genes fall below 0.05, which is impossible for the p-value that
  produced the arrhythmic call.

The real values were in the SI PDF's Table S2 the whole time: canonical *kaiB*
(`alr2885`) has rhythmicity p = **0.046**, amplitude 0.246, peak 21 h. The long-form
paralog `all3328` has no published rhythmicity p at all, because Table S2 lists only the
78 cycling genes and it is arrhythmic.

**This was a misattribution, not a fabrication — the distinction matters.** Every value
involved is a real cell, correctly read and correctly rounded: the 0.018 quoted in
`runningNotes.md` is literally `(0.0175558176693055)` at row 3546 col 7, and the 0.399
is `(0.399104984899883)` at row 1431 col 7. Nothing was invented and nothing was
misparsed. What was wrong was the *name* the column travelled under, and therefore the
quantity every later reader believed it represented. Note also that this register's own
2026-07-27 entry already recorded that Data Set S1 "has no amplitude/P columns at all,
only raw time series" — so the correct fact was on record, and the notebook's *naming*
contradicted it a day later without the conflict being noticed. The log was right
throughout; only the label was wrong.

**What I did:** Corrected the kaiB paralog table in `runningNotes.md`, which had been
citing 0.018 and 0.399 as rhythmicity p-values for the two paralogs; both are
enrichment p-values. Confirmed blast radius: `corr_p_anabaena` is carried in
`data/processed/all_orthologs_rbh.csv` but was never used as a filter or a statistic in
any headline result — the 565-pair set was cut on peak-time availability, and the
Rayleigh, delta-phi and amplitude analyses use tier calls and amplitude ranks. So the
main Version B numbers stand. The column needs renaming to `ri_ttest_p_anabaena`, and a
real `corr_p_anabaena` can only be populated for the 78 genes in Table S2.

**Transfers?:** Both cases, and it is nastier than it looks. Merged and multi-row headers
are near-universal in supplementary spreadsheets, and the failure is silent in the worst
way: no exception, values in the right range, and a variable name that actively asserts
the wrong meaning to every later reader — including the person who wrote it. Two
generalisable defences: (1) after parsing any supplementary table, assert that each
column behaves as its name claims — a rhythmicity p must be below the stated threshold
for every gene the paper calls rhythmic, and that one assertion would have caught this
immediately; (2) treat every `Unnamed: N` column pandas produces as an unresolved
header, never as data, until its provenance is traced back to the raw cell layout. For
the metagenomics pipeline the same applies to any externally-supplied metadata or
annotation table.

**Wider point about where to aim scrutiny.** The instinct on finding a wrong number in a
log is to ask whether the values were invented. Here they were not, and checking took
one command — read the literal cells and compare. That check is cheap and should be the
first move every time, because it partitions the problem: if the values are real, the
bug is in naming, joining or interpretation, and no amount of re-verifying the source
data will find it. Sources of error in this project so far have been overwhelmingly of
that second kind — a sign-flipped acrophase, an inferred sample map, a mislabelled
column — all cases where the numbers were real and the meaning attached to them was not.
Verification effort should be weighted accordingly: check provenance once, then spend
the rest of the effort on assertions about what each quantity is supposed to mean.

Moved to `runningNotes.md` (2026-07-27) to keep this file focused on the hurdle
register itself. Open questions, side findings (kai genes / long-form Anabaena kaiB
paralog, etc.), and build notes for next time now live there.

### 2026-07-30 — H11 fix executed: rename, real corr_p_anabaena, post-parse assertion

**What bit:** Executing the fix itself (PART D0) surfaced one thing the discovery entry
above didn't: the tier-a/b assertion it recommended ("every tier a gene must have
p < 0.01, every tier b < 0.05") failed on first run — one gene, `alr0611` (nrtD),
tier a, `corr_p_anabaena = 0.010`. Not a second bug: Table S2 reports p to 3 decimal
places, so a true value of e.g. 0.0098 displays as "0.010" and trips a strict `< 0.01`.
Fixed by adding a ±0.0005 tolerance (half the last reported digit) to the boundary
check, not by loosening the rule itself.

**What I did:**
1. Renamed the Ri-derived column to `ri_ttest_p_anabaena` throughout
   `VersionB_AllOrthologs.ipynb` §3 and `data/processed/all_orthologs_rbh.csv`.
2. Populated a real `corr_p_anabaena` by joining Table S2's 78-gene parse
   (`data/interim/kushige_rhythmic.csv`, from `src/01_parse_kushige.py`, re-run fresh
   to confirm it still reproduces 78 rows) onto the full 5,336-gene table on locus tag —
   null for every gene outside the 78, not reconstructed from the spreadsheet.
3. Added the post-parse assertion as new cells (§3b): checked it PASSES on the corrected
   column (78/78, with the rounding tolerance) and FAILS on the old mislabelled column
   (46/78 tier a/b genes, 59%, violate the bound) — the failure is the proof the rename
   was necessary, not just cosmetic.
4. Updated §8 (kaiB paralogs) to print "not published (arrhythmic, absent from Table S2)"
   for `all3328` instead of a bare `nan`, now that `corr_p_anabaena` is genuinely null
   there rather than silently wrong.
5. Re-ran the full notebook end to end (disabled the auto-push cell first so a
   verification run can't also commit). Every headline number reproduced exactly:
   Rayleigh p=1.373e-07, R-bar=0.1667 (n=565); |Δφ| vs. amplitude rho=+0.014, p=0.743;
   narM/gap1/kaiB spot-checks all matched their previously-recorded percentiles. The
   565-pair set is unchanged, as expected — it was always cut on peak-time availability,
   never on `corr_p_anabaena`.

**Transfers?:** The rounding-boundary edge case is the generalizable part, and it's the
same lesson H11 itself already drew, one level down: an assertion meant to catch a wrong
*meaning* can itself misfire on a wrong *tolerance* if it doesn't account for the
source's reporting precision. Any post-parse assertion built against a paper's own
published thresholds needs a tolerance sized to that paper's stated decimal precision,
not the raw float — otherwise a legitimate boundary case reads as a second bug and
either gets wrongly "fixed" by loosening the real rule, or wastes time being
re-investigated as if it were novel. For the metagenomics pipeline: when asserting a
downstream call against an upstream publication's own thresholds, check how many
significant figures that publication reports before picking an equality tolerance.

### 2026-07-30 — PART A: RefSeq WP_ protein-sharing and a genuine duplicate row both broke the "clean 1:1 ortholog table" assumption

**What bit:** Building `08_amplitude_by_function.ipynb` (amplitude-vs-amplitude
correlation on the full ~1,787-pair set, no rhythmicity filter), the spec's required
"assert no duplicated locus_tag" check did NOT pass cleanly, and tracing why surfaced two
independent, previously-unlogged issues:

1. **RefSeq `WP_` accessions are a non-redundant-protein namespace, not a 1:1 gene key.**
   18 of the 1,846 RBH pairs (11 on the 7942 side, 7 on the Anabaena side) share a
   `protein_id` with another pair -- e.g. `Synpcc7942_0893` and `Synpcc7942_1389` both
   map to `WP_011242480.1`. Confirmed at the raw GFF-derived locus<->protein map, not an
   RBH-step artifact: RefSeq deliberately assigns ONE `WP_` accession to every genomic
   locus encoding a byte-identical protein (tandem duplicates, here). RBH computed at the
   protein level therefore fans out to every locus tag sharing that protein when joined
   back to locus tags on either side.
2. **Ito2009's own spreadsheet lists `Synpcc7942_2452` twice**, with two different
   measured amplitudes (0.123 vs 0.157) -- verified as a real duplicate row in
   `Ito2009_sd1.xls`, not a parsing artifact.

Neither is large (combined <2% of pairs) and neither moved the headline numbers --
rho=+0.056, n=1787 reproduces exactly with both included. But (2) caused a real bug
while building the notebook: an early draft computed amplitude percentiles by merging
`ito_full`/`kushige_full` onto the pair table a SECOND time (once to build the pairs,
once again to attach percentile ranks). Because both the pair table and `ito_full` still
carried the `Synpcc7942_2452` duplicate key, the second merge squared the fan-out instead
of carrying it through once -- n silently went from 1787 to 1789, a difference of exactly
the affected gene's duplicate multiplicity. No exception, no obviously wrong number,
just two rows where the reproduction target expected two others.

**What I did:** Fixed the immediate bug by computing percentile ranks on `ito_full` /
`kushige_full` BEFORE the single join, not by re-merging afterward -- eliminates the
double-merge fan-out risk structurally rather than special-casing this one gene. Kept
both duplicate-row types in the analysis (reported, not silently dropped or
deduplicated), since removing them would change n away from the spec's own validated
1,787 without instruction to do so. Documented both directly in the notebook (§1b) with
root cause, not just the symptom.

**Transfers?:** Both cases, and (1) especially generalizes past this project. Any
ID-mapping pipeline that joins through a "non-redundant" accession namespace (RefSeq
`WP_`, UniRef cluster IDs, vsearch/CD-HIT cluster representatives -- all common in
metagenomic annotation) inherits the same many-to-one hazard: the accession is 1:1 with
a *sequence*, not with a *genomic locus* or a *MAG bin*, and duplicated/tandem genes
across bins will fan out silently on any join keyed by it. The squared-fan-out mechanism
generalizes further: re-merging a lookup table onto an already-joined table, when either
side still carries a duplicate key, multiplies rather than adds the duplication --
general lesson is to compute per-key derived columns (ranks, in this case) BEFORE the
join that consumes them, not after, whenever the join key's uniqueness isn't guaranteed.

---

### 2026-07-30 — PART A ad hoc keyword overlay: short gene-name abbreviations collide with ordinary English words in free-text annotation

**What bit:** The spec's suggested keyword list for the nitrogen-fixation flag was
`nif, fdxH, hetR, hetN, nitrogenase`, matched as a bare case-insensitive substring
against the `annotation` free text on either side. It returned 3 pairs. The user
caught the problem immediately on seeing the output: *S. elongatus* PCC 7942 is not a
diazotroph and has no nitrogen-fixation genes at all -- so a nonzero count was itself
the tell, before any number needed to be doubted on statistical grounds.

Pulling the three matched rows showed the mechanism: 7942's annotations read "cysteine
desulfurase **NifS**", "nitrogen regulation protein **NifR3** homolog", "**NifU**-like
protein". These gene names genuinely derive from the nitrogenase gene cluster
historically, but the protein families they name -- general Fe-S cluster biogenesis
(NifS/NifU-type cysteine desulfurase and scaffold, also called IscS/IscU or SufS/SufU in
non-diazotrophs), and tRNA-dihydrouridine synthase (NifR3) -- are broadly conserved
housekeeping machinery, not restricted to nitrogen fixers. The match was on legacy
nomenclature, not on pathway membership. Confirmed exhaustively: of the 18 genes in
Anabaena's own spreadsheet with a genuine `nif*` gene symbol (the real nitrogenase
structural/FeMo-cofactor genes -- nifHDKENB, nifWXTZ), **zero** have a 7942 RBH
ortholog. The only two that do (`nifJ2` = pyruvate:flavodoxin oxidoreductase, a second
`nifS` paralog = class-V aminotransferase) are the same generic-homolog trap.

**Checking further, unprompted, because the failure mode is generic and not
gene-specific:** the other two ad hoc flags in the same section used the identical
bare-substring approach and had the identical problem, worse. `"pec"` (intended for
phycoerythrin genes pecA/B/C) matched almost entirely inside the word "**spec**ific"
(occurs constantly in annotation prose: "heterocyst-specific", "pyrimidine-specific",
"site-specific", ...). `"sod"` (intended for superoxide dismutase) matched almost
entirely inside "**sod**ium" (sodium-dependent transporters, sodium ATP synthases,
sodium symporters -- zero of them superoxide dismutase). Checked directly against the
raw annotation columns with a word-boundary regex: `nif`, `sod`, `kat`, and bare `pec`
as literal standalone tokens **occur zero times** in either paper's annotation column.
The real genes are there (1-2 "superoxide" mentions, 1 "catalase" mention, plenty of
real nitrogenase genes on the Anabaena side) but the annotations spell the terms out in
full; nobody abbreviates them the way the keyword list assumed.

**What I did:** Replaced the single bare-substring approach with two different fixes
depending on what actually distinguishes signal from noise for that keyword:
- `psa`/`psb`/`apc`/`cpc`: kept as short forms, but matched only when they appear as a
  **capitalized gene-name token** (`PsbE`, `CpcF`) in the ORIGINAL-case annotation text,
  not a lowercased substring search -- this is what separates a real gene symbol
  embedded in prose from an English word sharing the same three letters. Verified by
  hand against all matches: 18 (Ito) + 15 (Kushige), 100% genuine, zero collisions.
  Applying the same capitalized-token method to `nif` did NOT fix it -- NifS/NifU/NifR3
  are themselves genuinely capitalized gene-name tokens, just not diazotrophy-specific
  ones, so capitalization alone can't distinguish this case.
- `sod` -> `superoxide`, `kat` -> `catalase`: match the spelled-out biological term
  instead of an abbreviation the data never uses.
- `nif` -> dropped entirely; kept as the phrase `nitrogenase` / `nitrogen fixation`,
  which is what actually separates the real nitrogenase-pathway genes from the
  NifS/NifU/NifR3 false positives (none of which are annotated with either phrase).

Result after the fix: nitrogen_fixation n=0 (correct -- reported explicitly as a
finding, not left as a silent blank: the category structurally cannot be populated in a
7942-anchored RBH comparison, because 7942 has nothing for real nitrogenase genes to be
orthologous to). photoprotection n=3, all three now genuine (two photolyase genes, plus
`sodB` -- a real superoxide dismutase gene that the ORIGINAL "sod" abbreviation search
never found even once, because it's spelled out in full and only the "sodium" collisions
were ever showing up). photosynthesis/phycobilisome n=14 (down from a false n=33), all
verified genuine.

**Transfers?:** Directly, and this is probably the highest-transfer entry in the
register alongside the H11 mislabelled-column entry, for the same underlying reason:
short abbreviations searched against free-text annotation collide with ordinary
language at a rate that isn't obvious until someone actually reads the matches. Two
separable lessons for the lab's metagenomic annotation pipeline: (1) a 3-4 letter
keyword needs either a word-boundary check or a length/frequency sanity check
(`\bKEYWORD\b` occurring zero times in the full annotation corpus, as happened here for
sod/kat/pec, is itself informative -- it means the abbreviation isn't how the data
actually writes that gene, not that the gene is absent) before trusting any hit count;
(2) even a "clean" gene-name-style match is not the same as pathway membership when
gene nomenclature carries historical baggage -- Nif, but also Rec, Mot, Fla and many
other bacterial gene-name prefixes were coined in one pathway's discovery context and
later found in paralogous families serving unrelated general functions. A biologically
plausible zero (no diazotrophy genes in a non-diazotroph) is a stronger signal to
double-check a keyword search than an biologically implausible nonzero would have
been -- and in this case it was the user, not the code, that caught it, which is itself
worth recording: a domain-knowledge sanity check ("does this organism even have these
genes") caught what no statistical or code-level check would have.

---

### 2026-07-30 — PART B: one cosinor caller on both organisms turns three DESCRIBED hurdles into DEMONSTRATED ones

**What happened:** `src/09_recompare_with_vijayan.py` fits the identical cosinor
(imported from `07_fetch_vijayan_markson.py`, not re-derived) to log2-transformed
Kushige N+ intensities, giving both organisms a genuinely commensurable amplitude and
peak time for the first time in this project. Executing it (network access to NCBI
confirmed working first, contra the earlier "sandbox has no outbound access" note --
that was an environment-specific limitation of an earlier session, not a permanent
constraint) produced concrete before/after numbers for three hurdles that were
previously only described:

- **H3 (CV as mean-dependent, cross-platform amplitude)** -- concretely demonstrated,
  not just avoided. kaiB (`alr2885`) sat at the 27th amplitude percentile in Anabaena
  under Kushige's raw CV; under the new log2-cosinor amplitude it sits at the **90th**.
  Verified this is not a join bug: the new fit's peak time (21.1h) matches Kushige's own
  published Table S2 value (21h) almost exactly, r2=0.51, p=0.0006. The old CV number
  was penalizing kaiB for being highly-expressed (CV is SD/mean, and high mean deflates
  CV even at constant absolute swing) -- exactly the H3 failure mode named at project
  start, now visible in a single named gene rather than only in an aggregate median
  comparison.
- **H4 (attenuation by measurement noise)** -- executed for the first time (status was
  "permanently blocked" as recently as 2026-07-27, then "upgraded to executable" on
  2026-07-29 once the Markson replicate was identified; this is that upgrade cashed in).
  7942-side reliability (Vijayan vs Markson, gene-wise rho=+0.611) and Anabaena-side
  reliability (Kushige's 1st vs 2nd N+ replicate series, rho=+0.197) combine to an
  attenuation factor of 0.347 -- raw amplitude-vs-amplitude rho=+0.132 disattenuates to
  +0.381 [CI +0.249, +0.511]. Both reported; the raw number is not replaced by the
  corrected one, and the CI is wider after disattenuation, not narrower, as the method
  predicts (a corrected point estimate is not a more certain one).
- **H1 (rhythmicity-threshold selection bias)** -- the clearest single demonstration in
  the project so far. The circular phase-offset Rayleigh statistic on the OLD,
  Ito-peak-time-filtered set (n=565, implicitly "cycling in 7942" by construction, H8)
  gave R-bar=0.167. The NEW set, unfiltered because Vijayan's cosinor gives every gene a
  peak time regardless of significance (n=1811, essentially the full RBH table), gives
  R-bar=0.092 -- still highly significant (p=2.5e-07, n is much larger) but the apparent
  clustering strength drops by nearly half once the implicit rhythmicity filter is
  removed. This is H1's predicted mechanism (selection on a noisy variable inflates
  apparent structure) shown with a real before/after number on the same data, not just
  argued for.

**What I did:** log2(0) is undefined and 684/128,064 Kushige N+ readings (0.53%,
concentrated in 464/5336 genes, almost always 1-2 of 24 points) are exactly zero -- a
real detection floor, not a rounding artifact (checked: all exactly 0.0, no small
negatives). Dropped per-gene before fitting rather than pseudocounted, since an
arbitrary pseudocount would arbitrarily inflate that gene's fitted amplitude; the
cosinor's own n>=4 guard confirms no gene is left underdetermined (worst case 12/24
points remain). Documented in the script's module docstring rather than left as an
undocumented `vals[vals>0]` filter.

**Transfers?:** All three, directly, and this entry is really the payoff of several
earlier ones (H1's original prediction, H3's original demonstration, H4's status
upgrade) landing together once the better dataset made a real "before vs. after" run
possible instead of only an argued one. For the lab's own pipeline: (1) when a paper's
own summary index (CV, a fold-change call, a normalized score) is swapped for a
raw-value refit, expect specific genes to move a lot, not just the aggregate correlation
to shift a little -- check named genes, not only summary statistics, before trusting a
methods change hasn't silently reversed a conclusion. (2) A disattenuated correlation
is a legitimate number to report but must ship with its (wider) CI attached, never as a
bare point estimate -- otherwise "correcting for noise" reads as "found a stronger
result" rather than what it actually is. (3) H1's mechanism is not specific to this
dataset: any comparison that implicitly restricts to "significant in condition A" before
joining to condition B will show inflated apparent structure relative to the unbiased
join, by a magnitude worth actually measuring (here, roughly halved R-bar) rather than
assuming is small.

---

### 2026-07-30 — PART C: getting real external data out from behind three different bot walls, and what showed up once it was in hand

**What bit:** Every one of the three external sources PART C needs (Arbel-Goren 2024's
motif table, Markson 2013's ChIP target list, Mitschke 2011's TSS annotations) is
gated by a DIFFERENT automated-access defense: the ASM journal site is Cloudflare
JS-challenge-protected; PMC's own new download endpoint requires solving a
proof-of-work challenge; PNAS is also Cloudflare-gated. Plain `curl` cleared none of
them -- this generalizes the H7 finding (SI PDFs needed manual download) to 2026:
publisher sites now actively defend against exactly this kind of automated retrieval,
not just passively lack an API.

**What I did:** Found the actual working route is Europe PMC's `supplementaryFiles`
REST endpoint, which mirrors the raw files for any article in the PMC **open-access**
subset, bypassing both bot walls entirely (`https://www.ebi.ac.uk/europepmc/webservices/
rest/{PMCID}/supplementaryFiles`). Worked for Arbel-Goren 2024 (CC BY, in the OA
subset). Did NOT work for Markson 2013 or Mitschke 2011 -- neither is in the OA
subset (checked via the same API, which reports this explicitly rather than silently
failing). For those two, used NCBI's `oa.fcgi` service to confirm the same thing
independently, then found Markson's raw ChIP-seq deposit on GEO (GSE51093) -- but it
turned out to contain only per-timepoint `.wig` signal tracks, no called peaks or
gene assignments. Re-deriving a target gene list from raw signal would mean
reimplementing their whole peak-calling and gene-assignment pipeline -- a categorically
larger and more error-prone undertaking than downloading a table, and exactly the kind
of "re-run the analysis as though it were open" the spec explicitly says not to do for
the FIMO scan. Stopped there rather than doing it, and used the partial 58-gene subset
already present in Arbel-Goren's own Table S2 as an explicitly-labelled stand-in
instead (documented precisely in `data/raw/rpaA/SOURCES.txt`, including which
direction its incompleteness biases the result). Same treatment for the missing
Mitschke background: used Arbel-Goren's own scanned-gene set as a lower-bound proxy,
paired with a second "full array" background as an upper-bound sensitivity check, so
the result is reported as a bracket rather than a single number resting on an
unverifiable denominator.

**What the real data actually showed, once obtained (not simulated, not assumed):**
- Reproduced the paper's "81 genes" exactly from the raw spreadsheet (81 rows at FIMO
  q<0.05) -- but only 55 of those 81 carry an assigned CyanoBase gene ID; the other 26
  are TSS/intergenic windows with no gene to test for overlap with anything. The
  paper's own headline number silently includes rows that cannot enter any gene-level
  analysis.
- **pecB — the gene the paper names as its headline example (P=1.5e-5) — does not
  actually clear their own q<0.05 threshold** (q=0.189). Found by pulling the real row,
  not by re-deriving anything: both the exact p-value and the caption text check out
  against the published PDF, so this isn't a parsing error on my end, it's the paper
  citing a raw p-value for a gene its own multiple-testing-corrected criterion excludes.
- The hypergeometric overlap (81/55 motif genes x Kushige's 78 rhythmic genes) is
  **null on both background definitions** (p=0.30 proxy, p=0.19 sensitivity; observed
  2 overlaps against an expected ~0.8-1.1) -- this directly tests
  `journal_club_kushige_2013.md`'s CONJECTURE REGISTER item 4 ("Prediction: the 81x78
  hypergeometric will show enrichment. Untested.") and the prediction does not hold.
  Cross-validated against the paper's own "Is in Kushige list" annotation column,
  independent of my own join logic: exact match (2 genes both ways), which is real
  confidence the locus-tag join itself isn't the problem.
- The two overlap genes (`alr4077`, `asr1667`) peak within 0.6h of each other (3.4h,
  4.0h) -- n=2 is nowhere near enough to test, but worth stating plainly rather than
  omitting since it's the kind of descriptive detail C3 asks be reported individually.
- C1, by contrast, comes back strong and clean: 7942's highest-amplitude genes (Vijayan
  cosinor fit, Part B) really are enriched for RpaA ChIP-target status even on the
  partial 58-gene proxy (median 88th vs 50th percentile, 95% CI [+15,+44] points,
  Mann-Whitney p=1.2e-8). Mike's premise holds; the follow-up relationship to the
  rhythmic gene SET in Anabaena does not, at least not by this test with this
  background. Two different questions, two different answers -- exactly why C1 and C3
  had to be kept separate rather than treated as one result.
- GC content, flagged "verify" in the journal club notes: computed directly from the
  actual RefSeq assemblies rather than trusted from memory -- 55.43% (7942) / 41.27%
  (Anabaena), confirming the "~55%/~41%" figures used throughout.

**Transfers?:** The access-barrier finding transfers directly and is worth its own
line: as of this session, three major publisher/repository platforms independently
chose to block plain HTTP retrieval of their own supplementary data with
JS-execution-requiring challenges, while an EU-funded aggregator (Europe PMC) mirrors
the same files with no such gate, for anything in the open-access subset. For the
lab's own literature-mining needs: check Europe PMC's supplementaryFiles API before
assuming a paywalled or bot-gated source is unreachable, but check `isOpenAccess`
first rather than discovering the gate case-by-case. The methodological finding
transfers too, and is the sharper one: a paper's own headline number ("81 genes") can
silently include entries that don't survive contact with a gene-level analysis (26 of
them here), and a paper's own cited example (pecB) can fail the paper's own stated
significance threshold -- both are only visible by pulling the real underlying rows,
never by trusting the prose summary of them. That is the same lesson as H11
(mislabelled column) and the nitrogen-fixation keyword bug, a third time: check
provenance once, at the row level, before trusting what a table is asserted to contain.

---

### 2026-07-30 — PART C follow-up: quantifying the RpaA/RpaB ambiguity, and a real multiple-testing bug caught mid-build

**What bit:** Following up on C4's qualitative RpaA/RpaB objection (a FIMO scan with the
RpaA matrix can't distinguish RpaA sites from RpaB sites, since RpaB binds a different,
known motif -- HLR1 -- that's documented to overlap RpaA sites at the *kaiBC* promoter),
a real HLR1 position weight matrix was built from Riediger et al. 2019's published,
curated motif list (959 real instances, obtained the same way as the other RpaA data --
Europe PMC's open-access mirror) and used to scan the same promoter windows Arbel-Goren
scanned for RpaA. The very first version of this scan reported a **100% "hit rate"** --
every one of the 55 RpaA-motif genes' promoters also scored above the chosen HLR1
significance threshold. That number was wrong, and wrong in the classic silent way: no
error, a plausible-looking print statement, and a result dramatic enough that it should
have been the first thing distrusted rather than the first thing reported.

**Root cause:** the significance threshold had been built from the distribution of a
**single position's** PWM score on random genomic sequence (20,000 random 18-mers), but
it was being compared against a **window-best** statistic -- the maximum score over
~1,064 positions (both strands) in each 551bp promoter window. Taking the best of ~1,000
draws from a distribution will blow past that distribution's own single-draw 99th
percentile almost regardless of whether there's real signal, purely from the number of
attempts. This is an unannotated multiple-testing problem: the "null" wasn't null for
the statistic actually being thresholded.

**What I did:** built the *correct* null -- the distribution of the window-best score
over many random 551bp windows (same procedure, same window length, just random
genomic position instead of a real promoter) -- and used its 95th/99th percentiles
instead. Kept the wrong version's number in the notebook explicitly, printed alongside
the corrected one ("For contrast, the WRONG single-position 99th pct was X -- using
that would have flagged nearly every window as a 'hit'"), rather than quietly replacing
it, so the failure mode is visible to a future reader who builds a similar scan-and-max
statistic. Also validated the PWM itself before either threshold was trusted: 100% of
the 959 real training instances score above the single-position null's 99th percentile
-- confirming the scoring function discriminates correctly, which is a separate check
from "is the threshold itself calibrated for what's being compared against it."

**Corrected result, once the null was fixed:** ~40% of the 55 named RpaA-motif genes
(95th-pct threshold, 8x the 5% expected by chance) also carry a plausible HLR1/RpaB site
in the same promoter window; only ~5% clear a stricter 99th-pct bar, close to the 1%
chance rate. The group as a whole is significantly shifted above the corrected null
(Mann-Whitney p=1.4e-15) -- a real, broad enrichment, not noise, but a minority-to-
substantial-minority one, not "all of them" and not "none of them". Both genes
responsible for C3's small hypergeometric overlap (`alr4077`, `asr1667`) fall in the
ambiguous set at the 95th-pct threshold. This narrows C4's objection from a blanket
qualitative caveat to a quantified fraction, without overturning C3's null conclusion --
if anything it reinforces it, since the two genes driving that overlap turn out to be
exactly the kind of site the objection warns about.

**Transfers?:** The substantive finding is Anabaena/RpaA-specific, but the bug
transfers directly and is worth its own line in the register: **whenever a threshold is
built from scanning a window and taking the best score, the null distribution must be
built the same way (scan-a-window-and-take-the-best on random sequence), not from a
single-position score.** This is not a rare mistake -- it is the default failure mode of
any "does this window contain a match" analysis (motif scanning here, but the identical
shape recurs in peak calling, anomaly detection over a sliding window, or any
best-of-k-candidates test), and it fails silently: the output is a plausible number, not
an error, exactly like the acrophase sign-flip bug (2026-07-29) and the sodium/pec
false positives (2026-07-30). The general defence, again: before trusting a threshold
derived from one procedure, confirm the null was generated by the *same* procedure
applied to the *same* statistic, not a related-but-simpler one.

---

### 2026-07-30 — a diagnostic that only prints, doesn't filter: 103 non-gene probes rode silently through Part B and C

**What bit:** Investigating raw data for the new Part F (global amplitude comparison)
required inspecting `data/interim/7942_ll_timeseries_long.csv` directly for the first
time since Part B was built, and it surfaced 53 rows with `locus_tag_7942` values like
`"640711016"` -- not a real locus tag. Traced to `07_fetch_vijayan_markson.py`'s
`read_platform_map()`: it computes `good = tbl["locus_tag_7942"].str.match(...)`,
*prints* the match percentage as a diagnostic ("GPL9534: 2715 probes, 2612 (96.2%)
match..."), warns if the rate is below 90% -- and then returns `tbl` completely
unfiltered. The diagnostic looked healthy (96.2%, above its own 90% warning bar) and
was read as "checked and fine" when it was actually "checked, printed, and ignored."
103 non-gene probes flowed through every downstream step: the cosinor fit
(`7942_ll_cosinor.csv`, reported as "2715 genes" when only 2612-2662 were real), the
percentile ranking (`amp_pct_vijayan`, denominator inflated by up to 103), and both
Part B's and Part C's read of that file.

Traced the bad rows to their source: GPL9534's `SPOT_ID` column reads `"JGI_IG: <id>"`
for all of them -- JGI intergenic-region tiling probes. Agilent/GEO leaves their `ORF`
field as a self-referential echo of the probe's own ID when there's no gene to name,
which is correct platform behavior, not a data error. The bug was entirely in this
project's own code not enforcing the check it already computed.

**Second-order catch, worth its own line:** the first fix (`^Synpcc7942_\d{4}[a-z]?$`,
requiring exactly 4 bare digits) dropped the 103 junk probes correctly but ALSO
silently dropped 50 genuinely real genes -- `Synpcc7942_B####`, the plasmid-encoded
genes (confirmed against every `old_locus_tag` in the RefSeq GFF: 48 `B`-prefixed, 53
`R`-prefixed, both real replicons). This was caught immediately, not accepted, because
of a habit from earlier in this project: after any change to a join key's filter,
recheck the RBH match count. It moved from 1846/1846 to 1834/1846 -- a regression that
had no business happening from removing rows that never matched anything -- and that
was the signal to investigate rather than ship. Corrected regex:
`^Synpcc7942_[A-Z]?\d{3,4}[a-z]?$`, verified against literally every `old_locus_tag` in
the genome with no exceptions, restored 1846/1846.

**Impact on already-reported results, checked not assumed:** re-ran Part B and Part C
end to end on the corrected data. Every headline number is unchanged at the reported
precision (rho=+0.132, R-bar=0.092, kaiB 99th/90th, C1's amplitude gap +38.3->+38.5
points) -- the junk probes never carried a real locus tag, so they never entered any
RBH-joined or gene-symbol-matched analysis. What DOES change: the "2715 genes" /
"~2715 probes" population-size figures used in prose throughout Part B and C, which
should read ~2,662 real genes -- relevant now because Part F explicitly needs a clean
full-population gene count for both organisms.

**Transfers?:** Directly, and it's a distinct failure mode from the others logged here
(H11's mislabelled column, the nitrogen-fixation keyword bug, the HLR1 null
miscalibration) -- those were all *wrong values under a plausible name*. This one is
*a check that exists in the code, executes, prints a correct diagnostic number, and
still doesn't do anything*. The general lesson: a `print()` after a boolean mask is not
a filter, no matter how good the code around it looks -- if a row shouldn't be trusted,
the code needs an actual `.loc[good]`/`.drop()`, not a warning message a human might
not read closely enough to act on. For the lab's own annotation pipelines, this is
exactly the shape of bug most likely to hide in an ID-mapping step: a validity check
computed for logging purposes that a later refactor (or the original author, under
deadline) never wired up to actually gate the data. And the second-order lesson stacks
on top of the recurring one in this project: any fix to a filter needs the SAME kind of
before/after sanity check (here, the RBH match count) that caught the original bug --
tightening a regex is exactly as capable of introducing a false negative as the
original diagnostic-only code was of missing a false positive.

---

### 2026-07-30 — PART F: the raw global-amplitude comparison pointed the wrong way, and the noise floor was the whole answer

**What bit:** Part F asks whether Anabaena and 7942 differ in overall oscillation
amplitude, motivated by needing to know whether Anabaena *kaiB*'s Part-B jump (27th to
90th percentile once CV was replaced by a cosinor fit) reflects real absolute signal or
just a high rank inside a globally quieter organism. The RAW answer, computed first and
taken at face value initially, said the opposite of both hypotheses on the table:
**Anabaena's raw fitted amplitude is higher than 7942's**, not lower (matched-pair
median 0.192 vs 0.141, n=1811, Wilcoxon p=6.5e-26; same direction in the full
populations, p=2.0e-57). Taken naively, that would have read as "Anabaena is the louder
oscillator" -- backwards from Uzumaki et al. 2004's independent report of lower
amplitude for the Anabaena Kai system, and worth being suspicious of for that reason
alone before trusting it.

**What I did:** built the noise floor the spec called for rather than treating it as
optional -- permuted time labels within each gene and refit the identical cosinor
(500 randomly sampled genes x 1,000 permutations per dataset, ~500,000 null fits each,
~25s per dataset at ~15-20k fits/sec). Anabaena's noise floor came back **2.42x** higher
than 7942's (0.170 vs 0.070 median null amplitude) -- larger than the 1.37x raw
amplitude gap itself. Expressing each gene's amplitude as a ratio to its own dataset's
noise floor and re-running the paired comparison **flipped the sign**: 7942 sits at a
median 2.0x its own noise floor, Anabaena only 1.1x its own (paired Wilcoxon
p=8.7e-119, opposite direction from the raw comparison). The raw comparison was
measuring which assay is noisier, not which organism oscillates more.

*kaiB* specifically: Anabaena *kaiB* is 90th percentile within Anabaena but only ~2.5x
Anabaena's noise floor; its 7942 ortholog is 99th percentile within 7942 and ~13x
7942's noise floor -- a 5x gap in noise-relative terms that the raw percentile alone
cannot show. Kushige's "kai genes barely cycle in Anabaena" headline stands.

**Why this is trusted and not just a plausible-looking flip:** (1) the one identified
confound that could rescue the raw reading -- two-colour ratio-array compression
understating 7942's amplitude -- points the WRONG way to explain the flip (it would
suppress 7942's number, not inflate it), so it can't be the artifact producing this
result; naming a confound and then checking whether it actually helps the alternative
explanation, rather than just listing it, is what makes this differ from `hurdles.md`'s
other sign-flip catches. (2) the magnitude is not marginal -- p=8.7e-119 is not a
borderline call decided by which correction was chosen.

**Transfers?:** Directly, and it's the cleanest demonstration in this project of why
F2's framing ("the noise floor is the analysis, not a caveat") is correct rather than
rhetorical. Any cross-platform or cross-organism amplitude/effect-size comparison
built on differently-sampled data (different n, different replication, different assay
type) needs its OWN per-dataset null before the raw comparison is trusted even
directionally -- not just for calibrating a p-value, but because the raw SIGN can be
wrong. For the lab's metagenomics work this is the single highest-transfer item of the
three amplitude-related hurdles logged today (platform-map filtering, HLR1 null
miscalibration, this one): comparing rhythmicity strength across samples, taxa, or
sequencing depths without a matched noise floor per comparison arm risks reporting the
assay's noise characteristics as if they were biology, and the direction of the error
is not predictable in advance -- it has to be measured, the way it was here.

---

### 2026-07-27 — Bulk transcriptomics averages away minority cell types
**What bit:** Kushige's bulk N− array is a whole-filament average, but heterocysts
are only ~1 per 10–15 cells (~5–10% of the population). So ~90% of the signal in the
N− array is vegetative-cell RNA and heterocyst-specific expression is diluted roughly
10-fold in the mix. This is why Kushige could NOT read heterocyst behavior off the
bulk N− data and had to run a separate heterocyst-enrichment experiment (80%-heterocyst
fraction) plus a hetR-null comparison to isolate genuinely heterocyst-specific genes
(all1427, hesAB). A minority-cell-type signal is only visible after enrichment or via
single-cell methods; in a bulk average it can be invisible even when it's strong per-cell.
**What I did / implication:** Do not treat the bulk N− array as reporting heterocyst
metabolism — it reports the vegetative majority. Any heterocyst-specific question must
go to the enrichment dataset (a different data type; see DEFERRED thread).
**Transfers?:** DIRECTLY and severely to the metagenomics work. Every environmental
sample is a mixed community; a metagenome/metatranscriptome is a population average
weighted by each taxon's abundance. A rhythm (or any expression pattern) in a 5–10%
member is diluted ~10–20x in bulk signal and can be undetectable — while an abundant
member's pattern dominates and can look like "the community's" behavior when it's one
taxon's. Quantify each taxon's fractional abundance BEFORE interpreting any bulk
pattern, and treat "signal absent in bulk" as "absent OR diluted below detection,"
never as "absent." This is the mixed-population version of the CV/normalization problem:
low apparent amplitude can be a dilution artifact, not biology. Enrichment, binning to
per-MAG resolution, or single-cell/single-taxon methods are the analogs of Kushige's
heterocyst fraction.

---

### 2026-07-30 — Reannotating "hypothetical" genes: three hurdles in one pipeline, all caught before shipping

**What bit (1): "manual lookup" in a spec is not automatable just because everything
else is.** Mike's suggestion was HHpred via the MPI Bioinformatics Toolkit web server,
by hand, for the 28 rhythmic genes still labelled hypothetical -- exactly what he did
for `all0232`. A `WebSearch` for a REST/API path came back negative: the Toolkit is
interactive-only for HHpred specifically (a different, unrelated tool on the same site,
COMER, does have an API -- easy to conflate the two and wrongly conclude "the Toolkit
has an API"). **What I did:** did not attempt to script a browser-form submission to
approximate automation, and did not fabricate placeholder results. Instead ran the
cheap, real, automatable step first (HMMER vs Pfam-A) and used it to *shrink* the manual
step: 11 of the 28 already had a `confident` Pfam call (including `all0232` itself,
independently landing on the same NYN domain Mike found by hand -- a real cross-method
agreement, not copied from his result), and 6 more had no current RefSeq sequence to
search at all (see hurdle 3). Only 11 of the original 28 actually needed a human at the
web form; those 11 sequences were staged as ready-to-paste FASTA files rather than left
as a vague "go do HHpred" instruction. **Transfers?:** Directly. "Automate everything"
and "sequence the work so the expensive/manual step might shrink or vanish" are
different instincts, and the second one is usually right when a tool is genuinely
interactive-only. Before scripting around a web tool's absence of an API, check whether
the *need* for that tool can be reduced first -- cheaper automatable methods often
already answer most of what the manual step was going to be asked.

**What bit (2): naive lexical matching between two naming vocabularies manufactures
disagreement that isn't there.** Comparing Pfam family descriptions to RefSeq/PGAP
product names via stopword-stripped word overlap flagged `"Major Facilitator
Superfamily"` vs. `"MFS transporter"`, `"Radical SAM superfamily"` vs. `"...
methyltransferase RlmN"`, and `"AhpC/TSA family"` vs. `"thioredoxin family protein"` as
disagreements. All three are the same call: MFS is the literal acronym, RlmN is a
textbook radical-SAM enzyme, AhpC/TSA is a thioredoxin-fold peroxiredoxin family. Pfam
names domains by family history; PGAP/NCBIfam names proteins by function -- literal
word overlap can't see through that gap, and naively reporting the raw overlap rate as
"the agreement rate" would have been wrong in the misleading direction (undersells
agreement, oversells disagreement). **What I did:** caught this by spot-checking the
"disagreement" list before writing it into the notebook, not after -- reported the
lexical-overlap number explicitly as a documented underestimate/floor rather than a
corrected true rate, and restricted the actual "disagreements to review" table to the
`confident` evidence tier only, so the list a human would actually read is the smallest,
highest-quality slice rather than 700+ rows dominated by vocabulary-gap noise.
**Transfers?:** Directly, and probably the single most reusable lesson of this session
for the metagenomics pipeline: any automated cross-database annotation comparison
(Pfam vs. KEGG vs. eggNOG vs. a MAG-caller's own product names) will have this same
vocabulary-mismatch failure mode. An "agreement rate" computed by string matching alone
is a floor, not a measurement, unless the naming conventions are known to be aligned --
state which one, on every table this shape produces.

**What bit (3): current RefSeq is not a superset of the 2013/2009 gene calls.** 406 of
the 3,934 hypothetical genes checked (400 Anabaena, 6 Synechococcus) have no protein_id
in the current RefSeq annotation at all -- confirmed by direct search of the current
GFF (zero occurrences of the locus tag anywhere in the file, not just missing from the
old_locus_tag cross-reference), not a parsing artifact. These genes cannot be searched
by HMMER, HHpred, or anything else run against current sequences, because there is no
current sequence. Six of the 28 rhythmic hypotheticals fall in this bucket. **What I
did:** did not skip silently or drop these rows -- flagged them explicitly
(`in_current_refseq = False`) in `reannotation_summary.csv`, reported the count as its
own line in the unknown-register breakdown rather than folding it into "no Pfam hit,"
and did not attempt to recover sequences from the original 2013 CyanoBase gene models
(a real, different data-recovery task, out of scope here and flagged rather than
quietly worked around). **Transfers?:** Directly. Any comparison between an old gene
catalog and a current reference annotation needs an explicit "still exists in the
current annotation at all" check before any downstream absence is interpreted as
"unknown function" rather than "gene call retired." Silently treating a missing
protein_id as just another `unknown` would have hidden a genuinely different failure
mode (annotation churn) inside a bucket meant for "we looked and found nothing."

---

### 2026-07-31 — HHpred: neither hit count nor top-hit probability is the signal, convergence is

**What bit:** 11 rhythmic-hypothetical genes went to manual HHpred after the Pfam pass.
Several came back with 150-250 hits above 90% probability, which looks at first glance
like an ambiguity problem (too many high-confidence answers to pick from). It isn't --
and treating "top probability" as the confidence measure would have produced both a
false negative and a false positive in this same batch. `all3516` (737 residues): the
top 8 hits are named completely different proteins from unrelated organisms and
pathways (a bacteriophage aspartate phosphatase, a human G-protein-signalling
modulator, a fly cell-polarity protein, a trypanosome flagellar motor protein) -- but
nearly every description contains "Tetratricopeptide repeat," all ~99.7% probability,
all covering the same ~400-residue region. That's not ambiguity, it's convergent
evidence for one recurring fold, and the top-ranked hit's specific protein identity
("Rap105") is not the answer -- treating hit #1 as *the* result would have produced a
confident, specific, and wrong label. Conversely, `all4578`'s top hit alone was 88.5%
probability -- comfortably above any of this pipeline's own Pfam confidence thresholds
-- but the next several hits were a lipoprotein, a zinc-resistance protein, a malaria
antigen, and a pilus protein, no shared theme, matching only a ~20-residue N-terminal
fragment: a single number that looked trustworthy in isolation was actually the
short-fragment false-positive regime.

**What I did:** built `summarize_convergence()` (`src/parse_hhr.py`) instead of reading
off probability or hit count directly -- checks whether independent hits *agree*: same
region of the query (>=50% reciprocal range overlap with the top hit) and a shared
description theme across >=40% of those corroborating hits. Caught a second, smaller
bug while building it: an early version's keyword extraction let `{Homo sapiens}`
organism tags and bare "SCOP" mentions (present on many unrelated hits purely as
citation boilerplate, not as a shared classification) masquerade as the "consensus"
keyword for two genes, including manufacturing a false consensus for `alr4939` that a
manual read of the same 15 hits didn't support -- fixed by stripping organism tags,
ligand codes, crystal resolution, and raw SCOP-code text before keyword counting, same
category of fix as the RefSeq lexical-overlap hurdle earlier in this file (metadata
words masquerading as content words). `alr4939` itself is the cleanest illustration of
why this matters: the single best hit (91.7%, mannose-6-phosphate isomerase) has almost
no corroboration, while a different, internally consistent 6-hit cluster (51-65%,
Type-III-secretion/flagellar-motor-switch family) competes for the same region -- two
real candidate identities, correctly left unresolved rather than reporting the higher
single number as the winner.

**Transfers?:** Directly, and it generalizes past HHpred specifically to any tool that
returns a ranked list of hits against a redundant database (BLAST against nr, Foldseek
against PDB, any structure/sequence search where the same fold or family is represented
by dozens of independently deposited entries). The redundancy of the target database
itself produces exactly this pattern -- many hits, high scores, all "agreeing" with each
other in substance while disagreeing in the specific identity attached to hit #1 -- and
a pipeline that reports "best hit" without checking whether the runners-up corroborate
it will silently convert a fold-level convergence signal into a false specific-identity
claim. The fix is the same shape every time: define corroboration in terms of the
things that are hard to fake (same region, independently repeated theme after metadata
is stripped out), and treat a lone high-scoring hit as weaker evidence than several
consistent moderate-scoring ones, not stronger.