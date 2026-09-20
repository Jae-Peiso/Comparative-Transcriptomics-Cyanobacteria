# Analysis specification — six linked pieces of work

Pre-declared specification for the next round of analyses. Written before the work was
carried out, so the statistical choices below (primary vs. secondary measures, tiering,
expected coverage) are pre-registration, not post-hoc description — `src/amplitude_stats.py`
cites this file's Tier 1 for exactly that reason.

Facts throughout — column names, row counts, category coverage, expected statistics — were
verified against the files on disk at the time of writing.

| Order | Part | What | New file |
|---|---|---|---|
| 1 | **0** | **Correction** — a mislabelled column. Run first; everything inherits it | patch existing |
| 2 | **A** | Amplitude correlation coloured by functional category | `src/08_amplitude_by_function.py` |
| 3 | **B** | Fetch the better 7942 dataset, redo the comparison with one caller | `src/07_*` (exists, unrun) + `src/09_recompare_with_vijayan.py` |
| 4 | **C** | **RpaA — questions 3 and 4.** Must follow B | `src/12_rpaa_overlap.py` |
| 5 | **D** | Heterocyst-enrichment layer | `src/10_heterocyst_layer.py` |
| 6 | **E** | The *kai* paralog question | `src/11_kai_paralogs.py` |
| 7 | **F** | Do the two organisms differ in overall amplitude? | `src/13_global_amplitude.py` |

**Two hard ordering constraints:**

- **Part 0 before everything.** It corrects data already on disk.
- **Part C immediately after Part B.** The RpaA premise is a claim about *the
  highest-amplitude 7942 genes*, and Part B is what produces properly recomputed
  amplitudes for them. Running C on Ito's published amplitudes instead would be testing
  his premise against the very numbers Part B exists to replace.

Parts D and E are independent of each other and of C; run them whenever.

Part A is listed before B because the current-data version of it is worth keeping — but
B changes A's inputs, so if you do both, run A first, keep the result, then redo it
after B and compare. That comparison is itself a deliverable (see B6).

---

# PART 0 — CORRECTION: a mislabelled column. Run this before anything else.

`VersionB_AllOrthologs.ipynb` §3 builds `corr_p_anabaena` from
`kushige_raw["Unnamed: 7"]` and treats it as Kushige's rhythmicity p-value.
**It is the t-test p-value for Ri, the heterocyst-enrichment ratio.** The cause is a
merged header: row 5 labels column 6 `Ri value (t-test: p-value)`, and that one label
spans two cells — col 6 is the ratio, col 7 the parenthesised p. Full write-up is H11 in
`hurdles.md`.

Do all of the following:

1. Rename the column to **`ri_ttest_p_anabaena`** everywhere it appears — the notebook,
   `data/processed/all_orthologs_rbh.csv`, and anything derived from them.
2. **Kushige's Data Set S1 contains no rhythmicity p-value at all.** If a real
   `corr_p_anabaena` is wanted it has to be parsed from the SI PDF's Table S2, which
   lists only the 78 cycling genes — so it will be null for every other gene. Populate
   it that way or leave it absent; do not reconstruct it from anything in the
   spreadsheet.
3. Add a **post-parse assertion** that would have caught this, and apply the same idea
   to every other parsed column: *every gene the paper calls tier a must have
   rhythmicity p < 0.01, and every tier b < 0.05.* Fail loudly if not. Confirm the
   mislabelled column fails this and a correctly-parsed one passes.
4. Re-run the Version B headline numbers and confirm they are unchanged — the 565-pair
   set was cut on peak-time availability, not on p, so the Rayleigh result
   (p = 1.4e-7, R-bar = 0.167), the amplitude/|delta-phi| Spearman (rho = +0.014) and the
   amplitude-vs-amplitude result (rho = +0.056, n = 1787) should all reproduce exactly.
   **If any of them move, that's a stop-and-check signal** — it would mean the blast radius is larger
   than assessed.

---

# PART A — amplitude correlation by functional category

## Task

Answer two of the original questions with a traceable, reproducible artifact:

1. **"What is the correlation between amplitude of rhythmic genes in 7942 and amplitude
   of homologous Anabaena genes?"** — currently answered only as a single number
   (Spearman rho = +0.056, n = 1787) computed ad hoc and not committed anywhere.
   Make it a script, a figure, and a CSV.
2. **"If little correlation on a gene-by-gene basis, is there correlation on a pathway
   level?"** - this was conditoned on (1) being weak. It is weak. So this is
   now a required analysis, not a stretch goal.

Write `src/08_amplitude_by_function.py`. Follow the conventions in the existing `src/`
scripts and in `src/07_fetch_vijayan_markson.py` (module docstring explaining WHY, not
just what; every non-obvious decision commented with its reason).

## Inputs — all already on disk, do not download anything

| File | What to take from it |
|---|---|
| `data/processed/rbh_orthologs.csv` | 1,846 reciprocal-best-BLAST ortholog pairs. Keys: `locus_tag_7942`, `locus_tag_anabaena` |
| `Ito2009/Ito2009_sd1.xls`, sheet `Table S1`, `header=5` | 7942 side. Columns: `7942ID`, `Gene Name`, `Annotation`, `Amplitude*`, `Peak time*`, `Correlation p-value*`, `Cycling property**` |
| `Kushige2013/Kushige2013_so2.xlsx`, sheet `Sheet1`, `header=5` | Anabaena side. Columns: `ORF No.`, `Gene name`, `Annotation`, **`Category`**, `Cycling property (N+, n=2)*`, and per-timepoint `N+_1st_LL{t}` / `N+_2nd_LL{t}` for t in 4..48 step 4 |

**`Category` is the key column and has not been used anywhere in this project yet.**
It is Kushige's own curated functional assignment, 16 categories, non-null for all
5,336 genes — the scheme behind their Figure S1. 1,810 of the 1,846 RBH pairs carry
one; 1,105 are non-hypothetical.

## Method — specific requirements

**Amplitude.** Ito publishes `Amplitude*` directly. Kushige does not publish an
amplitude, so compute CV = SD/mean across all 24 N+ columns (both replicate series
pooled), exactly as `VersionB_AllOrthologs.ipynb` §3 does. Reuse that code path rather
than reimplementing it; if the two disagree, that needs resolving before moving on.

**Rank within each platform's own FULL distribution before comparing — not within the
ortholog subset, and never compare raw CVs.** Kushige's raw CVs run ~2.7x higher than
Ito's across the board (median 0.341 vs 0.124) because they are different array
platforms with different noise floors. This is logged as hurdle H3 in `hurdles.md`.
Ranking against the subset instead of the full distribution silently changes what a
percentile means.

**Use all pairs with an amplitude on both sides (expect n ≈ 1,787), not the 565.**
The 565 figure comes from requiring a defined *peak time*, which Ito only reports for
genes it calls cycling. Amplitude has no such restriction, so this analysis gets the
larger and less filtered set. Report both n's so the difference is visible.

**Assign categories from the Anabaena side** (that is where `Category` lives) and say so
explicitly in the figure caption and the CSV header — the pair inherits one organism's
annotation, which is a real limitation, not a detail to smooth over.

## Statistics — read this whole section before writing any test

The default framing of this analysis invites a false positive: 14 categories at
alpha = 0.05 produces a spurious hit roughly half the time. The structure below is
deliberately conservative and stays fixed while the analysis runs — no renegotiating it partway through.
It is organised as three tiers, and each tier is gated on the one above it.

### Tier 0 — bound the global effect. Report this even if nothing else runs.

Do not report the global result as "significant but small". Report it as a bound:

- Spearman rho with a **95% CI via Fisher z**. Expected: rho = +0.056, CI
  [+0.010, +0.102], n = 1787. Reproduce these before continuing.
- State the upper bound in plain language — at rho = 0.10, shared function explains
  about **1% of amplitude variance**. Phrase the conclusion as *exclusion*
  ("we can rule out any correlation above ~0.10"), not as absence of evidence.
  This is an equivalence-style claim: more conservative than a significance test AND
  more informative, because it answers "how big could it be" rather than "is it zero".
- **Robustness to non-independence:** neighbouring genes have correlated amplitude
  (lag-1 rho approximately +0.20 in both organisms), so genes are not independent
  observations. Check it with a moving-block bootstrap over genomic position (sort by
  the numeric part of `locus_tag_7942`; block sizes 1, 5, 10, 20; 3,000 resamples).
  **This has already been run and the CI is stable** — [+0.008, +0.103] at block 1
  versus [+0.013, +0.102] at block 20, because the ortholog join breaks the local
  correlation structure. Reproduce it and report it as a passed check. If the numbers
  disagree materially, track down why before continuing.

### Tier 1 — ONE omnibus test. This is the gate.

Ask a single question: does functional category explain any of the amplitude
relationship at all?

- Permutation test: shuffle `category` labels across pairs 10,000 times, recomputing
  the per-category rho each time. Test statistic = the **variance of per-category rho**
  (or equivalently max |rho|; report both, pre-declare variance as primary).
- One p-value. That is the whole inferential claim about pathway structure.
- **If this test is null, make no per-category claims whatsoever.** Report the Tier 2
  estimates as description only and say explicitly in the output that the omnibus test
  did not license per-category inference. Do not soften this into "suggestive trends".

### Tier 2 — estimation, not testing.

For every category with **n >= 40** (not 20 — at n = 24 the CI on rho spans roughly
[-0.40, +0.40] and cannot support a conclusion in either direction; state the dropped
categories and their n):

- Report rho with a **bias-corrected bootstrap 95% CI**.
- **Report no p-values here.** Effect size plus interval is more informative than a
  significance verdict and does not accumulate multiplicity.
- Sort by effect size, plot with a reference line at zero, and overlay the 2.5-97.5
  percentile band from the Tier 1 permuted null so a reader can see what noise looks
  like on this dataset.

### Tier 3 — two pre-specified confirmatory hypotheses, kept separate.

`PROJECT_BRIEF_rhythmic_gene_comparison.md` §4 states a prior in advance:
phycobilisome genes (apc/pec/cpc) rhythmic in Anabaena but not 7942, and kai genes the
reverse. These are genuinely pre-specified, so they are confirmatory:

- Test these two and only these two at alpha = 0.05, **Holm-corrected across the two**.
- Use Holm rather than Benjamini-Hochberg anywhere p-values appear in this analysis.
  BH controls the expected proportion of false discoveries; Holm controls the
  probability of making any. The goal here is not to claim something false.
- Label them CONFIRMATORY in the output and everything in Tier 2 EXPLORATORY. Never
  merge the two sets in a table or a figure.

### Fixed before running, not adjustable after

- Drop `Hypothetical`, `Other categories` and `-` from all tested sets — they are not
  coherent functional groups. Keep them in the scatter as grey background points so the
  reader sees how much of the genome is unclassified (741 + 145 pairs, the largest bins).
- No dropping, merging or renaming categories after seeing any correlation.
- If I find myself wanting to change a threshold mid-analysis, that's a sign to stop and reconsider rather than push through.

## Outputs

1. `data/processed/amplitude_by_category.csv` — one row per ortholog pair:
   `locus_tag_7942, locus_tag_anabaena, gene_symbol_7942, gene_symbol_anabaena,
   annotation_7942, annotation_anabaena, category, amplitude_cv_7942,
   amplitude_cv_anabaena, amp_pct_7942, amp_pct_anabaena, cycling_7942, cycling_anabaena`

2. `data/processed/amplitude_category_stats.csv` — one row per tested category:
   `category, n, spearman_rho, p_value, q_value_bh, median_amp_pct_7942,
   median_amp_pct_anabaena`

3. `figures/07_amplitude_scatter_by_category.png` — percentile vs percentile scatter,
   coloured by category, excluded categories in grey. Include the global rho, p and n
   in the title. Do not add a regression line to the pooled data; at rho = 0.06 it
   implies a relationship that isn't there. Faceting into small multiples per category
   is likely to read better than one crowded panel — use judgement here, and keep both versions
   if it is close.

4. `figures/08_category_effect_sizes.png` — per-category rho with confidence intervals,
   sorted, with a reference line at zero. This is the actual answer to the pathway
   question and should be readable on its own.

## Also do this — the specific interests that don't fit the scheme

We care about nitrogen fixation and photoprotection, and neither is a top-level
category here (*nif* genes land under Energy metabolism or Other). Add a keyword-derived
boolean overlay from the annotation strings, kept **separate** from the main category
analysis and clearly labelled as ad hoc:

- nitrogen fixation: `nif`, `fdxH`, `hetR`, `hetN`, `nitrogenase`
- photosynthesis / phycobilisome: `psa`, `psb`, `apc`, `cpc`, `pec`, `phycobili`
- photoprotection / UV / oxidative stress: `uvr`, `phr`, `sod`, `kat`, `ocp`, `photolyase`

Report n and median amplitude percentile per organism for each flag. Note in the code
that ~36% of annotation strings are hypothetical/unknown, so recall on these flags is
poor and absence of a flag does not mean absence of the function.

## Validation — do not skip

- Assert the ortholog table has no duplicated `locus_tag_7942` or `locus_tag_anabaena`.
- Recompute the global Spearman and confirm it reproduces **rho = +0.056, p = 0.018,
  n = 1787**. If it does not, track down what differs before going further.
- Confirm the 200 highest-amplitude 7942 genes have a **median Anabaena amplitude
  percentile near the 46th** (i.e. no enrichment). This is the sanity check that the
  join is not scrambled.
- Spot-check three genes against `VersionB_AllOrthologs.ipynb` §7: `kaiB`
  (Synpcc7942_1217 <-> alr2885, 98th/27th percentile), `gap1` (Synpcc7942_0245 <->
  all2566, 100th/87th), `narM` (Synpcc7942_0933 <-> alr0614, 44th/92nd).

Do not touch `VersionB_AllOrthologs.ipynb`. New script, new outputs.

---

# PART B — fetch the better 7942 dataset, redo the comparison with one caller

## Why

The whole 7942 side currently rests on Ito 2009, which publishes only summary indices —
no per-timepoint data — and gives no peak time at all for the 1,715 of 2,515 genes it
does not call cycling. That reimposes a rhythmicity filter on one side of the
comparison, from what got published rather than from any choice made here
(`hurdles.md`, H8). The dataset decision is recorded in `runningNotes.md`.

## B1 — run the fetch script

`src/07_fetch_vijayan_markson.py` is written and its `cosinor()` and
`verify_sample_map()` functions are unit-tested against synthetic data and against the
real GEO sample titles. **Everything else in it — the fetch, parse, join and reporting
path — has never been executed**, because the machine it was written on had no NCBI
access. Expect to fix things. Run:

```
python src/07_fetch_vijayan_markson.py
```

It fetches GSE18902 (Vijayan 2009, 16 timepoints every 4 h across 60 h of constant
light), GSE52486 (Markson 2013 WT arm — GEO's own summary calls it a replicate of
GSE18902, 7 timepoints) and the GPL9534 platform table, verifies the sample maps
against GEO's titles, fits the cosinor, and prints a replicate-reliability estimate.

Note what it prints, in particular how many of the 1,846 RBH pairs gain a peak time
compared with the 565 under the Ito join.

## B2 — the amplitude semantics trap. Do not skip this.

GSE18902 and GSE52486 are **two-colour** arrays: each timepoint is hybridised against a
pool of equal-mass RNA from all timepoints, so a value is log2(sample / time-averaged
reference), **not an intensity**. CV is meaningless on these.

So do not carry Kushige's CV across. Instead:

- log2-transform Kushige's N+ intensities,
- run the **same** `cosinor()` on them (12 distinct timepoints, two replicate series),
- compare **fitted log2 amplitudes** on both sides.

This is the first time the two organisms will have had a genuinely commensurable
amplitude measure. It is also H3 reappearing in a new costume — write that down.

## B3 — one caller, both sides

Extend into `src/09_recompare_with_vijayan.py`. Fit the identical cosinor to both
organisms from raw values and derive amplitude, peak time and a p-value yourself rather
than inheriting either paper's published calls. Keep the papers' own calls as columns
for comparison, never as filters.

Carry forward the H6 caveat already logged: the cosinor p-value assumes independent
residuals, and 4 h-sampled series are autocorrelated, so effective df is below the
timepoint count and p runs anticonservative. Use it for ranking, not as a calibrated FDR.

## B4 — redo the two comparisons

- **Amplitude vs amplitude**, the Part A analysis, on the new fitted amplitudes.
- **Circular phase offset**, repeating `VersionB_AllOrthologs.ipynb` §4-6: delta-phi
  wrapped to +/-12 h, 4 h bins, Rayleigh test, amplitude overlay.

## B5 — H4, which is now executable

Logged in `hurdles.md` as blocked by data. It is not — Vijayan and Markson-WT are
independent biological replicates on identical hardware, and Kushige is already n=2. So:

- Estimate reliability from the Vijayan-vs-Markson gene-wise correlation.
- Report the attenuation factor and the disattenuated amplitude correlation alongside
  the raw one. Present both; never replace the raw number with the corrected one.
- State plainly that disattenuation raises an estimate and does not make it more
  certain — the CI widens.

## B6 — the A/B comparison. This is a deliverable, not a by-product.

Write `data/processed/ab_comparison.csv` and a figure, structured exactly like this:

| Quantity | A: Ito 2009, published calls | B: Vijayan 2009, recomputed | What the difference is |
|---|---|---|---|
| Amplitude rho (95% CI) | +0.056 [+0.010, +0.102] | ? | |
| n, amplitude pairs | 1,787 | ? | |
| n, pairs with peak time both sides | 565 | ? | ← expect the biggest change |
| Phase distribution shape | unimodal, long tail | ? | |
| Rayleigh p, R-bar | 1.4e-7, 0.167 | ? | |
| Genes with no peak time on the 7942 side | 1,715 of 2,515 | ? | |
| Amplitude measure | CV, published | fitted log2 amplitude | not comparable — say so |
| Rhythmicity call | each paper's own | one caller, both sides | |

Rules for it:

- **Every row must be labelled with which quantities are and are not comparable.** The
  amplitude row in particular compares a published CV against a fitted log2 amplitude —
  those are different quantities, and the row exists to show the change in *conclusion*,
  not to equate the numbers. Say that in the table, not in a footnote.
- **If the two columns agree, that is a real and reportable result** — it would mean the
  Ito-based conclusions were robust despite the thin data, which is worth one clear line.
  Do not bury agreement because disagreement makes a better story.
- **Do not frame column A as a mistake.** Ito is Kushige's own reference 7 with
  deliberately matched index definitions; anchoring there was the defensible choice and
  the limitation is in what Ito published, not in the decision to use it. The framing is
  "here is what the matched published dataset could support, and here is what archived
  measurements add" — not "here is the wrong answer and the right one".
- Keep both columns reproducible from committed scripts, since this table is the thing
  most likely to be questioned in a meeting.

## B7 — sanity checks before trusting anything

- `kaiB` (Synpcc7942_1217) should still be high-amplitude in 7942 and middling in
  Anabaena. If it is not, something is wrong with the join.
- Vijayan's paper reports 1,748 of 2,724 ORFs oscillating at 22-26 h. Your cosinor at a
  fixed 24 h period will not reproduce that exactly — but if you get wildly fewer,
  suspect the parse rather than the biology.
- Peak times should be spread across the cycle, not clustered on 0 h or on bin edges.

---

# PART C — RpaA: questions 3 and 4

## What was actually asked, verbatim

> the highest amplitude genes in 7942 are direct binding targets of RpaA (see the
> Markson et al paper)
>
> can RpaA binding sites be detected informatically in Anabaena? any relationship to
> the rhythmic genes they report?

## C1 — the premise is now directly checkable, so check it

The highest-amplitude 7942 genes are direct RpaA targets. Markson's ChIP
target list makes that testable rather than assumed. Two things to report:

- Do 7942 amplitude ranks (recomputed in Part B, or Ito's if B has not run) differ
  between RpaA ChIP targets and non-targets? Effect size with CI, not just a test.
- One correction to carry, from `journal_club_kushige_2013.md`: Markson found ~110
  binding sites giving **134 target transcripts** (93 protein/tRNA, ~170 genes once
  operons are counted; 41 ncRNAs). Against ~856 cycling genes in 7942, **most of the
  cycling transcriptome is indirect**, relayed through an RpaA-driven sigma-factor
  cascade. So "the highest-amplitude genes are direct targets" is fine; "the cycling
  transcriptome is the RpaA regulon" is not. That distinction changes what a positive
  result in E3 would mean.

Also relevant: in 7942 the clock output is delivered by RpaA **phosphorylation**, not
abundance (Gutu & O'Shea 2013; SasA kinase, CikA phosphatase). So *rpaA* transcript
level is not the readout, and any transcript-level argument about RpaA itself is weak.

## C2 — clause one is already answered in the literature

Arbel-Goren et al. 2024 scanned the *Anabaena* genome with **FIMO (MEME Suite v5.4.1)**
using Markson's RpaA position-specific probability matrix, both strands, min match
P < 10⁻⁴, restricted to −500/+50 bp around TSSs from **Mitschke et al. 2011**. *Anabaena*
RpaA is **All0129**, ~96% aa similarity to the 7942 protein. They report **81 genes**
with putative sites at FIMO q < 0.05, including **two sites upstream of *pecB***
(P = 1.5 × 10⁻⁵) — the very gene Kushige built their reporter on.

So do not re-run the scan as though it were an open question. Reproducing it is optional;
obtaining their 81-gene list is required.

## C3 — clause two is the open question, and it is a specific test

Arbel-Goren crossed **motif hits × orthologs of the 7942 ChIP targets**. They never
crossed **motif hits × Kushige's 78 rhythmic genes**, and their treatment of the 81 is a
hand-picked list of regulators with **no enrichment analysis of any kind**.

The test still sitting unclaimed:

- **Hypergeometric overlap of the 81 motif-bearing genes with the 78 rhythmic genes.**
- **Background must be genes having both a Mitschke TSS annotation and a slot on the
  Ehira/Ohmori array — not all 5,336.** Getting this wrong inflates the enrichment, and
  it is the single most likely way to produce a wrong answer here. State the background
  size explicitly in the output.
- Under independence expect ~1–2 overlaps. Five or six would be a large odds ratio.
- Then ask whether any overlapping genes are **phase-concentrated**, as direct RpaA
  targets are in 7942, or scattered. Use the circular machinery from
  `VersionB_AllOrthologs.ipynb` §4–5. With n this small, report the phases individually
  rather than testing — a Rayleigh test on five genes is not meaningful.

## C4 — two objections that must appear in any output

These are not optional caveats; without them the analysis overstates itself.

1. **A FIMO scan with the RpaA matrix cannot distinguish RpaA from RpaB.** RpaB binds
   HLR1, a pair of imperfect 8-nt direct repeats (G/T)TTACA(T/A)(T/A) separated by two
   nucleotides. The HLR1 element in the *kaiBC* promoter **overlaps the reported RpaA
   site**, is bound by RpaB in vivo, and peaks ~12 h out of phase with RpaA (Hanaoka
   et al. 2012; Espinosa et al. 2015). Both are OmpR/PhoB-family; RpaB is essential and
   conserved in *Anabaena*. The honest label for the 81 is **"RpaA/RpaB-family sites."**
2. **GC-content mismatch.** *Anabaena* is ~41% GC, *S. elongatus* ~55%. A PWM trained on
   a GC-rich genome and applied to an AT-rich one inflates hit rates unless the
   background is an *Anabaena*-specific higher-order Markov model. FIMO's q-values fix
   multiple testing, not base-composition mismatch in the null. **The GC figures are
   marked "verify" in the journal club notes — confirm them from the assemblies in
   `data/raw/genomes/` before repeating them.**

## C5 — what has to be obtained first

None of this is on disk. Before any analysis:

- Arbel-Goren et al. 2024's **81-gene list** with FIMO q-values (supplementary).
- Markson et al. 2013's **ChIP target list** (~134 transcripts) and, if the scan is to be
  reproduced, the **RpaA PSPM**.
- **Mitschke et al. 2011** TSS annotations for *Anabaena*, needed for the background set.

Apply the Part B lesson before committing to any of them: check the experimental design
and what is actually tabulated, not just that a supplement exists. Log every source in
`data/raw/SOURCES.txt` with URL and retrieval date.

## C6 — scope discipline
Do not drift into the RpaA~P mechanism, SasA/CikA conservation,
or whether the transducer is functionally intact in *Anabaena* — nobody has shown
*Anabaena* SasA phosphorylates *Anabaena* RpaA, or measured RpaA~P there at all, and that
gap is a separate project. If the hypergeometric is unremarkable, that is a clean and
reportable answer to this question.
# PART D — heterocyst-enrichment layer

## What is actually available

Kushige's own heterocyst-enrichment experiment (an 80%-heterocyst fraction plus a
*hetR*-null comparison) is already **genome-wide in the spreadsheet you have**, not just
the 23-row Table S3 excerpt in the SI PDF. In `Kushige2013/Kushige2013_so2.xlsx`,
sheet `Sheet1`, `header=5`:

| Column | Content | Coverage |
|---|---|---|
| `Ri value (t-test: p-value)` | heterocyst-enrichment ratio | 5,263 / 5,336 |
| `Unnamed: 7` | its p-value, parenthesised — strip `()` then `to_numeric` | 5,264 |
| `Sum(N-)/Sum(N+)` | nitrogen-deprived vs replete total | 5,336 |
| `N-_1st_LL{t}`, t = 4..48 step 4 | full N- time series, **n = 1** | 5,336 |
| `Cycling property (N-, n=1)**` | Kushige's N- rhythmicity call | 5,336 |

Reference values already computed — reproduce them as a check: **359** genes cycling
(tier a or b) under N- versus 78 under N+, with **39** in both; **202** genes at
Ri > 2 and p < 0.05, of which only **17** cycle under N- and **2** under N+.

## D1 — parse and join

Add the columns to the tidy Anabaena table and join onto the existing RBH ortholog
table. No new downloads, no new orthology work — same locus tags throughout.

Compute a peak time for the N- series with the same cosinor used everywhere else, and
flag it clearly as n = 1: it gets no replicate-based reliability, so H4 stays blocked on
that specific series even after Part B.

## D2 — the confound to handle explicitly, not to discover later

**Dilution and measured amplitude are not independent.** Heterocysts are roughly 1 cell
in 10-15, so a transcript confined to heterocysts is diluted about 10-fold in a bulk
filament measurement. Its *measured* amplitude is suppressed by that dilution regardless
of how strongly it actually oscillates per cell.

So a negative relationship between `Ri` and amplitude is the **expected artifact**, not
a finding. Any analysis relating enrichment to rhythmicity must say so up front and
must not present such a correlation as biology. This is hurdle H9 and it is the closest
thing in this dataset to the mixed-community metagenomics problem — treat it as the
point of the exercise rather than an obstacle to it.

Where possible, compare *within* enrichment strata rather than across them, so dilution
is held roughly constant.

## D3 — comparisons worth making

1. **Amplitude and phase, N+ versus N-, for the same genes.** Both series come from the
   same platform and the same lab, so this is the cleanest within-organism comparison
   available anywhere in the project — no cross-platform amplitude problem at all.
   The 39 genes cycling in both conditions are the natural focus.
2. **Do heterocyst-enriched genes differ in phase** from the bulk population? Use the
   circular machinery from `VersionB_AllOrthologs.ipynb` §4-5. Phase is far less
   dilution-sensitive than amplitude, which makes it the safer axis here — say why.
3. **Where do the enriched genes sit in the Part A functional categories?** Expect
   nitrogen fixation and related categories to dominate; confirm rather than assume.

## D4 — figures

- `figures/09_Ri_vs_amplitude.png` — enrichment against N- amplitude percentile, with
  the dilution expectation annotated directly on the plot so it cannot be misread.
- `figures/10_phase_Nplus_vs_Nminus.png` — per-gene peak time under the two conditions,
  for genes cycling in both. Circular axes, and mark the identity line.
- `figures/11_heterocyst_category_composition.png` — functional composition of the
  Ri > 2, p < 0.05 set against the genome background.

## D5 — the honest framing

This layer will **not** give a second cross-species phase comparison. `Ri` is a spatial
enrichment ratio from a physically separated fraction: it says *where* a transcript is,
not *when* it oscillates in that cell type. What it does give is a worked demonstration
of signal dilution in a mixed population — which is the hurdle that transfers most
severely to the lab's metagenomic work, where every sample is a mixed community and an
abundant member's pattern can read as the whole community's.

Write the outputs so that framing is unmissable, in the figure captions as well as the
notes. Do not let it get presented as a phase result.

---

# PART E — the *kai* paralog question that will be asked

## E1 — the answer, and the two follow-ups

We want to know about the long-form *kaiB*. The answer is that it is not rhythmic, but the
two follow-ups that will be asked next are the interesting part, and both are already
answerable from data on disk. Build `src/11_kai_paralogs.py` to produce it properly
rather than from a notebook cell.

Ground truth, verified — reproduce these exactly as a check:

| Gene | ORF | Call N+ | Call N− | Raw CV (pct) | Ri (pct) | Sum(N−)/Sum(N+) |
|---|---|---|---|---|---|---|
| *kaiA* | alr2884 | AR | c | 0.193 (4th) | 1.092 (57th) | 1.67 |
| *kaiB* | alr2885 | **b** | c | 0.275 (27th) | 1.167 (68th) | 1.68 |
| *kaiC* | alr2886 | AR | AR | 0.396 (66th) | 1.141 (65th) | 1.04 |
| *kaiB* long | all3328 | AR | AR | 0.404 (68th) | **0.709 (6.7th)** | 1.11 |

From Table S2 in the SI PDF, canonical *kaiB* (`alr2885`) alone: rhythmicity p = 0.046,
amplitude 0.246, peak 21 h under N+; p = 0.065, amplitude 0.213, peak 16 h under N−.
`all3328` is absent from Table S2 because it is arrhythmic, so no published p exists —
say that explicitly rather than leaving a blank.

**The two pre-empting answers:**

- *Is the long form's flat call just heterocyst dilution?* **No.** Ri = 0.709, the 6.7th
  percentile — it is *depleted* in the heterocyst fraction, not enriched. A transcript
  confined to a ~1-in-10 cell type could look flat through dilution alone (H9), so this
  needed checking; it is ruled out.
- *Does it switch on under nitrogen deprivation?* **No.** Sum(N−)/Sum(N+) = 1.11 against
  1.68 for canonical *kaiB*, and it stays AR in N−.

## E2 — figures

- `figures/12_kai_locus_timeseries.png` — N+ and N− phase-folded profiles for all four
  genes on shared axes. Annotate each with its tier call. The point a reader should take
  is that *kaiC* and long *kaiB* have nearly identical, fairly high raw CVs (66th/68th
  percentile) and are both arrhythmic, while *kaiB* cycles at less than half their CV.
  That contrast inside one locus is the cleanest CV-without-a-test illustration in the
  project — make it the visual argument, not a caption.
- `figures/13_kai_in_context.png` — the four genes located on the genome-wide
  CV-versus-Ri plane, so their position relative to the whole transcriptome is visible.

## E3 — one thing genuinely worth checking

`all3328` has no 7942 ortholog at all — confirmed independently, not merely absent from
the RBH table. Verify that once more directly against the proteome (search for KaiB-like
hits by BLAST rather than trusting the RBH output, since RBH discards paralogs by
construction, which is exactly the case here). If 7942 really has no long-form KaiB, say
so as a positive finding with the method that establishes it; if a weak hit exists,
that materially changes the story and needs reporting.

---

---

# PART F — do the two organisms differ in overall amplitude?

## Why this is now askable, and why it wasn't before

Every amplitude comparison in this project so far has been **within-organism percentile
rank**. That was the deliberate H3 mitigation: raw CVs were not comparable across
platforms, so we ranked each gene inside its own distribution. The cost is that the
ranking is **blind to a global difference between organisms by construction**.

That matters right now because of the *kaiB* result. It sits at the 90th percentile of
amplitude in *Anabaena* under a cosinor fit, against the 27th under CV. It is tempting to
read that as overturning Kushige's "kai genes barely cycle in *Anabaena*" headline. **It
does not**, because "high relative to other *Anabaena* genes" and "*Anabaena* oscillates
less than 7942 across the board" are perfectly compatible statements. Independent support
for the second: Uzumaki et al. 2004 (*Nat Struct Mol Biol* 11:623, Kushige's ref 15) report
that *Anabaena kaiA* in *Synechococcus* gives a ~40 h period **and lower amplitude**.

Part B changed what is possible. A **fitted cosine amplitude in log2 units is a
fold-change** — dimensionless, and far more portable across platforms than CV ever was.
So the global comparison is now defensible in a way it was not. Not free, though: see F2.

## F1 — the comparison

Write `src/13_global_amplitude.py`. Compare the distribution of fitted log2 amplitude
between *Anabaena* (Kushige, log2 intensities) and 7942 (Vijayan, log2 ratios), at two
levels, and report both:

1. **Matched ortholog subset** — the ~1,811 RBH pairs. This is the apples-to-apples
   comparison: same genes, so gene-content differences cannot drive it. Paired, so use a
   paired test (Wilcoxon signed-rank on the per-pair difference) and report the median
   paired difference with a CI.
2. **Full gene populations** — all ~2,715 (7942) vs ~5,336 (*Anabaena*). This answers the
   population-level question but confounds amplitude with gene content, since *Anabaena*
   has roughly twice the genes and a large hypothetical fraction. Say that.

Report medians, quartiles and the full ECDFs, not just a test. The shape matters more than
the p-value — a uniform shift and a difference confined to the tail mean different things.

## F2 — the noise-floor problem. This is the analysis, not a caveat.

**A cosinor fit to pure noise returns a non-zero amplitude.** That inflation depends on
the number of timepoints and the noise level, and the two datasets differ in both:

| | *Anabaena* (Kushige) | 7942 (Vijayan) |
|---|---|---|
| Distinct timepoints | 12 | 16 |
| Replication | n = 2 | n = 1 |
| Values entering the fit | 24 | 16 |
| Data type | log2 intensity | log2 ratio vs pooled reference |

So a raw comparison of fitted amplitudes is **not** apples to apples — the dataset with
noisier, sparser data will show systematically higher fitted amplitude for reasons that
have nothing to do with biology. Do not report the raw comparison alone.

Build a per-dataset null and compare against it:

- For each dataset separately, **permute the time labels** within each gene and refit
  (≥1,000 permutations over a random sample of genes, or all genes if it runs). That gives
  the amplitude distribution expected from noise alone, *for that dataset's own sampling
  design and noise level*.
- Report the null median for each dataset. If they differ, that difference is the
  measurement artifact and it needs subtracting from, or dividing out of, the comparison.
- Then report the organism comparison **relative to each one's own null** — e.g. amplitude
  expressed as a z-score or ratio against that dataset's permuted distribution.
- Present raw and null-corrected side by side. If the sign of the difference flips between
  them, the raw comparison was measuring the assay, not the organism.

## F3 — close out the *kai* question

With absolute amplitudes in hand, put *kaiA* (alr2884), *kaiB* (alr2885), *kaiC* (alr2886)
and long-form *kaiB* (all3328) on the same axis as their 7942 counterparts. Give for each:
fitted log2 amplitude, within-organism percentile, and where it sits against that
dataset's noise null.

The specific question to answer in one sentence: **is *Anabaena kaiB*'s absolute amplitude
comparable to 7942's, or is it merely high within a globally quieter distribution?** Those
lead to opposite readings of Kushige's headline, and right now we cannot tell them apart.

## F4 — confounds to state explicitly, not to discover in a meeting

- **Ratio arrays compress fold-changes.** Cross-hybridisation and dye saturation pull
  two-colour log ratios toward zero, so 7942 amplitudes may be systematically understated
  relative to intensity-derived ones. This biases *against* finding *Anabaena* quieter, so
  if the result comes out that way anyway it is the more robust direction. Say so.
- **Different arrays, different dynamic range.** No way to correct this from the data we
  have; name it.
- **Kushige is n=2, Vijayan n=1.** Averaging replicates suppresses noise-driven amplitude,
  which pushes the opposite way from the point above. The two biases do not cancel in any
  principled way — do not claim they do.
- **Different growth conditions beyond LL.** Both entrained with 2 LD cycles then released
  into constant light, which is a good match; note anything else that differs.

## F5 — figures

- `figures/14_amplitude_distributions.png` — overlaid ECDFs for the two organisms, on the
  matched ortholog subset, with each dataset's permuted null shown as a dashed line.
- `figures/15_kai_absolute_amplitude.png` — the four *kai* genes and their orthologs on a
  shared absolute axis, with the two noise floors marked.

## F6 — what would actually change the *kaiB* claim

State the decision rule before running anything:

- If *Anabaena* amplitudes are **globally lower** and *kaiB*'s absolute amplitude is
  unremarkable against 7942's, then Kushige's headline stands and the 90th-percentile
  result is only about rank within a quiet genome.
- If the distributions are **comparable** and *kaiB* holds up in absolute terms, then the
  headline is at least partly an artifact of CV, and that is worth saying out loud.
- If the two datasets' **noise floors** differ enough to explain the gap, neither
  conclusion is available and the honest answer is that this cannot be settled with these
  two datasets.

Write down which of the three I land on. Don't let it stay ambiguous.
