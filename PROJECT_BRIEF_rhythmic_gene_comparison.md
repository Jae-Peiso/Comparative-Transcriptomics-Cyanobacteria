# Project Brief — Cross-species comparison of circadian rhythmic genes (Kushige 2013 vs. Ito 2009)

> Handoff document for continuing this work in VS Code. Specification and context
> primer, not code. Read top to bottom before writing anything.
>
> **SCOPE: this is a few-days task.** It has been deliberately scaled to fit that
> window. The finishable core is in §4–§6. Work that is scientifically better but
> does NOT fit the window is parked in §7 (Deferred) — parked, not deleted. Do not
> pull deferred items back in without renegotiating the timeline.

---

## 1. What this project actually is

A **methodological rehearsal**, not primarily a biology project. Mike asked for a
quantitative comparison of rhythmic genes in an *Anabaena* paper vs. a published
*S. elongatus* PCC 7942 dataset. His stated reason:

> "whatever data you collect will basically be like the data in this paper, so any
> analysis hurdles that come up will also apply to any new results."

**Confirmed with Mike: he wants the hurdles characterized, not a polished cross-species
result.** So the deliverable is:

1. **One honest, finishable quantitative comparison** (Version A, §4).
2. **A register of analysis hurdles** (`hurdles.md`, §6) — each either *demonstrated*
   by the work actually done, or *described* with "here's the trap and what you'd do
   about it." A described hurdle is a legitimate deliverable item; you do NOT have to
   execute every hurdle to characterize it.

The register is the real product. The Version A analysis is the worked example that
makes several of the hurdles concrete.

### Downstream application these hurdles must transfer to
The lab runs a **cyanobacteria metagenomics** project generating its own
rhythmic-expression data, comparing **across organisms** — both
**environment-vs-reference** AND **environment-vs-environment (MAG-to-MAG)**. This
matters for H2: the statistical hurdles (H1, H3, H5) transfer to both; the
orthology hurdle transfers cleanly only to the reference-anchored case. See §6/H2.

---

## 2. The two datasets

| | *Anabaena* sp. PCC 7120 | *S. elongatus* PCC 7942 |
|---|---|---|
| Paper | Kushige et al. 2013, *J Bacteriol* 195:1276–1284 | Ito et al. 2009, *PNAS* 106:14168–14173 |
| DOI | 10.1128/JB.02067-12 | 10.1073/pnas.0902587106 |
| PMCID | (JB, check article page) | PMC2729038 |
| Platform | Ehira & Ohmori oligo array, 5,336 of 5,368 ORFs | high-density oligo array, ~2,515 genes |
| "Rhythmic" set | **78 genes** (P<0.05 AND CV>10^-0.7) | **97 higher-amplitude genes** (87 sDusk + 10 sDawn) |
| Index defs | amplitude = CV (SD/mean); cosine-corr P; peak time | **same definitions** — Kushige adopted Ito's, deliberately |

**Anchor on Ito, not Markson.** Ito is Kushige's own ref 7 with identical index
definitions, so the two sides are methodologically matched for free. The "856 genes"
figure traces to Markson 2013 (different measurement, different thresholds) and drags
in the RpaA regulon questions — explicitly out of scope. The 7942 rhythmic set at
Ito's working threshold is **97**, not 856. Do not let 856 leak in.

---

## 3. Data-availability gate (TIMEBOXED — half a day, hard stop)

The two SIs will almost certainly need **manual download** (every automated route
reCAPTCHA/bot-blocked during scoping). A step-by-step manual procedure lives in
`DOWNLOADS.md` (see §8). Half a day, then decide — do not exceed it.

- [ ] **Ito 2009 SI** (PMC2729038 / PNAS DCSupplemental): does it tabulate per-ORF
  amplitude + peak time, or only the 97 rhythmic genes? Version A needs *at least*
  the 97 with their indices — almost certainly present. The full ~2,515 table is only
  needed for the deferred Version B.
- [ ] **Kushige 2013 Data Set S1** (*J Bacteriol* article page): full ~5,336-ORF table,
  or only the 78? Version A needs *at least* the 78 with indices. Raw arrays are in
  the **dead KEGG Expression Database** (ex0001892–ex001947) — do NOT rely on it.
- [ ] **Genome/proteome inputs** for orthologs (see H2 fallback before committing).

**Hard rule for this timeline:** if the full per-ORF tables are not sitting there in
usable form, **Version A only**, and "the unbiased full-transcriptome data was not
archived" becomes a one-line written finding in `hurdles.md` (that IS one of the
recurring lessons — see H7). No spelunking into dead databases on a few-days budget.

---

## 4. The analysis — Version A ONLY

Map orthologs among the genes called rhythmic in at least one dataset, then:

- **Amplitude:** compare in **ranks** (Spearman), NOT raw CV — two different array
  platforms, absolute CVs not comparable (H3).
- **Phase:** compare peak times as a **dawn/dusk 2×2 contingency table** (Fisher +
  odds ratio). Do NOT run Pearson on peak times — phase is circular (H5).

That is a real, defensible, finishable result and it maps onto exactly what Mike
pictured.

**Version A is deliberately the *biased* version** (it conditions on rhythmicity — H1).
That is a feature here: doing it and noting what an unbiased analysis would fix is how
you demonstrate the selection-bias trap without building the unbiased pipeline.

### Expected qualitative result (state as prior, then test)
Gene-level correlation likely weak — the paper itself shows the two most salient gene
sets going *opposite* directions (kai: top-amplitude in 7942, flat in *Anabaena*;
phycobilisome apc/pec/cpc: rhythmic in *Anabaena*, not in 7942).

---

## 5. Task decomposition (scaled to the window)

```
project/
  DOWNLOADS.md        # manual SI-download procedure (§8) — do this first
  data/
    raw/              # downloaded SIs, proteomes — NEVER edit in place
    interim/          # parsed per-gene tables, one tidy CSV per dataset
    processed/        # ortholog-joined table
  src/
    01_parse_kushige.py  # SI -> tidy CSV: gene_id, gene_symbol, amplitude_cv, corr_p, peak_time_h, rhythmic_flag
    02_parse_ito.py      # SI -> tidy CSV, same schema
    03_orthologs.py      # RBH between proteomes — OR name-match fallback (H2)
    04_join.py           # ortholog/name-keyed merge -> processed table
    05_stats.py          # Spearman on amplitude ranks; dawn/dusk 2x2 on phase
    06_figures.py        # amplitude scatter (ranked); phase contingency viz
  hurdles.md          # THE DELIVERABLE — living log, append as things bite
  README.md
```

- **Tidy schema, identical columns both datasets:**
  `gene_id, gene_symbol, amplitude_cv, corr_p, peak_time_h, rhythmic_flag`.
  Making the columns line up forces the definitional questions into the open.
- **Persist every intermediate as CSV** (see H7 — provenance is itself rehearsed).
- **Keep `hurdles.md` open while working**; append the moment something bites. Do not
  reconstruct it from memory at the end — the friction IS the data.

---

## 6. Hurdle register (the spine of the deliverable)

Each entry marked **[DEMONSTRATED]** (Version A makes it concrete) or **[DESCRIBED]**
(named + mitigation sketched, not executed on this timeline).

### H1 — Rhythmicity-threshold selection bias  [DEMONSTRATED]
Conditioning a correlation on "rhythmic in dataset X" is selection on a noisy variable
-> range restriction + regression to the mean -> correlation biased down, plus a
spurious "*Anabaena* is less rhythmic" story. **Version A is the biased version by
construction**; you demonstrate the trap by running it and noting that an
all-orthologs analysis (deferred Version B, §7) is what removes it. Recurs identically
the first time the lab compares rhythmic genes across its own conditions. HIGHEST
transfer value.

### H2 — Orthology mapping  [DEMONSTRATED, easy instance — with a split-transfer warning]
Here: 5,368 vs. ~2,700 genes, paralog families in *Anabaena*. **Reciprocal best BLAST
hits (RBH)** is the default: clean, 1:1, discards paralogs (state as caveat, don't
solve). This is orthology's most forgiving instance — two finished, curated genomes,
stable locus tags.

**FALLBACK IF RBH BECOMES A TIME SINK (protect the timeline):** both papers print real
gene symbols for many of the 78 and 97 (kaiABC, apcE, pecB, psbA, groEL, etc.).
Name-matching a few dozen rhythmic genes by symbol is *finishable in an afternoon*; a
full proteome RBH with ID reconciliation (array probe -> locus tag -> protein
accession, on two era-specific annotation schemes) is the classic multi-day swamp. For
Version A you only need orthologs *among the rhythmic sets*, so name-matching is a
legitimate primary method here, not just a fallback. Use RBH only if it's already
familiar and quick.

**Split transfer to the two downstream cases:**
- **Env-vs-reference — WELL rehearsed.** One curated anchor; ambiguity one-sided.
  Transfers almost directly.
- **Env-vs-env (MAG-to-MAG) — only PARTIALLY rehearsed, strictly harder.** No ground
  truth on either side; a missing hit = genuine absence OR unassembled region OR
  split/fused gene, on *both* sides. MAG incompleteness is non-random (GC/coverage
  skewed) -> "gene absent" correlates with genome features -> with gene class -> with
  rhythmicity class = a confound shaped like biology, invisible to RBH. Two finished
  genomes do not exercise this. **Do not read a clean Kushige-vs-Ito result as evidence
  the method survives two MAGs** — that case needs its own de-risking
  (completeness-aware ortholog calling, absence-as-missing-not-zero, checking
  presence/absence against coverage & GC before interpreting).

### H3 — CV as a mean-dependent, cross-platform amplitude measure  [DEMONSTRATED]
"Amplitude" = CV = SD/mean, inflates for low-expression genes; two different arrays,
so absolute CVs not comparable. **Mitigation used in Version A: compare in ranks.**
Recurs, likely worse, on coverage-dependent metagenomic data.

### H5 — Phase is circular  [DEMONSTRATED]
Peak time wraps at ~24 h. Do NOT Pearson peak times. **Version A uses a dawn/dusk 2x2
contingency table** (Fisher + odds ratio) — more robust and interpretable than a
circular correlation coefficient given 4-h sampling and bimodal phase. (Fuller circular
stats — Rayleigh, mean resultant length, ccc — are in Schmelling 2017 if wanted.)
Recurs identically.

### H4 — Attenuation of correlation by measurement noise  [DESCRIBED — not executed]
Both amplitudes noisy -> any r attenuated. In principle: estimate reliability from
replicate–replicate correlation (both studies n=2) and disattenuate. **Not done on this
timeline** — requires per-replicate values that may not be in the SIs, and it's a
mini-project. Named here as the correction a rigorous Version B would need, pending
replicate-level data. For Mike's purposes, naming it is most of the value.

### H6 — FDR / rhythmicity-call calibration  [DESCRIBED — not executed]
Both papers' rhythmicity = cosine-correlation P on a short, autocorrelated series ->
effective df below timepoint count -> P anticonservative; amplitude thresholds
stringent; gene counts (78, 97) are loose-test x tight-filter products. **Use the
papers' published calls as-is.** Re-calling with JTK_CYCLE / RAIN is a project — parked.
Named as the caller-choice decision the lab faces for its own data.

### H7 — Data provenance / recoverability  [DEMONSTRATED — it already bit]
Manifested during scoping: two of three primary sources unverifiable from an automated
context; Kushige's raw arrays in a dead database (KEGG Expression DB). Lesson for the
lab's pipeline: archive every intermediate per-gene table yourself, durable and
machine-readable — journal SIs and project databases are not reliable five-year stores.
The §3 gate and the "was-it-archived" finding are this hurdle in action.

---

## 7. DEFERRED — better science, does NOT fit a few-days window

Parked, not deleted. Pitch to Mike as a phase two.

- **Version B (all-orthologs, unbiased comparison).** The scientifically correct
  analysis and the fix for H1 — but needs (a) full per-ORF tables from BOTH SIs
  (unconfirmed, may not exist), and (b) whole-transcriptome ortholog mapping, not just
  among the rhythmic sets. If the §3 gate shows the data isn't archived, Version B is
  foreclosed and that becomes a written finding rather than a task.
- **Attenuation correction (H4).** Needs replicate-level data; mini-project.
- **Re-calling rhythmicity with JTK/RAIN (H6).**
- **Pathway-level comparison (Q2).** Aggregate by functional category and re-test —
  only if Version A finishes with time to spare. Stretch goal, not a commitment.

---

## 8. Manual download procedure -> put in DOWNLOADS.md

See the companion file `DOWNLOADS.md` (provided alongside this brief) for the
step-by-step. Summary of what it covers:
- Ito 2009 SI via the PNAS article page and the PMC2729038 supplementary section.
- Kushige 2013 Data Set S1 via the *J Bacteriol* / ASM article page.
- What to check inside each file the moment it opens (full table vs. rhythmic subset).
- Institutional-access notes and fallbacks.

---

## 9. Explicitly OUT of scope
- RpaA/RpaB regulon, motif scanning, the 81-gene Arbel-Goren list, the 78x81
  hypergeometric. Different project (Q3/Q4 layer).
- Any biological claim about *why* the programs differ. This characterizes whether/how
  much + the method hurdles.
- Proteomics. Everything is transcript-level.

---

## 10. First moves
1. `DOWNLOADS.md` -> grab both SIs (timeboxed, §3). Check full-table-vs-subset on open.
2. Parse both into the identical tidy schema (§5). Definitional questions surface here
   — start `hurdles.md` now.
3. Orthologs among the rhythmic sets — **name-match first** (H2), RBH only if quick.
4. Version A stats (§4): Spearman on amplitude ranks, dawn/dusk 2x2 on phase.
5. Figures + finalize `hurdles.md`.
