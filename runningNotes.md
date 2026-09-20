# Running Notes — Kushige 2013 vs. Ito 2009 rhythmic-gene comparison

> Scratch space for things that aren't methodology hurdles (that's `hurdles.md`, the
> actual deliverable) but are worth keeping around: open questions, side findings,
> things to follow up on, build notes for next time. Moved out of `hurdles.md` PART 3
> on 2026-07-27 to keep that file focused on the hurdle register itself.

---

## DECISION LOG

### 2026-07-31 — Anabaena-side dataset: searched for a Vijayan-equivalent, found none

Prompted by Part F (`13_global_amplitude.ipynb`, F2): Kushige 2013's noise floor is
markedly higher than Vijayan's (12 timepoints, n=2 replication, single-channel log2
intensity, vs. 16 timepoints for 7942) -- exactly the kind of dataset-quality gap that
got Ito swapped for Vijayan on the 7942 side (see the 2026-07-29 entry above). Checked
whether an equivalent upgrade exists for Anabaena.

**Searched for:** genome-wide (not single-gene), free-running/constant-light (not a
stimulus-response design), time-course (not a single timepoint or mutant-vs-WT
comparison) transcriptomics for Anabaena/Nostoc sp. PCC 7120, any technology, published
after 2013.

| Candidate | Type | Why it doesn't work |
|---|---|---|
| Kushige et al. 2013 (already in use) | Microarray, LL time course | -- |
| Zhu et al. 2021 (eLife) + 2023 follow-up (bioRxiv) | Single-cell fluorescence imaging (*pecB*-GFP, *hetR*-GFP) + RT-qPCR (*kaiABC*, *rpaA*, *pecB*) | Best circadian work on Anabaena since Kushige, but a handful of genes with single-cell/spatial resolution, not genome-wide -- a different data type, not a substitute. |
| GSE26633 (2011) | RNA-seq | Real RNA-seq, but a nitrogen-deprivation time course (0/6/12/21h post-removal) -- no free-running cycling design to fit a cosinor to. |
| PRJEB24960 | RNA-seq | WT vs. *hetZ*/*hetP* mutant comparison, not a time series. |

**Verdict: no substitute exists.** Kushige 2013 is still the only genome-wide,
free-running circadian dataset for Anabaena PCC 7120 -- a real ceiling on the Anabaena
side of every comparison in this project, not a gap in the search. Logged here so this
doesn't get silently re-searched later without knowing it was already checked.

**Transfers?:** The asymmetry itself is worth remembering going forward, not just the
negative result. This project's cross-organism comparisons are structurally lopsided in
data quality -- 7942 has a modern, better-designed comparator available (Vijayan);
Anabaena does not, and nothing published since 2013 fills that role. Any future finding
that differs between the two organisms should be checked against "is this biology, or
is this the Anabaena side running on the only (and noisier) data that exists" before
being reported as a real difference -- the same discipline Part F's noise-floor
correction already applies, but with an added wrinkle: that noise floor is set by the
*only* dataset available, not the *best* one, so there is no independent replicate to
sanity-check whether Kushige's 2013 noise characteristics are typical or unusual for a
well-run Anabaena circadian experiment.

### 2026-07-29 — 7942-side dataset: Vijayan primary, Markson replicate, Ito retained

Decided with Jae. Replacing the Ito 2009 side, keeping Ito as a secondary comparator,
and re-calling rhythmicity with one identical caller on both organisms.

| | Design | Timepoints | Verdict |
|---|---|---|---|
| **GSE18902** Vijayan, Zuzow & O'Shea 2009 | 2× LD → **free-running LL** | **16, every 4 h, T=24→84 h** (60 h, 2.5 cycles) | **PRIMARY** |
| **GSE52486** Markson 2013, WT arm only | same protocol, **stated replicate of GSE18902** | 7, every 4 h, T=36→64 h (T=52 dropped by submitters) | **REPLICATE — unblocks H4** |
| GSE59112 Pattanayak/Phong/Rust 2014 | 2× LD → LL, WT + *cikA*-null | 6, every 4 h, 24 h only | Too thin; *cikA* arm is scope creep |
| GSE104204 Piechura 2017 | simulated natural light, diel + pulses | 8 over 12 h since dawn | **DISQUALIFIED** — under one cycle (see hurdles H10) |
| GSE50908/50919/50920 Markson Δ*rpaA* arms | LL, Δ*rpaA* | — | Out of scope per brief §9 |
| Yuan 2024 PNAS compendium | 300 mixed samples, ICA/iModulons | — | Not a controlled time course |

**Why Vijayan wins on the merits, not just on availability:**

- Sampling interval is **identical to Kushige's** (4 h), duration is longer (60 h vs 48 h).
- **GPL9534's `ORF` column is already `Synpcc7942_XXXX`** — the exact key in
  `data/processed/rbh_orthologs.csv`. No probe→protein reconciliation at all. This is the
  single biggest cost difference between candidates and it is invisible from the paper.
- Published result: **1,748 of 2,724 ORFs oscillate** at 22–26 h, against Ito's 97 at his
  working threshold. Every gene gets a computed peak time, so H8's one-sided filter
  disappears rather than being mitigated.
- Processed per-gene values sit in the GEO Series Matrix — no image processing.

**Known cost:** we lose the "Kushige deliberately adopted Ito's index definitions"
matching that the brief anchored on (§2). That matching bought commensurable *definitions*
but cost us the *data*. Since Kushige publishes all raw timepoints too, running one caller
over both raw series recovers commensurability on a better footing — and collapses H3, H6
and H8 together, because we stop inheriting two papers' analysis choices and make one
choice ourselves, applied symmetrically. Ito is retained as the "published calls"
comparator; the gap between the two answers is itself a hurdle demonstration.

**Amplitude semantics — do not skip this.** GSE18902/GSE52486 are **two-colour** arrays:
each timepoint is hybridised against a pool of equal-mass RNA from all timepoints, so a
value is log2(sample / time-averaged reference), not an intensity. **CV is meaningless on
these.** The commensurable quantity is the fitted cosine amplitude in log2 units, which
means Kushige's N+ intensities must be log2-transformed and put through the same cosinor
rather than carrying their published CV across. This is H3 resurfacing in a new costume.

**Still no Ito raw arrays anywhere.** No GEO or ArrayExpress deposit found. Ito's
Affymetrix array was also designed against the **PCC 6301** genome rather than 7942, which
adds a probe→locus-tag step and may partly explain why only derived indices were ever
published. H7 stands.

**Next actions:** DONE — 2026-07-30, `src/09_recompare_with_vijayan.py`. All three
completed: (1) `07_fetch_vijayan_markson.py` ran cleanly (NCBI access confirmed working
this session — the earlier "sandbox has no outbound access" note was specific to a prior
environment, not permanent); every RBH pair now has a peak time (1811/1846, vs 565 under
Ito). (2) Kushige N+ log2-transformed and fit with the identical cosinor, both
amplitude-vs-amplitude and Δφ/Rayleigh rerun on the new joint table — full results and
the old-vs-new comparison table in `hurdles.md` (2026-07-30, "PART B" entry). (3) kaiB
sanity-checked directly against Kushige's own published peak time (21.1h fit vs 21h
published) — confirmed correct, though its Anabaena amplitude percentile moved a lot
(27th under old CV -> 90th under the new commensurable measure), a real finding not a
bug; see the hurdles.md entry for why. gap1 also checked: still high in both (99th/99th
under the new fit, vs 100th/87th under the old CV) -- consistent, no concerns.

---

## Open questions / follow-ups

- ~~RpaA / SasA / CikA / RpaB: PARKED (out of scope)~~ — **RESOLVED 2026-07-30, Part C**
  (`12_rpaa_overlap.ipynb`, `hurdles.md` "PART C" entries). Scope did get renegotiated:
  Mike's Slack Q3 ("highest-amplitude 7942 genes are RpaA ChIP targets") CONFIRMED,
  p=1.2e-8. Q4a (informatic RpaB/HLR1 sites in Anabaena) already done by Arbel-Goren
  2024, reproduced here. Q4b (overlap with Kushige's rhythmic genes) tested fresh: NULL,
  p=0.19-0.30. A real HLR1/RpaB motif model was also built to partially resolve the
  RpaA/RpaB identity ambiguity flagged here — ~40% of Arbel-Goren's 81 RpaA-motif genes
  also carry a plausible RpaB site, including both genes driving Q4b's small overlap.
- ~~pecBAC / phycobilisome side of the brief's §4 prior — not specifically checked
  yet~~ — **RESOLVED 2026-07-30, Part A** (`08_amplitude_by_function.ipynb`). Confirmed:
  phycobilisome genes are significantly higher-amplitude in Anabaena than 7942 after
  Holm correction (p=0.0034, n=13). kai does not reach significance at n=3 but all three
  kai genes go the predicted direction — floored by minimum-achievable-p at that sample
  size, not by a wrong-direction result. The broader exploratory sweep across all 14
  Kushige categories found nothing beyond chance (Tier 1 omnibus p=0.712) — these two
  pre-specified hypotheses are the only real signal, not a general pathway-level effect.
- RBH build notes for next time: genome sources are `data/raw/genomes/` (RefSeq
  GCF_000012525.1 / GCF_000009705.1, provenance in that dir's `SOURCES.txt`); pipeline
  is inline in `VersionB_AllOrthologs.ipynb` §1-3, not a standalone `src/` script yet —
  worth extracting if this needs to run again on updated data.
- Part D (heterocyst-gene enrichment, brief's remaining unstarted section) — not
  touched yet as of 2026-07-30.
- Reannotation of "hypothetical protein" genes (28 rhythmic + whole-genome background,
  `14_reannotate_hypotheticals.ipynb`) — done 2026-07-30, separate side thread from
  Mike's four questions, feeds forward into Part A's grey bucket. 11 of 28 rhythmic
  hypotheticals still need a manual HHpred lookup (FASTA staged in
  `data/interim/reannotation/hhpred_queries/`); not yet run.

---

## Findings (not hurdles, but worth keeping)

### kai genes (kaiA/B/C, + long-form kaiB in Anabaena)

Original plan (pre-Version-B): pull per-gene amplitude+phase for kai genes from both
SIs, put side by side. Expected: flat/low-amp in Anabaena, high-amp kaiBC in 7942
(Kushige Fig. 1).

**2026-07-27 — confirmed, plus an extra finding.** kaiB in 7942 (`Synpcc7942_1217`,
via RBH pair with `alr2885`): 98th-percentile amplitude in 7942, only 27th-percentile
in Anabaena — matches the brief's prior directly, found via the RBH ortholog table
without a targeted lookup.

Separately (this note's original "+ kaiBL" flag) — Anabaena turns out to have a SECOND
gene annotated "kaiB" by Kushige2013 (`all3328`), standalone elsewhere in the genome,
not the operonic one (`alr2885`, next to `kaiA`/`kaiC`). No 7942 ortholog exists for it
at all (confirmed independently, not just absent from the RBH table). Confirmed via
protein length it's a real "long-form" KaiB paralog: 254 aa vs. the canonical 108 aa.

Compared rhythmicity between the two directly (`VersionB_AllOrthologs.ipynb` §8).

> **CORRECTED 2026-07-30.** An earlier version of this table gave "P-value 0.018 /
> 0.399". Those are **heterocyst-enrichment t-test p-values**, not rhythmicity
> p-values — see the H11 entry in `hurdles.md`. Real rhythmicity values below are from
> the SI PDF's Table S2. The qualitative conclusion is unchanged, because it rests on
> Kushige's tier calls, which were always right.

| | canonical (`alr2885`) | long-form (`all3328`) |
|---|---|---|
| Length | 108 aa | 254 aa |
| Genomic context | in the *kaiABC* operon | standalone |
| Kushige's call, N+ | **b — cycling** | **AR — arrhythmic** |
| Rhythmicity p, N+ (Table S2) | **0.046** | not published — AR genes are absent from Table S2 |
| Published amplitude, N+ (Table S2) | 0.246 | not published |
| Peak time, N+ / N− (Table S2) | 21 h / 16 h | — |
| Kushige's call, N− | c — cycling (looser tier) | **AR — arrhythmic** |
| Raw CV recomputed from N+ series | 0.275 (27th pct) | **0.404 (68th pct)** |
| Heterocyst enrichment Ri | 1.167 (68th pct) | **0.709 (6.7th pct — depleted)** |
| Sum(N−)/Sum(N+) | 1.68 | 1.11 |

Counter to the naive hypothesis that the long form would be the rhythmic one: it's the
**canonical operonic KaiB that's rhythmic**, and the long-form paralog reads as flat/
arrhythmic despite a raw CV in the 68th percentile — well above canonical kaiB's 27th.
A clean illustration of why CV alone without a paired significance test is misleading
(same lesson as H3/H5, different gene). Consistent with the literature's general picture
of long-form KaiB paralogs being structurally/functionally distinct from the canonical
fold-switching KaiB that runs the core oscillator, though that specific mechanistic
claim wasn't independently verified here — this is an expression-level observation only.

**Two obvious follow-ups, both already answerable — have these ready for Mike:**

- *Is the flat call just heterocyst dilution?* **No.** `all3328` has Ri = 0.709, the
  6.7th percentile — it is *depleted* in the heterocyst-enriched fraction, not enriched.
  A transcript confined to a ~1-in-10 cell type would be diluted and could look flat
  (H9), but this one is vegetative-biased, so that explanation is ruled out.
- *Does it switch on under nitrogen deprivation?* **No.** Sum(N−)/Sum(N+) = 1.11,
  essentially unchanged, against 1.68 for canonical *kaiB*. It stays AR in N− too.

### The whole kai locus, for context

| Gene | ORF | Call N+ | Call N− | Raw CV (pct) |
|---|---|---|---|---|
| *kaiA* | alr2884 | AR | c | 0.193 (4th) |
| *kaiB* | alr2885 | **b** | c | 0.275 (27th) |
| *kaiC* | alr2886 | AR | AR | 0.396 (66th) |
| *kaiB* long | all3328 | AR | AR | 0.404 (68th) |

*kaiC* and long-form *kaiB* have nearly identical raw CVs in the 66th–68th percentile and
are both called arrhythmic, while *kaiB* cycles at less than half their CV. That contrast
inside one operon is the tidiest demonstration of the CV-without-a-test problem anywhere
in this project. It also matches the Kushige paper's own headline that the *kai* genes
barely cycle in *Anabaena*.
