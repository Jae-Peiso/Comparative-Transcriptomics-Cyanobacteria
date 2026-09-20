# Vijayan 2009 and Markson 2013 — papers and datasets

Reference sheet for the new 7942 side. Provenance is marked on every claim:
**[GEO]** = read verbatim from the GEO record · **[paper]** = from the abstract/summary ·
**[JC]** = from `journal_club_kushige_2013.md` · **[unverified]** = needs checking.

---

## Vijayan, Zuzow & O'Shea 2009 — the primary dataset

**Oscillations in supercoiling drive circadian gene expression in cyanobacteria.**
*PNAS* 106(52):22564–22568. doi:10.1073/pnas.0912673106. PMID 20018699.
GEO: **GSE18902**. Platform **GPL9534**. Public 30 Nov 2009. Contact: Vikram Vijayan (Harvard).

### What the paper is about **[GEO summary, verbatim]**

7942 shows 24-h oscillations in transcript abundance under constant light, and neither
the *cis* nor *trans* factors setting circadian phase had been identified. The paper
shows **the topological status of the chromosome is highly correlated with circadian
gene expression state** — that genes monotonically activated or repressed by chromosomal
relaxation across the cycle resemble supercoiling-responsive genes in *E. coli*, and that
perturbing superhelicity within the physiological range produces global expression
changes resembling the normal circadian cycle.

So: **a chromosome-topology paper, not a clock-architecture paper.** That matters for us
only in that its circadian time course was built to characterise expression genome-wide,
not to defend a particular rhythmicity threshold.

### Design **[GEO overall design, verbatim]**

> Synechococcus elongatus PCC 7942 was subjected to two consecutive light/dark cycles and
> released into continuous light (T = 0). Cells were sampled every 4 hours from T = 24 to
> T = 84 hours for microarray analysis to characterize circadian gene expression.

Plus a **separate novobiocin experiment** — sampled at T = 56 and 64 h, treated with
novobiocin (0.1 µg/ml), measured at 5, 10, 30, 90 and 150 min — to test whether forced
chromosome relaxation reproduces circadian-like expression change.

| | |
|---|---|
| Entrainment | 2× 12h:12h light/dark, then release into **constant light** |
| Circadian arm | **16 timepoints**, every 4 h, T = 24 → 84 h (**60 h, ~2.5 cycles**) |
| Accessions | GSM468463 (T=24) … GSM468478 (T=84), all "experiment 1" |
| Novobiocin arm | GSM468479–468485 — **7 samples, exclude these** |
| Replication | **n = 1** for the LL series |

### Reported result **[unverified — confirm against the paper before quoting]**

**1,748 of 2,724 ORFs (~64%) oscillate with a 22–26 h period.** This came from a search
summary, not from the paper itself. It is the number that makes Vijayan attractive
against Ito's 97, so it is worth 10 minutes to verify directly.

---

## Markson et al. 2013 — the replicate, and the RpaA paper

**Circadian control of global gene expression by the cyanobacterial master regulator RpaA.**
*Cell* 155:1396–1408. doi:10.1016/j.cell.2013.11.005. PMID 24315105.
SuperSeries **GSE50922**. Contact: Joseph Markson (O'Shea lab, Harvard).

### Two separate things live under this citation — do not conflate them

**1. The WT LL time course — GSE52486 — this is what we are using.**

Its GEO summary states its purpose outright **[GEO, verbatim]**:

> The goal of the experiment was to obtain **a replicate of the wild-type LL circadian
> timecourse published in Vijayan et al, PNAS 106: 22564-22568 (2009)**, in order to
> identify reproducible circadian genes in LL.

Design **[GEO, verbatim]**: turbidostat as in Vijayan except 3 L rather than 4.5 L
culture; two consecutive light/dark cycles then release into constant light at T = 0;
sampled every 4 h from T = 36 to T = 64 h inclusive; **T = 52 h omitted for poor data
quality**. Gene expression at each timepoint compared to time-averaged expression using
a pool of equal-mass RNA from all timepoints, on a **two-colour Agilent microarray**.

| | |
|---|---|
| Timepoints | **7** — GSM1267750 (T=36) … GSM1267756 (T=64), T=52 absent |
| Relationship to Vijayan | purpose-built biological replicate, same protocol, same platform |
| Why it matters | gives **H4** its reliability estimate — previously logged as permanently blocked |

**2. The RpaA regulon work — the ChIP-seq and the Δ*rpaA* arms.** This is what Mike
pointed at, and it is Part C, not Part B. Other subseries on the same platform:
GSE50908 (Δ*rpaA* timecourse), GSE50919 (Δ*rpaA* vs WT), GSE50920 (clock rescue strains).

### The RpaA findings **[JC]**

- ~**110 RpaA binding sites** → **134 target transcripts**: 93 encode proteins or tRNAs
  (~170 genes once operons are counted), 41 are ncRNAs.
- Against ~856 cycling genes in 7942, **most of the cycling transcriptome is indirect**,
  relayed through an RpaA-driven sigma-factor cascade (Fleming & O'Shea).
- Mechanism: SasA (kinase) and CikA (phosphatase) set RpaA~P (Gutu & O'Shea 2013). In
  7942 the clock output is delivered by RpaA **phosphorylation, not abundance** — so
  *rpaA* transcript level is not the readout.

**Consequence for Mike's premise:** "the highest-amplitude genes are direct RpaA targets"
is supportable. "The cycling transcriptome *is* the RpaA regulon" is not.

---

## GPL9534 — the shared platform

**Agilent-020846 *S. elongatus* PCC 7942 oligo microarray**, 8×15k, Agilent array ID
020846. **2,723 features.** Submitted by Vijayan, Nov 2009. **[GEO]**

Platform table columns: `ID` (= JGI ID), **`ORF` (= locus tag)**, `JGI_ID`, `Chromosome`,
`RANGE_GB`, `RANGE_START`, four `SEQUENCE` columns (probes 1–4, Agilent linker
`TATCCTACTATACGTATCACATAGC` stripped), `SPOT_ID`.

**The `ORF` column is already in `Synpcc7942_XXXX` form** — the exact key used by
`data/processed/rbh_orthologs.csv`, which takes it from GFF `old_locus_tag`. No
probe→protein reconciliation. This is the single largest cost difference between this
dataset and the alternatives, and it is invisible from either paper.

Eight series share this platform, which is why cross-study comparison is clean here:
GSE18902 (Vijayan), GSE50908/50919/50920/50922/52486 (Markson), GSE59112
(Pattanayak/Rust *cikA*), GSE102914 (ZnO toxicity, irrelevant).

---

## What actually matters for our analysis

**Both are two-colour arrays**, and every value is **log2(sample / time-averaged
reference)** — a ratio, not an intensity.

### What "equal-mass RNA" means, and why it matters

A two-colour array measures *two* samples at once, in different dye channels, and
reports their ratio. So you need something to compare each timepoint against. Vijayan
and Markson built that reference by **taking the same weight of total RNA from every
timepoint — say 500 ng from each of the 16 — and mixing them into one pooled sample.**
Equal *mass*, not equal volume and not equal number of cells: it doesn't matter that one
timepoint yielded more culture or more RNA overall, each contributes the same amount to
the pool.

Two consequences:

- **The pool approximates the average expression profile across the whole cycle.** A
  gene high at dawn and low at dusk contributes both to the mix, so it sits mid-range in
  the reference.
- **Every array in the series shares that identical denominator.** This is a "common
  reference" design, and it's what makes ratios comparable from one timepoint to the
  next — otherwise each array would be on its own arbitrary scale.

So a reported value of +1 means "twice the cycle average for this gene", −1 means half,
and 0 means at the average.

### Which part of this actually matters — the pooling, not the mass

Three separable design choices are bundled in that sentence. They are worth pulling
apart, because only the first two affect us:

| | Choice | What it causes | Matters to us? |
|---|---|---|---|
| 1 | **Two-colour at all** | you get a **ratio**, not an intensity | **Yes** — the primary difference from Kushige |
| 2 | **Reference pooled from all timepoints** | the ratio is **centred on the cycle mean**, so values sit near zero | **Yes** — this is what makes CV undefined |
| 3 | **Equal mass rather than equal volume** | that average is **unweighted** rather than skewed toward whichever timepoints yielded more RNA | Barely |

### The ratio format is not a handicap — be careful how this gets described

It is easy to write this up as though ratio data were a degraded version of intensity
data. It is not, and for our purpose it is arguably the better design:

- **Probe-specific systematics cancel.** Numerator and denominator pass through the same
  physical spot on the same array, so probe affinity, spot size and dye effects divide
  out. Single-channel intensities carry all of that into the number.
- **It is already the form we are going to construct from Kushige anyway.** Fitting a
  cosinor to Kushige's intensities means log2-transforming them; the fit then splits into
  a **mesor** (mean log expression) plus an **amplitude**. Vijayan's values are
  log2(sample / cycle mean) = log2(sample) − log2(mean) — the same quantity with the
  mesor already subtracted. They arrive one step further along the road we are taking,
  and since we discard the mesor either way, nothing we need is lost.

What a ratio genuinely costs is the **absolute expression level**: you cannot ask "is
this gene highly expressed?" from Vijayan alone. A real limitation, and irrelevant to
amplitude and phase.

So the accurate claim is narrow: **CV is what breaks, not the data.** And CV was a weak
amplitude measure on intensities in the first place — mean-dependent, inflating for
low-expression genes, which is the original H3 complaint. The ratio format does not
degrade anything; it takes a statistic that was already poorly chosen and makes it
undefined.

Choice 2 was not forced by choice 1. A two-colour experiment could use direct pairwise
or loop designs, or hybridise against a fixed external reference — in which case values
would still be ratios but would *not* be centred on zero. Building the common reference
out of every timepoint is the specific decision that produces the zero-centring.

**Choice 3 is nearly irrelevant to this analysis.** Uneven pooling would shift each
gene's **mesor** — the vertical offset of the fitted wave — by a constant. It would not
change the **amplitude** or the **peak time**, and those are the only two quantities we
take from the fit. So equal-mass pooling is good practice on their part rather than
something our analysis depends on.

Short version for explaining to someone else: *it's a ratio against a pooled reference*
is the load-bearing part; *equal-mass* is just how they made the pool fair.

### Why this kills CV

CV = SD / mean. These values are log ratios centred near zero by construction, because
the reference *is* the average. **A gene's mean across the series is therefore ≈ 0 —
and may be slightly negative.** Dividing by it gives a number that explodes toward
infinity, flips sign arbitrarily, and carries no information about how strongly the gene
oscillates. CV is not merely a worse choice here; it is undefined in practice.

The right measure is the **fitted cosine amplitude in log2 units** — the half-height of
the fitted wave, which is exactly what `cosinor()` in `src/07_fetch_vijayan_markson.py`
returns. Kushige's intensities have to be log2-transformed and put through the same fit
rather than carrying their published CV across, or the two sides are not comparable.

This is H3 reappearing in a new costume, and it is the easiest thing here to get silently
wrong: computing CV on log ratios produces numbers, not an error.

**Design match to Kushige**, which is why Vijayan won on merit rather than availability:

| | Kushige 2013 (*Anabaena*) | Vijayan 2009 (7942) | Markson WT (7942) |
|---|---|---|---|
| Condition | constant light | constant light | constant light |
| Sampling interval | **4 h** | **4 h** | **4 h** |
| Duration | 48 h | 60 h | 28 h |
| Timepoints | 12 distinct | 16 | 7 |
| Replication | n = 2 | n = 1 | n = 1 (replicate of Vijayan) |

**Not yet done:** neither dataset has been downloaded or parsed.
`src/07_fetch_vijayan_markson.py` is written and its `cosinor()` and
`verify_sample_map()` are unit-tested, but the fetch, parse and join path has never run.

---

## Citations

- Vijayan V, Zuzow R, O'Shea EK. 2009. Oscillations in supercoiling drive circadian gene
  expression in cyanobacteria. *PNAS* 106(52):22564–22568. doi:10.1073/pnas.0912673106
- Markson JS, Piechura JR, Puszynska AM, O'Shea EK. 2013. Circadian control of global
  gene expression by the cyanobacterial master regulator RpaA. *Cell* 155:1396–1408.
  doi:10.1016/j.cell.2013.11.005
- GEO: [GSE18902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18902) ·
  [GSE52486](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52486) ·
  [GSE50922](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50922) ·
  [GPL9534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL9534)
