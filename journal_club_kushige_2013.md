---
output:
  html_document: default
  pdf_document: default
---
# Journal Club — Kushige et al. 2013

**Paper under discussion:** Kushige H, Kugenuma H, Matsuoka M, Ehira S, Ohmori M, Iwasaki H. 2013. Genome-wide and heterocyst-specific circadian gene expression in the filamentous cyanobacterium *Anabaena* sp. strain PCC 7120. *J Bacteriol* 195:1276–1284. doi:10.1128/JB.02067-12

---

## How these notes are organized

The four questions do **not** all live at the same level, and blurring them will confuse the room. Three layers:

| Layer | What it is | Which questions live here |
|---|---|---|
| **0** | What the 2013 paper itself establishes | — |
| **1** | The paper + its own reference 7 (Ito 2009), which defines its circadian indices | **Q1, Q2** |
| **2** | Literature published *after* this paper | **Q3, Q4** |

**The hard boundary:** Markson et al. (RpaA ChIP-seq) appeared in *Cell* in **December 2013**. Kushige was **accepted 4 January 2013**, published March 2013. RpaA is **never mentioned** in the paper under discussion. Q3 and Q4 are not questions about this paper — they are questions the paper *provoked*, and answering them requires work that did not exist when it was written. Say that out loud early; it protects the discussion.

---

# LAYER 0 — What this paper actually establishes

Everything in this section is in the paper. Nothing here needs an outside citation.

### The clock is real
Three criteria, all met (Fig. 4):  
- **Free-running:** ~24 h under LL for ≥6 days, from *pecB* and *all3173* luciferase reporters.  
- **Temperature compensation:** Q₁₀ = 0.95 (+N) and 0.97 (−N) across 25/30/35 °C.  
- **Entrainment:** PRC to 5-h dark pulses, similar in shape to *Synechococcus*.  
- **Kai-dependent:** Δ*kaiABC* nullifies all tested rhythms (Fig. 5).  

### The headline: *kai* genes barely cycle
Northern + microarray show no reliably high-amplitude rhythm in *kaiA*, *kaiB*, or *kaiC* — in sharp contrast to the high-amplitude *kaiBC* rhythm in *Synechococcus*. *kaiB* squeaks in as a low-amplitude cycler; *kaiA* and *kaiC* do not pass their filter. The extra *kaiB*-like paralog (*kaiBL*, ~3× longer, different locus, function unknown) is flat.

### 78 rhythmic genes
- Array: Ehira & Ohmori oligo array, 5,336 of 5,368 ORFs. Continuous culture, LL after 2× 12:12 LD, sampled every 4 h to 48 h, two independent cultures.  
- Filter: cosine-correlation *P* < 0.05 **and** amplitude > 10⁻⁰·⁷, where **amplitude ≡ CV (SD/mean)**.  
- Result: **78 genes (1.5% of the genome)**. 18 at stringent thresholds; 192 at loose.  

### Phase structure
Bimodal, peaks near CT 8–12 and CT 20–24 — same two-peak structure as *Synechococcus* — but with **more genes at intermediate phases** (CT 4–8, 12–20) than in *Synechococcus*.  

### Functional categories (this is the paper's own answer to Q2, for *Anabaena*)
Rhythmic fraction above the 1.5% genome baseline in: photosynthesis/respiration (9.2%, 14/153), translation (3.7%, 7/187), amino acid biosynthesis (3.6%, 4/111), cellular processes (3.3%, 3/92), energy metabolism (3.1%, 3/98).  

- **Subjective day:** photosynthesis/respiration/cell envelope. Includes *apcEABC*, *pecBACEF*, the *cpcG* cluster.  
- **Subjective night:** translation, amino acid biosynthesis, energy metabolism. Includes the ribosomal cluster (*all4198*–*all4215*), *secE*, *groEL*.  

### The comparative observation buried in the Results
The paper states plainly that **the homologues of the rhythmic phycobilisome genes are *not* rhythmic in *Synechococcus***. Combined with the *kai* result (top-amplitude in 7942, flat in *Anabaena*), the paper has already told you the gene-level output program does not match. **This is Q1's answer, hiding in the text.**

### Δ*kaiABC* does something *Synechococcus* doesn't
Nullifying *kaiABC* not only abolishes the rhythms — it **reduces the expression level** of *pecBAC* and *all3173*. In *Synechococcus*, no rhythmic gene is dramatically downregulated by *kaiABC* nullification. The clock in *Anabaena* is doing more than modulating; it's contributing baseline drive.

### Total mRNA doesn't crash in the dark
Unlike *Synechococcus* (total mRNA → ~20% at dark onset), *Anabaena* shows no dramatic change across the dark-to-dark transition. Consistent with *Cyanothece*, not with *Synechococcus*.

### Nitrogen deprivation and heterocysts
- 39 of the 78 clock-controlled genes keep circadian profiles in BG-11₀.  
- **Heterocyst-specific cyclers:** *all1427* (high amplitude, peaks subjective dusk); *hesA*/*hesB* (*all1432*–*all1431*, peak subjective dawn). Estimated 27-, 15-, and 10-fold higher in a single heterocyst than a vegetative cell.  
- Not expressed in the *hetR*-null → genuinely heterocyst-specific.  
- Rhythms abolished in Δ*kaiABC* under −N.  
- *kaiABC* expression is roughly equal in heterocyst-enriched vs. control (1.09 ≤ R_kai ≤ 1.17).  
- Δ*kaiABC* still forms heterocysts at semiregular intervals → **the clock is not required for patterned heterocyst formation.**. 
- *nifHDK* and *fdxN* were **not** clearly rhythmic under LL in their hands.  

---

# LAYER 0 — What this paper *cannot* tell you

### Limitations the authors state themselves
1. **They cannot exclude that heterocyst rhythms are driven indirectly** by time-dependent intercellular signals from oscillators in neighboring vegetative cells. (They say this outright, and propose single-cell fluorescence reporters as the fix.). 
2. **Transcript rhythm ≠ protein rhythm.** They flag this twice. Nothing here is proteomic.  

### Limitations they don't flag
3. **The FDR on the 78 is soft.** Their "correlation *P* value" is a Pearson correlation of a 13-timepoint profile against a 24-h cosine. With autocorrelated noise, the effective df is well below 13, so the *P* value is anticonservative. Meanwhile the amplitude threshold (CV > 10⁻⁰·⁷ ≈ 0.20) is *stringent*. The 78 is the product of a loose test and a tight filter, and neither is obviously well-calibrated.  
4. **"Amplitude" = CV, which is mean-dependent.** Lowly-expressed genes get inflated CV. There is no expression-level correction anywhere.  
5. **Low *kai* amplitude is not evidence of a weak clock** — and *their own references establish this.* Tomita 2005 and Hosokawa 2011 (their refs 8, 9) show the transcription-translation loop is dispensable in *Synechococcus*; the post-translational oscillator is the core timer. The headline finding is far less alarming than it reads, and the authors know it. **Watch whether the room catches this.**  
6. **The heterocyst enrichment math rests on two assumptions.** They invert Rᵢ = n(0.8Hᵢ + 0.2Vᵢ)/[Hᵢ + (n−1)Vᵢ] to get the 27/15/10-fold figures, assuming (i) total RNA content of a heterocyst ≈ that of a vegetative cell, and (ii) the max-enrichment gene (*asl3656*) is heterocyst-exclusive. Assumption (i) is doing a lot of work for a cell type that has shut down PSII and rebuilt its envelope. Their validation — that the max ratio (9.28) slightly exceeds *nifK*'s (8.98) — is suggestive, not dispositive. **[Criticism is mine; the assumptions are theirs, stated on p. 1278.]**  
7. **They never actually do the cross-species comparison.** The phycobilisome observation is an aside. No ortholog mapping, no correlation, no systematic analysis. This is the gap Mike's Q1/Q2 sit in.

---

# LAYER 1 — Q1 and Q2

These are answerable **within the paper's own methodological frame**, because Kushige explicitly adopt the circadian index definitions (amplitude = CV, correlation *P*, peak time) from **Ito et al. 2009 (their ref 7)**, which is the equivalent genome-wide dataset for *S. elongatus*. Same definitions, different array platforms. That's why the comparison is even possible — and why it's tricky.

### Q1: gene-by-gene amplitude and phase correlation

**Short answer, from the paper alone:** weak to nonexistent. The two most salient gene sets go in *opposite* directions (*kai*: top-amplitude in 7942, flat in *Anabaena*; phycobilisome: rhythmic in *Anabaena*, not in 7942).

**Why the question is hard to answer properly — raise these:**

- **"Homologous" isn't operational.** *Anabaena* has 5,368 genes to 7942's ~2,700, with large paralog families. Need an explicit criterion (reciprocal best hits vs. orthogroups).
- **CV is platform-dependent.** Two different arrays, different dynamic ranges, different noise floors. Absolute amplitudes are not comparable — only ranks.
- **Selection on rhythmicity biases the answer.** Conditioning the sample on "rhythmic in 7942" is selection on a noisy variable → range restriction + regression to the mean in the *Anabaena* column. Guarantees a depressed correlation *and* a spurious "*Anabaena* amplitudes are lower" story.
- **Attenuation.** Both amplitudes are noisy, so any observed *r* is attenuated. Both datasets have n = 2 independent experiments → reliability is estimable → a disattenuated *r* with bounds is possible. Without it you cannot distinguish "not conserved" from "conserved but measured badly."
- **Phase is circular** and meaningless at low SNR. With 4-h sampling and roughly bimodal phase distributions, a contingency table on phase class (dawn vs. dusk) is more robust and more interpretable than a circular correlation coefficient.

**The elephant:** ~856 genes reported rhythmic in *S. elongatus* vs. 78 here — roughly 32% vs. 1.5% of the transcriptome. Before calling a weak correlation "biology," bound how much of that ~20× asymmetry is detection threshold. *(Note: the 856 figure is not in the paper under discussion — see Layer 2 and the Verify list.)*

### Q2: pathway-level correlation

**Why it's the right robustness check:** averaging phase across *n* genes in a category shrinks independent measurement noise by ~1/√n. A conserved module-level program can survive even when gene-level *r* is attenuated to ~0.

**Kushige gives you the *Anabaena* half for free** (see Functional categories above). The open half is *S. elongatus*'s pathway phases, from Ito's supplement.

**Two separate questions — don't conflate:** (a) is *rhythmicity* conserved by category? (b) is *phase* conserved by category? Null model must permute gene→category labels within species, since category sizes and family expansions differ.

---

# LAYER 2 — Q3 and Q4 (post-publication; *not* in this paper)

**Flag this transition explicitly in the room.** None of what follows was available to Kushige.

### Q3 — a premise, with one correction

Markson et al. 2013 found ~**110 RpaA binding sites** → **134 target transcripts**; 93 encode proteins/tRNAs (~170 genes in operons), 41 are ncRNAs. Against ~856 cycling genes in 7942, that means **most of the cycling transcriptome is *indirect***, relayed through an RpaA-driven sigma-factor cascade (Fleming & O'Shea).

"The highest-amplitude genes are direct RpaA targets" — fine. "The cycling transcriptome *is* the RpaA regulon" — no. The difference matters for Q4.

Mechanism: SasA (kinase) and CikA (phosphatase) set RpaA~P (Gutu & O'Shea 2013). In 7942 the output is delivered by RpaA **phosphorylation**, not abundance.

### Q4 — already done, by Arbel-Goren et al. 2024

They scanned the *Anabaena* genome with **FIMO (MEME Suite v5.4.1)** using Markson's RpaA position-specific probability matrix, both strands, min match *P* < 10⁻⁴, restricted to −500/+50 bp around TSSs from Mitschke et al. 2011. *Anabaena* RpaA = **All0129**, ~96% aa similarity to the 7942 protein.

- **81 genes** with putative RpaA sites at FIMO q < 0.05.
- **Two putative RpaA sites upstream of *pecB*** (*P* = 1.5 × 10⁻⁵) — i.e. upstream of the very gene Kushige built their reporter on.
- Regulators among the hits: FurC (Alr0957), Alr2325 (cAMP-binding), SigE, Alr3646, and three two-component regulators (All3822, Alr5272, Alr5069).
- They also BLASTP'd 89 reported 7942 ChIP targets against *Anabaena* proteins (E ≤ 0.005, >36% similarity).

### What they did NOT do — and this is the opening

They crossed **motif hits × orthologs of the 7942 ChIP targets**. They never crossed **motif hits × Kushige's 78 rhythmic genes**. Their treatment of the 81 is a hand-picked list of regulators with **no enrichment analysis of any kind**.

**The test that's still sitting there, unclaimed:** hypergeometric overlap of the 81 motif-bearing genes with the 78 rhythmic genes. Background = genes with *both* a Mitschke TSS annotation *and* a slot on the Ehira/Ohmori array (not all 5,336 — getting this wrong inflates the enrichment). Under independence you'd expect ~1–2 overlaps. Even five or six is a large odds ratio. Then ask whether the overlapping genes are **phase-concentrated** (as direct RpaA targets are in 7942) or scattered.

### Two objections to the motif scan

1. **RpaA and RpaB are not distinguishable by sequence.** RpaB binds HLR1 — a pair of imperfect 8-nt direct repeats, (G/T)TTACA(T/A)(T/A), separated by two nucleotides. The HLR1 element in the *kaiBC* promoter **overlaps the reported RpaA site**, is bound by RpaB *in vivo*, and peaks ~12 h out of phase with RpaA (Hanaoka et al. 2012; Espinosa et al. 2015). Both are OmpR/PhoB-family regulators; RpaB is essential and conserved in *Anabaena*. **A FIMO scan with the RpaA PSPM cannot tell you which one you found.** The honest label for the 81 is "RpaA/RpaB-family sites."
   - **Sharpest version:** Arbel-Goren cite Espinosa 2015 (their ref 40) — *but only* to note that SasA is regulated by RpaB. They had the reference in hand and never turned it on their own method. **[The components are sourced; applying it as a criticism of their scan is mine.]**
2. **GC-content mismatch.** *Anabaena* is ~41% GC; *S. elongatus* ~55%. A PWM trained on a GC-rich genome and applied to an AT-rich one inflates hit rates unless the background model is an *Anabaena*-specific higher-order Markov model. FIMO's q-values fix multiple testing, not base-composition mismatch in the null. **[Criticism is mine. GC figures — verify.]**

### One thing the later work *closed* — worth saying

Kushige explicitly flagged that they could not exclude indirect induction of heterocyst rhythms by neighboring vegetative cells. **Arbel-Goren 2024 closed it** with single-cell imaging: P*pecB*-*gfp* oscillates in individual heterocysts, with the phase inherited from the progenitor's sister cell. It took eleven years and a different technique to answer a limitation the authors named themselves. That's a good note to end the discussion on.

---

# Discrepancies between this paper and the later work

**These are real, and I can't resolve them from the texts.**

| | Kushige 2013 | Arbel-Goren 2021/2024 |
|---|---|---|
| **Free-running period** | ~24–25 h (bioluminescence; the entire temperature-compensation argument rests on it) | **20.8 ± 0.4 h** vegetative; **21.3 ± 0.6 h** heterocysts (single-cell GFP) |
| ***pecB* phase** | CT 8–12 — **subjective day**, approaching dusk | *rpaA* and *pecB* peak during **subjective night** |

Same organism, same reporter gene, same entrainment protocol (2× 12:12 LD → LL). ~15% gap on the period and a substantial phase disagreement.

**[CONJECTURE — mine]** Probably culture geometry: Kushige used agar plates under white fluorescent at 30 µmol m⁻² s⁻¹; Arbel-Goren used agarose pads under LED + tungsten halogen at ~25 µmol total. But I don't actually know. **The question worth asking the room: how much of the "*Anabaena* clock is weird" story is the clock, and how much is the culture?**

---

# Optional synthesis to offer

**[ALL OF THIS IS MY SYNTHESIS, NOT ANY PAPER'S CLAIM. Present it as your read.]**

> The oscillator is functionally conserved. The output transducer is conserved *by homology*. The amplifier is not.

Unpacked, with the evidence grade on each clause:

- **"Oscillator functionally conserved"** — *supported, with a caveat.* Kushige Fig. 4 gives a proper 24-h temperature-compensated entrainable Kai-dependent oscillator. But the **parts are divergent**: *Anabaena* KaiA is missing the N-terminal two-thirds (Dvornyk 2003), there's an extra *kaiBL*, and transplanting *Anabaena kaiA* into *Synechococcus* yields a **~40 h period** at lower amplitude (Uzumaki 2004). Say "functionally conserved as a 24-h temperature-compensated oscillator, with visibly divergent parts." Someone *will* hit you with the truncated KaiA.  
- **"Transducer conserved by homology"** — *thin.* The chain is: sequence homologs of RpaA/SasA/CikA exist (Schmelling 2017, i.e. BLAST); RpaA is 96% similar; a motif scan finds sites in plausible places. **Nobody has shown *Anabaena* SasA phosphorylates *Anabaena* RpaA, or that CikA dephosphorylates it, or measured RpaA~P in *Anabaena* at all.** Arbel-Goren lean on that gap in their own Discussion — their explanation for reduced *pecB* in heterocysts hedges twice in one sentence, and the load-bearing quantity (RpaA~P) is one nobody has measured. **Naming the missing experiment is more interesting than asserting conservation.**  
- **"Amplifier is not conserved"** — *partly supported.* *Anabaena* lacks the high-amplitude *kaiBC* TTFL (Kushige). *Anabaena* instead oscillates the ***rpaA* transcript itself** with high amplitude, kai-dependently (Arbel-Goren 2021, RT-qPCR) — an additional layer that 7942 does not use to relay the signal. So the gain of the master output node is time-varying in one organism and not the other.

---

# CONJECTURE REGISTER

Everything below is **mine**. None of it is a claim made by any of the papers. Do not attribute it to them.

1. **"*Anabaena*'s sigma-factor relay is different."** Pure inference, from the fact that 7942's cycling transcriptome runs through an RpaA-driven sigma cascade and *Anabaena*'s sigma complement differs. **Nobody has shown this.** Drop it or label it.
2. **The oscillator/transducer/amplifier framing.** My synthesis.  
3. **Prediction:** gene-level amplitude/phase correlation (Q1) will be weak; pathway-level (Q2) will survive. Untested.  
4. **Prediction:** the 81 × 78 hypergeometric will show enrichment. Untested.  
5. **Speculation** that the period/phase discrepancies between the 2013 and 2021 papers are culture-condition artifacts. I don't know.  
6. **The RpaB/HLR1 identifiability criticism** as applied to Arbel-Goren's scan. Components sourced (Hanaoka 2012; Espinosa 2015); the application is mine.  
7. **The GC-content criticism** of the motif scan. Mine.  
8. **The attenuation / selection-effect critique** of Q1. Standard statistics; my application here.  
9. **The critique of the heterocyst enrichment math.** Mine. (Their assumptions are stated; the objection is not.). 

---

# VERIFY BEFORE THE MEETING

- [ ] **Arbel-Goren 2024's supplement captions are wrong.** The published caption lists Table S1 as P*pecB*-*gfp* scatterplots; the body text cites Table S1 twice as the RpaA motif table (the 81 genes, and the *pecB* p-value). **Download the .xlsx and check what's actually in it** before promising anyone it has the motif hits.
- [ ] **The "856" figure.** Cited in Arbel-Goren 2024 to Markson 2013 + Markson & O'Shea 2009 (a review). The entire 78-vs-856 comparison — the elephant in Q1 — rests on this denominator. Trace it to the primary source.
- [ ] **Ito et al. 2009 supplement:** confirm it actually tabulates per-ORF amplitude / correlation *P* / peak time for all ~2,700 genes. Kushige's methods imply it, but confirm.
- [ ] **Kushige Data Set S1:** confirm it covers all ~5,336 ORFs, not just the 78.
- [ ] **GC contents** (~41% *Anabaena*, ~55% 7942). From memory. Check.
- [ ] **Raw arrays:** KEGG Expression Database accessions ex0001892–ex001947. KEGG's expression database is very likely long dead. If you need raw time series (to recompute phase with JTK/RAIN, or to estimate replicate reliability for the attenuation correction), budget on emailing Iwasaki (Waseda) or Ehira.
- [ ] Full citations for Hanaoka 2012, Riediger (HLR1 PWM), and Fleming & O'Shea — flagged below.

---

# REFERENCES

### Under discussion
- **Kushige H, Kugenuma H, Matsuoka M, Ehira S, Ohmori M, Iwasaki H. 2013.** *J Bacteriol* 195:1276–1284. doi:10.1128/JB.02067-12

### Cited *within* Kushige (useful for the room — these are the paper's own refs)
- **Ito H, et al. 2009.** *PNAS* 106:14168–14173. — [Kushige ref 7] **Source of the circadian index definitions and the 7942 comparison set. This is the dataset Q1/Q2 need.**
- **Tomita J, Nakajima M, Kondo T, Iwasaki H. 2005.** *Science* 307:251–254. — [ref 8] No TTFL required for the KaiC phosphorylation rhythm.
- **Hosokawa N, et al. 2011.** *PNAS* 108:15396–15401. — [ref 9] Circadian transcription without *de novo kai* expression.
- **Nakajima M, et al. 2005.** *Science* 308:414–415. — [ref 12] In vitro reconstitution.
- **Dvornyk V, Vinogradova O, Nevo E. 2003.** *PNAS* 100:2495–2500. — [ref 14] KaiA N-terminal truncation in filamentous species.
- **Uzumaki T, et al. 2004.** *Nat Struct Mol Biol* 11:623–631. — [ref 15] *Anabaena kaiA* → *Synechococcus*: ~40 h period, lower amplitude.
- **Garces RG, Wu N, Gillon W, Pai EF. 2004.** *EMBO J* 23:1688–1698. — [ref 16] *Anabaena* KaiA/KaiB structures.
- **Ehira S, Ohmori M. 2006.** *Mol Microbiol* 59:1692–1703. — [ref 26] The oligonucleotide array used here.
- **Schmitz O, Katayama M, Williams SB, Kondo T, Golden SS. 2000.** *Science* 289:765–768. — [ref 31] CikA; the 7942 PRC they compare against.

### Post-publication — NOT available to Kushige
- **Markson JS, Piechura JR, Puszynska AM, O'Shea EK. 2013.** Circadian control of global gene expression by the cyanobacterial master regulator RpaA. *Cell* 155:1396–1408. doi:10.1016/j.cell.2013.11.005 — **Published Dec 2013. Kushige was accepted 4 Jan 2013.**
- **Gutu A, O'Shea EK. 2013.** Two antagonistic clock-regulated histidine kinases time the activation of circadian gene expression. *Mol Cell* 50:288–294. doi:10.1016/j.molcel.2013.02.022
- **Arbel-Goren R, Buonfiglio V, Di Patti F, Camargo S, Zhitnitsky A, Valladares A, Flores E, Herrero A, Fanelli D, Stavans J. 2021.** Robust, coherent, and synchronized circadian clock-controlled oscillations along *Anabaena* filaments. *eLife* 10:e64348. doi:10.7554/eLife.64348 — **Source of the high-amplitude *rpaA* transcription result (RT-qPCR) and the ~20.8 h period.**
- **Arbel-Goren R, Dassa B, Zhitnitsky A, Valladares A, Herrero A, Flores E, Stavans J. 2024.** Spatio-temporal coherence of circadian clocks and temporal control of differentiation in *Anabaena* filaments. *mSystems* 9:e00700-23. doi:10.1128/msystems.00700-23 — **Source of the FIMO scan, the 81 genes, All0129, and the Δ*kaiABC* diazotrophy phenotype.**
- **Espinosa J, Boyd JS, Cantos R, Salinas P, Golden SS, Contreras A. 2015.** Cross-talk and regulatory interactions between the essential response regulator RpaB and cyanobacterial circadian clock output. *PNAS* 112:2198–2203. doi:10.1073/pnas.1424632112
- **Mitschke J, Vioque A, Haas F, Hess WR, Muro-Pastor AM. 2011.** *PNAS* 108:20130–20135. — *Anabaena* TSS annotations used by the FIMO scan.
- **Grant CE, Bailey TL, Noble WS. 2011.** FIMO. *Bioinformatics* 27:1017–1018.
- **Schmelling NM, et al. 2017.** Minimal tool set for a prokaryotic circadian clock. *BMC Evol Biol* 17:169. — The BLAST-based survey behind "*Anabaena* has RpaA/SasA/CikA."

### Cited above but **full citations not verified** — check before quoting
- **Hanaoka M, et al. 2012.** "RpaB, another response regulator operating circadian clock-dependent transcriptional regulation in *Synechococcus elongatus* PCC 7942." *J Biol Chem.* — HLR1 in P*kaiBC* overlaps the RpaA site.
- **Riediger M, et al.** HLR1 / RpaB regulon; source of the HLR1 position weight matrix built from ~90 aligned motifs. (bioRxiv 443713; a published version exists — find it.)
- **Fleming KE, O'Shea EK. 2018.** An RpaA-dependent sigma factor cascade sets the timing of circadian transcriptional rhythms in *S. elongatus.* — The indirect layer between RpaA and the ~856 cyclers.
