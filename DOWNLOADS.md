# DOWNLOADS.md — Manual retrieval of the two supplements

> Why manual: every automated fetch of these SIs hit a reCAPTCHA / bot wall during
> scoping. That is not a bug to route around — it is exactly the barrier a human
> clears in seconds. Do this by hand, from a browser, ideally on institutional
> network or VPN (both papers' full SIs are free, but the hosts gate scrapers hard).
>
> Time budget: ~30 min if access is clean. Hard stop at half a day (see brief §3).
> The GOAL is not just to download — it is to answer, for each file:
> **"Does this contain per-gene amplitude + peak time, or only the rhythmic subset?"**

---

## Before you start
- Be on your institution's network or VPN (for full-text/SI access if prompted).
- Have a spreadsheet app ready (Excel / LibreOffice / Numbers) — both SIs are
  spreadsheet or spreadsheet-embedded.
- Create `project/data/raw/` and save everything there with clear names. Do NOT edit
  raw files in place — copy to `interim/` before touching.

---

## FILE 1 — Ito et al. 2009 (the *S. elongatus* 7942 side)

**Citation:** Ito H et al. 2009, PNAS 106(33):14168-14173. doi:10.1073/pnas.0902587106
**PMCID:** PMC2729038

This is an old-style PNAS paper (2009): SI is a set of individual files
("Supporting Information" / "DCSupplemental"), not one merged appendix. Datasets are
usually separate .xls files labeled Dataset S1, S2, ...

### Route A — PNAS article page (try first)
1. Go to the DOI: **https://doi.org/10.1073/pnas.0902587106** (resolves to pnas.org).
2. If you hit a login/paywall gate on the article: the SI is still usually free —
   look for a **"Supporting Information"** or **"Figures & SI"** tab/link, not the
   full-text PDF.
3. On the SI page you'll see items like *SI Appendix*, *SI Methods*, and one or more
   **Dataset** files (.xls / .xlsx). Download **every Dataset file** — you don't yet
   know which holds the per-ORF index table.
4. If PNAS itself throws a CAPTCHA at you repeatedly: switch to Route B.

### Route B — PMC (free full text + SI mirror)
1. Search Google/PubMed for: `Ito 2009 Cyanobacterial daily life PMC2729038`.
2. Open the **PMC** article (pmc.ncbi.nlm.nih.gov/articles/PMC2729038/).
3. Scroll to the **"Supplementary Materials"** / **"Associated Data"** section near
   the bottom. Download the linked dataset files there.
   - If PMC shows the "Checking your browser" reCAPTCHA: solve it once in the browser
     (that's the whole point of doing this manually), then the SI links work.
4. PMC sometimes stores SIs under a "supplementary-material" sub-path; the section
   links handle that for you.

### Route C — last resort
- **PNAS SI direct path pattern** (2009-era):
  `https://www.pnas.org/doi/suppl/10.1073/pnas.0902587106` — browse the file list.
- Email the corresponding author (Hideo Iwasaki, Waseda) for the dataset if both
  routes fail. Only worth it if this becomes the actual blocker — for a few-days task,
  if you can't get it in half a day, fall back to Version A with just the 97 rhythmic
  genes (which are also reported in the paper's tables/figures).

### THE CHECK (do this the moment a file opens)
- Open each Dataset .xls. Look at the **row count** and column headers.
- **Full table** = thousands of rows (~2,515) with columns for amplitude/CV,
  correlation P (or FDR), and peak time / phase per gene. -> Version B is *possible*.
- **Subset only** = ~97 rows (or split 87 + 10). -> Version A only for the 7942 side;
  note in hurdles.md that the full per-ORF table wasn't archived (H7).
- Either way, Version A is fine — it only needs the 97 with their indices.

---

## FILE 2 — Kushige et al. 2013 (the *Anabaena* 7120 side)

**Citation:** Kushige H et al. 2013, J Bacteriol 195(6):1276-1284. doi:10.1128/JB.02067-12
**Publisher:** American Society for Microbiology (ASM), journals.asm.org

### Route A — ASM article page (try first)
1. Go to the DOI: **https://doi.org/10.1128/JB.02067-12** (resolves to journals.asm.org).
2. Find the **"Supplemental Material"** section/tab on the article page.
3. Kushige lists several supplements — you specifically want **Data Set S1** (the
   per-ORF microarray output: amplitude, correlation P, peak time). Also grab
   **Table S2** (the functional-category / gene-list table) — useful for the
   `functional_category` column and for sanity-checking the 78.
4. Download the .xls/.xlsx files to `data/raw/`.

### Route B — PMC (if ASM gates you)
1. Search: `Kushige 2013 genome-wide heterocyst circadian Anabaena PMC`.
2. This paper has a PMC record (JB deposits to PMC). Open it, go to
   **"Supplementary Materials" / "Associated Data,"** download the dataset files.
   Solve the one-time reCAPTCHA if shown.

### Route C — last resort
- The raw arrays are in the **KEGG Expression Database** (accessions ex0001892–
  ex001947). **This database is effectively dead — do not plan around it.** If Data
  Set S1 doesn't have what you need and ASM/PMC both fail, email the authors
  (Iwasaki, or Shigeki Ehira who made the arrays) rather than chasing KEGG.

### THE CHECK (the moment it opens)
- **Full table** = ~5,336 rows with per-ORF amplitude / P / peak time. -> Version B
  possible on the *Anabaena* side.
- **Subset only** = ~78 rows. -> Version A only; note in hurdles.md (H7).
- Version A needs only the 78 with their indices — so you're covered regardless.

---

## After both are down
1. Record, in `hurdles.md`, one line per file: full-table or subset? (This is the H7
   provenance finding, captured live.)
2. If BOTH have full tables -> tell Mike Version B is on the table as a phase two.
3. If EITHER is subset-only -> Version B is foreclosed; that's a written result, and
   Version A proceeds. Do not burn the timeline trying to reconstruct full tables.
4. Copy raw -> interim, and start `01_parse_*.py` against the interim copies.

## Provenance note to keep
For each downloaded file, jot in `data/raw/SOURCES.txt`: the URL you got it from, the
date, and the exact filename as published. Future-you (and the metagenomics pipeline
this is rehearsing) will want to know precisely where each number came from — that's
H7 in practice.
