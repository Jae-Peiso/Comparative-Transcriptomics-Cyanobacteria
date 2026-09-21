# Rhythmic gene comparison

Comparative transcriptomics workflow contrasting rhythmic gene expression in *Synechococcus elongatus* PCC 7942 (Ito et al. 2009) and *Anabaena* sp. PCC 7120 (Kushige et al. 2013).

## Project goal

Characterize analysis hurdles while producing one defensible quantitative comparison of rhythmic genes across the two studies.

## Repository structure

- `PROJECT_BRIEF_rhythmic_gene_comparison.md` — project scope, analysis plan, and hurdle register
- `ANALYSIS_SPEC.md` — pre-registered specification for the current round of analyses (Part 0 correction, then Parts A–F)
- `DATASETS_vijayan_markson.md` — notes on the Vijayan 2009 / Markson 2013 GEO datasets used in Parts B and C
- `DOWNLOADS.md` — manual download instructions for the supplementary information
- `hurdles.md` / `hurdles.pdf` — working log of analysis issues and lessons learned
- `journal_club_kushige_2013.md` — journal club notes on Kushige et al. 2013
- `runningNotes.md` — running notes on dataset and method decisions
- `git_workflow.py` — git automation helper for saving and pushing analysis updates
- `data/interim/` — parsed and joined intermediate tables (rhythmicity calls, percentiles, BLAST reciprocal-best-hit results)
- `data/manual/` — hand-downloaded inputs no script can regenerate (HHpred queries and results)
- `data/processed/` — final tidy tables and comparison outputs used in the figures
- `notebooks/` — analysis notebooks, including `ComparativeTranscriptomics.ipynb` (exploratory) and `VersionB_AllOrthologs.ipynb` (the main ortholog comparison)
- `src/` — parsing, ortholog-matching, and statistics scripts

## Workflow notes

1. Download the supplementary materials described in `DOWNLOADS.md`.
2. Parse the datasets into a common tidy schema (`src/01_parse_kushige.py`, `src/02_parse_ito.py`).
3. Compare rhythmic genes following the plan in `PROJECT_BRIEF_rhythmic_gene_comparison.md` and `ANALYSIS_SPEC.md`.
