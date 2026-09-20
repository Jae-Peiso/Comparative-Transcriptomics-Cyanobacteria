"""Ortholog matching among the two rhythmic-gene sets (H2).

Method: gene-symbol name-matching, per the brief's own H2 guidance -- "for
Version A you only need orthologs among the rhythmic sets, so name-matching is
a legitimate primary method here." This is NOT reciprocal-best-BLAST-hit
(RBH): there are no proteome FASTA files in this repo, and standing up a real
RBH pipeline (proteome retrieval + alignment + reconciliation) is its own
sub-project, not a quick prerequisite. Every row below carries
ortholog_method="name-match" and rbh_confirmed=False explicitly, rather than
letting that distinction get lost downstream.
"""
from pathlib import Path

import pandas as pd

repo_root = Path(__file__).resolve().parents[1]
kushige = pd.read_csv(repo_root / "data" / "interim" / "kushige_rhythmic.csv")
ito = pd.read_csv(repo_root / "data" / "interim" / "ito_rhythmic.csv")
out_path = repo_root / "data" / "processed" / "name_matched_orthologs.csv"

kushige = kushige.dropna(subset=["gene_symbol"]).copy()
ito = ito.dropna(subset=["gene_symbol"]).copy()
kushige["gene_symbol_norm"] = kushige["gene_symbol"].str.lower()
ito["gene_symbol_norm"] = ito["gene_symbol"].str.lower()

merged = ito.merge(
    kushige, on="gene_symbol_norm", how="inner", suffixes=("_7942", "_anabaena")
)
merged["ortholog_method"] = "name-match"
merged["rbh_confirmed"] = False

print(f"Kushige rhythmic set: {len(kushige)} named genes (of 78 total)")
print(f"Ito rhythmic set:     {len(ito)} named genes (of 95 total)")
print(f"Name-matched pairs:   {len(merged)}")

cols = [
    "gene_symbol_norm", "gene_id_7942", "locus_tag", "annotation_7942",
    "peak_time_h_7942", "corr_p_7942", "amplitude_cv_7942",
    "gene_id_anabaena", "annotation_anabaena",
    "peak_time_h_anabaena", "corr_p_anabaena", "amplitude_cv_anabaena",
    "ortholog_method", "rbh_confirmed",
]
merged = merged.rename(columns={"gene_symbol_norm": "gene_symbol"})[
    ["gene_symbol"] + [c for c in cols if c != "gene_symbol_norm"]
]
out_path.parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(out_path, index=False)
print(f"Wrote {len(merged)} rows to {out_path}")
print(sorted(merged["gene_symbol"].tolist()))
