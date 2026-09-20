"""Redefine "rhythmic" per-platform by amplitude PERCENTILE RANK within each
dataset's own distribution, instead of the papers' shared raw CV/amplitude
number (confirmed non-comparable: 10^-0.6 sits at the 90th percentile of
Ito2009's full amplitude distribution but only the 19th percentile of
Kushige2013's -- Kushige's raw CVs run ~2.7x higher across the board, a
cross-platform noise-floor difference, not a biological one; see hurdles.md H3).

This does NOT touch each platform's own P-value significance criterion --
only the amplitude side, since that's the part shown to be platform-scale-
dependent. Ito keeps P<0.001 (its own headline criterion). Kushige uses
P<0.05 (its own "standard" tier's P criterion) combined with the new
percentile-matched amplitude bar, computed from the FULL ~5336-gene table
(not just Table S2's 78, since a different cutoff surfaces different genes).
"""
from pathlib import Path

import numpy as np
import pandas as pd

repo_root = Path(__file__).resolve().parents[1]
PERCENTILE = 0.90

# --- Ito2009: full distribution already has per-gene amplitude/P/peak time ---
ito_df = pd.read_excel(repo_root / "Ito2009" / "Ito2009_sd1.xls", sheet_name="Table S1", header=5)
ito_amp_thresh = ito_df["Amplitude*"].quantile(PERCENTILE)
ito_rhythmic = ito_df[
    (ito_df["Amplitude*"] > ito_amp_thresh) & (ito_df["Correlation p-value*"] < 0.001)
].copy()
ito_rhythmic = ito_rhythmic.rename(columns={
    "orfID": "gene_id", "7942ID": "locus_tag", "Gene Name": "gene_symbol",
    "Annotation": "annotation", "Correlation p-value*": "corr_p",
    "Amplitude*": "amplitude_cv", "Peak time*": "peak_time_h",
})[["gene_id", "locus_tag", "gene_symbol", "annotation", "peak_time_h", "corr_p", "amplitude_cv"]]
ito_rhythmic["dataset"] = "Ito2009"

print(f"Ito2009 amplitude @ {PERCENTILE:.0%} percentile (own distribution): {ito_amp_thresh:.4f}")
print(f"Ito2009 revised rhythmic set: {len(ito_rhythmic)} genes (P<0.001, amplitude>{ito_amp_thresh:.3f})")

# --- Kushige2013: compute CV for the full raw table, find its own percentile ---
kushige_df = pd.read_excel(repo_root / "Kushige2013" / "Kushige2013_so2.xlsx", sheet_name="Sheet1", header=5)
ll_timepoints = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48]
phases = sorted({t % 24 for t in ll_timepoints})
nplus_cols = [f"N+_1st_LL{t}" for t in ll_timepoints] + [f"N+_2nd_LL{t}" for t in ll_timepoints]
vals = kushige_df[nplus_cols].astype(float)
kushige_df["amplitude_cv"] = vals.std(axis=1, ddof=1) / vals.mean(axis=1)
kushige_df["corr_p"] = pd.to_numeric(kushige_df["Unnamed: 7"].astype(str).str.strip("()"), errors="coerce")

kushige_amp_thresh = kushige_df["amplitude_cv"].quantile(PERCENTILE)
print(f"\nKushige2013 CV @ {PERCENTILE:.0%} percentile (own distribution): {kushige_amp_thresh:.4f}")

kushige_hits = kushige_df[
    (kushige_df["amplitude_cv"] > kushige_amp_thresh) & (kushige_df["corr_p"] < 0.05)
].copy()
print(f"Kushige2013 revised rhythmic set: {len(kushige_hits)} genes (P<0.05, CV>{kushige_amp_thresh:.3f})")

# Peak time only for the genes that actually pass the filter (not all 5336).
def peak_time_nplus(row):
    by_phase = {p: [] for p in phases}
    for t in ll_timepoints:
        p = t % 24
        by_phase[p].append(row[f"N+_1st_LL{t}"])
        by_phase[p].append(row[f"N+_2nd_LL{t}"])
    avg = {p: sum(v) / len(v) for p, v in by_phase.items()}
    return max(avg, key=avg.get)

kushige_hits["peak_time_h"] = kushige_hits.apply(peak_time_nplus, axis=1)
kushige_rhythmic = kushige_hits.rename(columns={
    "ORF No.": "gene_id", "Gene name": "gene_symbol", "Annotation": "annotation",
})[["gene_id", "gene_symbol", "annotation", "peak_time_h", "corr_p", "amplitude_cv"]]
kushige_rhythmic["dataset"] = "Kushige2013"

interim = repo_root / "data" / "interim"
interim.mkdir(parents=True, exist_ok=True)
ito_rhythmic.to_csv(interim / "ito_rhythmic_percentile.csv", index=False)
kushige_rhythmic.to_csv(interim / "kushige_rhythmic_percentile.csv", index=False)

# --- Name-match join ---
ito_named = ito_rhythmic.dropna(subset=["gene_symbol"]).copy()
kushige_named = kushige_rhythmic.dropna(subset=["gene_symbol"]).copy()
ito_named["sym"] = ito_named["gene_symbol"].str.lower()
kushige_named["sym"] = kushige_named["gene_symbol"].str.lower()

merged = ito_named.merge(kushige_named, on="sym", how="inner", suffixes=("_7942", "_anabaena"))
merged["ortholog_method"] = "name-match"
merged["rbh_confirmed"] = False

print(f"\nNamed genes: Ito {len(ito_named)}, Kushige {len(kushige_named)}")
print(f"Name-matched pairs: {len(merged)}")
print(sorted(merged["sym"].tolist()))
print("\nnarM in revised Ito set?", ito_df.loc[ito_df["Gene Name"].str.lower() == "narm", ["Correlation p-value*", "Amplitude*"]].to_string(index=False) if (ito_df["Gene Name"].str.lower() == "narm").any() else "present in raw table but check threshold below")
narm_row = ito_df[ito_df["Gene Name"].str.lower() == "narm"].iloc[0]
print(f"narM (Ito): P={narm_row['Correlation p-value*']}, amplitude={narm_row['Amplitude*']:.3f} vs required P<0.001 & amplitude>{ito_amp_thresh:.3f}")

out_path = repo_root / "data" / "processed" / "name_matched_orthologs_percentile.csv"
out_path.parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(out_path, index=False)
print(f"\nWrote {len(merged)} rows to {out_path}")
