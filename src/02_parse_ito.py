"""Filter Ito2009 Dataset S1 to the paper's own "97 higher-amplitude genes" set.

Threshold per the SI PDF text (Ito2009_supporting information.pdf, Table S2
footnote): P < 0.001 AND amplitude > 10^-0.6, in the wild-type strain. This is
NOT the same as the notebook's earlier "cycling_property != AR" flag, which
matches the SI's own LOOSER a/b-tier threshold and reproduces the 800-gene
population reported in the SI's Table S1 -- a different, larger set.

NOTE: this filter reproduces 95 genes, not the paper's stated 97. Close but
not exact -- likely a rounding/boundary-inclusivity difference in how the
threshold was applied, or the per-replicate P/amplitude combination formula
(SI Methods: p_value = max(1st, 2nd); amplitude = min(1st, 2nd)) not being
exactly reconstructable from the single combined columns Dataset S1 provides.
Logged as a hurdle rather than silently forcing the count to match.
"""
from pathlib import Path

import pandas as pd

repo_root = Path(__file__).resolve().parents[1]
xls_path = repo_root / "Ito2009" / "Ito2009_sd1.xls"
out_path = repo_root / "data" / "interim" / "ito_rhythmic.csv"

ito_df = pd.read_excel(xls_path, sheet_name="Table S1", header=5)

P_THRESH = 0.001
AMP_THRESH = 10**-0.6  # ~0.2512

high_amp = ito_df[
    (ito_df["Correlation p-value*"] < P_THRESH) & (ito_df["Amplitude*"] > AMP_THRESH)
].copy()

print(f"P<{P_THRESH} & amplitude>{AMP_THRESH:.4f} -> {len(high_amp)} genes (paper reports 97)")

high_amp = high_amp.rename(columns={
    "orfID": "gene_id",
    "7942ID": "locus_tag",
    "Gene Name": "gene_symbol",
    "Annotation": "annotation",
    "Correlation p-value*": "corr_p",
    "Amplitude*": "amplitude_cv",
    "Peak time*": "peak_time_h",
    "Cycling property**": "rhythmic_tier",
})
high_amp["rhythmic_flag"] = True
high_amp["dataset"] = "Ito2009"

cols = [
    "dataset", "gene_id", "locus_tag", "gene_symbol", "annotation",
    "peak_time_h", "corr_p", "amplitude_cv", "rhythmic_tier", "rhythmic_flag",
]
high_amp = high_amp[cols]
out_path.parent.mkdir(parents=True, exist_ok=True)
high_amp.to_csv(out_path, index=False)
print(f"Wrote {len(high_amp)} rows to {out_path}")
