"""Parse Kushige2013 Table S2 (the 78 'cycling' genes) into the tidy schema.

Table S2 lives in the SI PDF (Kushige2013/Kushige2013_Info.pdf), NOT in the
Data Set S1 spreadsheet (Kushige2013_so2.xlsx) -- that spreadsheet is the full
unfiltered ~5,336-ORF raw time-series table and has no amplitude/P-value
columns at all. Table S2 is the authors' own curated, already-scored list, so
this uses their reported peak time / P-value / amplitude directly rather than
recomputing CV from the raw arrays (which risks not matching their exact
per-replicate combination formula).

Source text extraction requires poppler's pdftotext (brew install poppler).
"""
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

repo_root = Path(__file__).resolve().parents[1]
pdf_path = repo_root / "Kushige2013" / "Kushige2013_Info.pdf"
out_path = repo_root / "data" / "interim" / "kushige_rhythmic.csv"

result = subprocess.run(
    ["pdftotext", "-layout", str(pdf_path), "-"], capture_output=True, check=True
)
# Split on '\n' only -- the PDF (PowerPoint-exported) contains stray '\r'
# characters from overstruck text that Python's universal-newline mode would
# otherwise split on too, desynchronizing line numbers from what grep/wc see.
lines = result.stdout.decode("utf-8", errors="replace").split("\n")

header_idx = next(i for i, l in enumerate(lines) if l.startswith("  ORF No.        Gene name"))
end_idx = next(i for i, l in enumerate(lines) if "TABLE S3." in l)
data_lines = [l for l in lines[header_idx + 1 : end_idx] if l.strip()]

pattern = re.compile(
    r"^\s*(?P<gene_id>\S+)\s+(?P<gene_symbol>\S+)\s+(?P<annotation>.*?)\s+"
    r"(?P<peak_time_h>[\d.]+)\s+(?P<corr_p>[\d.]+)\s+(?P<amplitude_cv>[\d.]+)\s+(?P<rhythmic_tier>[abc])\s+"
    r"(?P<peak_time_h_nminus>[\d.]+|-)\s+(?P<corr_p_nminus>[\d.]+|-)\s+(?P<amplitude_cv_nminus>[\d.]+|-)\s+(?P<rhythmic_tier_nminus>[abc]|AR)\s*$"
)

rows = []
for line in data_lines:
    m = pattern.match(line)
    if not m:
        print(f"WARNING: unparsed line: {line!r}", file=sys.stderr)
        continue
    rows.append(m.groupdict())

df = pd.DataFrame(rows)
if len(df) != 78:
    print(f"WARNING: parsed {len(df)} rows, expected 78 (Table S2's own headline count)", file=sys.stderr)

df["gene_symbol"] = df["gene_symbol"].replace("-", pd.NA)
for col in ["peak_time_h", "corr_p", "amplitude_cv"]:
    df[col] = pd.to_numeric(df[col])
# rhythmic_flag: Table S2's headline "78 cycling genes" = tier a (stringent) + b
# (standard) under N+ combined; tier c (lenient, +114 more genes) is a separate,
# looser category the paper does NOT fold into the "78" figure -- but every row
# in Table S2 is already restricted to a/b, so rhythmic_flag is True throughout.
df["rhythmic_flag"] = True
df["dataset"] = "Kushige2013"

cols = [
    "dataset", "gene_id", "gene_symbol", "annotation", "peak_time_h", "corr_p",
    "amplitude_cv", "rhythmic_tier", "rhythmic_flag",
    "peak_time_h_nminus", "corr_p_nminus", "amplitude_cv_nminus", "rhythmic_tier_nminus",
]
df = df[cols]
out_path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out_path, index=False)
print(f"Wrote {len(df)} rows to {out_path}")
