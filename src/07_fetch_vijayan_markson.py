#!/usr/bin/env python3
"""
07_fetch_vijayan_markson.py — replace the Ito 2009 side with archived LL time series.

WHY THIS EXISTS
---------------
Ito 2009's Dataset S1 publishes only summary indices (amplitude, corr P, cycling
call, peak time) and stores a literal "-" for the peak time of the 1,715 / 2,515
genes it does not call cycling. There is no raw time series to recompute from, so
joining on peak time reimposed a rhythmicity filter on the 7942 side (hurdles.md,
H8). This script pulls two archived LL time courses that fix that.

DATASETS
--------
GSE18902  Vijayan, Zuzow & O'Shea 2009, PNAS 106:22564   [PRIMARY]
          2x LD entrainment -> release into constant light (LL).
          Sampled every 4 h, T = 24 -> 84 h. 16 timepoints, 2.5 circadian cycles.
          (Samples GSM468479-485 are the novobiocin arm — EXCLUDED here.)

GSE52486  Markson et al. 2013, Cell 155:1396 (WT arm only)  [REPLICATE]
          GEO summary states verbatim that its goal was "to obtain a replicate of
          the wild-type LL circadian timecourse published in Vijayan et al". Same
          protocol, same platform. Sampled every 4 h, T = 36 -> 64 h, 7 timepoints
          (T = 52 h omitted by the submitters for poor data quality).
          NB: this is the WT SubSeries only. The dRpaA arms (GSE50908/50919/50920)
          are deliberately NOT pulled — that is the RpaA scope the brief closed.

GPL9534   Agilent-020846 S. elongatus PCC 7942 oligo array, 2,723 features.
          Its `ORF` column is already in `Synpcc7942_XXXX` form — the exact key
          used by data/processed/rbh_orthologs.csv. No probe->protein
          reconciliation needed. This is why integration is cheap.

TWO THINGS TO KNOW BEFORE TRUSTING THE OUTPUT
---------------------------------------------
1. These are TWO-COLOUR arrays. Each timepoint is hybridised against a pool of
   equal-mass RNA from all timepoints, so a value is log2(sample / time-averaged
   reference), NOT an intensity. CV is therefore meaningless here — the natural
   amplitude measure is the fitted cosine amplitude in log2 units.

   Consequence for H3: do NOT carry Kushige's CV across. Log2-transform Kushige's
   intensities and fit the SAME cosinor, then compare fitted log2 amplitudes.
   That is a genuinely commensurable quantity, which raw CV never was.

2. Vijayan is n=1 and Markson-WT is n=1, but they are independent biological
   replicates of each other on identical hardware. Their gene-wise correlation is
   a reliability estimate — which is exactly the input H4's disattenuation needs
   and which we previously recorded as permanently blocked.

USAGE
-----
    python src/07_fetch_vijayan_markson.py            # fetch + parse + fit
    python src/07_fetch_vijayan_markson.py --no-fetch # reuse existing downloads

Writes to data/raw/geo/, data/interim/, and data/processed/.
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RAW = REPO / "data" / "raw" / "geo"
INTERIM = REPO / "data" / "interim"
PROCESSED = REPO / "data" / "processed"

GEO_MATRIX = ("https://ftp.ncbi.nlm.nih.gov/geo/series/"
              "{stub}nnn/{acc}/matrix/{acc}_series_matrix.txt.gz")
GEO_PLATFORM = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
                "GPL9nnn/GPL9534/soft/GPL9534_family.soft.gz")

# GSE18902: keep only the continuous-light arm. The novobiocin samples share the
# series but are a different experiment entirely.
VIJAYAN_LL = {
    "GSM468463": 24, "GSM468464": 28, "GSM468465": 32, "GSM468466": 36,
    "GSM468467": 40, "GSM468468": 44, "GSM468469": 48, "GSM468470": 52,
    "GSM468471": 56, "GSM468472": 60, "GSM468473": 64, "GSM468474": 68,
    "GSM468475": 72, "GSM468476": 76, "GSM468477": 80, "GSM468478": 84,
}
MARKSON_LL = {
    "GSM1267750": 36, "GSM1267751": 40, "GSM1267752": 44, "GSM1267753": 48,
    "GSM1267754": 56, "GSM1267755": 60, "GSM1267756": 64,
}

PERIOD_H = 24.0


# --------------------------------------------------------------------------
# fetch
# --------------------------------------------------------------------------
def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  cached  {dest.name}")
        return dest
    print(f"  fetching {url}")
    with urllib.request.urlopen(url, timeout=120) as r:
        dest.write_bytes(r.read())
    print(f"  wrote   {dest.name}  ({dest.stat().st_size/1e6:.1f} MB)")
    return dest


def fetch_all() -> None:
    for acc in ("GSE18902", "GSE52486"):
        stub = acc[:-3]  # GSE18902 -> GSE18
        download(GEO_MATRIX.format(stub=stub, acc=acc), RAW / f"{acc}_series_matrix.txt.gz")
    download(GEO_PLATFORM, RAW / "GPL9534_family.soft.gz")

    (RAW / "SOURCES.txt").write_text(
        "GSE18902  Vijayan, Zuzow & O'Shea 2009 PNAS 106:22564 — LL circadian timecourse\n"
        "          https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18902\n"
        "GSE52486  Markson et al. 2013 Cell 155:1396, WT arm — stated replicate of GSE18902\n"
        "          https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52486\n"
        "GPL9534   Agilent-020846 S. elongatus PCC 7942 oligo array, 2,723 features\n"
        "          https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL9534\n"
        f"Retrieved {pd.Timestamp.today().date()} by src/07_fetch_vijayan_markson.py\n"
    )


# --------------------------------------------------------------------------
# parse
# --------------------------------------------------------------------------
def read_series_matrix(path: Path) -> pd.DataFrame:
    """Extract the expression block between the !series_matrix_table markers.

    Read as bytes and split on \\n only — the same \\r hazard that desynchronised
    the SI PDF parsing (hurdles.md, 2026-07-27) applies to any text GEO serves.
    """
    with gzip.open(path, "rb") as fh:
        raw = fh.read().decode("utf-8", errors="replace")
    lines = raw.split("\n")

    try:
        start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
        end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    except StopIteration:
        sys.exit(f"ERROR: no expression table found in {path.name}")

    block = "\n".join(lines[start + 1:end])
    df = pd.read_csv(io.StringIO(block), sep="\t", na_values=["", "null", "NA"])
    df = df.rename(columns={df.columns[0]: "probe_id"})
    df.columns = [c.strip('"') for c in df.columns]
    df["probe_id"] = df["probe_id"].astype(str).str.strip('"')
    return df


def verify_sample_map(path: Path, expected: dict[str, int], label: str) -> None:
    """Check the hardcoded GSM -> timepoint map against the titles GEO actually ships.

    The maps below were written from the GEO web records. That is an assumption, and a
    wrong sample-to-time assignment produces a perfectly plausible phase estimate with
    no error anywhere — the same silent-failure shape as the acrophase sign bug
    (hurdles.md, 2026-07-29). So assert it rather than trust it.

    Titles look like 'continuous light T = 24 hours experiment 1' (GSE18902) or
    'WT_LL_36h_replicate' (GSE52486); both yield the hour with one number-grab.
    """
    with gzip.open(path, "rb") as fh:
        head = fh.read(400_000).decode("utf-8", errors="replace")
    lines = head.split("\n")

    def row(tag):
        for l in lines:
            if l.startswith(tag):
                return [c.strip().strip('"') for c in l.split("\t")[1:]]
        return []

    accs, titles = row("!Sample_geo_accession"), row("!Sample_title")
    if not accs or len(accs) != len(titles):
        print(f"  {label}: could not read sample titles — MAP UNVERIFIED, check by hand")
        return

    seen = dict(zip(accs, titles))
    problems = []
    for gsm, hours in expected.items():
        title = seen.get(gsm)
        if title is None:
            problems.append(f"{gsm} absent from series")
            continue
        # First integer followed by h / hr / hours. Do NOT use \b after the 'h':
        # underscore is a word character, so \b never fires in 'WT_LL_36h_replicate'.
        # The negative lookahead also keeps 'T = 5 minutes' from matching.
        m = re.search(r"(\d+)\s*h(?:ours?|rs?)?(?![a-z])", title, flags=re.I)
        if not m:
            problems.append(f"{gsm} no hour in title {title!r}")
        elif int(m.group(1)) != hours:
            problems.append(f"{gsm} map says T={hours}h, GEO title says {title!r}")

    if problems:
        print(f"  {label}: SAMPLE MAP MISMATCH — refusing to continue")
        for p in problems:
            print(f"     - {p}")
        sys.exit(1)
    print(f"  {label}: sample map verified against GEO titles "
          f"({len(expected)} timepoints)")


def read_platform_map(path: Path) -> pd.DataFrame:
    """GPL9534 ID -> ORF (locus tag). The ORF column is already Synpcc7942_XXXX."""
    with gzip.open(path, "rb") as fh:
        raw = fh.read().decode("utf-8", errors="replace")
    lines = raw.split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin"))
    end = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    tbl = pd.read_csv(io.StringIO("\n".join(lines[start + 1:end])), sep="\t", dtype=str)

    tbl = tbl.rename(columns={"ID": "probe_id", "ORF": "locus_tag_7942"})
    tbl = tbl[["probe_id", "locus_tag_7942"]].dropna()
    tbl["probe_id"] = tbl["probe_id"].str.strip()
    tbl["locus_tag_7942"] = tbl["locus_tag_7942"].str.strip()

    # Sanity-check the key format before relying on it — this is the whole reason
    # this dataset was chosen over the alternatives. Real 7942 locus tags (confirmed
    # against every old_locus_tag in the RefSeq GFF) are Synpcc7942_ + an OPTIONAL
    # single uppercase prefix letter (B = the plasmid, 48 genes; R = a second replicon,
    # 53 genes) + 3-4 digits + an optional trailing lowercase letter (split/duplicate
    # loci, e.g. Synpcc7942_1912a). An earlier version of this regex required exactly 4
    # bare digits and silently rejected all 50 real B-prefixed plasmid genes -- caught
    # immediately by cross-checking RBH match counts before/after (1846 -> 1834,
    # investigated rather than accepted). Bare 9-digit numeric strings do NOT match
    # this pattern and are correctly excluded: 103 probes have ORF == their own probe
    # ID verbatim (e.g. "640711016"), which SPOT_ID reveals as "JGI_IG: <id>" -- JGI
    # intergenic-region tiling probes, not gene probes. GEO/Agilent leaves ORF
    # blank-equivalent for these by echoing the probe ID; they are correctly
    # locus-tag-less, not a parsing failure. Fixed here (2026-07-30) after they were
    # found riding silently through Part B/C's "full population" gene counts and
    # percentile rankings -- this function used to only PRINT the mismatch rate
    # without ever dropping the rows, so the diagnostic passed while the bad rows
    # shipped anyway.
    good = tbl["locus_tag_7942"].str.match(r"^Synpcc7942_[A-Z]?\d{3,4}[a-z]?$")
    n_before = len(tbl)
    print(f"  GPL9534: {n_before} probes, "
          f"{good.sum()} ({good.mean()*100:.1f}%) match Synpcc7942_####[letter]")
    if good.mean() < 0.9:
        print("  WARNING: locus-tag format is not what was expected — inspect before trusting.")
    tbl = tbl[good].copy()
    print(f"  dropped {n_before - len(tbl)} non-gene probes (JGI intergenic tiling probes, "
          f"ORF echoing probe ID) -- {len(tbl)} real gene probes remain")
    dup = tbl["locus_tag_7942"].duplicated().sum()
    print(f"  duplicated locus tags (ambiguous probe->gene): {dup}")
    return tbl


def to_long(matrix: pd.DataFrame, sample_times: dict[str, int],
            pmap: pd.DataFrame, label: str) -> pd.DataFrame:
    keep = [s for s in sample_times if s in matrix.columns]
    missing = sorted(set(sample_times) - set(keep))
    if missing:
        print(f"  {label}: WARNING missing samples {missing}")
    print(f"  {label}: {len(keep)} timepoints, {len(matrix)} probes")

    df = matrix[["probe_id"] + keep].merge(pmap, on="probe_id", how="inner")
    df = df.groupby("locus_tag_7942", as_index=False)[keep].mean()  # collapse replicate features
    long = df.melt(id_vars="locus_tag_7942", var_name="gsm", value_name="log2_ratio")
    long["time_h"] = long["gsm"].map(sample_times)
    long["dataset"] = label
    return long.dropna(subset=["log2_ratio"])


# --------------------------------------------------------------------------
# cosinor — ONE caller, applied identically to every dataset
# --------------------------------------------------------------------------
def cosinor(times: np.ndarray, values: np.ndarray, period: float = PERIOD_H) -> dict:
    """Least-squares fit of y = M + A*cos(2*pi*(t - phi)/period).

    Returns mesor, amplitude, acrophase (peak time in h, wrapped to [0, period)),
    R^2, and an F-test p-value. Deliberately the same model both papers describe,
    so results stay comparable to the published calls rather than replacing them
    with something incommensurable (hurdles.md, H6).

    NOTE: the p-value assumes independent residuals. On a 4 h-sampled series they
    are autocorrelated, so effective df is below n and p runs anticonservative —
    the exact caveat logged under H6. Use it for ranking, not as a calibrated FDR.
    """
    n = len(times)
    if n < 4 or np.all(np.isnan(values)):
        return dict(mesor=np.nan, amplitude=np.nan, peak_time_h=np.nan,
                    r2=np.nan, p_value=np.nan, n_points=n)

    w = 2 * np.pi / period
    X = np.column_stack([np.ones(n), np.cos(w * times), np.sin(w * times)])
    try:
        beta, *_ = np.linalg.lstsq(X, values, rcond=None)
    except np.linalg.LinAlgError:
        return dict(mesor=np.nan, amplitude=np.nan, peak_time_h=np.nan,
                    r2=np.nan, p_value=np.nan, n_points=n)

    mesor, b, g = beta
    amplitude = float(np.hypot(b, g))
    # y = M + A*cos(w*(t - phi)) expands to b*cos(wt) + g*sin(wt) with
    # b = A*cos(w*phi), g = A*sin(w*phi)  =>  phi = atan2(g, b) / w.
    # (Sign matters: atan2(-g, b) offsets every peak time by a constant, which is
    # invisible within one dataset and corrupts every cross-dataset delta-phi.)
    peak = float(np.arctan2(g, b) / w) % period
    # A peak at 0 h can come back as period - 1e-14 through the modulo; snap it,
    # otherwise a gene sitting on the bin edge lands in the wrong 4 h bin.
    if period - peak < 1e-9:
        peak = 0.0

    fitted = X @ beta
    ss_res = float(np.sum((values - fitted) ** 2))
    ss_tot = float(np.sum((values - values.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    p = np.nan
    # r2 can round to exactly 1.0 on a noiseless series, which would divide by
    # zero below; treat that as a perfect fit (p -> 0) rather than crashing.
    if n > 3 and np.isfinite(r2):
        if r2 >= 1.0:
            p = 0.0
        elif ss_res > 0:
            f_stat = (r2 / 2) / ((1 - r2) / (n - 3))
            try:
                from scipy import stats
                p = float(stats.f.sf(f_stat, 2, n - 3))
            except ImportError:
                print("  (scipy unavailable — p-values skipped)", file=sys.stderr)

    return dict(mesor=float(mesor), amplitude=amplitude, peak_time_h=float(peak),
                r2=r2, p_value=p, n_points=n)


def fit_dataset(long: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for locus, grp in long.groupby("locus_tag_7942"):
        g = grp.sort_values("time_h")
        res = cosinor(g["time_h"].to_numpy(float), g["log2_ratio"].to_numpy(float))
        res["locus_tag_7942"] = locus
        rows.append(res)
    out = pd.DataFrame(rows)
    out["amp_pct"] = out["amplitude"].rank(pct=True)   # H3: rank within platform
    out = out.add_suffix(f"_{label}").rename(columns={f"locus_tag_7942_{label}": "locus_tag_7942"})
    print(f"  {label}: fitted {len(out)} genes, "
          f"{out[f'p_value_{label}'].lt(0.05).sum()} with p < 0.05")
    return out


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fetch", action="store_true", help="reuse existing downloads")
    args = ap.parse_args()

    for d in (RAW, INTERIM, PROCESSED):
        d.mkdir(parents=True, exist_ok=True)

    if not args.no_fetch:
        print("[1/5] fetching from GEO")
        fetch_all()

    print("[2/5] platform map")
    pmap = read_platform_map(RAW / "GPL9534_family.soft.gz")

    print("[3/5] expression matrices")
    vij_path = RAW / "GSE18902_series_matrix.txt.gz"
    mar_path = RAW / "GSE52486_series_matrix.txt.gz"
    verify_sample_map(vij_path, VIJAYAN_LL, "vijayan")
    verify_sample_map(mar_path, MARKSON_LL, "markson")
    vij = to_long(read_series_matrix(vij_path), VIJAYAN_LL, pmap, "vijayan")
    mar = to_long(read_series_matrix(mar_path), MARKSON_LL, pmap, "markson")
    pd.concat([vij, mar]).to_csv(INTERIM / "7942_ll_timeseries_long.csv", index=False)

    print("[4/5] cosinor fits")
    fv = fit_dataset(vij, "vijayan")
    fm = fit_dataset(mar, "markson")
    merged = fv.merge(fm, on="locus_tag_7942", how="outer")

    # --- H4: the reliability estimate we recorded as permanently blocked -----
    both = merged.dropna(subset=["amplitude_vijayan", "amplitude_markson"])
    if len(both) > 10:
        r_amp = both["amplitude_vijayan"].corr(both["amplitude_markson"], method="spearman")
        dphi = ((both["peak_time_h_vijayan"] - both["peak_time_h_markson"] + 12) % 24) - 12
        print(f"\n  --- replicate reliability (n={len(both)}) ---")
        print(f"  amplitude, Spearman rho     : {r_amp:+.3f}")
        print(f"  peak time, median |delta|   : {dphi.abs().median():.2f} h")
        print(f"  within one 4 h sampling bin : {(dphi.abs() <= 4).mean()*100:.1f}% of genes")
        print(f"  attenuation factor sqrt(r)  : {np.sqrt(max(r_amp, 0)):.3f}"
              "   <- H4 disattenuation input")
        merged["repro_dphi_h"] = np.nan
        merged.loc[both.index, "repro_dphi_h"] = dphi

    merged.to_csv(PROCESSED / "7942_ll_cosinor.csv", index=False)
    print(f"\n  wrote {PROCESSED/'7942_ll_cosinor.csv'}  ({len(merged)} genes)")

    print("[5/5] coverage against the existing RBH ortholog table")
    rbh_path = PROCESSED / "rbh_orthologs.csv"
    if rbh_path.exists():
        rbh = pd.read_csv(rbh_path)
        j = rbh.merge(merged, on="locus_tag_7942", how="left")
        have = j["peak_time_h_vijayan"].notna()
        print(f"  RBH pairs                              : {len(rbh)}")
        print(f"  with a Vijayan peak time               : {have.sum()}")
        print(f"  ...vs 565 under the Ito join (H8)      : "
              f"{have.sum() - 565:+d} genes recovered")
        j.to_csv(PROCESSED / "rbh_with_7942_ll.csv", index=False)
        print(f"  wrote {PROCESSED/'rbh_with_7942_ll.csv'}")
    else:
        print(f"  SKIP — {rbh_path} not found; run the Version B notebook first.")

    print("\nNEXT: log2-transform Kushige's N+ intensities and run this same cosinor "
          "on them, so both sides carry a fitted log2 amplitude rather than two "
          "incommensurable CVs (see module docstring, note 1).")


if __name__ == "__main__":
    main()
