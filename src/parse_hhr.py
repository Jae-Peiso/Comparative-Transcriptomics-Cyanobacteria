"""Parse HHsuite .hhr output files (HHpred/HHsearch) into a tidy per-hit table.

One row per query-hit pair, matching the schema already used in
data/processed/homology_evidence.csv: target id/description, probability,
E-value, and BOTH query and template residue ranges (multidomain proteins get
different real hits on different regions -- collapsing to one row per gene
throws that away, same reasoning as the Pfam pass).

Parses the detailed per-hit blocks (">target_id description" + the
"Probab=... E-value=... Score=..." stats line + the Q/T alignment blocks),
not the fixed-width summary table at the top of the file -- the detailed
blocks use stable "key=value" and "Q <name> <start> ... <end>" patterns that
don't depend on exact column byte-offsets, which do vary across HHsuite
versions.
"""
import re
from collections import Counter
from pathlib import Path

import pandas as pd

HIT_HEADER_RE = re.compile(r"^No\s+(\d+)\s*$")
TARGET_RE = re.compile(r"^>(\S+)\s*(.*)$")
STATS_RE = re.compile(
    r"Probab=([\d.]+)\s+E-value=([\d.eE+-]+)\s+Score=([\d.]+)\s+"
    r"Aligned_cols=(\d+)\s+Identities=(\d+)%\s+Similarity=([\-\d.]+)\s+"
    r"Sum_probs=([\d.]+)(?:\s+Template_Neff=([\d.]+))?"
)
QLINE_RE = re.compile(r"^Q\s+(?!Consensus\b)(\S+)\s+(\d+)\s+[A-Za-z.\-~]+\s+(\d+)")
TLINE_RE = re.compile(r"^T\s+(?!Consensus\b)(\S+)\s+(\d+)\s+[A-Za-z.\-~]+\s+(\d+)")


def parse_hhr(path) -> pd.DataFrame:
    """Parse one .hhr file into one row per hit. Returns empty DataFrame (with the
    right columns) if no hits are found, rather than raising -- some genes may
    legitimately come back with zero results above HHpred's own reporting floor.
    """
    lines = Path(path).read_text().splitlines()
    query_name = None
    for line in lines:
        if line.startswith("Query"):
            query_name = line.split(None, 1)[1].strip()
            break

    hits = []
    i, n = 0, len(lines)
    while i < n:
        m = HIT_HEADER_RE.match(lines[i])
        if not m:
            i += 1
            continue
        hit_no = int(m.group(1))
        i += 1
        while i < n and not lines[i].startswith(">"):
            i += 1
        if i >= n:
            break
        tm = TARGET_RE.match(lines[i])
        target_id, target_desc = tm.group(1), tm.group(2).strip()
        i += 1

        stats = None
        while i < n and lines[i].strip() != "":
            sm = STATS_RE.search(lines[i])
            i += 1
            if sm:
                stats = sm
                break
        if stats is None:
            continue
        prob, evalue, score, cols, ident, sim, sumprobs, neff = stats.groups()

        q_starts, q_ends, t_starts, t_ends = [], [], [], []
        while i < n and not HIT_HEADER_RE.match(lines[i]):
            qm = QLINE_RE.match(lines[i])
            if qm:
                q_starts.append(int(qm.group(2)))
                q_ends.append(int(qm.group(3)))
            tmatch = TLINE_RE.match(lines[i])
            if tmatch:
                t_starts.append(int(tmatch.group(2)))
                t_ends.append(int(tmatch.group(3)))
            i += 1
        if not q_starts or not t_starts:
            continue

        hits.append(dict(
            query_name=query_name, hit_rank=hit_no, target_id=target_id,
            target_description=target_desc, probability=float(prob),
            evalue=float(evalue), score=float(score), aligned_cols=int(cols),
            identity_pct=int(ident), similarity=float(sim), sum_probs=float(sumprobs),
            query_residue_range=f"{min(q_starts)}-{max(q_ends)}",
            target_residue_range=f"{min(t_starts)}-{max(t_ends)}",
        ))
    cols_order = ["query_name", "hit_rank", "target_id", "target_description",
                  "probability", "evalue", "score", "aligned_cols", "identity_pct",
                  "similarity", "sum_probs", "query_residue_range", "target_residue_range"]
    return pd.DataFrame(hits, columns=cols_order)


def guess_target_database(target_id: str) -> str:
    """Best-guess database from target-ID format alone -- HHpred mixes hits from
    every database you checked into one ranked list and doesn't label each hit's
    source database in the .hhr file. Pattern-based, not authoritative: PDB IDs
    are 4 chars + underscore + chain (e.g. 6xyz_A), SCOPe domain IDs start with
    'd' (e.g. d1dlwa1), Pfam accessions start with PF. Flags anything else as
    'unknown (check manually)' rather than guessing wrong silently.
    """
    if re.match(r"^PF\d{5}", target_id):
        return "Pfam-A"
    if re.match(r"^[0-9][A-Za-z0-9]{3}_[A-Za-z0-9]+$", target_id):
        return "PDB_mmCIF70"
    if re.match(r"^d[0-9][A-Za-z0-9]{3}[A-Za-z0-9_]\d?$", target_id):
        return "SCOPe70"
    return "unknown (check manually)"


def parse_hhr_folder(folder) -> pd.DataFrame:
    """Parse every *.hhr in `folder` into a single combined table, one row per
    query-hit pair. locus_tag comes from the file's own internal "Query" line
    (query_name), not the filename -- robust to however the file happened to be
    named on download (e.g. a "hhpred_" prefix some browsers add automatically)."""
    folder = Path(folder)
    frames = []
    for f in sorted(folder.glob("*.hhr")):
        df = parse_hhr(f)
        df["locus_tag"] = df["query_name"]
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"no .hhr files found in {folder}")
    out = pd.concat(frames, ignore_index=True)
    out["target_database"] = out["target_id"].map(guess_target_database)
    return out


_DESC_STOPWORDS = {
    "protein", "domain", "family", "containing", "of", "unknown", "function", "putative",
    "like", "domains", "repeat", "system", "component", "subunit", "type", "fold", "conserved",
    "region", "the", "and", "with", "related", "complex", "binding", "structural", "genomics",
    "het", "chain", "uncharacterized", "structure", "proteins", "scop", "psi", "initiative",
    "genesis", "center", "joint", "nesg", "mcsg",
}


def _clean_description(text: str) -> str:
    text = re.sub(r"\{[^}]*\}", " ", text)              # strip {Organism name}
    text = re.sub(r"HET:\s*[A-Za-z0-9, ]+;?", " ", text)  # strip ligand codes
    text = re.sub(r"\d+\.?\d*\s*[AÅ]\b", " ", text)      # strip crystal resolution, e.g. "2.4A"
    text = re.sub(r"SCOP:\s*\S+", " ", text)             # strip SCOP codes (not a description word)
    return text


def _description_keywords(text: str) -> set:
    text = _clean_description(text)
    return {w for w in re.findall(r"[a-zA-Z][a-zA-Z0-9]{2,}", text.lower()) if w not in _DESC_STOPWORDS}


def _range_overlap_frac(r1: str, r2: str) -> float:
    a0, a1 = map(int, r1.split("-"))
    b0, b1 = map(int, r2.split("-"))
    inter = max(0, min(a1, b1) - max(a0, b0))
    shorter = min(a1 - a0, b1 - b0)
    return inter / shorter if shorter > 0 else 0.0


def summarize_convergence(evidence: pd.DataFrame, min_prob: float = 50.0, top_n: int = 15,
                           min_corroborating: int = 3, consensus_frac_threshold: float = 0.4) -> pd.DataFrame:
    """One row per gene: does the top hit have independent corroboration, or is it a
    lone outlier? This is the question that actually matters for HHpred output, NOT
    "how high is the single best probability" and NOT "how many hits came back" --
    a query can return 250 hits above 50% probability and still be genuinely
    unresolved if they disagree with each other, and a single 99% hit with nothing
    else agreeing with it is weaker evidence than five independent structures that
    all land in the same place and describe the same thing.

    Method: take the top `top_n` hits at or above `min_prob` probability, keep only
    those whose query_residue_range overlaps the single best hit's range by >=50%
    (reciprocal, on the shorter of the two spans) -- these are "corroborating" in the
    sense of addressing the same region of the query, not just scoring well anywhere
    in a big multidomain protein. Among those, extract description keywords (organism
    tags, ligand codes, resolution, and SCOP-code text stripped out first -- they are
    metadata, not evidence) and find the single most common keyword. If it's shared by
    >= consensus_frac_threshold of the corroborating hits AND there are at least
    min_corroborating of them: CONFIDENT CONSENSUS -- report the shared theme, not any
    one hit's specific protein identity (the individual top hit is very likely not an
    ortholog, just the highest-ranked example of a fold/family that recurs across many
    unrelated proteins -- exactly the all0232/Marf1 caveat, generalized). A single
    high-probability hit with no corroboration is flagged separately, not silently
    trusted just because its own probability is high. Anything else -- few
    corroborating hits, or corroborating hits that don't agree with each other --
    is UNRESOLVED, even if the raw top probability looks high: several genuinely
    different families competing for the same region is a real result (see alr4939
    below), not a parsing failure.
    """
    rows = []
    for gene, sub in evidence.groupby("locus_tag"):
        sub = sub.sort_values("probability", ascending=False).reset_index(drop=True)
        top = sub.iloc[0]
        topN = sub[sub["probability"] >= min_prob].head(top_n)
        overlap_mask = topN["query_residue_range"].apply(
            lambda r: _range_overlap_frac(r, top["query_residue_range"]) >= 0.5)
        corroborating = topN[overlap_mask]

        kwsets = [_description_keywords(d) for d in corroborating["target_description"]]
        counts = Counter(w for s in kwsets for w in s)
        if counts and len(corroborating) >= 2:
            best_kw, best_n = counts.most_common(1)[0]
            consensus_frac = best_n / len(corroborating)
        else:
            best_kw, consensus_frac = None, 0.0

        if len(corroborating) >= min_corroborating and consensus_frac >= consensus_frac_threshold:
            verdict = "CONFIDENT CONSENSUS"
        elif len(corroborating) <= 2 and top["probability"] >= 90:
            verdict = "SINGLE STRONG HIT (uncorroborated)"
        else:
            verdict = "UNRESOLVED (no convergent theme)"

        rows.append(dict(
            locus_tag=gene, n_hits_total=len(sub), top_probability=round(float(top["probability"]), 1),
            top_hit_id=top["target_id"], top_hit_description=top["target_description"],
            n_corroborating=len(corroborating), consensus_keyword=best_kw,
            consensus_fraction=round(consensus_frac, 2), verdict=verdict,
        ))
    return pd.DataFrame(rows).sort_values("top_probability", ascending=False).reset_index(drop=True)


if __name__ == "__main__":
    # self-test against a synthetic .hhr snippet built from the documented HHsuite
    # format, so a parsing bug is caught here rather than on the first real file.
    import tempfile

    synthetic = """Query         all0232
Match_columns 320
No_of_seqs    1 out of 1
Neff          1.0
Searched_HMMs 92000
Date          Thu Jul 30 00:00:00 2026
Command       hhsearch -i all0232.a3m -d pdb70

 No Hit                               Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 6xyz_A Marf1 NYN domain            98.2 3.1E-21 2.6E-25  120.5  85.3  126    15-140     3-128 (130)

No 1
>6xyz_A Marf1 NYN domain, X-ray structure
Probab=98.24  E-value=3.1e-21  Score=120.54  Aligned_cols=126  Identities=25%  Similarity=0.456  Sum_probs=110.2  Template_Neff=8.500

Q all0232               15 MSTKVAAAGGWETC   80 (320)
Q Consensus              15 ~~~~~~~~~~~~~~   80 (320)
                             |||+.+++++|.|
T Consensus               3 ~~~~~~~~~~~~~~   68 (130)
T 6xyz_A                   3 msTkVaaagGwetc   68 (130)


Q all0232               81 DEFGHIKLMNPQRST  140 (320)
Q Consensus              81 ~~~~~~~~~~~~~~~  140 (320)
                             |||+.+++++|.|.
T Consensus               69 ~~~~~~~~~~~~~~~  128 (130)
T 6xyz_A                   69 defghiklmnpqrst  128 (130)

"""
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "all0232.hhr"
        p.write_text(synthetic)
        result = parse_hhr(p)
        assert len(result) == 1, f"expected 1 hit, got {len(result)}"
        row = result.iloc[0]
        assert row["target_id"] == "6xyz_A"
        assert abs(row["probability"] - 98.24) < 1e-6
        assert abs(row["evalue"] - 3.1e-21) < 1e-30
        assert row["query_residue_range"] == "15-140", row["query_residue_range"]
        assert row["target_residue_range"] == "3-128", row["target_residue_range"]
        assert guess_target_database("6xyz_A") == "PDB_mmCIF70"
        assert guess_target_database("PF01936") == "Pfam-A"
        assert guess_target_database("d1dlwa1") == "SCOPe70"
        print("self-test passed:")
        print(result.to_string(index=False))
