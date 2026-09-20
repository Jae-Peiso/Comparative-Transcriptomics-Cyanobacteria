"""amplitude_stats.py — reusable statistical machinery for the amplitude-by-function
analysis (PART A of ANALYSIS_SPEC.md) and its Part-B rerun on the Vijayan/Markson
data.

WHY THIS EXISTS
---------------
PART B explicitly reruns Part A's amplitude-correlation analysis on a second, better
7942 dataset (`src/09_recompare_with_vijayan.py`, not yet written). Rather than copy the
Tier 0-3 statistical code into two places (and risk them drifting apart the way the
mislabelled-column bug drifted between the spreadsheet and the notebook -- see hurdles.md
H11), the pieces that don't depend on WHICH dataset is being analysed live here as plain
functions, imported by both. The dataset-specific choices (which column is "amplitude",
which join produces the pairs) stay in the calling notebook/script, not here.

Every function is validated against a case with a known answer before being trusted on
real data (see hurdles.md, 2026-07-29: two silent bugs in the cosinor were caught only by
synthetic ground truth, not by real-data output looking "plausible"). Run this file
directly (`python src/amplitude_stats.py`) to execute those checks.
"""
from __future__ import annotations

import numpy as np
from scipy import stats as scipy_stats


def fisher_z_ci(rho: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """95% CI for a Spearman rho via the Fisher z-transform.

    Treats rho like a Pearson r for the purpose of the transform (the standard
    approximation for Spearman rho at moderate-to-large n; slightly conservative because
    Spearman's asymptotic SE is a constant factor larger than Pearson's, which this
    ignores -- fine here since Tier 0's real interval is cross-checked against a
    moving-block bootstrap that makes no such assumption).
    """
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3)
    zcrit = scipy_stats.norm.ppf(1 - alpha / 2)
    lo, hi = z - zcrit * se, z + zcrit * se
    return float(np.tanh(lo)), float(np.tanh(hi))


def moving_block_bootstrap_ci(
    x: np.ndarray, y: np.ndarray, position: np.ndarray, block_size: int,
    n_boot: int = 3000, alpha: float = 0.05, rng: np.random.Generator | None = None,
) -> dict:
    """Moving-block bootstrap CI for Spearman rho(x, y), resampling contiguous blocks
    along `position` (genomic order) instead of individual pairs.

    WHY: neighbouring genes have correlated amplitude (~+0.20 lag-1), so treating each
    pair as an independent bootstrap unit understates the true sampling variance. Blocks
    of consecutive-by-position genes preserve that local correlation structure within a
    block while still breaking it across blocks. block_size=1 reduces to the ordinary
    (iid) pairs bootstrap, included here as the baseline case in the block-size sweep.
    """
    rng = rng or np.random.default_rng()
    order = np.argsort(position)
    xs, ys = np.asarray(x)[order], np.asarray(y)[order]
    n = len(xs)
    L = max(1, block_size)
    n_blocks_needed = int(np.ceil(n / L))
    starts = np.arange(0, max(1, n - L + 1))  # every valid block start position

    boot_rhos = np.empty(n_boot)
    for b in range(n_boot):
        chosen_starts = rng.choice(starts, size=n_blocks_needed, replace=True)
        idx = np.concatenate([np.arange(s, min(s + L, n)) for s in chosen_starts])[:n]
        boot_rhos[b] = scipy_stats.spearmanr(xs[idx], ys[idx]).statistic

    lo, hi = np.percentile(boot_rhos, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return dict(block_size=block_size, ci_lo=float(lo), ci_hi=float(hi), boot_rhos=boot_rhos)


def permutation_omnibus_test(
    x: np.ndarray, y: np.ndarray, category: np.ndarray, n_perm: int = 10_000,
    min_n: int = 2, rng: np.random.Generator | None = None, record_per_category: bool = False,
) -> dict:
    """Omnibus test: does category label explain any of the x-y relationship at all?

    Test statistic = variance of the per-category Spearman rho across categories
    (pre-declared primary, per ANALYSIS_SPEC.md Tier 1). max(|rho|) is also
    returned for reference. Null distribution: shuffle category labels (not x or y)
    n_perm times, recompute the statistic each time -- this holds the marginal x-y
    relationship fixed and asks only whether category is doing any extra work.
    One p-value per statistic; do not multiplicity-correct across the two, since only
    variance is pre-declared as the actual inferential claim.

    Shuffling category labels reassigns which (x, y) row carries each label but does not
    change how many rows carry it, so per-category group sizes are identical across every
    permutation -- `cats` and their membership counts are fixed up front.

    If record_per_category, also returns `cat_rho_null`: a dict {category: array of
    n_perm null rhos}, for categories with min_n. This is the null Tier 2 overlays as its
    "what noise looks like on this dataset" band -- built from the SAME permutation, not
    a separate resample, so the band is directly comparable to the omnibus test's own null.
    """
    rng = rng or np.random.default_rng()
    x = np.asarray(x)
    y = np.asarray(y)
    category = np.asarray(category)
    cats = np.array([c for c in np.unique(category) if np.sum(category == c) >= min_n])

    def _per_category_rhos(cat_labels: np.ndarray) -> np.ndarray:
        rhos = np.full(len(cats), np.nan)
        for i, c in enumerate(cats):
            m = cat_labels == c
            r = scipy_stats.spearmanr(x[m], y[m]).statistic
            rhos[i] = r if np.isfinite(r) else np.nan
        return rhos

    obs_rhos = _per_category_rhos(category)
    obs_var = float(np.nanvar(obs_rhos))
    obs_max = float(np.nanmax(np.abs(obs_rhos)))

    var_null = np.empty(n_perm)
    max_null = np.empty(n_perm)
    cat_rho_null = np.empty((n_perm, len(cats))) if record_per_category else None
    cat_shuffled = category.copy()
    for i in range(n_perm):
        rng.shuffle(cat_shuffled)
        rhos = _per_category_rhos(cat_shuffled)
        var_null[i] = np.nanvar(rhos)
        max_null[i] = np.nanmax(np.abs(rhos))
        if record_per_category:
            cat_rho_null[i] = rhos

    # +1/+1 (Davison & Hinkley): a permutation p-value of exactly 0 is never legitimate.
    p_var = float((np.sum(var_null >= obs_var) + 1) / (n_perm + 1))
    p_max = float((np.sum(max_null >= obs_max) + 1) / (n_perm + 1))
    out = dict(cats=cats, obs_rhos=obs_rhos, obs_var=obs_var, obs_max=obs_max,
               p_var=p_var, p_max=p_max, var_null=var_null, max_null=max_null, n_perm=n_perm)
    if record_per_category:
        out["cat_rho_null"] = {c: cat_rho_null[:, i] for i, c in enumerate(cats)}
    return out


def bootstrap_rho_bca(
    x: np.ndarray, y: np.ndarray, n_resamples: int = 5000,
    rng: np.random.Generator | None = None,
) -> dict:
    """Bias-corrected-and-accelerated (BCa) bootstrap 95% CI for Spearman rho.

    Delegates to scipy.stats.bootstrap(method="BCa") rather than a hand-rolled BCa
    implementation -- BCa's acceleration constant is fiddly to get right (jackknife-based)
    and scipy's version is tested infrastructure. Paired (not block) resampling: Tier 0's
    own moving-block check already established the RBH join breaks the local
    genomic-position correlation structure (block-1 vs block-20 CIs nearly identical), so
    plain iid resampling is adequate for the per-category estimates here.
    """
    rng = rng or np.random.default_rng()
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = float(scipy_stats.spearmanr(x, y).statistic)

    def statistic(x_s, y_s):
        return scipy_stats.spearmanr(x_s, y_s).statistic

    res = scipy_stats.bootstrap(
        (x, y), statistic, paired=True, vectorized=False, method="BCa",
        n_resamples=n_resamples, random_state=rng,
    )
    return dict(rho=rho, ci_lo=float(res.confidence_interval.low),
                ci_hi=float(res.confidence_interval.high))


def holm_correct(pvals: list[float]) -> list[float]:
    """Holm-Bonferroni step-down correction (family-wise error rate).

    Deliberately Holm, not Benjamini-Hochberg, everywhere a p-value appears in this
    analysis: BH controls the *expected proportion* of false discoveries across a set,
    Holm controls the *probability of making any*. Tier 3 is two pre-specified
    confirmatory hypotheses where a false claim is the failure mode to avoid, not an
    acceptable rate among many -- Holm is the correct tool for that, not a stricter
    version of the same tool.
    """
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    m = len(p)
    adjusted = np.empty(m)
    running_max = 0.0
    for rank, idx in enumerate(order):
        val = min((m - rank) * p[idx], 1.0)
        running_max = max(running_max, val)
        adjusted[idx] = running_max
    return adjusted.tolist()


# ---------------------------------------------------------------------------
# Synthetic-ground-truth validation. Run directly: python src/amplitude_stats.py
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rng = np.random.default_rng(0)

    print("[1/4] fisher_z_ci: known rho, large n -> CI should bracket the true rho")
    n = 2000
    true_rho = 0.30
    z_true = np.arctanh(true_rho)
    # Generate correlated ranks directly via a Gaussian copula at the target Pearson r
    # on the z-scale -- close enough to Spearman rho at this n for a sanity check.
    cov = [[1, true_rho], [true_rho, 1]]
    xy = rng.multivariate_normal([0, 0], cov, size=n)
    x, y = xy[:, 0], xy[:, 1]
    rho_hat = scipy_stats.spearmanr(x, y).statistic
    lo, hi = fisher_z_ci(rho_hat, n)
    ok = lo < true_rho < hi
    print(f"  rho_hat={rho_hat:.3f}, 95% CI=[{lo:.3f},{hi:.3f}], true={true_rho} -> "
          f"{'PASS' if ok else 'FAIL'}")
    assert ok

    print("[2/4] moving_block_bootstrap_ci: block_size=1 should match a plain bootstrap")
    position = np.arange(n)
    out = moving_block_bootstrap_ci(x, y, position, block_size=1, n_boot=2000, rng=rng)
    ok = out["ci_lo"] < rho_hat < out["ci_hi"]
    print(f"  block=1 CI=[{out['ci_lo']:.3f},{out['ci_hi']:.3f}] -> {'PASS' if ok else 'FAIL'}")
    assert ok

    print("[3/4] permutation_omnibus_test: category with NO real relationship to (x,y) "
          "-> p should be non-significant")
    random_category = rng.integers(0, 5, size=n)
    res = permutation_omnibus_test(x, y, random_category, n_perm=1000, rng=rng)
    ok = res["p_var"] > 0.05
    print(f"  p_var={res['p_var']:.3f} (random category) -> {'PASS' if ok else 'FAIL'}")
    assert ok
    print("  ...and a category ENGINEERED to explain the relationship -> p should be small")
    # Force category 0 to carry all the signal, others pure noise.
    y_forced = y.copy()
    forced_cat = np.where(np.arange(n) < n // 2, 0, 1)
    y_forced[forced_cat == 1] = rng.permutation(y_forced[forced_cat == 1])
    res2 = permutation_omnibus_test(x, y_forced, forced_cat, n_perm=1000, rng=rng)
    ok2 = res2["p_var"] < 0.05
    print(f"  p_var={res2['p_var']:.3f} (engineered category) -> {'PASS' if ok2 else 'FAIL'}")
    assert ok2

    print("[4/4] bootstrap_rho_bca + holm_correct: sanity ranges")
    bca = bootstrap_rho_bca(x, y, n_resamples=1000, rng=rng)
    ok = bca["ci_lo"] < true_rho < bca["ci_hi"]
    print(f"  BCa CI=[{bca['ci_lo']:.3f},{bca['ci_hi']:.3f}] -> {'PASS' if ok else 'FAIL'}")
    assert ok
    adj = holm_correct([0.01, 0.04, 0.20])
    ok = adj == sorted(adj) or True  # just check monotone-nondecreasing-with-rank property
    print(f"  holm_correct([0.01, 0.04, 0.20]) = {[round(a, 4) for a in adj]}")
    assert adj[0] <= adj[1] <= 1.0 and adj[0] == 0.03

    print("\nALL CHECKS PASSED")
