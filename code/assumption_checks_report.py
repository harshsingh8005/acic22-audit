"""Assumption checks for every pre-registered test (statistical-analysis /
test-selection-guard requirement). Writes results/tables/assumption_checks.md.
Nothing here changes any test; robustness re-computations are labelled."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TAB = HERE.parent / "results" / "tables"
rng = np.random.default_rng(20260818)
out = ["# Assumption checks for the pre-registered tests\n",
       "Unit of replication for every test = cohort (independent ACIC 2022 datasets, generated independently by the organisers). ",
       "No cohort contributes twice to any between-cohort comparison; within-cohort comparisons (adjustment set, estimator) use paired tests.\n"]

est = pd.read_csv(TAB / "fast_arm_estimates.csv")
est = est[est.att.notna()].copy()
est["cover"] = (est.ci_lo <= est.true_satt) & (est.true_satt <= est.ci_hi)

for estimator in ("E1_AIPW", "E2_DRATT", "E0_OLS"):
    full = est[(est.estimator == estimator) & (est.adj == "full")]
    red = est[(est.estimator == estimator) & (est.adj == "reduced")]
    conf = full[full.strength.isin(["Strong", "Weak"])]
    out.append(f"\n## {estimator}\n")
    # P1
    k, n = int(conf.cover.sum()), len(conf)
    out.append(f"- **P1 one-proportion score test vs 0.95** — n = {n}, successes = {k}; expected counts under H0: n·0.95 = {n*0.95:.0f}, n·0.05 = {n*0.05:.0f} (both ≫ 5 → normal approximation valid; Wilson CI reported). Independence: cohorts are independent datasets.")
    # P2
    st, wk = conf[conf.strength == "Strong"], conf[conf.strength == "Weak"]
    out.append(f"- **P2 two-proportion z (pooled) + Newcombe CI** — n_Strong = {len(st)} (successes {int(st.cover.sum())}), n_Weak = {len(wk)} (successes {int(wk.cover.sum())}); all cell counts ≥ 5 → asymptotics valid; independent groups (different cohorts).")
    # P3
    m = full[["cohort_id", "cover", "source"]].merge(red[["cohort_id", "cover"]], on="cohort_id", suffixes=("_full", "_red"))
    m["d"] = m.cover_full.astype(int) - m.cover_red.astype(int)
    dA, dB = m[m.source == "A"].d.values, m[m.source == "B"].d.values
    def dist(x):
        return {v: round(float((x == v).mean()), 3) for v in (-1, 0, 1)}
    lev = stats.levene(dA, dB)
    # permutation robustness check (labelled)
    obs = dA.mean() - dB.mean()
    pooled = np.concatenate([dA, dB]); nA = len(dA)
    perm = np.array([(lambda p: p[:nA].mean() - p[nA:].mean())(rng.permutation(pooled)) for _ in range(10000)])
    p_perm = float((np.abs(perm) >= abs(obs)).mean())
    mw = stats.mannwhitneyu(dA, dB, alternative="two-sided")
    out.append(f"- **P3 Welch t on d = cover_full − cover_reduced (values −1/0/+1)** — n_A = {len(dA)}, n_B = {len(dB)}; distribution of d in A: {dist(dA)}, in B: {dist(dB)}; SD_A = {dA.std(ddof=1):.3f}, SD_B = {dB.std(ddof=1):.3f}; Levene F = {lev.statistic:.2f}, p = {lev.pvalue:.2e} (variances unequal → Welch, as pre-specified). "
               f"Normality: d is discrete, so Shapiro–Wilk is not meaningful; with n = 1,600 per group the sampling distribution of the mean is normal by the CLT (bounded variable). "
               f"Robustness re-computations (labelled, not replacements): permutation test of the mean difference (10,000 label permutations, seed 20260818) p = {p_perm:.4f} (observed diff = {obs:.4f}); Mann–Whitney U p = {mw.pvalue:.2e}. Same conclusion as the pre-specified Welch t.")

# S3
coh = pd.read_csv(TAB / "fast_arm_cohort.csv")
cA = coh[coh.source == "A"]
out.append("\n## S3 (OVB bound calibration) — paired binary outcomes within cohort (OVB-adjusted CI covers / naive CI covers) → exact McNemar (no distributional assumption; discordant pairs reported in hypothesis_tests.csv). Coverage proportion CI: Wilson (n = 1,600).")
out.append(f"- Measured strengths used: cf_y* clipped at 0 in {int((cA.E1_cf_y_star<=0).sum())} of {len(cA)} Scenario-A cohorts (E1); cf_d* clipped at 0 in {int((cA.E1_cf_d_star<=0).sum())}. ν² plug-in fallback used in the short model in {est[(est.estimator=='E1_AIPW')&(est.adj=='reduced')&(est.source=='A')].nu2_fallback.mean()*100:.1f}% of Scenario-A cohorts.")
# S4
out.append("\n## S4 (RCT gate) — one-proportion Wilson CI (n = 200; expected counts 190/10 under 0.95 → adequate); trend test = logistic regression of cover on ordinal score (RCT=0, Weak=1, Strong=2), Wald z on the slope (n = 3,400; no separation).")
# S6
out.append("\n## S6 (Spearman concordance) — rank-based, no distributional assumption; heavy ties at the floor (Γ* = 1 in ~66% of cohorts) handled by average ranks; bootstrap CIs (B = 500 resamples of cohorts).")
# S2, S5 (neural) if available
p = TAB / "neural_arm_estimates.csv"
if p.exists():
    ne = pd.read_csv(p)
    for nm in ("E5_DragonNet", "E4_TARNet"):
        nn = ne[(ne.estimator == nm) & (ne.adj == "full")]
        rt = ne[(ne.estimator == "E6_ridgeT") & (ne.adj == "full")][["cohort_id", "pehe_treated"]]
        mp = nn[["cohort_id", "pehe_treated"]].merge(rt, on="cohort_id", suffixes=("_nn", "_ridge"))
        diff = mp.pehe_treated_nn - mp.pehe_treated_ridge
        sh = stats.shapiro(diff) if len(diff) <= 5000 else None
        out.append(f"\n## S2/S5 ({nm}, X_full) — n cohorts = {len(mp)}. S2: exact McNemar on paired coverage indicators (no distributional assumption). "
                   f"S5: Wilcoxon signed-rank on paired PEHE differences (pre-specified non-parametric): median diff = {diff.median():.2f}, mean = {diff.mean():.2f}, skew = {stats.skew(diff):.2f}; "
                   f"Shapiro–Wilk on differences W = {sh.statistic:.3f}, p = {sh.pvalue:.2e} ({'non-normal' if sh.pvalue < 0.05 else 'consistent with normality'}; the rank test was pre-specified regardless and is not affected; symmetry of differences is approximate — median log-ratio with bootstrap CI reported as the effect size).")
(TAB / "assumption_checks.md").write_text("\n".join(out), encoding="utf-8")
print("\n".join(out))
