"""Addendum analysis for the install-gated legs (frozen plan §4 E3/E8; §6 S1/S3 apply to E3;
S2/S5 replication with catenets). Reads results/tables/e3e8_*.csv and catenets_arm_estimates.csv;
writes hypothesis_tests_addendum.csv, summary_e3e8_by_stratum.csv, e3e8_crosscheck.csv, key_numbers_addendum.json.
Same tests / same Holm-within-family rules as analyze.py; nothing tuned."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import confint_proportions_2indep, proportion_confint, proportions_ztest

HERE = Path(__file__).resolve().parent
TAB = HERE.parent / "results" / "tables"
rng = np.random.default_rng(20260818)


def wilson(k, n):
    lo, hi = proportion_confint(k, n, method="wilson"); return float(k / n), float(lo), float(hi)


def welch_ci(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = a.mean() - b.mean(); se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    df = (a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) ** 2 / ((a.var(ddof=1) / len(a)) ** 2 / (len(a) - 1) + (b.var(ddof=1) / len(b)) ** 2 / (len(b) - 1))
    t = stats.t.ppf(0.975, df); return float(d), float(d - t * se), float(d + t * se)


def paired_boot_diff(a, b, B=2000):
    a, b = np.asarray(a, float), np.asarray(b, float); n = len(a)
    vals = np.array([(a[i] - b[i]).mean() for i in (rng.integers(0, n, n) for _ in range(B))])
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def boot_ci(x, fn, B=2000):
    x = np.asarray(x); vals = np.array([fn(x[rng.integers(0, len(x), len(x))]) for _ in range(B)])
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


tests, key = [], {}


def add(fam, hid, plan_ref, desc, stat, p, effect, lo, hi, direction_ok, n, extra=None):
    tests.append({"family": fam, "id": hid, "plan_ref": plan_ref, "test": desc, "statistic": stat, "p_raw": p, "effect": effect,
                  "ci_lo": lo, "ci_hi": hi, "direction_as_predicted": direction_ok, "n": n, **(extra or {})})


est = pd.read_csv(TAB / "e3e8_estimates.csv"); coh = pd.read_csv(TAB / "e3e8_cohort.csv")
est = est[est.att.notna()].copy()
est["cover"] = (est.ci_lo <= est.true_satt) & (est.true_satt <= est.ci_hi)
est["width"] = est.ci_hi - est.ci_lo; est["bias"] = est.att - est.true_satt
est["confounded"] = est.strength.isin(["Strong", "Weak"])
est["ovb03_cover"] = (est.ovb03_ci_lo <= est.true_satt) & (est.true_satt <= est.ovb03_ci_hi) if "ovb03_ci_lo" in est else np.nan
key["n_cohorts_e3e8"] = int(coh.cohort_id.nunique())
key["e8_excluded_rows"] = int(pd.read_csv(TAB / "e3e8_estimates.csv").att.isna().sum())

# summary by stratum
def summ(g):
    return pd.Series({"n": len(g), "coverage": g.cover.mean(), "cov_lo": wilson(g.cover.sum(), len(g))[1], "cov_hi": wilson(g.cover.sum(), len(g))[2],
                      "bias_mean": g.bias.mean(), "rmse": np.sqrt((g.bias ** 2).mean()), "width_mean": g.width.mean(),
                      "se_ratio": g.bias.std(ddof=1) / g.se.mean(),
                      "ovb03_coverage": g.ovb03_cover.mean() if g.ovb03_ci_lo.notna().any() else np.nan})
parts = []
for keys in (["estimator", "adj"], ["estimator", "adj", "strength"], ["estimator", "adj", "strength", "source"], ["estimator", "adj", "source"]):
    s = est.groupby(keys, dropna=False).apply(summ).reset_index(); s["grouping"] = "|".join(keys); parts.append(s)
pd.concat(parts, ignore_index=True).to_csv(TAB / "summary_e3e8_by_stratum.csv", index=False)

# S1-family tests for E3 (and E8 for completeness)
for estimator, fam in (("E3_DML_ATTE", "S1_E3"), ("E8_EconML_DRLearner", "S1_E8")):
    full = est[(est.estimator == estimator) & (est.adj == "full")]; red = est[(est.estimator == estimator) & (est.adj == "reduced")]
    conf = full[full.confounded]
    k, n = int(conf.cover.sum()), len(conf); p_hat, lo, hi = wilson(k, n); z1, p1 = proportions_ztest(k, n, value=0.95, prop_var=0.95)
    add(fam, "P1", "§6 P1/S1", f"{estimator} X_full coverage < 0.95", z1, p1, p_hat, lo, hi, p_hat < 0.95, n)
    st, wk = conf[conf.strength == "Strong"], conf[conf.strength == "Weak"]
    k1, n1, k2, n2 = int(st.cover.sum()), len(st), int(wk.cover.sum()), len(wk)
    z2, p2 = proportions_ztest([k1, k2], [n1, n2]); dlo, dhi = confint_proportions_2indep(k1, n1, k2, n2, method="newcomb", compare="diff")
    add(fam, "P2", "§6 P2/S1", f"{estimator} coverage(Strong) − coverage(Weak) < 0", z2, p2, k1 / n1 - k2 / n2, dlo, dhi, (k1 / n1 - k2 / n2) < 0, n1 + n2, {"cov_strong": k1 / n1, "cov_weak": k2 / n2})
    m = full[["cohort_id", "cover", "source"]].merge(red[["cohort_id", "cover"]], on="cohort_id", suffixes=("_full", "_red"))
    m["d"] = m.cover_full.astype(int) - m.cover_red.astype(int); dA, dB = m[m.source == "A"].d, m[m.source == "B"].d
    t3, p3 = stats.ttest_ind(dA, dB, equal_var=False); d, lo3, hi3 = welch_ci(dA, dB)
    add(fam, "P3", "§6 P3/S1", f"{estimator} mean(cover_full−cover_reduced | A) − (· | B) > 0", t3, p3, d, lo3, hi3, d > 0, len(m),
        {"mean_d_A": float(dA.mean()), "mean_d_B": float(dB.mean()), "cov_full_A": float(m[m.source=="A"].cover_full.mean()), "cov_red_A": float(m[m.source=="A"].cover_red.mean()),
         "cov_full_B": float(m[m.source=="B"].cover_full.mean()), "cov_red_B": float(m[m.source=="B"].cover_red.mean())})
    # RCT gate
    r = full[full.strength == "RCT"]; k, n = int(r.cover.sum()), len(r); p_hat, lo, hi = wilson(k, n)
    add("S4_addendum", f"S4_{estimator}", "§6 S4", f"{estimator} X_full RCT coverage (Wilson)", np.nan, np.nan, p_hat, lo, hi, hi < 0.95, n)

# S3 with DoubleML's own benchmark (E3): Scenario A, X_reduced, OVB CI at (cf_y_bm, cf_d_bm, rho=1)
cA = coh[(coh.source == "A") & coh.E3_red_ovbstar_ci_lo.notna()].copy()
cA["ovbstar_cover"] = (cA.E3_red_ovbstar_ci_lo <= cA.true_satt) & (cA.true_satt <= cA.E3_red_ovbstar_ci_hi)
naive = est[(est.estimator == "E3_DML_ATTE") & (est.adj == "reduced")][["cohort_id", "cover"]].rename(columns={"cover": "naive_cover"})
cA = cA.merge(naive, on="cohort_id")
k, n = int(cA.ovbstar_cover.sum()), len(cA); p_hat, lo, hi = wilson(k, n)
b_ = int((cA.ovbstar_cover & ~cA.naive_cover).sum()); c_ = int((~cA.ovbstar_cover & cA.naive_cover).sum()); mc = mcnemar([[0, b_], [c_, 0]], exact=True)
bm_zero = cA["E3_bm_zero"].fillna(False).astype(bool) if "E3_bm_zero" in cA else pd.Series(False, index=cA.index)
add("S3_addendum", "S3_E3_dml_benchmark", "§6 S3 (DoubleML-native)", "E3 X_reduced Scenario A: DoubleML sensitivity_analysis at its own sensitivity_benchmark(['pre_trend']) strengths (rho=1) covers truth; Wilson lower ≥ 0.93; McNemar vs naive",
    float(mc.statistic), float(mc.pvalue), p_hat, lo, hi, lo >= 0.93, n,
    {"naive_coverage_same_cohorts": float(cA.naive_cover.mean()), "n_ovb_only": b_, "n_naive_only": c_,
     "median_cf_y_bm": float(cA.E3_cf_y_bm.median()), "median_cf_d_bm": float(cA.E3_cf_d_bm.median()), "frac_bm_unusable(cf=0 or 1)": float(bm_zero.mean()),
     "ovbstar_cover_where_bm_usable": float(cA[~bm_zero].ovbstar_cover.mean()) if (~bm_zero).any() else np.nan,
     "n_bm_usable": int((~bm_zero).sum())})
key["E3_benchmark_by_stratum"] = coh.groupby(["strength", "source"], dropna=False)[["E3_cf_y_bm", "E3_cf_d_bm", "E3_rho_bm", "E3_delta_theta_bm"]].median().reset_index().to_dict("records")

# software cross-check: my OVB layer on DoubleML nuisances vs DoubleML's own RV/RVa/CI (X_full)
e3 = est[est.estimator == "E3_DML_ATTE"].copy()
cc = e3[["cohort_id", "adj", "rv", "my_rv_on_dml", "rva", "my_rva_on_dml", "ovb03_ci_lo", "my_ovb03_ci_lo", "ovb03_ci_hi", "my_ovb03_ci_hi", "dml_sigma2", "my_sigma2_on_dml", "dml_nu2", "my_nu2_on_dml", "my_nu2_fallback"]].dropna(subset=["my_rv_on_dml"])
cc["d_rv"] = (cc.rv - cc.my_rv_on_dml).abs(); cc["d_rva"] = (cc.rva - cc.my_rva_on_dml).abs()
cc["d_ci_lo"] = (cc.ovb03_ci_lo - cc.my_ovb03_ci_lo).abs(); cc["d_ci_hi"] = (cc.ovb03_ci_hi - cc.my_ovb03_ci_hi).abs()
cc["d_sigma2_rel"] = (cc.dml_sigma2 - cc.my_sigma2_on_dml).abs() / cc.dml_sigma2; cc["d_nu2_rel"] = (cc.dml_nu2 - cc.my_nu2_on_dml).abs() / cc.dml_nu2.abs()
cc.to_csv(TAB / "e3e8_crosscheck.csv", index=False)
key["crosscheck_vs_doubleml"] = {"n": int(len(cc)), "median_abs_dRV": float(cc.d_rv.median()), "p90_abs_dRV": float(cc.d_rv.quantile(.9)), "max_abs_dRV": float(cc.d_rv.max()),
                                 "median_abs_dRVa": float(cc.d_rva.median()), "median_abs_dCI_lo": float(cc.d_ci_lo.median()), "p90_abs_dCI_lo": float(cc.d_ci_lo.quantile(.9)),
                                 "median_rel_dsigma2": float(cc.d_sigma2_rel.median()), "median_rel_dnu2": float(cc.d_nu2_rel.median()), "p90_rel_dnu2": float(cc.d_nu2_rel.quantile(.9)),
                                 "frac_my_nu2_fallback": float(cc.my_nu2_fallback.mean()),
                                 "spearman_rv": float(stats.spearmanr(cc.rv, cc.my_rv_on_dml).statistic)}
# EconML sensitivity failure rate
e8 = est[est.estimator == "E8_EconML_DRLearner"]
key["E8_econml_sensitivity_error_rate"] = {adj: float(e8[e8.adj == adj]["econml_sens_error"].notna().mean()) if "econml_sens_error" in e8 else np.nan for adj in ("full", "reduced")}
key["E8_econml_rv_ate_gt0_rate(X_full)"] = float((e8[(e8.adj == "full")]["rv_econml_ate"] > 0).mean()) if "rv_econml_ate" in e8 else np.nan

# catenets replication (S2/S5), X_full
p = TAB / "catenets_arm_estimates.csv"
if p.exists():
    ne = pd.read_csv(p); ne = ne[ne.adj == "full"].copy()
    ne["dr_cover"] = (ne.dr_ci_lo <= ne.true_satt) & (ne.true_satt <= ne.dr_ci_hi); ne["dr_width"] = ne.dr_ci_hi - ne.dr_ci_lo
    ne["ens_cover"] = (ne.ens_ci_lo <= ne.true_satt) & (ne.true_satt <= ne.ens_ci_hi)
    e1 = ne[ne.estimator == "E1_AIPW"][["cohort_id", "dr_cover", "dr_width", "dr_att", "true_satt"]]
    ridge = pd.read_csv(TAB / "neural_arm_estimates.csv"); ridge = ridge[(ridge.estimator == "E6_ridgeT") & (ridge.adj == "full")][["cohort_id", "pehe_treated"]]
    pyt = pd.read_csv(TAB / "neural_arm_estimates.csv"); pyt = pyt[pyt.adj == "full"]
    key["n_cohorts_catenets"] = int(ne.cohort_id.nunique())
    for nm, pyname in (("C5_catenets_DragonNet", "E5_DragonNet"), ("C4_catenets_TARNet", "E4_TARNet")):
        nn = ne[ne.estimator == nm]; m = nn.merge(e1, on="cohort_id", suffixes=("_nn", "_e1"))
        b_ = int((m.dr_cover_e1 & ~m.dr_cover_nn).sum()); c_ = int((~m.dr_cover_e1 & m.dr_cover_nn).sum()); mc = mcnemar([[0, b_], [c_, 0]], exact=True)
        diff = float(m.dr_cover_nn.mean() - m.dr_cover_e1.mean()); lo, hi = paired_boot_diff(m.dr_cover_nn.astype(float), m.dr_cover_e1.astype(float))
        add("S2_catenets", f"S2_{nm}", "§6 S2 (catenets replication)", f"{nm} neural-DR coverage < E1, same cohorts, X_full (McNemar; paired bootstrap CI)", float(mc.statistic), float(mc.pvalue), diff, lo, hi, diff < 0, len(m),
            {"cov_neural_dr": float(m.dr_cover_nn.mean()), "cov_e1": float(m.dr_cover_e1.mean()), "width_neural_dr": float(m.dr_width_nn.mean()), "width_e1": float(m.dr_width_e1.mean()),
             "rmse_neural_plugin": float(np.sqrt(((m.att - m.true_satt_nn) ** 2).mean())), "rmse_neural_dr": float(np.sqrt(((m.dr_att_nn - m.true_satt_nn) ** 2).mean())), "rmse_e1": float(np.sqrt(((m.dr_att_e1 - m.true_satt_nn) ** 2).mean())),
             "cov_seed_ensemble_interval(exploratory)": float(m.ens_cover.mean()), "n_discordant": b_ + c_})
        mp = nn[["cohort_id", "pehe_treated"]].merge(ridge, on="cohort_id", suffixes=("_nn", "_ridge"))
        w = stats.wilcoxon(mp.pehe_treated_nn, mp.pehe_treated_ridge); lr = np.log(mp.pehe_treated_nn / mp.pehe_treated_ridge); lo, hi = boot_ci(lr.values, np.median)
        add("S5_catenets", f"S5_{nm}", "§6 S5 (catenets replication)", f"{nm} PEHE-on-treated vs ridge T-learner, X_full (Wilcoxon; median log-ratio)", float(w.statistic), float(w.pvalue), float(np.median(lr)), lo, hi, None, len(mp),
            {"median_pehe_neural": float(mp.pehe_treated_nn.median()), "median_pehe_ridge": float(mp.pehe_treated_ridge.median())})
        # agreement with the PyTorch implementation
        pp = pyt[pyt.estimator == pyname][["cohort_id", "att", "pehe_treated", "dr_att"]]
        mm = nn[["cohort_id", "att", "pehe_treated", "dr_att"]].merge(pp, on="cohort_id", suffixes=("_cat", "_pyt"))
        key[f"{nm}_vs_pytorch"] = {"n": int(len(mm)), "spearman_plugin_att": float(stats.spearmanr(mm.att_cat, mm.att_pyt).statistic),
                                   "median_pehe_catenets": float(mm.pehe_treated_cat.median()), "median_pehe_pytorch": float(mm.pehe_treated_pyt.median()),
                                   "mean_abs_diff_plugin_att": float((mm.att_cat - mm.att_pyt).abs().mean())}

tdf = pd.DataFrame(tests); tdf["p_holm"] = np.nan
for fam, idx in tdf.groupby("family").groups.items():
    tdf.loc[idx, "p_holm"] = multipletests(tdf.loc[idx, "p_raw"].fillna(1.0).values, method="holm")[1]
def verdict(r):
    if r.family.startswith("S3"): return "CONFIRMED (criterion met)" if r.direction_as_predicted else "NOT MET"
    if r.family.startswith("S4"): return "RCT under-covers" if r.direction_as_predicted else "RCT coverage CI reaches 0.95"
    if r.direction_as_predicted is None or (isinstance(r.direction_as_predicted, float) and np.isnan(r.direction_as_predicted)): return "reported (no direction pre-declared)"
    return "CONFIRMED" if (r.p_holm < 0.05 and r.direction_as_predicted) else "NOT CONFIRMED"
tdf["verdict"] = tdf.apply(verdict, axis=1)
tdf.to_csv(TAB / "hypothesis_tests_addendum.csv", index=False)
(TAB / "key_numbers_addendum.json").write_text(json.dumps(key, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o)))
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(tdf[["family", "id", "effect", "ci_lo", "ci_hi", "p_raw", "p_holm", "n", "verdict"]].to_string())
print(json.dumps(key, indent=1, default=str)[:3000])
