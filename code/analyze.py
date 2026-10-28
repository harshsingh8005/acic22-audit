"""Pre-specified confirmatory analysis (frozen plan §6–§9, §13).

Reads results/tables/{fast_arm,neural_arm}_{estimates,cohort}.csv and writes
  results/tables/summary_by_stratum.csv        coverage/bias/RMSE/width/SE-ratio per estimator x adj x stratum
  results/tables/hypothesis_tests.csv          every pre-registered test with raw + Holm p, effect + 95% CI, verdict
  results/tables/sensitivity_summary.csv       E-value / RV / Gamma* / OVB-CI coverage / fallback rates by stratum
  results/tables/calibration_curve.csv         empirical coverage at nominal 50/80/90/95/99
  results/tables/key_numbers.json              the numbers quoted in the reports
Nothing here is tuned; every test maps to a plan line (column `plan_ref`).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import (confint_proportions_2indep, proportion_confint,
                                          proportions_ztest)

HERE = Path(__file__).resolve().parent
TAB = HERE.parent / "results" / "tables"
Z = {0.50: 0.6744897501960817, 0.80: 1.2815515655446004, 0.90: 1.6448536269514722,
     0.95: 1.959963984540054, 0.99: 2.5758293035489004}
BOOT_SEED = 20260818
B = 2000


def wilson(k, n):
    lo, hi = proportion_confint(k, n, method="wilson")
    return float(k / n), float(lo), float(hi)


def welch_ci(a, b, level=0.95):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = a.mean() - b.mean()
    se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    df = (a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) ** 2 / (
        (a.var(ddof=1) / len(a)) ** 2 / (len(a) - 1) + (b.var(ddof=1) / len(b)) ** 2 / (len(b) - 1))
    t = stats.t.ppf(0.5 + level / 2, df)
    return float(d), float(d - t * se), float(d + t * se), float(se), float(df)


def boot_ci(x, fn, rng, B=B):
    x = np.asarray(x)
    vals = np.array([fn(x[rng.integers(0, len(x), len(x))]) for _ in range(B)])
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def paired_boot_diff(a, b, rng, B=B):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = len(a)
    vals = np.array([(a[i] - b[i]).mean() for i in (rng.integers(0, n, n) for _ in range(B))])
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main():
    est = pd.read_csv(TAB / "fast_arm_estimates.csv")
    coh = pd.read_csv(TAB / "fast_arm_cohort.csv")
    nest = pd.read_csv(TAB / "neural_arm_estimates.csv") if (TAB / "neural_arm_estimates.csv").exists() else None
    ncoh = pd.read_csv(TAB / "neural_arm_cohort.csv") if (TAB / "neural_arm_cohort.csv").exists() else None
    key = {}
    rng = np.random.default_rng(BOOT_SEED)

    est = est[est["att"].notna()].copy()
    est["cover"] = (est["ci_lo"] <= est["true_satt"]) & (est["true_satt"] <= est["ci_hi"])
    est["width"] = est["ci_hi"] - est["ci_lo"]
    est["bias"] = est["att"] - est["true_satt"]
    est["z"] = est["bias"] / est["se"]
    for lv, z in Z.items():
        est[f"cover_{int(lv*100)}"] = est["z"].abs() <= z
    est["confounded"] = est["strength"].isin(["Strong", "Weak"])
    est["ovb03_cover"] = (est["ovb03_ci_lo"] <= est["true_satt"]) & (est["true_satt"] <= est["ovb03_ci_hi"])
    key["n_cohorts_fast"] = int(coh["cohort_id"].nunique())
    key["n_rows_fast"] = int(len(est))
    excl = pd.read_csv(TAB / "fast_arm_estimates.csv")
    key["excluded_rows"] = excl[excl["att"].isna()][["cohort_id", "adj", "estimator", "excluded_reason"]].to_dict("records") if "excluded_reason" in excl else []

    # ------------------------------------------------ summary by stratum
    def summ(g):
        return pd.Series({
            "n": len(g), "coverage": g["cover"].mean(),
            "cov_lo": wilson(g["cover"].sum(), len(g))[1], "cov_hi": wilson(g["cover"].sum(), len(g))[2],
            "bias_mean": g["bias"].mean(), "abs_bias_mean": g["bias"].abs().mean(),
            "rmse": np.sqrt((g["bias"] ** 2).mean()), "width_mean": g["width"].mean(),
            "se_mean": g["se"].mean(), "sd_err": g["bias"].std(ddof=1),
            "se_ratio": g["bias"].std(ddof=1) / g["se"].mean(),
            "ovb03_coverage": g["ovb03_cover"].mean() if g["ovb03_ci_lo"].notna().any() else np.nan,
            "nu2_fallback_rate": g["nu2_fallback"].mean() if "nu2_fallback" in g and g["nu2_fallback"].notna().any() else np.nan,
            "n_trimmed_mean": g["n_trimmed"].mean() if "n_trimmed" in g else np.nan,
        })
    strata = []
    for keys in (["estimator", "adj"], ["estimator", "adj", "strength"], ["estimator", "adj", "strength", "source"],
                 ["estimator", "adj", "source"], ["estimator", "adj", "cell"]):
        s = est.groupby(keys, dropna=False).apply(summ).reset_index()
        s["grouping"] = "|".join(keys)
        strata.append(s)
    summary = pd.concat(strata, ignore_index=True)
    summary.to_csv(TAB / "summary_by_stratum.csv", index=False)

    # calibration curve
    cal = []
    for (e_, a_, s_), g in est.groupby(["estimator", "adj", "strength"]):
        for lv in Z:
            cal.append({"estimator": e_, "adj": a_, "strength": s_, "nominal": lv,
                        "empirical": g[f"cover_{int(lv*100)}"].mean(), "n": len(g)})
    pd.DataFrame(cal).to_csv(TAB / "calibration_curve.csv", index=False)

    tests = []

    def add(fam, hid, plan_ref, desc, stat, p, effect, lo, hi, direction_ok, n, extra=None):
        tests.append({"family": fam, "id": hid, "plan_ref": plan_ref, "test": desc, "statistic": stat, "p_raw": p,
                      "effect": effect, "ci_lo": lo, "ci_hi": hi, "direction_as_predicted": direction_ok, "n": n,
                      **(extra or {})})

    # ------------------------------------------------ primary + S1 (per estimator families)
    for estimator, fam in (("E1_AIPW", "PRIMARY"), ("E2_DRATT", "S1_E2"), ("E0_OLS", "S1_E0")):
        full = est[(est.estimator == estimator) & (est.adj == "full")]
        red = est[(est.estimator == estimator) & (est.adj == "reduced")]
        conf = full[full["confounded"]]
        # P1
        k, n = int(conf["cover"].sum()), len(conf)
        p_hat, lo, hi = wilson(k, n)
        z1, p1 = proportions_ztest(k, n, value=0.95, prop_var=0.95)
        add(fam, "P1", "§6 P1", f"{estimator} X_full 95% coverage < 0.95 (score test vs 0.95; Wilson CI)", z1, p1,
            p_hat, lo, hi, p_hat < 0.95, n)
        # P2
        st, wk = conf[conf.strength == "Strong"], conf[conf.strength == "Weak"]
        k1, n1, k2, n2 = int(st["cover"].sum()), len(st), int(wk["cover"].sum()), len(wk)
        z2, p2 = proportions_ztest([k1, k2], [n1, n2])
        dlo, dhi = confint_proportions_2indep(k1, n1, k2, n2, method="newcomb", compare="diff")
        add(fam, "P2", "§6 P2", f"{estimator} X_full coverage(Strong) − coverage(Weak) < 0 (two-proportion z; Newcombe CI)",
            z2, p2, k1 / n1 - k2 / n2, dlo, dhi, (k1 / n1 - k2 / n2) < 0, n1 + n2,
            {"cov_strong": k1 / n1, "cov_weak": k2 / n2})
        # P3
        m = full[["cohort_id", "cover", "source"]].merge(red[["cohort_id", "cover"]], on="cohort_id", suffixes=("_full", "_red"))
        m["d"] = m["cover_full"].astype(int) - m["cover_red"].astype(int)
        dA, dB = m[m.source == "A"]["d"], m[m.source == "B"]["d"]
        t3, p3 = stats.ttest_ind(dA, dB, equal_var=False)
        d, lo3, hi3, se3, df3 = welch_ci(dA, dB)
        add(fam, "P3", "§6 P3", f"{estimator} mean(cover_full−cover_reduced | A) − mean(· | B) > 0 (Welch t)",
            t3, p3, d, lo3, hi3, d > 0, len(m), {"mean_d_A": float(dA.mean()), "mean_d_B": float(dB.mean()),
                                                 "cov_full_A": float(m[m.source=="A"]["cover_full"].mean()), "cov_red_A": float(m[m.source=="A"]["cover_red"].mean()),
                                                 "cov_full_B": float(m[m.source=="B"]["cover_full"].mean()), "cov_red_B": float(m[m.source=="B"]["cover_red"].mean())})

    # ------------------------------------------------ S3 OVB bound calibration (Scenario A, X_reduced, E1 and E2)
    for estimator, pref in (("E1_AIPW", "E1"), ("E2_DRATT", "E2")):
        cA = coh[coh.source == "A"].copy()
        cA = cA[cA[f"{pref}_red_ovbstar_ci_lo"].notna()]
        cA["ovbstar_cover"] = (cA[f"{pref}_red_ovbstar_ci_lo"] <= cA["true_satt"]) & (cA["true_satt"] <= cA[f"{pref}_red_ovbstar_ci_hi"])
        naive = est[(est.estimator == estimator) & (est.adj == "reduced")][["cohort_id", "cover"]].rename(columns={"cover": "naive_cover"})
        cA = cA.merge(naive, on="cohort_id")
        k, n = int(cA["ovbstar_cover"].sum()), len(cA)
        p_hat, lo, hi = wilson(k, n)
        tbl = pd.crosstab(cA["ovbstar_cover"], cA["naive_cover"])
        b_ = int(((cA["ovbstar_cover"]) & (~cA["naive_cover"])).sum()); c_ = int((~cA["ovbstar_cover"] & cA["naive_cover"]).sum())
        mc = mcnemar([[0, b_], [c_, 0]], exact=True)
        add("S3", f"S3_{pref}", "§6 S3", f"{estimator} X_reduced Scenario A: OVB-adjusted CI at measured (cf_y*,cf_d*,rho=1) covers truth; Wilson lower ≥ 0.93; McNemar vs naive CI",
            float(mc.statistic), float(mc.pvalue), p_hat, lo, hi, lo >= 0.93, n,
            {"naive_coverage_same_cohorts": float(cA["naive_cover"].mean()), "n_ovb_only": b_, "n_naive_only": c_,
             "median_cf_y_star": float(cA[f"{pref}_cf_y_star"].median()), "median_cf_d_star": float(cA[f"{pref}_cf_d_star"].median()),
             "median_rho_true": float(cA[f"{pref}_rho_true"].median()) if f"{pref}_rho_true" in cA else np.nan})
        # same in Scenario B (negative control) descriptive
        cB = coh[coh.source == "B"].copy()
        cB = cB[cB[f"{pref}_red_ovbstar_ci_lo"].notna()]
        cB["ovbstar_cover"] = (cB[f"{pref}_red_ovbstar_ci_lo"] <= cB["true_satt"]) & (cB["true_satt"] <= cB[f"{pref}_red_ovbstar_ci_hi"])
        key[f"S3_{pref}_scenarioB_ovbstar_cover"] = float(cB["ovbstar_cover"].mean())
        key[f"S3_{pref}_scenarioB_median_cf_y_star"] = float(cB[f"{pref}_cf_y_star"].median())
        key[f"S3_{pref}_scenarioB_median_cf_d_star"] = float(cB[f"{pref}_cf_d_star"].median())

    # ------------------------------------------------ S4 RCT gate + trend (E1, E2, E0)
    for estimator in ("E1_AIPW", "E2_DRATT", "E0_OLS"):
        full = est[(est.estimator == estimator) & (est.adj == "full")]
        r = full[full.strength == "RCT"]
        k, n = int(r["cover"].sum()), len(r)
        p_hat, lo, hi = wilson(k, n)
        # Cochran–Armitage trend RCT(0) < Weak(1) < Strong(2): logistic slope z
        import statsmodels.api as sm
        score = full["strength"].map({"RCT": 0, "Weak": 1, "Strong": 2}).astype(float)
        fit = sm.Logit(full["cover"].astype(int), sm.add_constant(score)).fit(disp=0)
        add("S4", f"S4_{estimator}", "§6 S4", f"{estimator} X_full RCT coverage (Wilson CI; upper<0.95 ⇒ non-confounding under-coverage); trend z = logistic slope RCT→Weak→Strong",
            float(fit.tvalues.iloc[1]), float(fit.pvalues.iloc[1]), p_hat, lo, hi, hi < 0.95, n,
            {"trend_slope_logit": float(fit.params.iloc[1]), "cov_RCT": p_hat,
             "cov_weak": float(full[full.strength == "Weak"]["cover"].mean()), "cov_strong": float(full[full.strength == "Strong"]["cover"].mean())})

    # ------------------------------------------------ S2, S5 (neural sub-sample)
    if nest is not None:
        nest = nest.copy()
        nest["dr_cover"] = (nest["dr_ci_lo"] <= nest["true_satt"]) & (nest["true_satt"] <= nest["dr_ci_hi"])
        nest["dr_width"] = nest["dr_ci_hi"] - nest["dr_ci_lo"]
        nest["ens_cover"] = (nest["ens_ci_lo"] <= nest["true_satt"]) & (nest["true_satt"] <= nest["ens_ci_hi"])
        nest["bias"] = nest["att"] - nest["true_satt"]
        key["n_cohorts_neural"] = int(nest["cohort_id"].nunique())
        for adj in ("full", "reduced"):
            e1 = nest[(nest.estimator == "E1_AIPW") & (nest.adj == adj)][["cohort_id", "dr_cover", "dr_width", "dr_att", "true_satt"]]
            for nm in ("E5_DragonNet", "E4_TARNet"):
                nn = nest[(nest.estimator == nm) & (nest.adj == adj)]
                m = nn.merge(e1, on="cohort_id", suffixes=("_nn", "_e1"))
                b_ = int((m["dr_cover_e1"] & ~m["dr_cover_nn"]).sum()); c_ = int((~m["dr_cover_e1"] & m["dr_cover_nn"]).sum())
                mc = mcnemar([[0, b_], [c_, 0]], exact=True)
                diff = float(m["dr_cover_nn"].mean() - m["dr_cover_e1"].mean())
                lo, hi = paired_boot_diff(m["dr_cover_nn"].astype(float), m["dr_cover_e1"].astype(float), rng)
                add("S2" if adj == "full" else "S2_reduced(exploratory)", f"S2_{nm}_{adj}", "§6 S2",
                    f"{nm} neural-DR interval coverage < E1 AIPW coverage, same cohorts, X_{adj} (exact McNemar; paired bootstrap CI)",
                    float(mc.statistic), float(mc.pvalue), diff, lo, hi, diff < 0, len(m),
                    {"cov_neural_dr": float(m["dr_cover_nn"].mean()), "cov_e1": float(m["dr_cover_e1"].mean()),
                     "width_neural_dr": float(m["dr_width_nn"].mean()), "width_e1": float(m["dr_width_e1"].mean()),
                     "rmse_neural_plugin": float(np.sqrt(((m["att"] - m["true_satt_nn"]) ** 2).mean())),
                     "rmse_neural_dr": float(np.sqrt(((m["dr_att_nn"] - m["true_satt_nn"]) ** 2).mean())),
                     "rmse_e1": float(np.sqrt(((m["dr_att_e1"] - m["true_satt_nn"]) ** 2).mean())),
                     "cov_seed_ensemble_interval(exploratory)": float(m["ens_cover"].mean()),
                     "n_discordant": b_ + c_})
                # S5 PEHE
                rt = nest[(nest.estimator == "E6_ridgeT") & (nest.adj == adj)][["cohort_id", "pehe_treated"]]
                mp = nn[["cohort_id", "pehe_treated"]].merge(rt, on="cohort_id", suffixes=("_nn", "_ridge"))
                w = stats.wilcoxon(mp["pehe_treated_nn"], mp["pehe_treated_ridge"])
                lr = np.log(mp["pehe_treated_nn"] / mp["pehe_treated_ridge"])
                lo, hi = boot_ci(lr.values, np.median, rng)
                add("S5" if adj == "full" else "S5_reduced(exploratory)", f"S5_{nm}_{adj}", "§6 S5",
                    f"{nm} PEHE-on-treated vs ridge T-learner, X_{adj} (paired Wilcoxon; median log-ratio + bootstrap CI); no direction pre-declared",
                    float(w.statistic), float(w.pvalue), float(np.median(lr)), lo, hi, None, len(mp),
                    {"median_pehe_neural": float(mp["pehe_treated_nn"].median()), "median_pehe_ridge": float(mp["pehe_treated_ridge"].median())})
        # DoWhy cross-check
        if ncoh is not None and "dowhy_ols_att" in ncoh:
            e0 = est[(est.estimator == "E0_OLS") & (est.adj == "full")][["cohort_id", "att", "rv", "rva"]]
            mm = ncoh.merge(e0, on="cohort_id")
            key["E7_dowhy_vs_E0_max_abs_diff_att"] = float((mm["dowhy_ols_att"] - mm["att"]).abs().max())
            key["E7_dowhy_vs_E0_max_abs_diff_rv"] = float((mm["dowhy_rv"] - mm["rv"]).abs().max())
            key["E7_dowhy_vs_E0_max_abs_diff_rva"] = float((mm["dowhy_rva"] - mm["rva"]).abs().max())
            key["E7_n"] = int(len(mm))
            if "dowhy_ipw_att" in mm:
                e1f = est[(est.estimator == "E1_AIPW") & (est.adj == "full")][["cohort_id", "att", "true_satt"]]
                mi = ncoh.merge(e1f, on="cohort_id")
                key["E7_dowhy_ipw_rmse(exploratory)"] = float(np.sqrt(((mi["dowhy_ipw_att"] - mi["true_satt"]) ** 2).mean()))
                key["E7_E1_rmse_same_cohorts"] = float(np.sqrt(((mi["att"] - mi["true_satt"]) ** 2).mean()))

    # ------------------------------------------------ S6 concordance + floors (X_full)
    conc = []
    for estimator in ("E1_AIPW", "E2_DRATT", "E0_OLS"):
        g = est[(est.estimator == estimator) & (est.adj == "full")]
        cols = [c for c in ("rv", "gamma_star", "e_pt") if c in g and g[c].notna().any()]
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                gg = g[[cols[i], cols[j]]].dropna()
                rho = stats.spearmanr(gg[cols[i]], gg[cols[j]]).statistic
                idx = np.arange(len(gg))
                bs = [stats.spearmanr(gg.iloc[s][cols[i]], gg.iloc[s][cols[j]]).statistic for s in (rng.integers(0, len(gg), len(gg)) for _ in range(500))]
                conc.append({"estimator": estimator, "pair": f"{cols[i]}~{cols[j]}", "spearman": rho,
                             "ci_lo": float(np.nanpercentile(bs, 2.5)), "ci_hi": float(np.nanpercentile(bs, 97.5)), "n": len(gg)})
        floors = {"rv_floor(<0.001)": float((g["rv"] < 1e-3).mean()) if "rv" in g else np.nan,
                  "gamma_floor(=1)": float((g["gamma_star"] <= 1.0).mean()) if "gamma_star" in g else np.nan,
                  "e_ci_floor(=1)": float((g["e_ci"] <= 1.0).mean()) if "e_ci" in g else np.nan,
                  "sig_at_gamma1(p<0.05)": float((g["p_gamma1"] < 0.05).mean()) if "p_gamma1" in g else np.nan}
        key[f"S6_floors_{estimator}"] = floors
    pd.DataFrame(conc).to_csv(TAB / "concordance.csv", index=False)

    # ------------------------------------------------ sensitivity summary by stratum
    sens_rows = []
    for (e_, a_, s_, src), g in est.groupby(["estimator", "adj", "strength", "source"], dropna=False):
        sens_rows.append({"estimator": e_, "adj": a_, "strength": s_, "source": src, "n": len(g),
                          "median_rv": g["rv"].median() if "rv" in g else np.nan,
                          "median_rva": g["rva"].median() if "rva" in g else np.nan,
                          "median_gamma_star": g["gamma_star"].median() if "gamma_star" in g else np.nan,
                          "frac_gamma_gt1": (g["gamma_star"] > 1).mean() if "gamma_star" in g else np.nan,
                          "median_e_pt": g["e_pt"].median(), "median_e_ci": g["e_ci"].median(),
                          "frac_e_ci_gt1": (g["e_ci"] > 1).mean(),
                          "ovb03_coverage": g["ovb03_cover"].mean() if g["ovb03_ci_lo"].notna().any() else np.nan,
                          "naive_coverage": g["cover"].mean(),
                          "nu2_fallback_rate": g["nu2_fallback"].mean() if "nu2_fallback" in g and g["nu2_fallback"].notna().any() else np.nan,
                          "median_n_pairs": g["n_pairs"].median() if "n_pairs" in g else np.nan,
                          "median_proxy_gamma": g["proxy_gamma"].median()})
    pd.DataFrame(sens_rows).to_csv(TAB / "sensitivity_summary.csv", index=False)

    # measured strengths + placebo by stratum (cohort table)
    ms = coh.groupby(["strength", "source"], dropna=False).agg(
        n=("cohort_id", "size"),
        E1_cf_y_star_median=("E1_cf_y_star", "median"), E1_cf_d_star_median=("E1_cf_d_star", "median"),
        E1_rho_hat_star_median=("E1_rho_hat_star", "median"), E1_rho_true_median=("E1_rho_true", "median"),
        E1_B_star_median=("E1_B_star", "median"),
        placebo_att_mean=("placebo_att", "mean"), placebo_att_sd=("placebo_att", "std"),
        placebo_se_mean=("placebo_se", "mean"),
        E1_full_minus_red_att_mean=("E1_full_att", lambda x: np.nan),
    ).reset_index()
    ms["placebo_att_ci_lo"] = ms["placebo_att_mean"] - 1.96 * ms["placebo_att_sd"] / np.sqrt(ms["n"])
    ms["placebo_att_ci_hi"] = ms["placebo_att_mean"] + 1.96 * ms["placebo_att_sd"] / np.sqrt(ms["n"])
    ms["E1_full_minus_red_att_mean"] = coh.assign(dd=coh["E1_full_att"] - coh["E1_red_att"]).groupby(["strength", "source"], dropna=False)["dd"].mean().values
    ms.to_csv(TAB / "measured_strength_placebo_by_stratum.csv", index=False)

    # ------------------------------------------------ Holm within families
    tdf = pd.DataFrame(tests)
    tdf["p_holm"] = np.nan
    for fam, idx in tdf.groupby("family").groups.items():
        p = tdf.loc[idx, "p_raw"].values
        tdf.loc[idx, "p_holm"] = multipletests(p, method="holm")[1]
    def verdict(r):
        if r["family"].startswith("S3"):
            return "CONFIRMED (criterion met)" if r["direction_as_predicted"] else "NOT MET"
        if r["family"].startswith("S4"):
            return "RCT under-covers (non-confounding under-coverage present)" if r["direction_as_predicted"] else "RCT coverage CI reaches 0.95"
        if r["direction_as_predicted"] is None:
            return "reported (no direction pre-declared)"
        return "CONFIRMED" if (r["p_holm"] < 0.05 and r["direction_as_predicted"]) else "NOT CONFIRMED"
    tdf["verdict"] = tdf.apply(verdict, axis=1)
    tdf.to_csv(TAB / "hypothesis_tests.csv", index=False)
    key["hypothesis_tests"] = tdf.to_dict("records")
    (TAB / "key_numbers.json").write_text(json.dumps(key, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o)))
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(tdf[["family", "id", "test", "effect", "ci_lo", "ci_hi", "p_raw", "p_holm", "n", "verdict"]].to_string())
    return tdf, summary, key


if __name__ == "__main__":
    main()
