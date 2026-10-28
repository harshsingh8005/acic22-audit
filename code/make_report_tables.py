"""Render the markdown tables used in sensitivity_report.md / findings_summary.md
directly from results/tables/*.csv so no number is transcribed by hand.
Writes results/tables/report_tables.md (sections keyed by name)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TAB = HERE.parent / "results" / "tables"


def md(df: pd.DataFrame, floatfmt=3) -> str:
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f":
            vals = df[c].dropna()
            if len(vals) and np.all(np.equal(np.mod(vals, 1), 0)) and c in ("n", "n_discordant", "n_ovb_only", "n_naive_only", "dof"):
                df[c] = df[c].map(lambda v: "" if pd.isna(v) else f"{int(v)}")
            else:
                df[c] = df[c].map(lambda v: "" if pd.isna(v) else f"{v:.{floatfmt}f}")
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def main():
    parts = []
    summ = pd.read_csv(TAB / "summary_by_stratum.csv")
    tests = pd.read_csv(TAB / "hypothesis_tests.csv")
    sens = pd.read_csv(TAB / "sensitivity_summary.csv")
    ms = pd.read_csv(TAB / "measured_strength_placebo_by_stratum.csv")
    conc = pd.read_csv(TAB / "concordance.csv")

    # T1 overall by estimator x adj
    t1 = summ[summ.grouping == "estimator|adj"][["estimator", "adj", "n", "coverage", "cov_lo", "cov_hi", "bias_mean", "rmse", "width_mean", "se_ratio", "nu2_fallback_rate", "n_trimmed_mean"]]
    parts.append("### T1 — Overall performance by estimator and adjustment set (all 3,400 cohorts)\n\n" + md(t1))
    # T2 by strength x source
    t2 = summ[summ.grouping == "estimator|adj|strength|source"][["estimator", "adj", "strength", "source", "n", "coverage", "cov_lo", "cov_hi", "bias_mean", "rmse", "width_mean", "se_ratio"]]
    parts.append("### T2 — Coverage ladder by estimator × adjustment × strength × source\n\n" + md(t2))
    # T3 hypothesis tests
    t3 = tests[["family", "id", "test", "n", "effect", "ci_lo", "ci_hi", "statistic", "p_raw", "p_holm", "verdict"]].copy()
    t3["p_raw"] = t3["p_raw"].map(lambda v: f"{v:.2e}" if v < 1e-3 else f"{v:.3f}")
    t3["p_holm"] = t3["p_holm"].map(lambda v: f"{v:.2e}" if v < 1e-3 else f"{v:.3f}")
    parts.append("### T3 — Pre-registered hypothesis tests (Holm within family)\n\n" + md(t3))
    # T4 sensitivity summary (X_full and X_reduced, E1 and E2)
    t4 = sens[sens.estimator.isin(["E1_AIPW", "E2_DRATT", "E0_OLS"])][["estimator", "adj", "strength", "source", "n", "median_rv", "median_rva", "median_gamma_star", "frac_gamma_gt1", "median_e_pt", "median_e_ci", "frac_e_ci_gt1", "naive_coverage", "ovb03_coverage", "nu2_fallback_rate", "median_n_pairs", "median_proxy_gamma"]]
    parts.append("### T4 — Sensitivity diagnostics by stratum (medians; fractions above floor)\n\n" + md(t4))
    # T5 measured strength & placebo
    parts.append("### T5 — Measured strength of the withheld confounder (E1) and pre-period placebo effect, by stratum\n\n" + md(ms, 4))
    # T6 concordance
    parts.append("### T6 — Cross-family Spearman concordance (X_full)\n\n" + md(conc))
    # extra detail columns of tests
    extra_cols = [c for c in tests.columns if c not in t3.columns and c not in ("plan_ref", "direction_as_predicted")]
    t7 = tests[["id"] + extra_cols]
    parts.append("### T7 — Test detail columns\n\n" + md(t7))
    (TAB / "report_tables.md").write_text("\n\n".join(parts), encoding="utf-8")
    print("wrote report_tables.md")


if __name__ == "__main__":
    main()
