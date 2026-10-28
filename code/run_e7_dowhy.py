"""E7 — DoWhy (PyWhy-native) cross-check, single process (frozen plan §4 E7 / §10 V5):
backdoor.linear_regression estimate + add_unobserved_common_cause(linear-partial-R2)
robustness values on X_full, compared with E0 OLS + our sensemakr-form RV; plus DoWhy
IPW-ATT (exploratory). Runs on the neural sub-sample cohort list.
Output: results/tables/e7_dowhy_crosscheck.csv"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLBACKEND", "Agg")
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acic_common import load_cohort  # noqa: E402
from estimators import ols_att, aipw_att, ols_att  # noqa: E402
from sensitivity import rv_linear  # noqa: E402

logging.getLogger("dowhy").setLevel(logging.ERROR)
from dowhy import CausalModel  # noqa: E402

TAB = HERE.parent / "results" / "tables"


def main():
    ids = json.loads((TAB / "neural_arm_cohort_ids.json").read_text())["ids"]
    nmax = [a for a in sys.argv[1:] if a.isdigit()]
    if nmax:
        ids = ids[: int(nmax[0])]
    rows = []
    t0 = time.perf_counter()
    for cid in ids:
        c = load_cohort(cid)
        X = c.X_full; names = [f"w{i}" for i in range(X.shape[1])]
        df = pd.DataFrame(X, columns=names); df["T"] = c.T; df["Y"] = c.Y
        o = ols_att(X, c.T, c.Y)
        model = CausalModel(data=df, treatment="T", outcome="Y", common_causes=names)
        ident = model.identify_effect(proceed_when_unidentifiable=True)
        est = model.estimate_effect(ident, method_name="backdoor.linear_regression", test_significance=False)
        ref = model.refute_estimate(ident, est, method_name="add_unobserved_common_cause",
                                    simulation_method="linear-partial-R2", benchmark_common_causes=[names[-1]],
                                    effect_fraction_on_treatment=[1.0], plot=False)
        an = getattr(ref, "refutation_result", None)
        dw_rv = getattr(an, "RVq", None) if an is not None else None
        dw_rva = getattr(an, "RV_qalpha", None) if an is not None else None
        if dw_rv is None:
            for attr in ("stats", "sensitivity_stats"):
                st = getattr(ref, attr, None)
                if isinstance(st, dict) and "robustness_value" in st:
                    dw_rv = st["robustness_value"]; dw_rva = st.get("robustness_value_alpha"); break
        row = {"cohort_id": cid, "true_satt": c.true_satt, "e0_att": o.att, "dowhy_ols_att": float(est.value),
               "rv_mine": rv_linear(o.t_stat, o.dof, 1.0, 1.0), "rva_mine": rv_linear(o.t_stat, o.dof, 1.0, 0.05),
               "rv_dowhy": float(dw_rv) if dw_rv is not None else np.nan,
               "rva_dowhy": float(dw_rva) if dw_rva is not None else np.nan}
        if "--ipw" in sys.argv:
            try:
                est_ipw = model.estimate_effect(ident, method_name="backdoor.propensity_score_weighting", target_units="att",
                                                method_params={"weighting_scheme": "ips_weight"}, test_significance=False)
                row["dowhy_ipw_att"] = float(est_ipw.value)
            except Exception as exc:  # noqa: BLE001
                row["dowhy_ipw_att"] = np.nan; row["dowhy_ipw_error"] = str(exc)[:120]
        rows.append(row)
        if len(rows) % 10 == 0:
            print(f"  {len(rows)}/{len(ids)} in {time.perf_counter()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "e7_dowhy_crosscheck.csv", index=False)
    print(f"n={len(df)}  max|dAtt|={np.nanmax(np.abs(df.e0_att-df.dowhy_ols_att)):.2e}  max|dRV|={np.nanmax(np.abs(df.rv_mine-df.rv_dowhy)):.2e}  "
          f"max|dRVa|={np.nanmax(np.abs(df.rva_mine-df.rva_dowhy)):.2e}  ipw_ok={df.dowhy_ipw_att.notna().sum() if 'dowhy_ipw_att' in df else 'skipped'}  ({time.perf_counter()-t0:.0f}s)")


if __name__ == "__main__":
    main()
