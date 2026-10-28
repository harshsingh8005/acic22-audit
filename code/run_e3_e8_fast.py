"""E3 (DoubleML-IRM ATTE) and E8 (EconML DRLearner) — the install-gated legs of the
frozen plan (§4), run after approval, in the isolated venv (.venv).

Per cohort x adjustment set (same seed-42 folds as E1/E2):
  E3  DoubleMLIRM(score='ATTE', n_folds=5, n_rep=1, trimming_threshold=0.02,
      ml_g=Ridge(alpha=10), ml_m=StandardScaler+LogisticRegression(C=0.5)), set_sample_splitting(folds)
      -> coef, se, CI; sensitivity_analysis(0.03, 0.03, rho=1) -> RV, RVa, OVB CI;
      on X_full: sensitivity_benchmark(['pre_trend']) -> (cf_y_bm, cf_d_bm, rho_bm, delta_theta);
      on X_reduced: sensitivity_analysis at the X_full benchmark strengths (cf_y_bm, cf_d_bm, rho=1)
      -> DoubleML's own version of S3 (OVB CI at measured strength).
      Also: my OVB layer evaluated on DoubleML's nuisance predictions (software cross-check).
  E8  EconML DRLearner(model_propensity=..., model_regression=..., model_final=StatsModelsLinearRegression,
      cv=folds, min_propensity=0.02) -> ATT = ate_inference(X[T==1]) (mean, se, CI);
      robustness_value(T=1) and sensitivity_summary(T=1, c_y=0.03, c_t=0.03) values (EconML defines these
      for the ATE of treatment 1 vs 0; reported as such).
Outputs: results/tables/e3e8_estimates.csv, results/tables/e3e8_cohort.csv, results/tables/e3e8_errors.jsonl
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse
import json
import sys
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acic_common import all_cohorts, load_cohort, make_folds, stratified_cohort_sample  # noqa: E402
from estimators import make_classifier, make_regressor, Z95  # noqa: E402
from sensitivity import evalues_continuous, ovb_bounds, ovb_elements_att, robustness_values  # noqa: E402

RESULTS = HERE.parent / "results" / "tables"
CF = 0.03


def _first(x):
    return float(np.asarray(x).ravel()[0])


def process_cohort(cid: int) -> dict:
    warnings.filterwarnings("ignore")
    import doubleml as dml
    from econml.dr import DRLearner
    from econml.sklearn_extensions.linear_model import StatsModelsLinearRegression
    t0 = time.perf_counter()
    out = {"cohort_id": cid, "rows": [], "cohort": {}, "error": None}
    try:
        c = load_cohort(cid)
        folds = make_folds(c.T)
        smpls = [[(tr, te) for tr, te in folds]]
        T, Y = c.T, c.Y
        sd_y = float(Y.std(ddof=1))
        base = {"cohort_id": cid, "true_satt": c.true_satt, "n": c.n, "n_treated": int(T.sum()), "sd_y": sd_y,
                **{k: c.factors[k] for k in ("strength", "source", "heterogeneity", "idiosyncrasy", "cell")}}
        coh = {**base}
        bm = None
        for adj in ("full", "reduced"):
            X = c.X_full if adj == "full" else c.X_reduced
            names = list(c.feature_names_full if adj == "full" else c.feature_names_reduced)
            # ------------------------------------------------ E3 DoubleML
            df = pd.DataFrame(X, columns=names); df["Y"] = Y; df["T"] = T
            data = dml.DoubleMLData(df, y_col="Y", d_cols="T", x_cols=names)
            m = dml.DoubleMLIRM(data, ml_g=make_regressor(42), ml_m=make_classifier(), n_folds=5, n_rep=1,
                                score="ATTE", trimming_threshold=0.02, draw_sample_splitting=False)
            m.set_sample_splitting(smpls)
            m.fit()
            ci = m.confint()
            att, se = _first(m.coef), _first(m.se)
            m.sensitivity_analysis(cf_y=CF, cf_d=CF, rho=1.0)
            sp = m.sensitivity_params
            row = {**base, "adj": adj, "estimator": "E3_DML_ATTE", "att": att, "se": se,
                   "ci_lo": float(ci.iloc[0, 0]), "ci_hi": float(ci.iloc[0, 1]),
                   "rv": _first(sp["rv"]), "rva": _first(sp["rva"]),
                   "ovb03_ci_lo": _first(sp["ci"]["lower"]), "ovb03_ci_hi": _first(sp["ci"]["upper"]),
                   "ovb03_theta_lo": _first(sp["theta"]["lower"]), "ovb03_theta_hi": _first(sp["theta"]["upper"])}
            e_pt, e_ci = evalues_continuous(att, row["ci_lo"], row["ci_hi"], sd_y)
            row.update({"e_pt": e_pt, "e_ci": e_ci})
            # my OVB layer on DoubleML's nuisance predictions (software cross-check)
            try:
                g0 = m.predictions["ml_g0"][:, 0, 0]; g1 = m.predictions["ml_g1"][:, 0, 0]; mh = m.predictions["ml_m"][:, 0, 0]
                mh = np.clip(mh, 0.02, 0.98)
                se_el = m.sensitivity_elements
                row.update({"dml_sigma2": _first(se_el["sigma2"]), "dml_nu2": _first(se_el["nu2"])})
                # influence values: psi / (-mean(psi_a)) reconstructed from DoubleML psi
                psi = m.psi[:, 0, 0]; psi_a = m.psi_deriv[:, 0, 0] if hasattr(m, "psi_deriv") else None
                infl = psi / (-np.mean(psi_a)) if psi_a is not None else psi
                el = ovb_elements_att(T, Y, mh, g0, g1, att, infl - infl.mean())
                rv_my, rva_my = robustness_values(el)
                b_my = ovb_bounds(el, CF, CF, 1.0)
                row.update({"my_sigma2_on_dml": el.sigma2, "my_nu2_on_dml": el.nu2, "my_nu2_fallback": el.nu2_fallback,
                            "my_rv_on_dml": rv_my, "my_rva_on_dml": rva_my, "my_ovb03_ci_lo": b_my["ci_lower"], "my_ovb03_ci_hi": b_my["ci_upper"]})
            except Exception as exc:  # noqa: BLE001
                row["crosscheck_error"] = str(exc)[:120]
            if adj == "full":
                try:
                    bmdf = m.sensitivity_benchmark(benchmarking_set=["pre_trend"])
                    bm = {"cf_y_bm": float(bmdf["cf_y"].iloc[0]), "cf_d_bm": float(bmdf["cf_d"].iloc[0]),
                          "rho_bm": float(bmdf["rho"].iloc[0]), "delta_theta_bm": float(bmdf["delta_theta"].iloc[0])}
                    coh.update({f"E3_{k}": v for k, v in bm.items()})
                except Exception as exc:  # noqa: BLE001
                    coh["E3_benchmark_error"] = str(exc)[:120]
            else:
                if bm is not None and bm["cf_y_bm"] > 0 and 0 < bm["cf_d_bm"] < 1:
                    m.sensitivity_analysis(cf_y=bm["cf_y_bm"], cf_d=bm["cf_d_bm"], rho=1.0)
                    sp2 = m.sensitivity_params
                    coh.update({"E3_red_ovbstar_ci_lo": _first(sp2["ci"]["lower"]), "E3_red_ovbstar_ci_hi": _first(sp2["ci"]["upper"]),
                                "E3_red_ovbstar_theta_lo": _first(sp2["theta"]["lower"]), "E3_red_ovbstar_theta_hi": _first(sp2["theta"]["upper"])})
                else:
                    coh.update({"E3_red_ovbstar_ci_lo": row["ci_lo"], "E3_red_ovbstar_ci_hi": row["ci_hi"],
                                "E3_red_ovbstar_theta_lo": att, "E3_red_ovbstar_theta_hi": att, "E3_bm_zero": True})
                coh.update({"E3_red_att": att, "E3_red_se": se})
            if adj == "full":
                coh.update({"E3_full_att": att, "E3_full_se": se})
            out["rows"].append(row)
            # ------------------------------------------------ E8 EconML DRLearner
            try:
                dr = DRLearner(model_propensity=make_classifier(), model_regression=make_regressor(42),
                               model_final=StatsModelsLinearRegression(), cv=[(tr, te) for tr, te in folds],
                               random_state=42, min_propensity=0.02)
                dr.fit(Y, T, X=X)
                inf = dr.ate_inference(X[T == 1])
                lo, hi = inf.conf_int_mean()
                r8 = {**base, "adj": adj, "estimator": "E8_EconML_DRLearner", "att": float(inf.mean_point), "se": float(inf.stderr_mean),
                      "ci_lo": float(lo), "ci_hi": float(hi)}
                try:
                    r8["rv_econml_ate"] = float(dr.robustness_value(T=1, alpha=0.05, interval_type="ci"))
                    r8["rv_econml_ate_theta"] = float(dr.robustness_value(T=1, alpha=0.05, interval_type="theta"))
                    si = dr.sensitivity_interval(T=1, alpha=0.05, c_y=CF, c_t=CF, rho=1.0, interval_type="ci")
                    r8["ovb03_ci_lo_econml_ate"], r8["ovb03_ci_hi_econml_ate"] = float(si[0]), float(si[1])
                    ate = dr.ate_inference(X)
                    r8["ate_econml"], r8["ate_se_econml"] = float(ate.mean_point), float(ate.stderr_mean)
                except Exception as exc:  # noqa: BLE001
                    r8["econml_sens_error"] = str(exc)[:120]
                e_pt, e_ci = evalues_continuous(r8["att"], r8["ci_lo"], r8["ci_hi"], sd_y)
                r8.update({"e_pt": e_pt, "e_ci": e_ci})
                out["rows"].append(r8)
            except Exception as exc:  # noqa: BLE001
                out["rows"].append({**base, "adj": adj, "estimator": "E8_EconML_DRLearner", "att": np.nan, "se": np.nan,
                                    "excluded_reason": str(exc)[:200]})
        coh["wall_s"] = time.perf_counter() - t0
        out["cohort"] = coh
    except Exception as exc:  # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", default="all")
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--chunk", type=int, default=200)
    ap.add_argument("--pilot", type=int, default=0)
    ap.add_argument("--tag", default="e3e8")
    args = ap.parse_args(argv)
    if args.pilot:
        ids = list(range(1, args.pilot + 1))
        t0 = time.perf_counter()
        outs = Parallel(n_jobs=min(args.n_jobs, len(ids)))(delayed(process_cohort)(i) for i in ids)
        errs = [o["error"] for o in outs if o["error"]]
        per = [o["cohort"].get("wall_s", np.nan) for o in outs]
        print(f"PILOT n={len(ids)} wall={time.perf_counter()-t0:.1f}s per-cohort={np.nanmean(per):.2f}s errors={len(errs)}")
        if errs:
            print(errs[0][-1500:])
        else:
            r = outs[0]["rows"]; print({k: v for k, v in r[0].items() if k in ("estimator", "att", "se", "rv", "rva")}, {k: v for k, v in r[1].items() if k in ("estimator", "att", "se", "rv_econml_ate")})
            print({k: v for k, v in outs[0]["cohort"].items() if k.startswith("E3_")})
        print(f"PROJECTED full at n_jobs={args.n_jobs}: ~{np.nanmean(per)*3400/args.n_jobs/60:.1f} min")
        return 0
    ids = all_cohorts() if args.cohorts == "all" else (stratified_cohort_sample(int(args.cohorts.replace("sample", ""))) if args.cohorts.startswith("sample") else [int(x) for x in args.cohorts.split(",")])
    RESULTS.mkdir(parents=True, exist_ok=True)
    est_path = RESULTS / f"{args.tag}_estimates.csv"; coh_path = RESULTS / f"{args.tag}_cohort.csv"; err_path = RESULTS / f"{args.tag}_errors.jsonl"
    done = set(pd.read_csv(coh_path, usecols=["cohort_id"])["cohort_id"].tolist()) if coh_path.exists() else set()
    todo = [i for i in ids if i not in done]
    print(f"{len(ids)} cohorts, {len(done)} done, {len(todo)} to run; n_jobs={args.n_jobs}", flush=True)
    t0 = time.perf_counter()
    for s in range(0, len(todo), args.chunk):
        chunk = todo[s:s + args.chunk]
        outs = Parallel(n_jobs=args.n_jobs)(delayed(process_cohort)(i) for i in chunk)
        rows = [r for o in outs if not o["error"] for r in o["rows"]]
        cohs = [o["cohort"] for o in outs if not o["error"]]
        errs = [o for o in outs if o["error"]]
        if rows:
            pd.DataFrame(rows).to_csv(est_path, mode="a", header=not est_path.exists(), index=False)
        if cohs:
            pd.DataFrame(cohs).to_csv(coh_path, mode="a", header=not coh_path.exists(), index=False)
        if errs:
            with err_path.open("a") as fh:
                for o in errs:
                    fh.write(json.dumps({"cohort_id": o["cohort_id"], "error": o["error"]}) + "\n")
        print(f"  {s + len(chunk)}/{len(todo)} done in {(time.perf_counter()-t0)/60:.1f} min ({len(errs)} errors)", flush=True)
    print(f"FINISHED {len(todo)} cohorts in {(time.perf_counter()-t0)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
