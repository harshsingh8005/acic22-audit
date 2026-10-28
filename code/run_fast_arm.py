"""Fast arm (frozen plan §4/§5/§8): E0 OLS-ATT, E1 AIPW-ATT, E2 DR-ATT on every
cohort x {X_full, X_reduced}, with the full sensitivity layer, the measured
omitted-confounder strength, and the pre-period negative-control outcome.

Usage
  python run_fast_arm.py --pilot 5                 # timing only, nothing written
  python run_fast_arm.py --cohorts all --n-jobs 8  # the pre-registered run
  python run_fast_arm.py --cohorts sample100       # fallback: 100/cell (plan §8)
Outputs (results/tables/):
  fast_arm_estimates.csv   one row per (cohort, adjustment set, estimator)
  fast_arm_cohort.csv      one row per cohort (factors, measured strengths, placebo, timing)
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
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acic_common import (all_cohorts, load_cohort, make_folds,  # noqa: E402
                         stratified_cohort_sample)
from estimators import aipw_att, dr_att, ols_att, Z95  # noqa: E402
from sensitivity import (evalues_continuous, measured_confounder_strength,  # noqa: E402
                         nn_match_pairs, ovb_bounds, ovb_elements_att,
                         robustness_values, rosenbaum_gamma_star, rv_linear,
                         widening_bound_proxy)

RESULTS = HERE.parent / "results" / "tables"
CF_DEFAULT = 0.03


def _sens_common(att, ci_lo, ci_hi, sd_y):
    e_pt, e_ci = evalues_continuous(att, ci_lo, ci_hi, sd_y)
    return {"e_pt": e_pt, "e_ci": e_ci}


def _rosenbaum(e_hat, T, Y, direction, keep=None):
    if keep is not None:
        idx = np.flatnonzero(keep)
        pairs = nn_match_pairs(e_hat[idx], T[idx])
        pairs = [(int(idx[i]), int(idx[j])) for i, j in pairs]
    else:
        pairs = nn_match_pairs(e_hat, T)
    if len(pairs) < 20:
        return {"gamma_star": np.nan, "p_gamma1": np.nan, "n_pairs": len(pairs), "pair_mean_diff": np.nan}
    D = np.array([Y[i] - Y[j] for i, j in pairs])
    g, p1, nnz = rosenbaum_gamma_star(D, direction)
    return {"gamma_star": g, "p_gamma1": p1, "n_pairs": len(pairs), "pair_mean_diff": float(D.mean())}


def process_cohort(cid: int) -> dict:
    t0 = time.perf_counter()
    out = {"cohort_id": cid, "rows": [], "cohort": {}, "error": None}
    try:
        c = load_cohort(cid)
        folds = make_folds(c.T)
        T, Y = c.T, c.Y
        sd_y = float(Y.std(ddof=1))
        base = {"cohort_id": cid, "true_satt": c.true_satt, "n": c.n, "n_treated": int(T.sum()),
                "sd_y": sd_y, **{k: c.factors[k] for k in ("strength", "source", "heterogeneity", "idiosyncrasy", "cell")}}
        el_store = {}
        for adj in ("full", "reduced"):
            X = c.X_full if adj == "full" else c.X_reduced
            # ---- E0 OLS ------------------------------------------------------
            o = ols_att(X, T, Y)
            row = {**base, "adj": adj, "estimator": "E0_OLS", "att": o.att, "se": o.se,
                   "ci_lo": o.ci_lo, "ci_hi": o.ci_hi, "n_trimmed": 0,
                   "rv": rv_linear(o.t_stat, o.dof, 1.0, 1.0),
                   "rva": rv_linear(o.t_stat, o.dof, 1.0, 0.05),
                   "t_stat": o.t_stat, "dof": o.dof, **_sens_common(o.att, o.ci_lo, o.ci_hi, sd_y),
                   "proxy_gamma": widening_bound_proxy(o.att, o.se)}
            out["rows"].append(row)
            # ---- E1 AIPW -----------------------------------------------------
            r1, nu1 = aipw_att(X, T, Y, folds)
            el1 = ovb_elements_att(T, Y, nu1.e_hat, nu1.mu0_hat, nu1.mu1_hat, r1.att, r1.infl)
            rv1, rva1 = robustness_values(el1)
            b1 = ovb_bounds(el1, CF_DEFAULT, CF_DEFAULT, 1.0)
            ros1 = _rosenbaum(nu1.e_hat, T, Y, np.sign(r1.att) or 1)
            row = {**base, "adj": adj, "estimator": "E1_AIPW", "att": r1.att, "se": r1.se,
                   "ci_lo": r1.ci_lo, "ci_hi": r1.ci_hi, "n_trimmed": 0,
                   "prop_min": r1.prop_min, "prop_max": r1.prop_max,
                   "e_raw_min": float(nu1.e_raw.min()), "e_raw_max": float(nu1.e_raw.max()),
                   "sigma2": el1.sigma2, "nu2": el1.nu2, "nu2_dr": el1.nu2_dr, "nu2_fallback": el1.nu2_fallback, "rv": rv1, "rva": rva1,
                   "ovb03_ci_lo": b1["ci_lower"], "ovb03_ci_hi": b1["ci_upper"], "ovb03_bias_bound": b1["bias_bound"],
                   **_sens_common(r1.att, r1.ci_lo, r1.ci_hi, sd_y),
                   **ros1, "proxy_gamma": widening_bound_proxy(r1.att, r1.se)}
            out["rows"].append(row)
            el_store[("E1", adj)] = el1
            # ---- E2 DR-ATT (Crump trim) --------------------------------------
            try:
                r2, nu2 = dr_att(X, T, Y, folds)
                el2 = ovb_elements_att(T, Y, nu2.e_hat, nu2.mu0_hat, nu2.mu1_hat, r2.att, r2.infl, keep=r2.keep)
                rv2, rva2 = robustness_values(el2)
                b2 = ovb_bounds(el2, CF_DEFAULT, CF_DEFAULT, 1.0)
                ros2 = _rosenbaum(nu2.e_hat, T, Y, np.sign(r2.att) or 1, keep=r2.keep)
                row = {**base, "adj": adj, "estimator": "E2_DRATT", "att": r2.att, "se": r2.se,
                       "ci_lo": r2.ci_lo, "ci_hi": r2.ci_hi, "n_trimmed": r2.n_trimmed,
                       "prop_min": r2.prop_min, "prop_max": r2.prop_max,
                       "sigma2": el2.sigma2, "nu2": el2.nu2, "nu2_dr": el2.nu2_dr, "nu2_fallback": el2.nu2_fallback, "rv": rv2, "rva": rva2,
                       "ovb03_ci_lo": b2["ci_lower"], "ovb03_ci_hi": b2["ci_upper"], "ovb03_bias_bound": b2["bias_bound"],
                       **_sens_common(r2.att, r2.ci_lo, r2.ci_hi, sd_y),
                       **ros2, "proxy_gamma": widening_bound_proxy(r2.att, r2.se)}
                out["rows"].append(row)
                el_store[("E2", adj)] = el2
            except ValueError as exc:  # pre-specified exclusion (plan §11)
                out["rows"].append({**base, "adj": adj, "estimator": "E2_DRATT", "att": np.nan,
                                    "se": np.nan, "excluded_reason": str(exc)})
        # ---- cross-adjustment quantities (plan §5.4, S3) ------------------------
        coh = {**base}
        for est in ("E1", "E2"):
            if (est, "full") in el_store and (est, "reduced") in el_store:
                el_l, el_s = el_store[(est, "full")], el_store[(est, "reduced")]
                cf_y, cf_d, rho_hat, b_star = measured_confounder_strength(el_l, el_s)
                bs = ovb_bounds(el_s, cf_y, cf_d, 1.0)
                coh.update({f"{est}_cf_y_star": cf_y, f"{est}_cf_d_star": cf_d,
                            f"{est}_rho_hat_star": rho_hat, f"{est}_B_star": b_star,
                            f"{est}_rho_true": ((el_s.theta - c.true_satt) / b_star) if b_star > 0 else np.nan,
                            f"{est}_red_att": el_s.theta, f"{est}_full_att": el_l.theta,
                            f"{est}_red_ovbstar_ci_lo": bs["ci_lower"], f"{est}_red_ovbstar_ci_hi": bs["ci_upper"],
                            f"{est}_red_ovbstar_theta_lo": bs["theta_lower"], f"{est}_red_ovbstar_theta_hi": bs["theta_upper"]})
        # ---- negative-control (pre-period placebo) outcome under X_reduced -----
        rp, _ = aipw_att(c.X_reduced, T, c.pre_trend, folds)
        coh.update({"placebo_att": rp.att, "placebo_se": rp.se})
        coh["wall_s"] = time.perf_counter() - t0
        out["cohort"] = coh
    except Exception as exc:  # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", default="all", help="all | sample100 | comma-separated ids")
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--pilot", type=int, default=0, help="timing pilot on N cohorts (nothing written)")
    ap.add_argument("--chunk", type=int, default=200)
    ap.add_argument("--tag", default="fast_arm")
    args = ap.parse_args(argv)

    if args.pilot:
        ids = list(range(1, args.pilot + 1))
        t0 = time.perf_counter()
        outs = Parallel(n_jobs=min(args.n_jobs, len(ids)))(delayed(process_cohort)(i) for i in ids)
        wall = time.perf_counter() - t0
        errs = [o["error"] for o in outs if o["error"]]
        per = [o["cohort"].get("wall_s", np.nan) for o in outs]
        print(f"PILOT n={len(ids)} wall={wall:.1f}s per-cohort(mean worker time)={np.nanmean(per):.2f}s "
              f"errors={len(errs)}")
        if errs:
            print(errs[0])
        proj = np.nanmean(per) * 3400 / args.n_jobs / 60
        print(f"PROJECTED full-arm wall at n_jobs={args.n_jobs}: ~{proj:.1f} min (+load overhead)")
        return 0

    if args.cohorts == "all":
        ids = all_cohorts()
    elif args.cohorts.startswith("sample"):
        ids = stratified_cohort_sample(int(args.cohorts.replace("sample", "")))
    else:
        ids = [int(x) for x in args.cohorts.split(",")]

    RESULTS.mkdir(parents=True, exist_ok=True)
    est_path = RESULTS / f"{args.tag}_estimates.csv"
    coh_path = RESULTS / f"{args.tag}_cohort.csv"
    err_path = RESULTS / f"{args.tag}_errors.jsonl"
    done = set()
    if coh_path.exists():
        done = set(pd.read_csv(coh_path, usecols=["cohort_id"])["cohort_id"].tolist())
    todo = [i for i in ids if i not in done]
    print(f"{len(ids)} cohorts requested, {len(done)} already done, {len(todo)} to run; n_jobs={args.n_jobs}")
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
        el = time.perf_counter() - t0
        print(f"  {s + len(chunk)}/{len(todo)} done in {el/60:.1f} min ({len(errs)} errors this chunk)", flush=True)
    print(f"FINISHED {len(todo)} cohorts in {(time.perf_counter()-t0)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
