"""Neural arm (frozen plan §4 E4/E5/E6, §5, §8): TARNet + DragonNet (PyTorch,
K seeds) and the ridge T-learner on the stratified cohort sub-sample, X_full and
X_reduced; PEHE-on-treated; the pre-registered neural-DR interval with its
sensitivity layer; the exploratory seed-ensemble interval; DoWhy PyWhy-native
cross-check (E7) on X_full.

Usage
  python run_neural_arm.py --pilot 2                 # timing only
  python run_neural_arm.py --per-cell 25 --seeds 5   # the pre-registered run
Outputs (results/tables/): neural_arm_estimates.csv, neural_arm_cohort.csv, neural_arm_errors.jsonl
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

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
from acic_common import load_cohort, make_folds, stratified_cohort_sample  # noqa: E402
from estimators import (aipw_att, aipw_att_from_nuisances, ridge_tlearner_cate,  # noqa: E402
                        Z95, CLIP_AIPW)
from sensitivity import (evalues_continuous, nn_match_pairs, ovb_bounds,  # noqa: E402
                         ovb_elements_att, robustness_values, rosenbaum_gamma_star,
                         widening_bound_proxy)

RESULTS = HERE.parent / "results" / "tables"
CF_DEFAULT = 0.03


def pehe(tau_hat, tau_true):
    m = np.isfinite(tau_true)
    return float(np.sqrt(np.mean((tau_hat[m] - tau_true[m]) ** 2)))


def _rosenbaum(e_hat, T, Y, direction):
    pairs = nn_match_pairs(e_hat, T)
    if len(pairs) < 20:
        return {"gamma_star": np.nan, "p_gamma1": np.nan, "n_pairs": len(pairs)}
    D = np.array([Y[i] - Y[j] for i, j in pairs])
    g, p1, _ = rosenbaum_gamma_star(D, direction)
    return {"gamma_star": g, "p_gamma1": p1, "n_pairs": len(pairs)}


def process_cohort(cid: int, seeds: list[int], do_dowhy: bool = True) -> dict:
    import torch  # noqa: F401  (import inside worker)
    from neural import fit_predict
    t0 = time.perf_counter()
    out = {"cohort_id": cid, "rows": [], "cohort": {}, "error": None}
    try:
        c = load_cohort(cid)
        folds = make_folds(c.T)
        T, Y = c.T, c.Y
        sd_y = float(Y.std(ddof=1))
        treated = T == 1
        base = {"cohort_id": cid, "true_satt": c.true_satt, "n": c.n, "n_treated": int(T.sum()), "sd_y": sd_y,
                **{k: c.factors[k] for k in ("strength", "source", "heterogeneity", "idiosyncrasy", "cell")}}
        coh = {**base}
        for adj in ("full", "reduced"):
            X = c.X_full if adj == "full" else c.X_reduced
            # E1 nuisances (for TARNet's propensity and as the paired comparator on identical cohorts)
            r1, nu1 = aipw_att(X, T, Y, folds)
            # E6 ridge T-learner
            m0r, m1r = ridge_tlearner_cate(X, T, Y)
            out["rows"].append({**base, "adj": adj, "estimator": "E6_ridgeT", "att": float((m1r - m0r)[treated].mean()),
                                "pehe_treated": pehe((m1r - m0r)[treated], c.true_catt_treated)})
            for name, dragon in (("E4_TARNet", False), ("E5_DragonNet", True)):
                mu0s, mu1s, es, taus, pehes, infos = [], [], [], [], [], []
                for s in seeds:
                    mu0, mu1, e_net, info = fit_predict(X, T, Y, seed=s, dragon=dragon)
                    mu0s.append(mu0); mu1s.append(mu1); es.append(e_net); infos.append(info)
                    tau_k = mu1 - mu0
                    taus.append(float(tau_k[treated].mean()))
                    pehes.append(pehe(tau_k[treated], c.true_catt_treated))
                mu0m, mu1m = np.mean(mu0s, 0), np.mean(mu1s, 0)
                tau_ens = mu1m - mu0m
                att_plug = float(tau_ens[treated].mean())
                # pre-registered neural-DR interval
                if dragon:
                    e_hat = np.clip(np.mean(es, 0), CLIP_AIPW[0], CLIP_AIPW[1])
                    e_src = "network"
                else:
                    e_hat = nu1.e_hat; e_src = "E1_crossfit"
                rd = aipw_att_from_nuisances(T, Y, e_hat, mu0m)
                el = ovb_elements_att(T, Y, e_hat, mu0m, mu1m, rd.att, rd.infl)
                rv, rva = robustness_values(el)
                b = ovb_bounds(el, CF_DEFAULT, CF_DEFAULT, 1.0)
                e_pt, e_ci = evalues_continuous(rd.att, rd.ci_lo, rd.ci_hi, sd_y)
                ros = _rosenbaum(e_hat, T, Y, np.sign(rd.att) or 1)
                sd_seed = float(np.std(taus, ddof=1)) if len(taus) > 1 else np.nan
                out["rows"].append({**base, "adj": adj, "estimator": name,
                                    "att": att_plug, "att_seed_sd": sd_seed,
                                    "ens_ci_lo": att_plug - Z95 * sd_seed, "ens_ci_hi": att_plug + Z95 * sd_seed,
                                    "pehe_treated": pehe(tau_ens[treated], c.true_catt_treated),
                                    "pehe_treated_seed_mean": float(np.mean(pehes)),
                                    "dr_att": rd.att, "dr_se": rd.se, "dr_ci_lo": rd.ci_lo, "dr_ci_hi": rd.ci_hi,
                                    "e_source": e_src, "prop_min": float(e_hat.min()), "prop_max": float(e_hat.max()),
                                    "sigma2": el.sigma2, "nu2": el.nu2, "nu2_dr": el.nu2_dr, "nu2_fallback": el.nu2_fallback,
                                    "rv": rv, "rva": rva, "ovb03_ci_lo": b["ci_lower"], "ovb03_ci_hi": b["ci_upper"],
                                    "e_pt": e_pt, "e_ci": e_ci, **ros, "proxy_gamma": widening_bound_proxy(rd.att, rd.se),
                                    "epochs_mean": float(np.mean([i["epochs"] for i in infos])),
                                    "best_epoch_mean": float(np.mean([i["best_epoch"] for i in infos]))})
            # E1 row on the same cohort/adjustment for the paired comparison
            out["rows"].append({**base, "adj": adj, "estimator": "E1_AIPW", "att": r1.att, "dr_att": r1.att,
                                "dr_se": r1.se, "dr_ci_lo": r1.ci_lo, "dr_ci_hi": r1.ci_hi})
        # E7 DoWhy PyWhy-native cross-check (X_full): OLS estimate + linear-partial-R2 RV; IPW-ATT exploratory
        if do_dowhy:
            import logging
            logging.getLogger("dowhy").setLevel(logging.ERROR)
            from dowhy import CausalModel
            X = c.X_full; names = [f"w{i}" for i in range(X.shape[1])]
            df = pd.DataFrame(X, columns=names); df["T"] = T; df["Y"] = Y
            model = CausalModel(data=df, treatment="T", outcome="Y", common_causes=names)
            ident = model.identify_effect(proceed_when_unidentifiable=True)
            est = model.estimate_effect(ident, method_name="backdoor.linear_regression", test_significance=False)
            ref = model.refute_estimate(ident, est, method_name="add_unobserved_common_cause",
                                        simulation_method="linear-partial-R2", benchmark_common_causes=[names[-1]],
                                        effect_fraction_on_treatment=[1.0], plot=False)
            an = getattr(ref, "refutation_result", None)
            coh.update({"dowhy_ols_att": float(est.value),
                        "dowhy_rv": float(getattr(an, "RVq", np.nan)) if an is not None else np.nan,
                        "dowhy_rva": float(getattr(an, "RV_qalpha", np.nan)) if an is not None else np.nan})
            try:
                est_ipw = model.estimate_effect(ident, method_name="backdoor.propensity_score_weighting",
                                                target_units="att", method_params={"weighting_scheme": "ips_weight"},
                                                test_significance=False)
                coh["dowhy_ipw_att"] = float(est_ipw.value)
            except Exception as exc:  # noqa: BLE001
                coh["dowhy_ipw_att"] = np.nan; coh["dowhy_ipw_error"] = str(exc)[:200]
        coh["wall_s"] = time.perf_counter() - t0
        out["cohort"] = coh
    except Exception as exc:  # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-cell", type=int, default=25)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--pilot", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=48)
    ap.add_argument("--tag", default="neural_arm")
    ap.add_argument("--no-dowhy", action="store_true")
    args = ap.parse_args(argv)
    seeds = list(range(args.seeds))

    if args.pilot:
        ids = stratified_cohort_sample(1)[: args.pilot]
        t0 = time.perf_counter()
        outs = Parallel(n_jobs=min(args.n_jobs, len(ids)))(delayed(process_cohort)(i, seeds, not args.no_dowhy) for i in ids)
        wall = time.perf_counter() - t0
        per = [o["cohort"].get("wall_s", np.nan) for o in outs]
        errs = [o["error"] for o in outs if o["error"]]
        print(f"PILOT n={len(ids)} seeds={len(seeds)} wall={wall:.1f}s per-cohort worker time={np.nanmean(per):.1f}s errors={len(errs)}")
        if errs:
            print(errs[0])
        for pc in (25, 10, 6):
            print(f"  projected {pc}/cell ({17*pc} cohorts) at n_jobs={args.n_jobs}: ~{np.nanmean(per)*17*pc/args.n_jobs/60:.0f} min")
        return 0

    ids = stratified_cohort_sample(args.per_cell)
    RESULTS.mkdir(parents=True, exist_ok=True)
    est_path = RESULTS / f"{args.tag}_estimates.csv"
    coh_path = RESULTS / f"{args.tag}_cohort.csv"
    err_path = RESULTS / f"{args.tag}_errors.jsonl"
    done = set(pd.read_csv(coh_path, usecols=["cohort_id"])["cohort_id"].tolist()) if coh_path.exists() else set()
    todo = [i for i in ids if i not in done]
    print(f"{len(ids)} cohorts ({args.per_cell}/cell), {len(done)} done, {len(todo)} to run; seeds={seeds}; n_jobs={args.n_jobs}")
    (RESULTS / f"{args.tag}_cohort_ids.json").write_text(json.dumps({"per_cell": args.per_cell, "seed": 20260818, "ids": ids, "seeds": seeds}))
    t0 = time.perf_counter()
    for s in range(0, len(todo), args.chunk):
        chunk = todo[s:s + args.chunk]
        outs = Parallel(n_jobs=args.n_jobs)(delayed(process_cohort)(i, seeds, not args.no_dowhy) for i in chunk)
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
