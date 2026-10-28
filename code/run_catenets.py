"""catenets replication of E4/E5 (frozen plan §4 E8, install-gated leg, run after approval).
catenets.models.torch TARNet / DragonNet with library DEFAULTS (lr 1e-4, weight_decay 1e-4,
n_iter 10000 iterations, batch 100, val_split 0.3, patience 10 x n_iter_print 50, n_iter_min 200,
batch_norm True, elu), seeds 0..K-1, on the same 170-cohort sub-sample and identical cohorts as
run_neural_arm.py.  Same outputs as the PyTorch nets: plug-in SATT, PEHE-on-treated, the
neural-DR interval (network mu0 + own propensity head for DragonNet / E1's cross-fitted e for TARNet),
seed-ensemble interval (exploratory).  Default: X_full (the pre-registered S2/S5 comparisons);
--adj both adds X_reduced.
Output: results/tables/catenets_arm_estimates.csv, catenets_arm_cohort.csv, catenets_arm_errors.jsonl
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
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acic_common import load_cohort, make_folds  # noqa: E402
from estimators import aipw_att, aipw_att_from_nuisances, Z95, CLIP_AIPW  # noqa: E402
from sensitivity import (evalues_continuous, nn_match_pairs, ovb_bounds, ovb_elements_att,  # noqa: E402
                         robustness_values, rosenbaum_gamma_star)

RESULTS = HERE.parent / "results" / "tables"
CF = 0.03


def pehe(tau_hat, tau_true):
    m = np.isfinite(tau_true)
    return float(np.sqrt(np.mean((tau_hat[m] - tau_true[m]) ** 2)))


def _arr(t):
    import torch
    return t.detach().cpu().numpy().ravel() if isinstance(t, torch.Tensor) else np.asarray(t, dtype=float).ravel()


def _dragon_propensity(net, X):
    """DragonNet's own propensity from its representation + propensity head; None if not extractable."""
    import torch
    try:
        net.eval()
        with torch.no_grad():
            Xt = net._check_tensor(X).float()
            rep = net._repr_estimator(Xt)
            out = net._propensity_estimator(rep)
            out = out.detach().cpu().numpy()
        if out.ndim == 2 and out.shape[1] == 2:
            z = out - out.max(axis=1, keepdims=True); p = np.exp(z); p = p / p.sum(axis=1, keepdims=True)
            e = p[:, 1]
        else:
            v = out.ravel()
            e = v if (v.min() >= 0 and v.max() <= 1) else 1.0 / (1.0 + np.exp(-v))
        return np.clip(e, CLIP_AIPW[0], CLIP_AIPW[1])
    except Exception:
        return None


def _rosenbaum(e_hat, T, Y, direction):
    pairs = nn_match_pairs(e_hat, T)
    if len(pairs) < 20:
        return {"gamma_star": np.nan, "p_gamma1": np.nan, "n_pairs": len(pairs)}
    D = np.array([Y[i] - Y[j] for i, j in pairs])
    g, p1, _ = rosenbaum_gamma_star(D, direction)
    return {"gamma_star": g, "p_gamma1": p1, "n_pairs": len(pairs)}


def process_cohort(cid: int, seeds: list[int], adjs: tuple) -> dict:
    warnings.filterwarnings("ignore")
    import torch
    torch.set_num_threads(1)
    from catenets.models.torch import TARNet, DragonNet
    t0 = time.perf_counter()
    out = {"cohort_id": cid, "rows": [], "cohort": {}, "error": None}
    try:
        c = load_cohort(cid)
        folds = make_folds(c.T)
        T, Y = c.T, c.Y
        sd_y = float(Y.std(ddof=1)); treated = T == 1
        base = {"cohort_id": cid, "true_satt": c.true_satt, "n": c.n, "n_treated": int(T.sum()), "sd_y": sd_y,
                **{k: c.factors[k] for k in ("strength", "source", "heterogeneity", "idiosyncrasy", "cell")}}
        for adj in adjs:
            X = c.X_full if adj == "full" else c.X_reduced
            r1, nu1 = aipw_att(X, T, Y, folds)
            out["rows"].append({**base, "adj": adj, "estimator": "E1_AIPW", "att": r1.att, "dr_att": r1.att,
                                "dr_se": r1.se, "dr_ci_lo": r1.ci_lo, "dr_ci_hi": r1.ci_hi})
            for name, cls in (("C4_catenets_TARNet", TARNet), ("C5_catenets_DragonNet", DragonNet)):
                mu0s, mu1s, es, taus, pehes = [], [], [], [], []
                for s in seeds:
                    net = cls(n_unit_in=X.shape[1], seed=s)
                    net.fit(X, Y, T)
                    cate, y0, y1 = net.predict(X, return_po=True)
                    mu0, mu1 = _arr(y0), _arr(y1)
                    mu0s.append(mu0); mu1s.append(mu1)
                    tau_k = mu1 - mu0
                    taus.append(float(tau_k[treated].mean())); pehes.append(pehe(tau_k[treated], c.true_catt_treated))
                    if cls is DragonNet:
                        e_k = _dragon_propensity(net, X)
                        if e_k is not None:
                            es.append(e_k)
                mu0m, mu1m = np.mean(mu0s, 0), np.mean(mu1s, 0)
                tau_ens = mu1m - mu0m
                att_plug = float(tau_ens[treated].mean())
                if cls is DragonNet and es:
                    e_hat = np.clip(np.mean(es, 0), CLIP_AIPW[0], CLIP_AIPW[1]); e_src = "network"
                else:
                    e_hat = nu1.e_hat; e_src = "E1_crossfit"
                rd = aipw_att_from_nuisances(T, Y, e_hat, mu0m)
                el = ovb_elements_att(T, Y, e_hat, mu0m, mu1m, rd.att, rd.infl)
                rv, rva = robustness_values(el)
                b = ovb_bounds(el, CF, CF, 1.0)
                e_pt, e_ci = evalues_continuous(rd.att, rd.ci_lo, rd.ci_hi, sd_y)
                ros = _rosenbaum(e_hat, T, Y, np.sign(rd.att) or 1)
                sd_seed = float(np.std(taus, ddof=1)) if len(taus) > 1 else np.nan
                out["rows"].append({**base, "adj": adj, "estimator": name, "att": att_plug, "att_seed_sd": sd_seed,
                                    "ens_ci_lo": att_plug - Z95 * sd_seed, "ens_ci_hi": att_plug + Z95 * sd_seed,
                                    "pehe_treated": pehe(tau_ens[treated], c.true_catt_treated), "pehe_treated_seed_mean": float(np.mean(pehes)),
                                    "dr_att": rd.att, "dr_se": rd.se, "dr_ci_lo": rd.ci_lo, "dr_ci_hi": rd.ci_hi, "e_source": e_src,
                                    "prop_min": float(e_hat.min()), "prop_max": float(e_hat.max()),
                                    "sigma2": el.sigma2, "nu2": el.nu2, "nu2_fallback": el.nu2_fallback, "rv": rv, "rva": rva,
                                    "ovb03_ci_lo": b["ci_lower"], "ovb03_ci_hi": b["ci_upper"], "e_pt": e_pt, "e_ci": e_ci, **ros})
        out["cohort"] = {**base, "wall_s": time.perf_counter() - t0}
    except Exception as exc:  # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--chunk", type=int, default=16)
    ap.add_argument("--adj", default="full", choices=["full", "both"])
    ap.add_argument("--pilot", type=int, default=0)
    ap.add_argument("--tag", default="catenets_arm")
    args = ap.parse_args(argv)
    seeds = list(range(args.seeds)); adjs = ("full",) if args.adj == "full" else ("full", "reduced")
    ids = json.loads((RESULTS / "neural_arm_cohort_ids.json").read_text())["ids"]   # identical 170-cohort sample
    if args.pilot:
        ids_p = ids[: args.pilot]; t0 = time.perf_counter()
        outs = Parallel(n_jobs=min(args.n_jobs, len(ids_p)))(delayed(process_cohort)(i, seeds, adjs) for i in ids_p)
        errs = [o["error"] for o in outs if o["error"]]; per = [o["cohort"].get("wall_s", np.nan) for o in outs]
        print(f"PILOT n={len(ids_p)} wall={time.perf_counter()-t0:.1f}s per-cohort={np.nanmean(per):.1f}s errors={len(errs)}")
        if errs: print(errs[0][-1500:])
        else: print([{k: r.get(k) for k in ("estimator", "att", "dr_att", "dr_se", "pehe_treated", "e_source")} for r in outs[0]["rows"]])
        print(f"PROJECTED 170 cohorts at n_jobs={args.n_jobs}: ~{np.nanmean(per)*170/args.n_jobs/60:.0f} min")
        return 0
    RESULTS.mkdir(parents=True, exist_ok=True)
    est_path = RESULTS / f"{args.tag}_estimates.csv"; coh_path = RESULTS / f"{args.tag}_cohort.csv"; err_path = RESULTS / f"{args.tag}_errors.jsonl"
    done = set(pd.read_csv(coh_path, usecols=["cohort_id"])["cohort_id"].tolist()) if coh_path.exists() else set()
    todo = [i for i in ids if i not in done]
    print(f"{len(ids)} cohorts, {len(done)} done, {len(todo)} to run; seeds={seeds}; adj={adjs}; n_jobs={args.n_jobs}", flush=True)
    t0 = time.perf_counter()
    for s in range(0, len(todo), args.chunk):
        chunk = todo[s:s + args.chunk]
        outs = Parallel(n_jobs=args.n_jobs)(delayed(process_cohort)(i, seeds, adjs) for i in chunk)
        rows = [r for o in outs if not o["error"] for r in o["rows"]]; cohs = [o["cohort"] for o in outs if not o["error"]]; errs = [o for o in outs if o["error"]]
        if rows: pd.DataFrame(rows).to_csv(est_path, mode="a", header=not est_path.exists(), index=False)
        if cohs: pd.DataFrame(cohs).to_csv(coh_path, mode="a", header=not coh_path.exists(), index=False)
        if errs:
            with err_path.open("a") as fh:
                for o in errs: fh.write(json.dumps({"cohort_id": o["cohort_id"], "error": o["error"]}) + "\n")
        print(f"  {s + len(chunk)}/{len(todo)} done in {(time.perf_counter()-t0)/60:.1f} min ({len(errs)} errors)", flush=True)
    print(f"FINISHED {len(todo)} cohorts in {(time.perf_counter()-t0)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
