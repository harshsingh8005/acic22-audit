"""Software validation checks V1–V6 (frozen plan §10). Not hypotheses.
Writes results/tables/software_validation.json and prints a summary."""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acic_common import REPO_TRACK2, load_cohort, make_folds, N_COHORTS  # noqa: E402
from estimators import aipw_att, dr_att, ols_att, aipw_att_from_nuisances  # noqa: E402
from sensitivity import (evalue_from_rr, ovb_bounds, ovb_elements_att, robustness_values,  # noqa: E402
                         rosenbaum_p_upper, rv_linear)

OUT = HERE.parent / "results" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
report = {}

# V1 and V2 compare this code against a separate reference implementation of the
# two baseline estimators.  That implementation is not part of this release, so
# these two checks run only when ACIC22_REPRODUCER_ROOT points at it; otherwise
# they are skipped.  V3 to V6 are self-contained and always run.  See README.
HAVE_REPRODUCER = False
if REPO_TRACK2 is not None and REPO_TRACK2.exists():
    sys.path.insert(0, str(REPO_TRACK2))
    try:
        from src import (aipw_att as repo_aipw, dr_att as repo_dr,  # noqa: E402
                         load_track2_cohort as repo_load)
        HAVE_REPRODUCER = True
    except ImportError as exc:
        print(f"V1/V2: reference implementation found but not importable ({exc}); skipping.")
else:
    print("V1/V2: ACIC22_REPRODUCER_ROOT not set; skipping (V3-V6 still run).")

# ---------------------------------------------------------------- V1 identity
if HAVE_REPRODUCER:
    print("V1: re-implementation vs reference implementation (X_reduced) ...")
    v1 = []
    for cid in (1, 500, 1000, 2000, 3400):
        c = load_cohort(cid); rc = repo_load(cid)
        folds = make_folds(c.T)
        same_X = c.X_reduced.shape == rc.X.shape and np.allclose(c.X_reduced, rc.X) and tuple(c.feature_names_reduced) == tuple(rc.feature_names)
        r1, _ = aipw_att(c.X_reduced, c.T, c.Y, folds); a = repo_aipw(rc.X, rc.T, rc.Y, random_state=42)
        r2, _ = dr_att(c.X_reduced, c.T, c.Y, folds); d = repo_dr(rc.X, rc.T, rc.Y, random_state=42)
        v1.append({"cohort": cid, "same_X": bool(same_X), "same_truth": abs(c.true_satt - rc.true_ate) < 1e-12,
                   "d_aipw_att": abs(r1.att - a.att), "d_aipw_se": abs(r1.se - a.se),
                   "d_dr_att": abs(r2.att - d.att), "d_dr_se": abs(r2.se - d.se), "same_ntrim": r2.n_trimmed == d.n_trimmed})
    v1df = pd.DataFrame(v1)
    report["V1"] = {"pass": bool(v1df["same_X"].all() and v1df["same_truth"].all() and
                                 (v1df[["d_aipw_att", "d_aipw_se", "d_dr_att", "d_dr_se"]].max().max() < 1e-9) and v1df["same_ntrim"].all()),
                    "detail": v1df.to_dict("records")}
    print(v1df.to_string()); print("V1 pass:", report["V1"]["pass"])
else:
    report["V1"] = {"skipped": "ACIC22_REPRODUCER_ROOT not available"}

# --------------------------------------------- V2 canonical 100-cohort contract
if HAVE_REPRODUCER:
    print("V2: canonical 100-cohort seed-42 sample vs expected_results.json ...")
    exp = json.loads((REPO_TRACK2 / "expected_results.json").read_text())["reproducer_targets"]
    all_ids = np.arange(1, N_COHORTS + 1)
    canon = sorted(int(i) for i in np.random.default_rng(42).permutation(all_ids)[:100])
    rows = []
    t0 = time.perf_counter()
    for cid in canon:
        c = load_cohort(cid); folds = make_folds(c.T)
        r1, _ = aipw_att(c.X_reduced, c.T, c.Y, folds)
        r2, _ = dr_att(c.X_reduced, c.T, c.Y, folds)
        rows.append({"cid": cid, "true": c.true_satt, "a": r1.att, "alo": r1.ci_lo, "ahi": r1.ci_hi,
                     "d": r2.att, "dlo": r2.ci_lo, "dhi": r2.ci_hi})
    df = pd.DataFrame(rows)
    def agg(e, lo, hi, tr):
        err = e - tr
        return {"bias_mean": float(err.mean()), "rmse": float(np.sqrt((err ** 2).mean())),
                "coverage_95": float(((tr >= lo) & (tr <= hi)).mean()), "ci_width_mean": float((hi - lo).mean())}
    v2 = {"aipw": agg(df.a, df.alo, df.ahi, df.true), "dr_att": agg(df.d, df.dlo, df.dhi, df.true), "wall_s": time.perf_counter() - t0}
    ok = True
    for m in ("aipw", "dr_att"):
        for k, spec in exp[m].items():
            within = abs(v2[m][k] - spec["target"]) <= spec["tol"]
            ok &= within
            print(f"  {m}.{k}: observed {v2[m][k]:.3f} target {spec['target']} +/- {spec['tol']} -> {'PASS' if within else 'FAIL'}")
    report["V2"] = {"pass": bool(ok), "observed": v2, "targets": exp}
    print("V2 pass:", ok, f"({v2['wall_s']:.0f}s for 100 cohorts x 2 estimators, single process)")
else:
    report["V2"] = {"skipped": "ACIC22_REPRODUCER_ROOT not available"}

# ---------------------------------------------------------------- V3 E-value
e = evalue_from_rr(3.9)
report["V3"] = {"pass": abs(e - 7.26) < 0.01, "evalue_rr3.9": e}
print(f"V3 E-value(RR=3.9) = {e:.3f} (expected 7.26) pass:", report["V3"]["pass"])

# --------------------------------------------------------------- V4 Rosenbaum
rng = np.random.default_rng(0)
D = rng.normal(0.4, 1.0, 120)
D = np.round(D, 3)
p_mine = rosenbaum_p_upper(D, 1.0)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    p_scipy = wilcoxon(D, alternative="greater", method="approx", correction=False).pvalue
gammas = [1.0, 1.2, 1.5, 2.0, 3.0]
ps = [rosenbaum_p_upper(D, g) for g in gammas]
mono = all(ps[i] <= ps[i + 1] + 1e-12 for i in range(len(ps) - 1))
Dpos = np.arange(1, 11, dtype=float)      # all positive, no ties
r = np.arange(1, 11); hand = float(__import__("scipy").stats.norm.sf(r.sum() / np.sqrt((r ** 2).sum())))
report["V4"] = {"pass": bool(abs(p_mine - p_scipy) < 1e-6 and mono and abs(rosenbaum_p_upper(Dpos, 1.0) - hand) < 1e-12),
                "p_gamma1_mine": p_mine, "p_gamma1_scipy": float(p_scipy), "p_by_gamma": dict(zip(gammas, ps)), "hand_case": hand}
print(f"V4 p(G=1) mine={p_mine:.6f} scipy={p_scipy:.6f}; monotone={mono}; pass:", report["V4"]["pass"])

# --------------------------------------------------- V5 RV: bound hits 0 at RV; DoWhy
c = load_cohort(7); folds = make_folds(c.T)
r1, nu = aipw_att(c.X_full, c.T, c.Y, folds)
el = ovb_elements_att(c.T, c.Y, nu.e_hat, nu.mu0_hat, nu.mu1_hat, r1.att, r1.infl)
rv, rva = robustness_values(el)
b = ovb_bounds(el, rv, rv, 1.0)
hit0 = b["theta_lower"] if r1.att > 0 else b["theta_upper"]
b2 = ovb_bounds(el, rva, rva, 1.0) if 0 < rva < 1 else None
ci0 = (b2["ci_lower"] if r1.att > 0 else b2["ci_upper"]) if b2 else 0.0
print(f"V5a cohort 7 E1: att={r1.att:.3f} se={r1.se:.3f} RV={rv:.4f} RVa={rva:.4f}; bound at RV -> {hit0:.2e}; CI at RVa -> {ci0:.2e}")
v5a = abs(hit0) < 1e-6 and abs(ci0) < 1e-5

# DoWhy linear-partial-R2 cross-check on 10 cohorts (X_reduced, first 12 features to keep DoWhy fast)
import logging
logging.getLogger("dowhy").setLevel(logging.ERROR)
from dowhy import CausalModel  # noqa: E402
v5b = []
for cid in range(1, 11):
    c = load_cohort(cid)
    X = c.X_reduced; names = [f"w{i}" for i in range(X.shape[1])]
    df = pd.DataFrame(X, columns=names); df["T"] = c.T; df["Y"] = c.Y
    o = ols_att(X, c.T, c.Y)
    mine_rv, mine_rva = rv_linear(o.t_stat, o.dof, 1.0, 1.0), rv_linear(o.t_stat, o.dof, 1.0, 0.05)
    model = CausalModel(data=df, treatment="T", outcome="Y", common_causes=names)
    ident = model.identify_effect(proceed_when_unidentifiable=True)
    est = model.estimate_effect(ident, method_name="backdoor.linear_regression", test_significance=False)
    ref = model.refute_estimate(ident, est, method_name="add_unobserved_common_cause",
                                simulation_method="linear-partial-R2", benchmark_common_causes=[names[0]],
                                effect_fraction_on_treatment=[1.0], plot=False)
    # the analyzer object carries the RVs
    an = getattr(ref, "refutation_result", None)
    dw_rv = getattr(an, "RVq", None) if an is not None else None
    dw_rva = None  # bound unconditionally: the fallback loop below may not match
    if dw_rv is None:  # fall back: search attributes
        for attr in ("stats", "sensitivity_stats"):
            s = getattr(ref, attr, None)
            if isinstance(s, dict) and "robustness_value" in s:
                dw_rv = s["robustness_value"]; dw_rva = s.get("robustness_value_alpha")
                break
    else:
        dw_rva = getattr(an, "RV_qalpha", None)
    v5b.append({"cohort": cid, "att_ols": o.att, "dowhy_est": float(est.value), "rv_mine": mine_rv, "rv_dowhy": dw_rv,
                "rva_mine": mine_rva, "rva_dowhy": dw_rva})
v5df = pd.DataFrame(v5b)
print(v5df.to_string())
d_rv = (v5df["rv_mine"] - v5df["rv_dowhy"].astype(float)).abs().max() if v5df["rv_dowhy"].notna().all() else np.inf
d_rva = (v5df["rva_mine"] - v5df["rva_dowhy"].astype(float)).abs().max() if v5df["rva_dowhy"].notna().all() else np.inf
d_est = (v5df["att_ols"] - v5df["dowhy_est"]).abs().max()
report["V5"] = {"pass": bool(v5a and d_rv < 1e-6 and d_rva < 1e-6 and d_est < 1e-6),
                "bound_at_rv": hit0, "ci_at_rva": ci0, "max_dRV_vs_dowhy": float(d_rv), "max_dRVa_vs_dowhy": float(d_rva),
                "max_dEst_vs_dowhy": float(d_est), "detail": v5df.to_dict("records")}
print("V5 pass:", report["V5"]["pass"], f"(max |dRV|={d_rv:.2e}, |dRVa|={d_rva:.2e}, |dEst|={d_est:.2e})")

# ------------------------------------------------------------ V6 neural smoke
from neural import fit_predict  # noqa: E402
rng = np.random.default_rng(1)
n, p = 500, 10
Xs = rng.normal(size=(n, p))
ps_ = 1 / (1 + np.exp(-(0.5 * Xs[:, 0] - 0.5 * Xs[:, 1])))
Ts = rng.binomial(1, ps_)
tau = 2.0
Ys = 1.0 + Xs[:, 0] + 0.5 * Xs[:, 1] + tau * Ts + rng.normal(scale=1.0, size=n)
t0 = time.perf_counter()
mu0, mu1, _, info_t = fit_predict(Xs, Ts, Ys, seed=0, dragon=False)
t_tar = time.perf_counter() - t0
t0 = time.perf_counter()
mu0d, mu1d, ed, info_d = fit_predict(Xs, Ts, Ys, seed=0, dragon=True)
t_dra = time.perf_counter() - t0
att_t = float((mu1 - mu0)[Ts == 1].mean()); att_d = float((mu1d - mu0d)[Ts == 1].mean())
res_d = aipw_att_from_nuisances(Ts, Ys, np.clip(ed, 0.02, 0.98), mu0d)
report["V6"] = {"pass": bool(abs(att_t - tau) < 0.5 and abs(att_d - tau) < 0.5 and abs(res_d.att - tau) < 2 * max(res_d.se, 0.1)),
                "tarnet_att": att_t, "dragonnet_att": att_d, "dragon_dr_att": res_d.att, "dragon_dr_se": res_d.se,
                "true": tau, "sec_tarnet": t_tar, "sec_dragonnet": t_dra, "epochs": [info_t, info_d]}
print(f"V6 TARNet att={att_t:.3f} ({t_tar:.1f}s, {info_t}); DragonNet att={att_d:.3f} ({t_dra:.1f}s, {info_d}); "
      f"Dragon-DR {res_d.att:.3f}+/-{res_d.se:.3f}; true 2.0; pass:", report["V6"]["pass"])

(OUT / "software_validation.json").write_text(json.dumps(report, indent=2, default=float))
_ran = {k: v for k, v in report.items() if "pass" in v}
_skipped = sorted(k for k, v in report.items() if "pass" not in v)
print(f"ALL PASS: {all(v['pass'] for v in _ran.values())} ({len(_ran)} of {len(report)} checks run)"
      + (f"; skipped {', '.join(_skipped)}" if _skipped else ""))
