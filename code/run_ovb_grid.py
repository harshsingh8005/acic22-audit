"""OVB coverage grid for figure F4 (frozen plan §14): for E1 AIPW-ATT on X_reduced,
evaluate whether the OVB-adjusted 95% CI at each (cf_y, cf_d) on a fixed grid
(rho = 1) covers the true SATT, for every cohort; also the smallest diagonal
cf = cf_y = cf_d at which the adjusted CI covers the truth ("cf needed").
Pure recomputation of pre-specified quantities (no new estimator, no tuning).
Output: results/tables/ovb_grid_E1_reduced.csv (one row per cohort, one column per grid cell)
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acic_common import all_cohorts, load_cohort, make_folds  # noqa: E402
from estimators import aipw_att  # noqa: E402
from sensitivity import ovb_bounds, ovb_elements_att  # noqa: E402

GRID = [0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50]
TAB = HERE.parent / "results" / "tables"


def one(cid, adj="reduced"):
    c = load_cohort(cid); folds = make_folds(c.T)
    X = c.X_reduced if adj == "reduced" else c.X_full
    r, nu = aipw_att(X, c.T, c.Y, folds)
    el = ovb_elements_att(c.T, c.Y, nu.e_hat, nu.mu0_hat, nu.mu1_hat, r.att, r.infl)
    row = {"cohort_id": cid, "adj": adj, "true_satt": c.true_satt, "att": r.att, "se": r.se,
           **{k: c.factors[k] for k in ("strength", "source", "heterogeneity", "idiosyncrasy")}}
    for cy in GRID:
        for cd in GRID:
            b = ovb_bounds(el, cy, min(cd, 0.999), 1.0)
            row[f"cov_{cy}_{cd}"] = bool(b["ci_lower"] <= c.true_satt <= b["ci_upper"])
    # smallest diagonal cf at which the adjusted CI covers the truth
    def gap(cf):
        b = ovb_bounds(el, cf, cf, 1.0)
        return min(c.true_satt - b["ci_lower"], b["ci_upper"] - c.true_satt)  # >0 when covered
    if gap(0.0) >= 0:
        row["cf_needed_diag"] = 0.0
    elif gap(0.999) < 0:
        row["cf_needed_diag"] = np.nan  # not covered even at cf ~ 1
    else:
        row["cf_needed_diag"] = float(brentq(gap, 0.0, 0.999, xtol=1e-6))
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--n-jobs", type=int, default=8); ap.add_argument("--adj", default="reduced")
    args = ap.parse_args(argv)
    t0 = time.perf_counter()
    rows = Parallel(n_jobs=args.n_jobs)(delayed(one)(i, args.adj) for i in all_cohorts())
    df = pd.DataFrame(rows)
    df.to_csv(TAB / f"ovb_grid_E1_{args.adj}.csv", index=False)
    print(f"wrote {len(df)} rows in {(time.perf_counter()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
