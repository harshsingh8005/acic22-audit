"""Shared data access for the confirmatory analysis (frozen plan §1–§2).

Loads one ACIC 2022 Track-2 cohort and reduces it to (X, T, Y) with the SAME
difference-in-differences change score as the separate reference implementation
used by software checks V1 and V2 (see README, "Software checks"), plus:

* ``pre_trend`` = Y(year 2) − Y(year 1) at the practice level — the designed
  Scenario-A confounder ("pre-intervention change in outcome").  X_full includes
  it; X_reduced (the reproducer's original set) does not.
* practice-level SATT truths for the treated practices (PEHE-on-treated).
* the DGP factor labels for every cohort (RCT cohorts have NaN labels).

Nothing here estimates anything.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

_HERE = Path(__file__).resolve().parent

# Where the ACIC 2022 Track-2 release is unpacked.  Override with the
# ACIC22_DATA_ROOT environment variable; see README ("Data").
DATA_ROOT = Path(os.environ.get("ACIC22_DATA_ROOT", _HERE.parent / "data"))
TRACK2 = DATA_ROOT / "track2"
TRUTHS = DATA_ROOT / "ACIC_estimand_truths.csv"

# Optional.  A separate reference implementation of the two baseline estimators,
# used by software checks V1 and V2 only and not part of this release.  Point
# ACIC22_REPRODUCER_ROOT at it to enable those two checks; left unset, V1 and V2
# are skipped and V3 to V6 still run.  See README ("Software checks").
_reproducer = os.environ.get("ACIC22_REPRODUCER_ROOT")
REPO_TRACK2 = Path(_reproducer) if _reproducer else None

N_COHORTS = 3400
FOLD_SEED = 42          # frozen plan §2
N_FOLDS = 5
SAMPLE_SEED = 20260818  # frozen plan §8


@dataclass(frozen=True)
class Cohort:
    cohort_id: int
    X_full: np.ndarray
    X_reduced: np.ndarray
    feature_names_full: tuple
    feature_names_reduced: tuple
    T: np.ndarray
    Y: np.ndarray                 # DiD change score (post mean − pre mean)
    pre_trend: np.ndarray         # Y2 − Y1 (also the negative-control outcome)
    practice_ids: np.ndarray
    true_satt: float
    true_catt_treated: np.ndarray  # per treated practice, aligned to practice_ids[T==1]
    factors: dict

    @property
    def n(self) -> int:
        return int(self.T.shape[0])


@lru_cache(maxsize=1)
def _truth_frame() -> pd.DataFrame:
    df = pd.read_csv(TRUTHS, dtype={"dataset.num": str})
    df["cohort_id"] = df["dataset.num"].astype(int)
    return df


@lru_cache(maxsize=1)
def factor_table() -> pd.DataFrame:
    """One row per cohort: design factors (NaN → 'RCT')."""
    df = _truth_frame().drop_duplicates("cohort_id")[
        ["cohort_id", "Confounding Strength", "Confounding Source",
         "Impact Heterogeneity", "Idiosyncrasy of Impacts"]
    ].copy()
    df = df.rename(columns={
        "Confounding Strength": "strength", "Confounding Source": "source",
        "Impact Heterogeneity": "heterogeneity", "Idiosyncrasy of Impacts": "idiosyncrasy",
    })
    df["strength"] = df["strength"].fillna("RCT")
    df["source"] = df["source"].fillna("RCT")
    df["source"] = df["source"].str.replace("Scenario ", "", regex=False)
    df["cell"] = (df["strength"] + "|" + df["source"] + "|" +
                  df["heterogeneity"] + "|" + df["idiosyncrasy"])
    return df.set_index("cohort_id").sort_index()


@lru_cache(maxsize=1)
def _overall_truth() -> dict:
    df = _truth_frame()
    ov = df[df["variable"] == "Overall"]
    return dict(zip(ov["cohort_id"].astype(int), ov["SATT"].astype(float)))


def _practice_truth(cohort_id: int) -> pd.Series:
    df = _truth_frame()
    sub = df[(df["cohort_id"] == cohort_id) & df["id.practice"].notna()]
    return pd.Series(sub["SATT"].astype(float).values, index=sub["id.practice"].astype(int).values)


def _wmean(v: np.ndarray, w: np.ndarray) -> float:
    tot = float(w.sum())
    return float(np.dot(v, w) / tot) if tot > 0 else float(v.mean())


def load_cohort(cohort_id: int) -> Cohort:
    tag = f"{cohort_id:04d}"
    practice = pd.read_csv(TRACK2 / "practice" / f"acic_practice_{tag}.csv")
    panel = pd.read_csv(TRACK2 / "practice_year" / f"acic_practice_year_{tag}.csv")

    mix_cols = [c for c in panel.columns if c.endswith("_avg")]
    recs = []
    for pid, g in panel.groupby("id.practice", sort=True):
        pre = g[g["post"] == 0]
        post = g[g["post"] == 1]
        if len(pre) == 0 or len(post) == 0:
            continue
        w_pre = pre["n.patients"].to_numpy(dtype=np.float64)
        w_post = post["n.patients"].to_numpy(dtype=np.float64)
        y_pre = _wmean(pre["Y"].to_numpy(dtype=np.float64), w_pre)
        y_post = _wmean(post["Y"].to_numpy(dtype=np.float64), w_post)
        rec = {
            "id.practice": int(pid),
            "Z": int(g["Z"].iloc[0]),
            "dY": y_post - y_pre,
            "pre_Y": y_pre,
            "pre_n_patients": float(np.mean(w_pre)),
        }
        for c in mix_cols:
            rec[f"pre_{c}"] = _wmean(pre[c].to_numpy(dtype=np.float64), w_pre)
        # designed Scenario-A confounder: pre-intervention change in outcome
        y_by_year = dict(zip(pre["year"].astype(int), pre["Y"].astype(float)))
        rec["pre_trend"] = float(y_by_year[2] - y_by_year[1])
        recs.append(rec)

    reduced = pd.DataFrame.from_records(recs)
    merged = reduced.merge(practice, on="id.practice", how="inner")

    T = merged["Z"].to_numpy(dtype=np.int_)
    Y = merged["dY"].to_numpy(dtype=np.float64)
    pre_trend = merged["pre_trend"].to_numpy(dtype=np.float64)
    pids = merged["id.practice"].to_numpy(dtype=np.int_)

    # --- X_reduced: exactly the reproducer's covariate construction ---------
    cov_cols_reduced = [c for c in merged.columns if c not in ("id.practice", "Z", "dY", "pre_trend")]
    Xr_df = pd.get_dummies(merged[cov_cols_reduced], drop_first=True, dtype=np.float64)
    # --- X_full = X_reduced + pre_trend --------------------------------------
    Xf_df = Xr_df.copy()
    Xf_df["pre_trend"] = pre_trend

    truth = _overall_truth()[cohort_id]
    ptruth = _practice_truth(cohort_id)
    treated_ids = pids[T == 1]
    catt = ptruth.reindex(treated_ids).to_numpy(dtype=np.float64)

    fac = factor_table().loc[cohort_id].to_dict()
    return Cohort(
        cohort_id=cohort_id,
        X_full=Xf_df.to_numpy(dtype=np.float64),
        X_reduced=Xr_df.to_numpy(dtype=np.float64),
        feature_names_full=tuple(Xf_df.columns),
        feature_names_reduced=tuple(Xr_df.columns),
        T=T, Y=Y, pre_trend=pre_trend, practice_ids=pids,
        true_satt=float(truth), true_catt_treated=catt, factors=fac,
    )


def make_folds(T: np.ndarray, seed: int = FOLD_SEED, n_folds: int = N_FOLDS):
    """The single shared fold assignment (frozen plan §2). Depends only on T."""
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    return [(tr, te) for tr, te in skf.split(np.zeros((len(T), 1)), T)]


def stratified_cohort_sample(per_cell: int, seed: int = SAMPLE_SEED) -> list[int]:
    """Nested stratified sample: `per_cell` cohorts from each of the 17 design cells."""
    ft = factor_table()
    rng = np.random.default_rng(seed)
    out = []
    for cell, grp in ft.groupby("cell", sort=True):
        ids = np.sort(grp.index.to_numpy())
        perm = rng.permutation(ids)
        out.extend(int(i) for i in perm[:per_cell])
    return sorted(out)


def all_cohorts() -> list[int]:
    return list(range(1, N_COHORTS + 1))
