"""Pre-specified estimators (frozen plan §4): OLS-ATT (E0), AIPW-ATT (E1),
DR-ATT with Crump trim (E2), ridge T-learner CATE (E6).

E1/E2 are re-implementations of a separate reference implementation's ``aipw_att`` /
``dr_att`` that additionally return the cross-fitted nuisances (needed by the
sensitivity layer and matching).  They must agree with the reproducer to 1e-9
(software check V1 in validate_software.py).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

Z95 = 1.959963984540054
CLIP_AIPW = (0.02, 0.98)
CLIP_DR = (0.01, 0.99)
TRIM_ALPHA = 0.10


def make_classifier():
    return Pipeline([("scale", StandardScaler()),
                     ("clf", LogisticRegression(C=0.5, max_iter=1000))])


def make_regressor(random_state: int = 42):
    return Ridge(alpha=10.0, random_state=random_state)


@dataclass
class Nuisances:
    e_hat: np.ndarray      # clipped cross-fitted propensity
    e_raw: np.ndarray      # unclipped
    mu0_hat: np.ndarray    # E[Y|X,T=0] cross-fitted
    mu1_hat: np.ndarray    # E[Y|X,T=1] cross-fitted (extra; not used by the point estimate)


def cross_fit_nuisances(X, T, Y, folds, clip, random_state: int = 42) -> Nuisances:
    n = X.shape[0]
    e = np.empty(n); mu0 = np.empty(n); mu1 = np.empty(n)
    for tr, te in folds:
        clf = make_classifier(); clf.fit(X[tr], T[tr])
        e[te] = clf.predict_proba(X[te])[:, 1]
        ctrl = tr[T[tr] == 0]
        reg0 = make_regressor(random_state); reg0.fit(X[ctrl], Y[ctrl])
        mu0[te] = reg0.predict(X[te])
        trt = tr[T[tr] == 1]
        reg1 = make_regressor(random_state); reg1.fit(X[trt], Y[trt])
        mu1[te] = reg1.predict(X[te])
    e_raw = e.copy()
    np.clip(e, clip[0], clip[1], out=e)
    return Nuisances(e_hat=e, e_raw=e_raw, mu0_hat=mu0, mu1_hat=mu1)


@dataclass
class AttResult:
    att: float
    se: float
    ci_lo: float
    ci_hi: float
    infl: np.ndarray = field(repr=False)     # per-unit influence values (mean≈0), se = sd/sqrt(n)
    keep: np.ndarray = field(repr=False)     # boolean mask of units used (trimming)
    n_treated: int = 0
    n_trimmed: int = 0
    prop_min: float = float("nan")
    prop_max: float = float("nan")


def aipw_att_from_nuisances(T, Y, e_hat, mu0_hat, keep=None) -> AttResult:
    """Hájek AIPW-ATT and EIF standard error — identical arithmetic to the
    reproducer's aipw_att / dr_att given the (possibly trimmed) sample."""
    if keep is None:
        keep = np.ones(len(T), dtype=bool)
    Tt, Yt, et, m0 = T[keep], Y[keep], e_hat[keep], mu0_hat[keep]
    m = Tt.shape[0]
    n_treated = int((Tt == 1).sum())
    p_hat = n_treated / m
    resid = Yt - m0
    w = et / (1.0 - et)
    is_ctrl = Tt == 0
    treated_term = float(resid[Tt == 1].mean())
    control_term = float(np.dot(w[is_ctrl], resid[is_ctrl]) / w[is_ctrl].sum())
    att = treated_term - control_term
    infl = np.where(Tt == 1, resid, -w * resid) / p_hat - (Tt / p_hat) * att
    se = float(infl.std(ddof=1) / np.sqrt(m))
    return AttResult(att=att, se=se, ci_lo=att - Z95 * se, ci_hi=att + Z95 * se,
                     infl=infl, keep=keep, n_treated=n_treated,
                     n_trimmed=int((~keep).sum()),
                     prop_min=float(e_hat.min()), prop_max=float(e_hat.max()))


def aipw_att(X, T, Y, folds, random_state: int = 42):
    """E1: cross-fitted AIPW-ATT (reproducer-identical) + nuisances."""
    if (T == 1).sum() < len(folds) or (T == 0).sum() < len(folds):
        raise ValueError("each arm must have >= k_folds units")
    nu = cross_fit_nuisances(X, T, Y, folds, CLIP_AIPW, random_state)
    res = aipw_att_from_nuisances(T, Y, nu.e_hat, nu.mu0_hat)
    return res, nu


def dr_att(X, T, Y, folds, random_state: int = 42, trim_alpha: float = TRIM_ALPHA):
    """E2: DR-ATT with Crump overlap trim (reproducer-identical) + nuisances."""
    if (T == 1).sum() < len(folds) or (T == 0).sum() < len(folds):
        raise ValueError("each arm must have >= k_folds units")
    nu = cross_fit_nuisances(X, T, Y, folds, CLIP_DR, random_state)
    keep = (nu.e_hat >= trim_alpha) & (nu.e_hat <= 1.0 - trim_alpha)
    Tt = T[keep]
    if (Tt == 1).sum() < 5 or (Tt == 0).sum() < 5:
        raise ValueError("overlap trim left too few units")
    res = aipw_att_from_nuisances(T, Y, nu.e_hat, nu.mu0_hat, keep=keep)
    return res, nu


@dataclass
class OlsResult:
    att: float
    se: float
    ci_lo: float
    ci_hi: float
    t_stat: float
    dof: int
    resid_var: float


def ols_att(X, T, Y) -> OlsResult:
    """E0: OLS of Y on [1, T, X] with HC1 SE; ATT = coefficient on T.
    t_stat/dof (classical, homoskedastic) are what the linear RV uses (sensemakr)."""
    D = np.column_stack([np.ones(len(T)), T.astype(float), X])
    fit_hc1 = sm.OLS(Y, D).fit(cov_type="HC1")
    fit_ols = sm.OLS(Y, D).fit()
    att = float(fit_hc1.params[1]); se = float(fit_hc1.bse[1])
    return OlsResult(att=att, se=se, ci_lo=att - Z95 * se, ci_hi=att + Z95 * se,
                     t_stat=float(fit_ols.tvalues[1]), dof=int(fit_ols.df_resid),
                     resid_var=float(fit_ols.mse_resid))


def ridge_tlearner_cate(X, T, Y, random_state: int = 42):
    """E6: in-sample ridge T-learner: mu0 fit on controls, mu1 on treated.
    Uses exactly E1's outcome-model family (Ridge alpha=10, unscaled), so the
    linear comparator is E1's own nuisance model used as a T-learner."""
    r0 = make_regressor(random_state).fit(X[T == 0], Y[T == 0])
    r1 = make_regressor(random_state).fit(X[T == 1], Y[T == 1])
    return r0.predict(X), r1.predict(X)
