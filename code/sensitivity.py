"""Unified sensitivity layer (frozen plan §5).

Three families, all computed from an estimate + its nuisances:

5.1  E-values (VanderWeele & Ding 2017), continuous outcome via d -> RR = exp(0.91 d).
5.2  Omitted-variable-bias bounds / robustness values in the Chernozhukov, Cinelli,
     Newey, Sharma & Syrgkanis (2022) "long story short" parametrisation as used
     by DoubleML: |bias| <= rho * sqrt(cf_y * cf_d/(1-cf_d)) * sqrt(sigma2 * nu2).
     Bound SEs by influence functions (delta method); RV, RVa; linear
     (sensemakr) forms for OLS.
5.3  Rosenbaum Gamma bounds on 1:1 propensity-matched pairs (Wilcoxon signed
     rank, rbounds::psens logic) -> Gamma* = sup{Gamma : p_upper(Gamma) < 0.05}.
5.4  Measured strength of a withheld confounder from long/short nuisances.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm, rankdata, t as student_t

Z95 = 1.959963984540054


# ---------------------------------------------------------------- 5.1 E-values
def evalue_from_rr(rr: float) -> float:
    rr = float(rr)
    if rr < 1.0:
        rr = 1.0 / rr
    return rr + np.sqrt(rr * (rr - 1.0))


def evalues_continuous(att: float, ci_lo: float, ci_hi: float, sd_y: float):
    """E_pt and E_CI for a continuous outcome (VanderWeele & Ding 2017:
    d = mean difference / SD; RR ~ exp(0.91 d)); E_CI = 1 if the CI covers 0."""
    d = att / sd_y
    e_pt = evalue_from_rr(np.exp(0.91 * abs(d)))
    if ci_lo <= 0.0 <= ci_hi:
        e_ci = 1.0
    else:
        lim = ci_lo if att > 0 else ci_hi          # limit closest to the null
        e_ci = evalue_from_rr(np.exp(0.91 * abs(lim / sd_y)))
    return float(e_pt), float(e_ci)


# ------------------------------------------------------- 5.2 OVB (nonparametric)
@dataclass
class OvbElements:
    theta: float
    n: int
    sigma2: float
    nu2: float
    psi_theta: np.ndarray   # influence values for theta (mean ~ 0)
    psi_sigma2: np.ndarray
    psi_nu2: np.ndarray
    nu2_dr: float = float("nan")      # the debiased (orthogonal) estimate before any fallback
    nu2_fallback: bool = False        # True if nu2_dr <= 0 and the plug-in mean(alpha^2) was used


def ovb_elements_att(T, Y, e_hat, mu0_hat, mu1_hat, theta, infl, keep=None) -> OvbElements:
    """sigma2 = E[(Y - g(T,X))^2] with g = T*mu1 + (1-T)*mu0;
    nu2 = 2 E[m(W;alpha)] - E[alpha^2] with the ATT Riesz representer
    alpha = (T - (1-T) e/(1-e))/p and m(W;alpha) = T/(p^2 (1-e)).
    If the debiased nu2 is not positive (extreme propensities), fall back to the
    non-orthogonal plug-in mean(alpha^2) — the same rule as DoubleML's reference
    implementation ("The estimated nu2 is not positive. Re-estimation based on
    riesz representer (non-orthogonal).") — and flag it."""
    if keep is None:
        keep = np.ones(len(T), dtype=bool)
    Tt, Yt, et = T[keep], Y[keep], e_hat[keep]
    m0, m1 = mu0_hat[keep], mu1_hat[keep]
    n = len(Tt)
    p = Tt.mean()
    g = Tt * m1 + (1 - Tt) * m0
    s2_el = (Yt - g) ** 2
    sigma2 = float(s2_el.mean())
    alpha = (Tt - (1 - Tt) * et / (1 - et)) / p
    m_alpha = Tt / (p ** 2 * (1 - et))
    nu2_el = 2.0 * m_alpha - alpha ** 2
    nu2_dr = float(nu2_el.mean())
    fallback = False
    if nu2_dr <= 0.0:
        fallback = True
        nu2_el = alpha ** 2
    nu2 = float(nu2_el.mean())
    return OvbElements(theta=float(theta), n=n, sigma2=sigma2, nu2=nu2,
                       psi_theta=np.asarray(infl, dtype=float),
                       psi_sigma2=s2_el - sigma2, psi_nu2=nu2_el - nu2,
                       nu2_dr=nu2_dr, nu2_fallback=fallback)


def ovb_bounds(el: OvbElements, cf_y: float, cf_d: float, rho: float = 1.0, level: float = 0.95):
    """Bias bound and OVB-adjusted CI at scenario (cf_y, cf_d, rho)."""
    z = norm.ppf(0.5 + level / 2)
    k = abs(rho) * np.sqrt(cf_y * cf_d / (1.0 - cf_d)) if (cf_y > 0 and 0 < cf_d < 1) else 0.0
    S = np.sqrt(max(el.sigma2 * el.nu2, 0.0))
    B = k * S
    if S > 0 and k > 0:
        psi_S2 = el.sigma2 * el.psi_nu2 + el.nu2 * el.psi_sigma2
        psi_bias = (k / (2.0 * S)) * psi_S2
    else:
        psi_bias = np.zeros_like(el.psi_theta)
    psi_lower = el.psi_theta - psi_bias
    psi_upper = el.psi_theta + psi_bias
    se_lower = float(psi_lower.std(ddof=1) / np.sqrt(el.n))
    se_upper = float(psi_upper.std(ddof=1) / np.sqrt(el.n))
    theta_lower = el.theta - B
    theta_upper = el.theta + B
    return {
        "bias_bound": float(B),
        "theta_lower": float(theta_lower), "theta_upper": float(theta_upper),
        "ci_lower": float(theta_lower - z * se_lower),
        "ci_upper": float(theta_upper + z * se_upper),
        "se_lower": se_lower, "se_upper": se_upper,
    }


def robustness_values(el: OvbElements, level: float = 0.95):
    """RV: cf_y = cf_d = rv, rho = 1 drives the bias bound to |theta|.
    RVa: smallest rv at which the OVB-adjusted CI reaches 0 (bisection)."""
    S = np.sqrt(max(el.sigma2 * el.nu2, 0.0))
    if S <= 0:
        return float("nan"), float("nan")
    f = abs(el.theta) / S
    rv = 0.5 * (np.sqrt(f ** 4 + 4 * f ** 2) - f ** 2)

    def ci_null_side(rv_):
        b = ovb_bounds(el, rv_, rv_, 1.0, level)
        return b["ci_lower"] if el.theta > 0 else -b["ci_upper"]

    if ci_null_side(0.0) <= 0.0:
        rva = 0.0                       # CI already covers 0 without any confounding
    elif ci_null_side(1 - 1e-9) > 0.0:
        rva = 1.0
    else:
        rva = float(brentq(ci_null_side, 0.0, 1 - 1e-9, xtol=1e-10))
    return float(rv), float(rva)


# ---------------------------------------------------- 5.2 linear (sensemakr) RV
def rv_linear(t_stat: float, dof: int, q: float = 1.0, alpha: float = 1.0) -> float:
    """sensemakr::robustness_value. alpha=1 -> RV_q ; alpha=0.05 -> RV_{q,alpha}."""
    fq = q * abs(t_stat) / np.sqrt(dof)
    f_crit = 0.0 if alpha >= 1.0 else abs(student_t.ppf(alpha / 2, dof - 1)) / np.sqrt(dof - 1)
    fqa = fq - f_crit
    if fqa < 0:
        return 0.0
    rv = 0.5 * (np.sqrt(fqa ** 4 + 4 * fqa ** 2) - fqa ** 2)
    if f_crit > 0 and fqa > 0 and fq > 1.0 / f_crit:      # sensemakr "extreme" branch
        rv = (fq ** 2 - f_crit ** 2) / (1.0 + fq ** 2)
    return float(rv)


# --------------------------------------------- 5.3 Rosenbaum bounds (matched pairs)
def nn_match_pairs(e_hat, T, seed: int = 42, caliper_sd: float = 0.2):
    """Greedy 1:1 nearest-neighbour matching without replacement on logit(e),
    caliper = caliper_sd * SD(logit e), treated processed in seeded random order."""
    lg = np.log(e_hat / (1 - e_hat))
    cal = caliper_sd * lg.std(ddof=1)
    rng = np.random.default_rng(seed)
    treated = np.flatnonzero(T == 1)
    controls = np.flatnonzero(T == 0)
    rng.shuffle(treated)
    ctrl_lg = lg[controls]
    used = np.zeros(len(controls), dtype=bool)
    pairs = []
    for i in treated:
        d = np.abs(ctrl_lg - lg[i])
        d[used] = np.inf
        j = int(np.argmin(d))
        if d[j] <= cal:
            used[j] = True
            pairs.append((int(i), int(controls[j])))
    return pairs


def rosenbaum_p_upper(D: np.ndarray, gamma: float) -> float:
    """Worst-case (upper bound) one-sided p-value for the Wilcoxon signed-rank
    statistic T+ = sum of ranks of positive differences under Gamma (Rosenbaum
    2002 §4.3.3; rbounds::psens): E = Gamma/(1+Gamma) * sum r, V = Gamma/(1+Gamma)^2 * sum r^2."""
    D = np.asarray(D, dtype=float)
    D = D[D != 0]
    if len(D) == 0:
        return 1.0
    r = rankdata(np.abs(D))
    Tplus = float(r[D > 0].sum())
    pplus = gamma / (1.0 + gamma)
    E = pplus * r.sum()
    V = pplus * (1 - pplus) * (r ** 2).sum()
    if V <= 0:
        return 1.0
    return float(norm.sf((Tplus - E) / np.sqrt(V)))


def rosenbaum_gamma_star(D: np.ndarray, direction: int, alpha: float = 0.05,
                         gamma_max: float = 20.0):
    """Gamma* = sup{Gamma >= 1 : p_upper(Gamma) < alpha} for the effect in the
    given direction (+1 / -1).  Returns (gamma_star, p_at_gamma1, n_pairs_nonzero)."""
    Dd = np.asarray(D, dtype=float) * (1 if direction >= 0 else -1)
    p1 = rosenbaum_p_upper(Dd, 1.0)
    n_nz = int((Dd != 0).sum())
    if p1 >= alpha:
        return 1.0, p1, n_nz
    if rosenbaum_p_upper(Dd, gamma_max) < alpha:
        return gamma_max, p1, n_nz
    g = brentq(lambda gm: rosenbaum_p_upper(Dd, gm) - alpha, 1.0, gamma_max, xtol=1e-6)
    return float(g), p1, n_nz


def widening_bound_proxy(att: float, se: float, alpha: float = 0.05) -> float:
    """The public repo's 'Rosenbaum' surrogate (NOT a Rosenbaum bound): |att|/(z_alpha se), clamped at 1."""
    z = norm.ppf(1 - alpha)
    if se <= 0:
        return 1.0 if att == 0 else float("inf")
    c = abs(att) / (z * se)
    return 1.0 if c <= 1.0 else float(c)


# ------------------------------------------- 5.4 measured strength of a withheld X
def measured_confounder_strength(el_long: OvbElements, el_short: OvbElements):
    """cf_y* = (sigma2_short - sigma2_long)/sigma2_short  (partial R^2 of the withheld
    variable with Y given T, X_reduced); cf_d* = 1 - nu2_short/nu2_long
    (so C_D^2 = cf_d/(1-cf_d) = (nu2_long - nu2_short)/nu2_short); rho_hat* = realised
    change in estimate divided by the rho=1 bound at (cf_y*, cf_d*)."""
    cf_y = max(0.0, (el_short.sigma2 - el_long.sigma2) / el_short.sigma2)
    cf_d = max(0.0, 1.0 - el_short.nu2 / el_long.nu2) if el_long.nu2 > 0 else 0.0
    cf_d = min(cf_d, 1 - 1e-9)
    b = ovb_bounds(el_short, cf_y, cf_d, 1.0)["bias_bound"]
    rho_hat = (el_short.theta - el_long.theta) / b if b > 0 else float("nan")
    return float(cf_y), float(cf_d), float(np.clip(rho_hat, -1, 1)) if np.isfinite(rho_hat) else float("nan"), float(b)
