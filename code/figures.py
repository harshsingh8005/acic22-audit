"""Pre-specified figures F1–F6 (frozen plan §14). Reads results/tables/*, writes
results/figures/F*.png (300 dpi) + .pdf and a provenance JSON per figure.
Encodings: position on common scales, Wilson 95% CIs (unit of replication =
cohort), colour + marker redundancy (Okabe–Ito), no truncated bar baselines,
opaque white background, constrained layout, no bbox_inches='tight'."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.proportion import proportion_confint

HERE = Path(__file__).resolve().parent
TAB = HERE.parent / "results" / "tables"
FIG = HERE.parent / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
BOOT_SEED = 20260818

OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
      "yellow": "#F0E442", "blue": "#0072B2", "verm": "#D55E00", "purple": "#CC79A7"}
EST_LABEL = {"E0_OLS": "E0 OLS-ATT", "E1_AIPW": "E1 AIPW-ATT (primary)", "E2_DRATT": "E2 DR-ATT (Crump trim)",
             "E4_TARNet": "TARNet", "E5_DragonNet": "DragonNet", "E6_ridgeT": "ridge T-learner"}
LADDER = [("RCT", "RCT"), ("Weak", "B"), ("Weak", "A"), ("Strong", "B"), ("Strong", "A")]
LADDER_LAB = ["RCT\n(no conf.)", "Weak\nScen. B", "Weak\nScen. A", "Strong\nScen. B", "Strong\nScen. A"]

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white", "pdf.fonttype": 42})


def wilson(k, n):
    lo, hi = proportion_confint(k, n, method="wilson")
    return k / n, lo, hi


def qa_checks(fig):
    """alterlab-figure-qa: text-bbox collisions, clipping outside the figure, font floor."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    texts = []
    hidden = set()
    for ax in fig.axes:
        x0, x1 = sorted(ax.get_xlim()); y0, y1 = sorted(ax.get_ylim())
        for tk in ax.xaxis.get_major_ticks() + ax.xaxis.get_minor_ticks():
            loc = tk.get_loc()
            if not (x0 - 1e-9 <= loc <= x1 + 1e-9):
                hidden.update({id(tk.label1), id(tk.label2)})
        for tk in ax.yaxis.get_major_ticks() + ax.yaxis.get_minor_ticks():
            loc = tk.get_loc()
            if not (y0 - 1e-9 <= loc <= y1 + 1e-9):
                hidden.update({id(tk.label1), id(tk.label2)})
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text().strip() or id(t) in hidden:
            continue
        try:
            bb = t.get_window_extent(r)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        axb = t.axes.bbox if t.axes is not None else None
        if axb is not None and not bb.overlaps(axb.expanded(1.6, 1.6)):
            continue  # tick label beyond the view limits (not rendered)
        texts.append((t, bb))
    collisions, clipped = [], []
    for i in range(len(texts)):
        ti, bi = texts[i]
        if bi.x0 < -0.5 or bi.y0 < -0.5 or bi.x1 > W + 0.5 or bi.y1 > H + 0.5:
            clipped.append(ti.get_text()[:40])
        for j in range(i + 1, len(texts)):
            tj, bj = texts[j]
            if ti.axes is not None and tj.axes is not None and ti.axes is not tj.axes:
                pass
            ov = bi.overlaps(bj)
            if ov:
                # ignore tick labels of the same axis that merely touch by <1px
                inter_w = min(bi.x1, bj.x1) - max(bi.x0, bj.x0); inter_h = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
                if inter_w > 1.0 and inter_h > 1.0:
                    collisions.append((ti.get_text()[:30], tj.get_text()[:30]))
    # legend over data markers: check legend bbox vs scatter/line data points bboxes
    min_font = min(t.get_fontsize() for t, _ in texts) if texts else None
    return {"n_text": len(texts), "collisions": collisions[:20], "n_collisions": len(collisions),
            "clipped": clipped[:20], "n_clipped": len(clipped), "min_font_pt": min_font,
            "figure_size_mm": [round(fig.get_figwidth() * 25.4, 1), round(fig.get_figheight() * 25.4, 1)]}


def save(fig, name, provenance):
    png = FIG / f"{name}.png"; pdf = FIG / f"{name}.pdf"
    qa = qa_checks(fig)
    fig.savefig(png, dpi=300); fig.savefig(pdf)
    from PIL import Image
    im = Image.open(png); qa["png_pixels"] = list(im.size); qa["png_dpi"] = list(im.info.get("dpi", (None, None)))
    provenance = {**provenance, "qa": qa}
    (FIG / f"{name}_provenance.json").write_text(json.dumps(provenance, indent=2))
    plt.close(fig)
    print("wrote", png.name, pdf.name, "| QA:", {k: qa[k] for k in ("n_collisions", "n_clipped", "min_font_pt", "png_dpi", "png_pixels")})
    if qa["n_collisions"] or qa["n_clipped"]:
        print("   collisions:", qa["collisions"], "clipped:", qa["clipped"])


def load_est():
    est = pd.read_csv(TAB / "fast_arm_estimates.csv")
    est = est[est["att"].notna()].copy()
    est["cover"] = (est["ci_lo"] <= est["true_satt"]) & (est["true_satt"] <= est["ci_hi"])
    est["z"] = (est["att"] - est["true_satt"]) / est["se"]
    return est


# ------------------------------------------------------------------ F1 ladder
def f1(est):
    ests = ["E0_OLS", "E1_AIPW", "E2_DRATT"]
    p3 = TAB / "e3e8_estimates.csv"
    if p3.exists():
        e3 = pd.read_csv(p3); e3 = e3[(e3.estimator == "E3_DML_ATTE") & e3.att.notna()].copy()
        e3["cover"] = (e3.ci_lo <= e3.true_satt) & (e3.true_satt <= e3.ci_hi); e3["z"] = (e3.att - e3.true_satt) / e3.se
        est = pd.concat([est, e3[[c for c in est.columns if c in e3.columns]]], ignore_index=True)
        ests.append("E3_DML_ATTE"); EST_LABEL["E3_DML_ATTE"] = "E3 DoubleML-IRM ATTE"
    fig, axes = plt.subplots(1, len(ests), figsize=(180 / 25.4, 62 / 25.4), layout="constrained", sharey=True)
    for ax, estr in zip(axes, ests):
        for adj, col, mk, off in (("full", OI["blue"], "o", -0.08), ("reduced", OI["verm"], "s", 0.08)):
            ys, los, his = [], [], []
            for s_, src in LADDER:
                g = est[(est.estimator == estr) & (est.adj == adj) & (est.strength == s_) & (est.source == src)]
                p, lo, hi = wilson(int(g["cover"].sum()), len(g)); ys.append(p); los.append(lo); his.append(hi)
            x = np.arange(len(LADDER)) + off
            ax.errorbar(x, ys, yerr=[np.array(ys) - np.array(los), np.array(his) - np.array(ys)], fmt=mk, color=col,
                        capsize=2, ms=4, lw=1, label=f"X_{adj}" + (" (+pre_trend)" if adj == "full" else " (pre_trend withheld)"))
        ax.axhline(0.95, color="grey", ls="--", lw=0.8); ax.text(4.45, 0.952, "nominal 0.95", ha="right", va="bottom", fontsize=6, color="grey")
        ax.set_xticks(range(len(LADDER))); ax.set_xticklabels(["RCT", "Wk-B", "Wk-A", "St-B", "St-A"], fontsize=6.5); ax.set_title(EST_LABEL[estr], fontsize=8)
        ax.set_xlabel("regime", fontsize=6.5)
        ax.set_ylim(0, 1.0); ax.grid(axis="y", lw=0.3, alpha=0.5)
    axes[0].set_ylabel("coverage of true SATT (Wilson 95% CI)", fontsize=7)
    axes[0].legend(loc="lower left", frameon=False)
    fig.suptitle("F1 — Coverage ladder by regime (Wk/St = Weak/Strong, A/B = Scenario) and adjustment set; 3,400 ACIC 2022 cohorts", fontsize=8)
    save(fig, "F1_coverage_ladder", {"source": "results/tables/fast_arm_estimates.csv", "n_per_stratum": "200 (RCT) or 800 per Strength×Source cell",
                                    "uncertainty": "Wilson 95% CI on the coverage proportion", "script": "code/figures.py::f1"})


# ------------------------------------------------------------ F2 calibration
def f2(est):
    levels = [0.50, 0.80, 0.90, 0.95, 0.99]
    zs = [stats.norm.ppf(0.5 + l / 2) for l in levels]
    fig, axes = plt.subplots(1, 3, figsize=(180 / 25.4, 62 / 25.4), layout="constrained", sharey=True)
    cols = {"RCT": OI["green"], "Weak": OI["blue"], "Strong": OI["verm"]}; mks = {"RCT": "^", "Weak": "o", "Strong": "s"}
    for ax, estr in zip(axes, ("E0_OLS", "E1_AIPW", "E2_DRATT")):
        xpos = list(range(len(levels)))   # categorical spacing of the nominal levels (labels carry the values)
        ax.plot(xpos, levels, color="grey", ls="--", lw=0.8, label="perfect calibration")
        for s_ in ("RCT", "Weak", "Strong"):
            g = est[(est.estimator == estr) & (est.adj == "full") & (est.strength == s_)]
            emp = [(g["z"].abs() <= z).mean() for z in zs]
            ax.plot(xpos, emp, marker=mks[s_], color=cols[s_], ms=4, lw=1, label=f"{s_} (n={len(g)})")
        ax.set_title(EST_LABEL[estr] + ", X_full"); ax.set_xlabel("Nominal coverage"); ax.set_ylim(0.2, 1.0)
        ax.set_xticks(xpos); ax.set_xticklabels([f"{l:.2f}" for l in levels]); ax.grid(lw=0.3, alpha=0.5)
    axes[0].set_ylabel("Empirical coverage"); axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle("F2 — Calibration curves: empirical vs nominal Wald-interval coverage, by confounding strength", fontsize=9)
    save(fig, "F2_calibration_curves", {"source": "results/tables/fast_arm_estimates.csv", "definition": "empirical = mean(|z| <= z_nominal), z=(att-truth)/se",
                                       "script": "code/figures.py::f2"})


# ---------------------------------------------------------------- F3 PEHE
def f3():
    p = TAB / "neural_arm_estimates.csv"
    if not p.exists():
        print("F3 skipped: neural table missing"); return
    ne = pd.read_csv(p)
    ne = ne[ne.adj == "full"]
    pc = TAB / "catenets_arm_estimates.csv"
    ests = [("E6_ridgeT", OI["black"], "o"), ("E4_TARNet", OI["blue"], "s"), ("E5_DragonNet", OI["verm"], "^")]
    if pc.exists():
        cn = pd.read_csv(pc); cn = cn[(cn.adj == "full") & cn.estimator.str.startswith("C")]
        ne = pd.concat([ne, cn[[c for c in ne.columns if c in cn.columns]]], ignore_index=True)
        ests += [("C4_catenets_TARNet", OI["sky"], "D"), ("C5_catenets_DragonNet", OI["orange"], "v")]
        EST_LABEL["C4_catenets_TARNet"] = "TARNet (catenets)"; EST_LABEL["C5_catenets_DragonNet"] = "DragonNet (catenets)"
        EST_LABEL["E4_TARNet"] = "TARNet (PyTorch)"; EST_LABEL["E5_DragonNet"] = "DragonNet (PyTorch)"
    ne["regime"] = ne["heterogeneity"].fillna("Large") + " het. /\n" + ne["idiosyncrasy"].fillna("Small") + " idio."
    regimes = ["Large het. /\nLarge idio.", "Large het. /\nSmall idio.", "Small het. /\nLarge idio.", "Small het. /\nSmall idio."]
    fig, ax = plt.subplots(figsize=(180 / 25.4, 78 / 25.4), layout="constrained")
    w = 0.8 / len(ests)
    for k, (e_, col, mk) in enumerate(ests):
        for i, r in enumerate(regimes):
            g = ne[(ne.estimator == e_) & (ne.regime == r)]["pehe_treated"].dropna()
            if len(g) == 0:
                continue
            x = i + (k - (len(ests) - 1) / 2) * w
            ax.boxplot(g.values, positions=[x], widths=w * 0.8, showfliers=False,
                       medianprops={"color": col, "lw": 1.5}, boxprops={"color": col}, whiskerprops={"color": col}, capprops={"color": col})
            rng = np.random.default_rng(1)
            ax.scatter(x + rng.uniform(-w * 0.25, w * 0.25, len(g)), g.values, s=5, color=col, alpha=0.5, marker=mk, lw=0)
        ax.scatter([], [], color=col, marker=mk, label=f"{EST_LABEL[e_]}")
    ax.set_xticks(range(len(regimes))); ax.set_xticklabels(regimes)
    ax.set_ylabel("PEHE on treated practices ($)"); ax.set_ylim(bottom=0); ax.grid(axis="y", lw=0.3, alpha=0.5)
    fig.legend(frameon=False, loc="outside lower center", fontsize=6.5, ncol=len(ests))
    n_c = ne["cohort_id"].nunique()
    ax.set_title(f"F3 — PEHE-on-treated by heterogeneity regime, X_full ({n_c} cohorts; boxes = median/IQR, points = cohorts)", fontsize=8)
    save(fig, "F3_pehe_by_regime", {"source": "results/tables/neural_arm_estimates.csv", "n_cohorts": int(n_c),
                                   "definition": "PEHE = sqrt(mean over treated practices of (tau_hat - true practice SATT)^2); neural = 5-seed ensemble",
                                   "script": "code/figures.py::f3"})


# ---------------------------------------------------------- F4 OVB contour
def f4():
    p = TAB / "ovb_grid_E1_reduced.csv"
    if not p.exists():
        print("F4 skipped: ovb grid missing"); return
    g = pd.read_csv(p)
    coh = pd.read_csv(TAB / "fast_arm_cohort.csv")
    grid = sorted({float(c.split("_")[1]) for c in g.columns if c.startswith("cov_")})
    fig = plt.figure(figsize=(180 / 25.4, 72 / 25.4), layout="constrained")
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.06, 1.1])
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1], sharey=axA); cax = fig.add_subplot(gs[0, 2]); axC = fig.add_subplot(gs[0, 3])
    show = {0.0, 0.02, 0.05, 0.1, 0.2, 0.5}
    lab = [f"{v:g}" if v in show else "" for v in grid]
    for ax, src, title in ((axA, "A", "A: pre_trend IS a confounder"), (axB, "B", "B: negative control")):
        sub = g[g.source == src]
        M = np.array([[sub[f"cov_{cy}_{cd}"].mean() for cd in grid] for cy in grid])
        im = ax.imshow(M, origin="lower", vmin=0.5, vmax=1.0, cmap="viridis", aspect="auto")
        cs = ax.contour(M, levels=[0.95], colors="white", linewidths=1)
        ax.clabel(cs, fmt="0.95", fontsize=6)
        ax.set_xticks(range(len(grid))); ax.set_xticklabels(lab); ax.set_yticks(range(len(grid))); ax.set_yticklabels(lab)
        ax.set_xlabel("posited cf_d (treatment share)"); ax.set_title(title, fontsize=8)
        cs_ = coh[coh.source == src]
        xi = np.interp(cs_["E1_cf_d_star"], grid, range(len(grid))); yi = np.interp(cs_["E1_cf_y_star"], grid, range(len(grid)))
        ax.scatter(xi, yi, s=3, color=OI["orange"], alpha=0.35, lw=0, label="measured (cf_y*, cf_d*)")
        ax.legend(frameon=False, fontsize=6, loc="upper right")
    axA.set_ylabel("posited cf_y (outcome share)")
    plt.setp(axB.get_yticklabels(), visible=False)
    fig.colorbar(im, cax=cax, label="coverage of OVB-adjusted 95% CI" + chr(10) + "(E1 AIPW-ATT, X_reduced)")
    ax = axC
    for src, col in (("A", OI["verm"]), ("B", OI["blue"])):
        v = g[g.source == src]["cf_needed_diag"]
        v_plot = v.fillna(1.0).values
        xs = np.sort(v_plot); ax.step(xs, np.arange(1, len(xs) + 1) / len(xs), color=col, lw=1.2, label=f"Scenario {src} (n={len(v)})")
        cs_ = coh[coh.source == src]
        ax.axvline(cs_["E1_cf_y_star"].median(), color=col, ls=":", lw=1)
    ax.set_xscale("symlog", linthresh=0.01); ax.set_xlim(0, 1)
    ax.set_xlabel("smallest diagonal cf covering the truth", fontsize=7)
    ax.set_ylabel("cumulative fraction of cohorts"); ax.set_title("cf needed to cover the truth" + chr(10) + "(dotted: median measured cf_y*)", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="lower right"); ax.grid(lw=0.3, alpha=0.5)
    fig.suptitle("F4 — OVB sensitivity contour: coverage of the OVB-adjusted CI over posited confounding strength (rho = 1)", fontsize=9)
    save(fig, "F4_ovb_sensitivity_contour", {"source": ["results/tables/ovb_grid_E1_reduced.csv", "results/tables/fast_arm_cohort.csv"],
                                            "grid": grid, "script": "code/figures.py::f4",
                                            "note": "left/middle: mean over cohorts of 1[truth in OVB-adjusted CI at (cf_y, cf_d)], white contour = 0.95; right: ECDF of the diagonal cf needed"})


# ------------------------------------------- F5 coverage vs measured strength
def f5(est):
    coh = pd.read_csv(TAB / "fast_arm_cohort.csv")
    e1 = est[(est.estimator == "E1_AIPW")][["cohort_id", "adj", "cover"]].pivot(index="cohort_id", columns="adj", values="cover")
    d = coh.merge(e1, left_on="cohort_id", right_index=True)
    d = d[d.strength.isin(["Strong", "Weak"])].copy()
    d["prod"] = d["E1_cf_y_star"] * d["E1_cf_d_star"]
    # zero-inflated product: one bin for measured strength == 0, quintiles of the positive values
    pos = d[d["prod"] > 0].copy()
    pos["bin"] = pd.qcut(pos["prod"].rank(method="first"), 5, labels=False) + 1
    d["bin"] = 0
    d.loc[pos.index, "bin"] = pos["bin"]
    rng = np.random.default_rng(BOOT_SEED)
    fig, axes = plt.subplots(1, 2, figsize=(180 / 25.4, 65 / 25.4), layout="constrained")
    ax = axes[0]
    xs_all, labels = [], []
    for adj, col, mk, lab in (("reduced", OI["verm"], "s", "X_reduced (pre_trend withheld)"), ("full", OI["blue"], "o", "X_full (pre_trend adjusted)")):
        xs, ys, los, his = [], [], [], []
        for k, gb in d.groupby("bin"):
            p, lo, hi = wilson(int(gb[adj].sum()), len(gb)); xs.append(k); ys.append(p); los.append(lo); his.append(hi)
        ax.errorbar(xs, ys, yerr=[np.array(ys) - np.array(los), np.array(his) - np.array(ys)], fmt=mk + "-", color=col, ms=4, lw=1, capsize=2, label=lab)
    ax.axhline(0.95, color="grey", ls="--", lw=0.8)
    bins = sorted(d["bin"].unique())
    tick = []
    for k in bins:
        gb = d[d.bin == k]
        NL = chr(10)
        tick.append(("= 0" + NL + "(n=%d)") % len(gb) if k == 0 else ("%.1e" + NL + "(n=%d)") % (gb["prod"].median(), len(gb)))
    ax.set_xticks(bins); ax.set_xticklabels(tick, fontsize=6)
    ax.set_xlabel("cf_y* × cf_d* of withheld pre_trend (bin median; n)", fontsize=7)
    ax.set_ylabel("E1 95% CI coverage (Wilson CI)"); ax.set_ylim(0, 1); ax.grid(lw=0.3, alpha=0.5); ax.legend(frameon=False, fontsize=6, loc="lower left")
    ax.set_title("coverage by measured strength", fontsize=8)
    ax = axes[1]
    ys, los, his = [], [], []
    for k, gb in d.groupby("bin"):
        diff = gb["full"].astype(float) - gb["reduced"].astype(float)
        ys.append(diff.mean())
        bs = [diff.values[rng.integers(0, len(diff), len(diff))].mean() for _ in range(1000)]
        los.append(np.percentile(bs, 2.5)); his.append(np.percentile(bs, 97.5))
    ax.errorbar(bins, ys, yerr=[np.array(ys) - np.array(los), np.array(his) - np.array(ys)], fmt="D-", color=OI["green"], ms=4, lw=1, capsize=2)
    ax.axhline(0, color="grey", ls="--", lw=0.8)
    ax.set_xticks(bins); ax.set_xticklabels(tick, fontsize=6)
    ax.set_xlabel("same bins as left (bootstrap 95% CI)", fontsize=7); ax.set_ylabel("coverage gain, X_full − X_reduced", fontsize=7)
    ax.set_title("benefit of adjusting vs measured strength", fontsize=8); ax.grid(lw=0.3, alpha=0.5)
    fig.suptitle("F5 — Coverage vs measured strength of the withheld confounder (E1 AIPW-ATT)", fontsize=9)
    save(fig, "F5_coverage_vs_measured_strength", {"source": ["results/tables/fast_arm_cohort.csv", "results/tables/fast_arm_estimates.csv"],
                                                  "binning": "cf_y*·cf_d* == 0 as one bin (zero-inflated), quintiles of positive values (rank-based); x = bin median",
                                                  "note": "plan §14 specified deciles of the product; 6/10 deciles were degenerate at 0, so binning was changed to 0 + quintiles(>0) (presentational; disclosed in RUN_LOG)",
                                                  "uncertainty": "Wilson CI (left); cohort-bootstrap 95% CI, B=1000, seed 20260818 (right)", "script": "code/figures.py::f5"})


# ------------------------------------------------------------ F6 concordance
def f6(est):
    g = est[(est.estimator == "E1_AIPW") & (est.adj == "full")]
    fig, axes = plt.subplots(1, 3, figsize=(180 / 25.4, 60 / 25.4), layout="constrained")
    pairs = [("rv", "gamma_star", "RV (OVB robustness value)", "Rosenbaum Γ*"), ("rv", "e_pt", "RV", "E-value (point)"), ("gamma_star", "e_pt", "Rosenbaum Γ*", "E-value (point)")]
    for ax, (a, b, la, lb) in zip(axes, pairs):
        gg = g[[a, b]].dropna()
        rho = stats.spearmanr(gg[a], gg[b]).statistic
        ax.scatter(gg[a], gg[b], s=4, alpha=0.4, color=OI["blue"], lw=0)
        note = ""
        if "gamma_star" in (a, b):
            note = f"; {100*(gg['gamma_star']<=1).mean():.0f}% at Γ*=1"
        ax.set_xlabel(la); ax.set_ylabel(lb); ax.set_title(f"Spearman ρ = {rho:.2f} (n={len(gg)}{note})", fontsize=7); ax.grid(lw=0.3, alpha=0.5)
    fig.suptitle("F6 — Cross-family concordance of sensitivity diagnostics (E1 AIPW-ATT, X_full)", fontsize=9)
    save(fig, "F6_concordance", {"source": "results/tables/fast_arm_estimates.csv", "script": "code/figures.py::f6"})


if __name__ == "__main__":
    est = load_est()
    which = sys.argv[1:] or ["f1", "f2", "f3", "f4", "f5", "f6"]
    for w in which:
        {"f1": lambda: f1(est), "f2": lambda: f2(est), "f3": f3, "f4": f4, "f5": lambda: f5(est), "f6": lambda: f6(est)}[w]()
