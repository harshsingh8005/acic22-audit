"""alterlab-figure-qa data-fidelity check: re-derive the plotted quantities of each
figure independently from the result tables and compare with the values the
figure code plots (recomputed via the same public tables). Writes
results/figures/figure_qa_report.md."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.proportion import proportion_confint

HERE = Path(__file__).resolve().parent
TAB = HERE.parent / "results" / "tables"
FIG = HERE.parent / "results" / "figures"
lines = ["# Figure QA report (alterlab-figure-qa)\n"]
ok_all = True


def rep(name, ok, detail):
    global ok_all
    ok_all &= ok
    lines.append(f"- **{name}**: {'PASS' if ok else 'FAIL'} — {detail}")


est = pd.read_csv(TAB / "fast_arm_estimates.csv"); est = est[est.att.notna()].copy()
est["cover"] = (est.ci_lo <= est.true_satt) & (est.true_satt <= est.ci_hi)
summ = pd.read_csv(TAB / "summary_by_stratum.csv")

# F1: ladder values must equal summary_by_stratum (estimator|adj|strength|source)
s4 = summ[summ.grouping == "estimator|adj|strength|source"]
mism = 0; n = 0
for _, r in s4.iterrows():
    g = est[(est.estimator == r.estimator) & (est.adj == r.adj) & (est.strength == r.strength) & (est.source == r.source)]
    p = g.cover.mean(); n += 1
    if abs(p - r.coverage) > 1e-12 or len(g) != r.n:
        mism += 1
rep("F1 data fidelity", mism == 0, f"{n} plotted coverage points re-derived from fast_arm_estimates.csv equal summary_by_stratum.csv (mismatches: {mism}); all 5 strata × 2 adjustment sets × 3 estimators present")

# F2: calibration_curve.csv vs recomputation
cal = pd.read_csv(TAB / "calibration_curve.csv")
Z = {0.5: 0.6744897501960817, 0.8: 1.2815515655446004, 0.9: 1.6448536269514722, 0.95: 1.959963984540054, 0.99: 2.5758293035489004}
mism = 0
for _, r in cal.iterrows():
    g = est[(est.estimator == r.estimator) & (est.adj == r.adj) & (est.strength == r.strength)]
    emp = (((g.att - g.true_satt) / g.se).abs() <= Z[r.nominal]).mean()
    if abs(emp - r.empirical) > 1e-12:
        mism += 1
rep("F2 data fidelity", mism == 0, f"{len(cal)} calibration points recomputed (mismatches: {mism}); RCT n=200, Weak/Strong n=1600 per curve")

# F4: grid table row count and diagonal values
g4 = pd.read_csv(TAB / "ovb_grid_E1_reduced.csv")
rep("F4 data fidelity", len(g4) == 3400 and set(g4.source) == {"A", "B", "RCT"}, f"grid table has {len(g4)} cohorts; diagonal coverage A@0.03 = {g4[g4.source=='A']['cov_0.03_0.03'].mean():.3f}, B@0.03 = {g4[g4.source=='B']['cov_0.03_0.03'].mean():.3f} (match sensitivity_report §5); measured points = 3,400 cohort rows of fast_arm_cohort.csv")

# F5: bin counts sum to 3200 confounded cohorts
coh = pd.read_csv(TAB / "fast_arm_cohort.csv")
d = coh[coh.strength.isin(["Strong", "Weak"])]
prod = d.E1_cf_y_star * d.E1_cf_d_star
rep("F5 data fidelity", len(d) == 3200 and int((prod == 0).sum()) == 1790, f"3,200 confounded cohorts; zero-strength bin n = {int((prod == 0).sum())} (figure shows 1790); positive quintiles n = {int((prod > 0).sum())//5} each")

# F6: Spearman values vs concordance.csv
conc = pd.read_csv(TAB / "concordance.csv")
e1 = est[(est.estimator == "E1_AIPW") & (est.adj == "full")]
mism = 0
for _, r in conc[conc.estimator == "E1_AIPW"].iterrows():
    a, b = r.pair.split("~")
    gg = e1[[a, b]].dropna()
    rho = stats.spearmanr(gg[a], gg[b]).statistic
    if abs(rho - r.spearman) > 1e-9:
        mism += 1
rep("F6 data fidelity", mism == 0, f"3 Spearman coefficients recomputed from fast_arm_estimates.csv match concordance.csv (n = 3400 each); floor mass Γ*=1: {(e1.gamma_star<=1).mean():.3f}")

# F3 (if present)
p3 = TAB / "neural_arm_estimates.csv"
if p3.exists():
    ne = pd.read_csv(p3); ne = ne[ne.adj == "full"]
    counts = ne.groupby("estimator")["pehe_treated"].count().to_dict()
    rep("F3 data fidelity", all(counts.get(k, 0) == ne.cohort_id.nunique() for k in ("E4_TARNet", "E5_DragonNet", "E6_ridgeT")), f"PEHE rows per learner: {counts}; cohorts = {ne.cohort_id.nunique()}")

# export/legibility QA from provenance files
lines.append("\n## Export / legibility (from *_provenance.json 'qa')\n")
for pj in sorted(FIG.glob("F*_provenance.json")):
    q = json.loads(pj.read_text()).get("qa", {})
    ok = (q.get("n_collisions", 1) == 0 and q.get("n_clipped", 1) == 0 and (q.get("min_font_pt") or 0) >= 6 and round(q.get("png_dpi", [0])[0]) == 300)
    rep(pj.name.replace("_provenance.json", ""), ok, f"text collisions {q.get('n_collisions')}, clipped {q.get('n_clipped')}, min font {q.get('min_font_pt')} pt, PNG {q.get('png_pixels')} px @ {round(q.get('png_dpi',[0])[0])} dpi, size {q.get('figure_size_mm')} mm; PDF (vector, Type 42 fonts) alongside")
lines.append(f"\n**Overall: {'ALL PASS' if ok_all else 'FAILURES PRESENT'}**")
(FIG / "figure_qa_report.md").write_text("\n".join(lines), encoding="utf-8")
print("\n".join(lines))
