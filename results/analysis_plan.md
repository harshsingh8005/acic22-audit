> **REDACTED READING COPY.** This is not the hash-bearing pre-registration.
> The frozen plan's SHA-256 is `961cc354c62e078a49467527b35a72becabcd2020f2c03f809f69d5b3cb0b953`
> and that value refers to the verbatim original, which is available from the author on
> request. This copy differs from it only by the removals listed in `REDACTIONS.md`:
> absolute local filesystem paths and one email address.
> No scientific content, hypothesis, criterion, estimator specification or number was
> changed. Verify this copy against the hash recorded for it in `REDACTIONS.md`.

# Pre-registered Analysis Plan — FROZEN

**Status:** FROZEN at 2026-08-18T20:50Z (UTC), before any estimate, coverage figure, or sensitivity number was computed or inspected. The SHA-256 of this file is recorded in `results/RUN_LOG.md` immediately after writing; any later edit is a *deviation* and must be disclosed under the Deviation Disclosure Gate (§13).
**Prepared by:** Claude Fable 5 (autonomous confirmatory-analysis run) on behalf of the author.
**Discipline skills applied:** `alterlab-preregistration-discipline` (Iron Law: no analysis before a frozen plan), `experimental-design` (unit/replicate/blocking/negative-control logic), `statistical-power` (MDE-based sample-size justification, §9).

---

## 0. Paper and claim under confirmation

**Paper (Phase-2-validated candidate, "GO"):** *A unified sensitivity-analysis layer — Rosenbaum Γ-bounds, VanderWeele–Ding E-values, and Chernozhukov–Cinelli omitted-variable-bias (OVB) robustness values — for PyWhy-ecosystem estimators, validated against ACIC 2022 ground truth.* (Phase-2 audit `compass_artifact_wf-b5113edf…`: the empty "column 8" — no existing sensitivity tool, R or Python, is validated on ACIC-style ground truth — is the differentiator this analysis supplies.)

**Empirical finding the manuscript needs (the "validated claim" this run confirms or refutes):**
> Nominal 95% intervals from standard doubly-robust ATT estimators (and neural CATE learners) under-cover the true SATT on ACIC 2022; under-coverage worsens as confounding strengthens; and *unobserved* confounding — created by a controlled omitted-confounder manipulation — degrades coverage further exactly where the withheld variable is a designed confounder. The unified sensitivity layer computed on those same estimates behaves as its theory dictates (bounds evaluated at the measured omitted-confounder strength recover the truth), which is the ground-truth validation the software paper claims.

### 0.1 Premise correction (recorded before any analysis)
The prompt's example claim reads "…coverage degrades as **unobserved** confounding strengthens." The ACIC 2022 organisers' design slides (Thal & Finucane, May 24 2022, `ACIC_data_challenge.pdf` p.4) state explicitly: *"Voluntary participation → confounding — But not unobserved confounding!"* All designed confounders are observed. The design knobs are:
- **Confounding strength** = % of treatment selection due to confounders: **Strong = 75 %, Weak = 50 %**.
- **Confounding source**: **Scenario A** = pre-intervention *change* in outcome **and** a complex varying function of covariates; **Scenario B** = the complex covariate function only.
- **Impact heterogeneity** (SD of practice-level CATTs: Large $42 / Small $17); **Idiosyncrasy of impacts** (80 % / 20 % idiosyncratic).
- 2×2×2×2 = 16 confounded cells × 200 datasets = 3,200, **plus 200 no-confounding (RCT) datasets** (large heterogeneity, small idiosyncrasy) — these are the 200 datasets with missing factor labels in `ACIC_estimand_truths.csv`.

Therefore: (i) the confirmatory *strength* factor is **observed** confounding strength (Strong vs Weak); (ii) genuine **unobserved** confounding is created by a **pre-registered omitted-confounder manipulation** — withholding the pre-intervention outcome trend from the adjustment set — which is a designed confounder in **Scenario A only**, making **Scenario B a built-in negative control** and the **RCT cohorts a no-confounding gate**. The manuscript must not phrase the ACIC 2022 factor as "unobserved confounding".

### 0.2 Second premise correction (existing "Rosenbaum" code)
The Rosenbaum module in the author's public benchmark repo (`ROSENBOUND-BENCHMARKS/benchmarks/acic22_track2/src/rosenbaum.py`) is a documented *widening-bound simplification* (`Γ* = |θ̂|/(z·se)`, a t-ratio surrogate), **not** a Rosenbaum randomisation-inference bound. This plan implements a genuine matched-pair Rosenbaum bound (§5.3) and reports the proxy only as a labelled comparator. The manuscript must not describe the proxy as a Rosenbaum bound.

---

## 1. Data and assets (paths resolved from `[PATH]` placeholders)

| Asset | Path | Notes |
|---|---|---|
| ACIC 2022 Track-2 (practice-level) canonical, 3,400 cohorts | `$ACIC22_DATA_ROOT/track2/{practice,practice_year}/` | verified: 3,400 + 3,400 files |
| Ground truth (SATT: Overall, subgroups X1–X5, Yearly, per treated practice) + DGP factors | `…\acic22\ACIC_estimand_truths.csv` (731,170 rows) | per-treated-practice SATTs enable PEHE-on-treated |
| Official competition scorecards (58 submissions) | `…\acic22\ACIC_results.csv` | external reference only (e.g. `aipwboot*` overall coverage 0.46, `dml_1` 0.71, `DiConf-v3*` 0.85) |
| Existing public pipeline (loader, AIPW, DR-ATT, widening-bound proxy) | `ROSENBOUND-BENCHMARKS/benchmarks/acic22_track2/src/` | reused verbatim for the data reduction and as the regression oracle for the AIPW/DR-ATT re-implementation |
| MIMIC-IV pipeline / count-isolation finding | `ROSENBOUND-BENCHMARKS/benchmarks/mimic_iv_w1_mortality/` | **not used in this paper** (belongs to the ML4H count-isolation paper); listed only because the prompt template carries it forward |
| Project (read/write scope) | `[project working directory]` (`results/`, `code/`, `logs/`) | data read from the data root above (read-only) |

Software (frozen versions): Python 3.11.0; numpy 1.26.4; pandas 2.3.3; scipy 1.15.3; scikit-learn 1.8.0; statsmodels 0.14.6; torch 2.11.0+cpu; dowhy 0.14; matplotlib 3.10.8; lightgbm 4.6.0. **Not installed:** econml, doubleml, catenets/jax — the corresponding legs (§4, E3/E8) are *pending user approval* per the stop condition "ask before any dependency install"; everything else runs now. Hardware: i7-1260P (16 threads), 16 GB RAM, no GPU.

---

## 2. Units, estimand, and data reduction (identical for every estimator)

- **Unit** = practice (≈500 per cohort); **replicate** = cohort (dataset); the cohort is the unit at which the design factors are randomised, so all inference about coverage/PEHE is at the cohort level (no pseudoreplication across practices or years).
- **Treatment** T = practice-level Z. **Outcome** Y = difference-in-differences change score = patient-weighted mean Y in post years (3–4) − patient-weighted mean Y in pre years (1–2) — exactly the existing loader's reduction (`data_loader.py`).
- **Estimand** = SATT (overall) — the challenge's scored estimand; truth taken from the `variable == "Overall"` row. Practice-level SATTs (rows with `id.practice`) are the truth for PEHE-on-treated.
- **Adjustment sets** (the pre-registered manipulation):
  - **X_reduced** = the existing loader's covariates: X1, X3, X5 (binary), one-hot X2, X4 (`drop_first`), X6–X9, `pre_Y` (pre-period weighted mean level), `pre_n_patients`, pre-period patient-mix means `pre_V1_avg … pre_V5_C_avg`.
  - **X_full** = X_reduced ∪ {**`pre_trend`** = Y(year 2) − Y(year 1) at the practice level} — the designed Scenario-A confounder ("pre-intervention change in outcome").
  - Withholding `pre_trend` (running on X_reduced) is the *omitted-confounder* condition; its true strength is measured from the data (§5.4) because the variable is in fact observed.
- **Splits**: one `StratifiedKFold(n_splits=5, shuffle=True, random_state=42)` fold assignment per cohort (stratified on T), computed once and shared by every cross-fitted estimator and by the sensitivity layer. Neural learners use a seeded 30 % early-stopping validation split (§4).

---

## 3. Design summary (experimental-design skill)

- Factors (between-cohort, fully crossed by the organisers): Strength {Strong, Weak}, Source {A, B}, Heterogeneity {Large, Small}, Idiosyncrasy {Large, Small}; plus the RCT arm.
- Our within-cohort manipulation: adjustment set {X_full, X_reduced} (paired within cohort → paired analyses).
- Negative controls: Scenario B (trend is not a designed confounder) and the RCT arm (no confounding).
- Blocking for the neural sub-sample: stratified by the 17 design cells (equal allocation) with a fixed seed.

---

## 4. Estimators (all on identical cohorts and folds; nuisance learners fixed a priori — no tuning on ACIC)

| ID | Estimator | Spec | Interval | Status |
|---|---|---|---|---|
| E0 | OLS-ATT (regression adjustment) | `statsmodels` OLS of Y on [T, X]; ATT = coefficient on T; HC1 SE | Wald 95 % | runs now |
| E1 | **AIPW-ATT (primary)** | cross-fitted (5 folds, seed 42), `logit_ridge` nuisances (StandardScaler+LogisticRegression C=0.5; Ridge α=10 for μ̂0 on controls and μ̂1 on treated), propensity clip [0.02, 0.98], Hájek control weights — a re-implementation of the repo's `aipw_att` that additionally returns ê, μ̂0, μ̂1 (verified identical, §10 V1) | efficient-influence-function SE, Wald 95 % | runs now |
| E2 | DR-ATT, Crump-trimmed | repo `dr_att` (trim α=0.10, clip [0.01, 0.99]) | IF SE, Wald 95 % | runs now |
| E3 | DoubleML-IRM ATTE | `DoubleMLIRM(score='ATTE', n_folds=5, n_rep=1, trimming_threshold=0.02, ml_g=Pipeline(scaler, Ridge(10)), ml_m=Pipeline(scaler, LogisticRegression(C=0.5)))`, `set_sample_splitting` = the same folds; `.confint()`, `.sensitivity_analysis(cf_y=0.03, cf_d=0.03, rho=1)`, `.sensitivity_params` (RV, RVa), `.sensitivity_benchmark(['pre_trend'])` on X_full | DoubleML CI | **pending install approval** (`pip install doubleml`) |
| E4 | TARNet (PyTorch) | Shalit et al. 2017 architecture, catenets `SNet1` sizes: representation 3×200 ELU, two outcome heads 2×100 ELU; inputs/outcome z-scored per cohort; loss = factual MSE + 1e-4·L2; Adam lr 1e-3, batch 100, ≤500 epochs, early stopping (patience 30 epochs, min 50) on a seeded 30 % validation split; **K = 5 seeds (0–4)**; in-sample μ̂0(x), μ̂1(x) | see "neural intervals" below | runs now (torch installed) |
| E5 | DragonNet (PyTorch) | E4 + linear-sigmoid propensity head on the representation (Shi et al. 2019 / catenets `SNet2`), loss = factual MSE + 1.0·BCE(propensity) + 1e-4·L2; no targeted regularisation; same training regime and seeds | see below | runs now |
| E6 | Ridge T-learner CATE | Ridge(α=10) per arm, in-sample μ̂0, μ̂1 → τ̂(x) — linear comparator for PEHE only | — | runs now |
| E7 | DoWhy PyWhy-native cross-check | `CausalModel` → `backdoor.linear_regression` (= E0) + refuter `add_unobserved_common_cause(simulation_method="linear-partial-R2")` → RV; run on the neural sub-sample only (software cross-check, §10 V5) | — | runs now (dowhy 0.14) |
| E8 | EconML DRLearner + `robustness_value`; catenets TARNet/DragonNet (JAX) | econml `DRLearner(cv=5, same folds; propensity LogisticRegression C=0.5; regression Ridge α=10; final StatsModelsLinearRegression)`, ATT = `.ate_inference(X[T==1])`, `.sensitivity_analysis`; catenets `TARNet`, `DragonNet` defaults, same seeds & sub-sample — a replication of E4/E5 | — | **pending install approval** (`pip install econml catenets jax jaxlib`) |

**Neural point estimates & intervals (pre-specified):**
- Plug-in SATT: τ̂ = mean_k mean_{T=1}[μ̂1_k(x) − μ̂0_k(x)] (used for bias/RMSE and, per unit, for PEHE-on-treated).
- **Neural-DR interval (the pre-registered neural interval for S2):** the E1 AIPW-ATT formula evaluated with the seed-averaged *in-sample* network nuisances — μ̂0 = mean_k μ̂0_k; ê = mean_k ê_k from DragonNet's propensity head (clipped [0.02, 0.98]) or, for TARNet, the same cross-fitted logistic ê used by E1 — with the IF SE. This is the Shi et al. (2019) "network nuisances → AIPW" pipeline; its known weakness (in-sample nuisance overfitting) is precisely the over-confidence mechanism under test.
- Exploratory only: seed-ensemble interval τ̂ ± 1.96·SD_k(τ̂_k).

---

## 5. Sensitivity layer (computed on E0, E1, E2, [E3], and the neural-DR estimates; both adjustment sets)

### 5.1 E-values (VanderWeele & Ding 2017)
Continuous outcome → standardised difference d = θ̂ / SD_pooled(Y) (SD over all practices in the cohort), RR = exp(0.91·|d|), **E_pt = RR + √(RR(RR−1))**; **E_CI** uses the CI limit closest to the null (E_CI = 1 if the CI contains 0).

### 5.2 OVB robustness values and bias bounds (Chernozhukov, Cinelli, Newey, Sharma & Syrgkanis 2022; DoubleML parametrisation)
With the short (observed-covariate) nuisances: σ̂² = mean[(Y − ĝ(T,X))²], ĝ = T·μ̂1 + (1−T)·μ̂0; ATT Riesz representer α̂ = (T − (1−T)·ê/(1−ê))/p̂; ν̂² = 2·mean[T/(p̂²(1−ê))] − mean[α̂²].
Bias bound at scenario (cf_y, cf_d, ρ): **B = ρ·√(cf_y · cf_d/(1−cf_d)) · √(σ̂²ν̂²)**; θ_lower/upper = θ̂ ∓ B; SEs of the bounds by influence functions (delta method over θ̂, σ̂², ν̂²); OVB-adjusted (1−α) CI = [θ_lower − z·se_lower, θ_upper + z·se_upper].
**RV** = value cf_y = cf_d = rv (ρ = 1) at which the bound reaches 0: rv = ½(√(f⁴+4f²) − f²), f = |θ̂|/√(σ̂²ν̂²). **RVa** = smallest rv at which the OVB-adjusted 95 % CI reaches 0 (bisection). For E0 (OLS) the linear Cinelli–Hazlett forms are used (f_q = |t|/√dof; RV_q=1 and RV_{q=1,α=.05}) and cross-checked against DoWhy (E7).
Benchmarks reported: (a) DoubleML default (0.03, 0.03, ρ=1); (b) **measured omitted-confounder strength** (§5.4).

### 5.3 Rosenbaum Γ (genuine, matched pairs — `rbounds::psens` logic)
1:1 nearest-neighbour matching without replacement on logit(ê) (E1's cross-fitted ê for the same adjustment set), caliper 0.2·SD(logit ê), greedy in seeded random order (seed 42); paired differences D = Y_treated − Y_control; Wilcoxon signed-rank statistic (average ranks; zeros dropped); for Γ ≥ 1 the worst-case one-sided p-value bound p_upper(Γ) uses E_Γ = Γ/(1+Γ)·Σr, V_Γ = Γ/(1+Γ)²·Σr² (normal approximation). Direction = sign of the E1 point estimate. **Γ\*** = sup{Γ ≥ 1 : p_upper(Γ) < 0.05} by bisection on [1, 20]; Γ\* = 1 when p_upper(1) ≥ 0.05 ("not significant even without hidden bias"). Also reported: number of pairs, matched-pair mean difference, and the repo's widening-bound proxy (labelled "t-ratio proxy").

### 5.4 Measured strength of the withheld confounder (ground-truth benchmark, unique to this design)
For each cohort, from the E1 nuisances fitted with X_full (long) and X_reduced (short):
**cf_y\*** = max(0, (σ̂²_short − σ̂²_long)/σ̂²_short) (partial R² of `pre_trend` with Y given T, X_reduced); **cf_d\*** = max(0, 1 − ν̂²_short/ν̂²_long) (so C_D² = cf_d/(1−cf_d) = (ν̂²_long − ν̂²_short)/ν̂²_short); **ρ̂\*** = (θ̂_short − θ̂_long)/B(cf_y\*, cf_d\*, ρ=1) clipped to [−1, 1]; and the truth-anchored **ρ_true** = (θ̂_short − θ_0)/B(cf_y\*, cf_d\*, ρ=1). If DoubleML is later installed, its `sensitivity_benchmark(['pre_trend'])` output is reported alongside with any definitional differences documented.

### 5.5 Negative-control outcome
The DiD change score computed on the **pre-period only** (Y_year2 − Y_year1) cannot be affected by treatment (it precedes it). Under X_reduced (which does not contain it) the E1 "effect" on this pre-period placebo outcome should be non-zero in Scenario A (selection on the trend) and ≈0 in Scenario B/RCT — a design-validity check of the manipulation, reported descriptively (mean placebo effect with 95 % CI per stratum; not a hypothesis).

---

## 6. Hypotheses (falsification criteria fixed in advance)

All tests two-sided at α = 0.05 with the predicted direction stated; a hypothesis is **confirmed** only if the Holm-adjusted p < 0.05 **and** the sign matches the prediction. Coverage = 1[θ_0 ∈ 95 % CI].

**Primary family (Holm, m = 3) — all on E1 AIPW-ATT:**
- **P1 (under-coverage):** empirical 95 % coverage over the 3,200 confounded cohorts (X_full) < 0.95. Test: one-proportion score test vs 0.95; Wilson CI. *Falsified if coverage ≥ 0.95 or CI includes 0.95.*
- **P2 (observed strength):** coverage(Strong) − coverage(Weak) < 0 (X_full; 1,600 vs 1,600). Test: two-proportion z; Newcombe hybrid-score CI. *Falsified if difference ≥ 0 or CI includes 0.*
- **P3 (induced unobserved confounding, negative-controlled):** with d_c = cover_full,c − cover_reduced,c ∈ {−1,0,1}, mean(d | Scenario A) − mean(d | Scenario B) > 0 (1,600 vs 1,600). Test: Welch t; 95 % CI. *Falsified if difference ≤ 0 or CI includes 0.*

**Secondary families (Holm within each family):**
- **S1** P1–P3 repeated for E2 DR-ATT and E0 OLS-ATT (and E3 if installed) — one family of 3 per estimator.
- **S2 (neural over-confidence; neural sub-sample; paired):** coverage of the neural-DR interval < coverage of E1 on the same cohorts (X_full), separately for DragonNet and TARNet (m = 2). Test: exact McNemar; report the paired difference with CI, mean widths, RMSE. Prediction: neural lower.
- **S3 (OVB bound calibration under induced confounding):** in Scenario-A cohorts on X_reduced (n = 1,600), the OVB-adjusted 95 % CI at (cf_y\*, cf_d\*, ρ = 1) covers θ_0 with Wilson lower bound ≥ 0.93 (criterion) and higher than the naive CI (McNemar). *Falsified if the lower bound < 0.93.*
- **S4 (RCT gate — interpretive):** E1 coverage in the 200 RCT cohorts (X_full), Wilson CI. If the upper bound < 0.95, under-coverage exists without any confounding, and P1 must be attributed mainly to variance mis-estimation/DGP structure, with P2/P3 as the confounding-specific increments. Trend RCT → Weak → Strong: Cochran–Armitage (secondary).
- **S5 (PEHE-on-treated):** TARNet and DragonNet vs ridge T-learner (paired Wilcoxon signed-rank on per-cohort PEHE, m = 2; median paired log-ratio with 2,000-bootstrap CI). No direction pre-declared.
- **S6 (cross-family concordance, descriptive):** Spearman ρ among {RV, Γ\*, E_pt} per estimator (X_full), bootstrap CIs; proportion of cohorts at the floor (RV = 0 / Γ\* = 1 / E_CI = 1) — pre-declared expectation: most cohorts are at the floor because per-cohort effects are small relative to noise (SATT SD ≈ $20 across cohorts vs DR CI widths ≈ $70–90), so the informative validation is S3/F4, not the explain-away statistics.

**Exploratory (labelled, no confirmatory p-values, replication caveat):** seed-ensemble neural intervals; `hist_gbm` nuisances on the neural sub-sample; DoWhy IPW; coverage as a function of *measured* omitted-confounder strength (deciles of cf_y\*·cf_d\*); RV vs realised |bias|; subgroup/yearly SATTs; anything not listed above.

---

## 7. Metrics (definitions)
Per cohort × estimator × adjustment set: θ̂, se, CI, coverage indicator, width, bias θ̂ − θ_0, |z| = |θ̂ − θ_0|/se, coverage at nominal 50/80/90/95/99 % (calibration curve), E_pt, E_CI, RV, RVa, OVB CI at (0.03, 0.03) and at measured strength, Γ\*, n_pairs, proxy Γ, cf_y\*, cf_d\*, ρ̂\*, ρ_true, placebo effect; neural: PEHE-on-treated = √mean_{T=1}(τ̂(x_i) − τ_i)², per seed and ensemble; runtime. Aggregates: proportions with Wilson CIs, RMSE, mean width, SE-ratio = SD(θ̂ − θ_0)/mean(se), Spearman ρ.

---

## 8. Sample sizes and replicate counts
- **Fast arm** (E0, E1, E2, [E3], full sensitivity layer): **all 3,400 cohorts × 2 adjustment sets** (6,800 estimator runs each), provided the timing pilot (5 cohorts, timing only, results discarded) projects ≤ 60 min at 8 parallel workers; otherwise a stratified random sub-sample of **100 cohorts per design cell (1,700)**, seed 20260818 (nested).
- **Neural arm** (E4, E5, E6, neural-DR, E7): stratified **25 cohorts per cell × 17 = 425** (nested; seed 20260818), K = 5 seeds; if the timing pilot (2 cohorts) projects > 90 min at 8 workers → 10 per cell (170); if still > 90 min → K = 3. Any downgrade is logged before the run starts.

## 9. Precision / power justification (statistical-power skill; α = 0.05 two-sided, power 0.80; computed with statsmodels 0.14.6 before any data)
| Test | n | MDE |
|---|---|---|
| P1 one-proportion vs 0.95 | 3,200 | 1.2 pp below 0.95 |
| P2 Strong vs Weak (base 0.70) | 1,600/arm | **5.0 pp** (7.0 pp at 800/arm) |
| P3 Welch on d (SD 0.5) | 1,600/scenario | **5.0 pp** coverage-loss difference (7.5 pp at 800) |
| S2 McNemar (30 % discordant) | 425 pairs | 8.4 pp paired coverage difference (12 pp at 170) |
| S3 one-sided vs 0.95 | 1,600 | 1.5 pp |
| S5 paired PEHE (SD log-ratio 0.3) | 425 | 4.6 % PEHE difference |
| S4 RCT Wilson half-width | 200 | ±4.2 pp at p = 0.90 |
Effect-size floor for "practically meaningful": 5 pp of coverage (≈ the gap between the ACIC top entries, 0.85, and `dml_1`, 0.71, is 14 pp — a 5 pp MDE is well inside the range that changes a practitioner's conclusion).

## 10. Software validation checks (run first; results in RUN_LOG; not hypotheses)
V1 AIPW/DR-ATT re-implementation equals the repo functions on 5 cohorts (|Δ| < 1e-9). V2 canonical 100-cohort seed-42 sample reproduces the repo's `expected_results.json` tolerance bands. V3 E-value: RR = 3.9 → 7.26 (VanderWeele & Ding worked example). V4 Rosenbaum: p_upper(Γ=1) equals the one-sided normal-approximation Wilcoxon p; monotone in Γ; hand-computed all-positive case. V5 RV: at cf = RV (ρ = 1) the bound hits 0 exactly; linear RV equals DoWhy's `linear-partial-R2` RV on ≥ 10 cohorts (|Δ| < 1e-6). V6 TARNet/DragonNet recover a synthetic linear ATT within 2 SE (smoke test).

## 11. Exclusions and failure handling
A cohort is excluded **for a given estimator only** if that estimator raises (DR-ATT trim leaving < 5 units per arm; < 20 matched pairs; NaN in neural training) — counts reported; no cohort excluded for E1 is expected (all 3,400 have ≥ 5 units per arm). No outlier removal on estimates or coverage indicators. Missing DGP labels (the 200 RCT cohorts) are a stratum, not an exclusion.

## 12. Stopping rule
One pass over the pre-specified sample; no interim looks; no re-runs except to fix a code error (logged); no hyperparameter tuning on ACIC; runtime caps as in §8. If a cap is hit, the run stops and the partial sample is reported as such (its cohort list logged) — it is not extended after inspection.

## 13. Multiplicity, effect sizes, deviations, exploratory labelling
Holm within each family (§6). All effects with 95 % CIs (Newcombe for proportion differences; Wilson for proportions; BCa bootstrap B = 2,000, seed 20260818, for medians/Spearman). Any deviation from this plan → STATE original spec, STATE deviation + reason, CLASSIFY affected result as exploratory unless provably outcome-independent, DISCLOSE in `findings_summary.md`. Contradicting results are reported, not rescued.

## 14. Pre-specified figures (300 dpi PNG + PDF; QA per `alterlab-figure-qa`)
F1 coverage ladder: RCT / Weak / Strong × Scenario A/B, X_full vs X_reduced, per estimator (Wilson CIs). F2 calibration curves (empirical vs nominal 50–99 %) per estimator by strength. F3 PEHE-on-treated by regime (heterogeneity × idiosyncrasy), neural vs ridge T-learner. F4 OVB sensitivity contour: coverage of the OVB-adjusted CI over a (cf_y, cf_d) grid (ρ = 1) for Scenario-A vs B on X_reduced, with measured-strength points overlaid. F5 coverage vs measured omitted-confounder strength (deciles). F6 cross-family concordance (RV, Γ\*, E_pt).

## 15. Interpretation rules (fixed now)
"The finding holds" ⇔ P1, P2, P3 all confirmed. Otherwise report exactly which hold and the effect sizes; if S4 shows RCT under-coverage, the headline is re-worded to "under-coverage is largely non-confounding, with a confounding-specific increment of X pp (Strong vs Weak) and Y pp (induced unobserved, Scenario A vs B)". S3 failing means the OVB implementation or its assumptions do not hold on this benchmark and the software claim must be weakened accordingly — it will be reported either way.

## 16. Run log requirement
`results/RUN_LOG.md` maps every number, table and figure to the script, seed, cohort list, package versions, wall time, and output file; every analysis attempted (including failures and pending legs) is listed in `findings_summary.md` (results-transparency gate).
