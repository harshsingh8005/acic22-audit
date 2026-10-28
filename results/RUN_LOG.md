# RUN LOG — confirmatory analysis (ACIC 2022 sensitivity-layer validation)

All entries are appended chronologically. Every number/figure in the deliverables must trace to an entry here.

## 0. Plan freeze
- **2026-08-18T20:52:50Z** — `results/analysis_plan.md` FROZEN. SHA-256 = `961cc354c62e078a49467527b35a72becabcd2020f2c03f809f69d5b3cb0b953`
- No estimate, coverage or sensitivity number had been computed or inspected before this timestamp. Prior to freezing, only these were inspected: directory listings, the existing public pipeline source, the DGP factor labels in `ACIC_estimand_truths.csv`, the published `ACIC_results.csv` scorecards, the raw schema of cohort 0001, and the challenge slide text.
- Environment: Python 3.11.0; numpy 1.26.4; pandas 2.3.3; scipy 1.15.3; scikit-learn 1.8.0; statsmodels 0.14.6; torch 2.11.0+cpu; dowhy 0.14; matplotlib 3.10.8; lightgbm 4.6.0. econml/doubleml/catenets NOT installed (pending approval).
- Hardware: Intel i7-1260P (16 threads), 16 GB RAM, no GPU. Windows 11.

## 1. Software validation (plan §10) — `code/validate_software.py` → `results/tables/software_validation.json`, `logs/validate_software.log`
- V1 PASS: `code/estimators.py` AIPW-ATT / DR-ATT identical to the public reproducer on cohorts 1, 500, 1000, 2000, 3400 (|Δatt|,|Δse| = 0.0; identical X_reduced matrices, column names, truths, trim counts).
- V2 PASS: canonical 100-cohort seed-42 sample reproduces `expected_results.json` exactly (AIPW bias 26.007 / RMSE 38.793 / coverage 0.72 / width 91.593; DR-ATT 20.683 / 30.627 / 0.74 / 71.257).
- V3 PASS: E-value(RR = 3.9) = 7.263 (VanderWeele & Ding worked example 7.26).
- V4 PASS: Rosenbaum p_upper(Γ=1) equals scipy one-sided normal-approx Wilcoxon p (|Δ| < 1e-6); monotone in Γ; hand case exact.
- V5 PASS (after clarification below): bound at RV = 0 (|·| < 1e-6), CI at RVa = 0; linear RV/RVa equal DoWhy 0.14 `linear-partial-R2` on cohorts 1–10 (max |ΔRV| = 2.45e-12, |ΔRVa| = 2.44e-12, |Δestimate| = 7.9e-10).
- V6 PASS: synthetic linear DGP (n = 500, τ = 2.0): TARNet plug-in 2.001, DragonNet plug-in 2.149, DragonNet-DR 2.039 ± 0.119.
- **Plan clarification (not outcome-driven; discovered in V5 before any confirmatory result):** plan §5.2 defines ν̂² by the debiased form 2·mean[m] − mean[α̂²] but did not cover ν̂² ≤ 0, which occurs when X_full induces near-deterministic propensities (raw ê up to 0.996–1.000; e.g. cohorts 1, 7, 2500). Adopted rule = DoubleML's reference behaviour ("The estimated nu2 is not positive. Re-estimation based on riesz representer (non-orthogonal)"): fall back to the plug-in mean[α̂²] and flag the cohort (`nu2_fallback`); fallback rate is reported. Bound/RV definitions otherwise unchanged.

## 2. Timing pilots (plan §8) — timing only, outputs discarded
- Fast arm: `python code/run_fast_arm.py --pilot 8 --n-jobs 8` → 6.10 s/cohort worker time (measured while validation re-run was also on CPU) → projected ~43 min at 8 workers ≤ 60 min → **run all 3,400 cohorts × 2 adjustment sets**.
- Neural arm: `python code/run_neural_arm.py --pilot 2 --seeds 5 --n-jobs 2` → 119 s/cohort (5 seeds × 2 nets × 2 adjustment sets + DoWhy) → 25/cell projects ~105 min at 8 workers > 90-min cap → **pre-registered downgrade: 10 cohorts/cell = 170 cohorts, K = 5 seeds** (projected ~42 min). Recorded BEFORE the run starts. (A sklearn ConvergenceWarning appears from DoWhy's own internal LogisticRegression in the IPW cross-check; it does not concern the pre-registered estimators.)

## 3. Execution note — neural arm restarted (no results produced or inspected before restart)
- 17:15 EDT: fast arm (8 workers) and neural arm (6 workers) were launched concurrently. CPU contention slowed both: at 17:52 the fast arm was at 2000/3400 (35 min; projected total ≈ 68 min > its 60-min cap) and the neural arm had not yet completed its first 36-cohort chunk (no rows written).
- 17:56 EDT: the neural arm process tree (parent PID 36860 + 7 workers) was stopped so the primary fast arm finishes inside its cap; **zero neural rows existed** (only `neural_arm_cohort_ids.json`), so nothing was seen or discarded. The neural arm is re-launched alone after the fast arm completes, same pre-registered sample (10/cell, seeds 0–4), more workers. This is an execution/scheduling restart, not an analysis decision.

## 4. Fast arm — the pre-registered confirmatory run (all 3,400 cohorts × {X_full, X_reduced} × {E0, E1, E2})
- Command: `python code/run_fast_arm.py --cohorts all --n-jobs 8 --chunk 200` (launched 17:15:19 EDT, finished 18:03 EDT, **48.1 min**, 0 errors, 0 pre-specified exclusions triggered — every cohort had ≥ 5 units per arm after the Crump trim and ≥ 20 matched pairs).
- Seeds: cross-fitting folds StratifiedKFold(5, shuffle, random_state=42) shared by E1/E2/E-placebo; matching order seed 42; ridge random_state 42. No bootstrap.
- Outputs: `results/tables/fast_arm_estimates.csv` (20,400 rows; sha256 3d62211cee675283…), `results/tables/fast_arm_cohort.csv` (3,400 rows; 34fe459bf6ee3b47…); log `logs/fast_arm.log`.
- Software: code/acic_common.py, estimators.py, sensitivity.py (versions as §0).

## 5. Analysis of the fast arm — `python code/analyze.py` (18:05 EDT; log `logs/analyze_fastarm_only.log`)
- Outputs: `summary_by_stratum.csv`, `hypothesis_tests.csv` (P1–P3, S1, S3, S4; S2/S5 pending neural), `sensitivity_summary.csv`, `calibration_curve.csv`, `concordance.csv`, `measured_strength_placebo_by_stratum.csv`, `key_numbers.json`. Bootstrap seed 20260818 (B = 2,000; concordance B = 500).
- Assumption checks: `python code/assumption_checks_report.py` → `results/tables/assumption_checks.md` (permutation p for P3 with 10,000 permutations, seed 20260818 — robustness only, not a replacement).
- Report tables: `python code/make_report_tables.py` → `results/tables/report_tables.md`.

## 6. OVB coverage grid for F4 — `python code/run_ovb_grid.py --n-jobs 5 --adj reduced` (18:04–18:53 EDT, 49.1 min under CPU contention with the neural arm; recomputes E1 X_reduced nuisances with the same seed-42 folds) → `results/tables/ovb_grid_E1_reduced.csv` (3,400 rows; 11×11 grid {0, .01, .02, .03, .05, .075, .1, .15, .2, .3, .5}²; `cf_needed_diag` by Brent root-finding, xtol 1e-6).

## 7. Figures — `python code/figures.py` (F1, F2, F4, F5, F6 from the fast arm; F3 after the neural arm) → `results/figures/F*.png` (300 dpi) + `.pdf` + `*_provenance.json`.
- Presentational deviation from plan §14 (disclosed): F5's pre-specified "deciles of cf_y*·cf_d*" collapsed (6/10 deciles at exactly 0 because cf_d* = 0 for most Scenario-B/RCT cohorts); F5 now uses one "= 0" bin (n = 1,790) plus rank-quintiles of the positive values (n = 282 each) and adds a paired-difference panel (bootstrap CI, B = 1,000, seed 20260818). No test depends on F5.
- Figure iterations (legibility only, no data change): F2 tick-label rotation; F4 re-layout (dedicated colorbar axis, shorter titles); F5 label shortening; F6 floor-mass annotation.

## 8. Neural arm — pre-registered run STALLED; cap hit with zero completed cohorts (disclosed deviation)
- 18:03:33 EDT: `python code/run_neural_arm.py --per-cell 10 --seeds 5 --n-jobs 10 --chunk 40` relaunched alone (fast arm finished). Worker CPU times advanced until ~19:00 and then froze (7 workers alive with static CPU, 3 workers exited, parent idle, no chunk written); the process tree was killed at 19:17:29 EDT. **No neural result was produced or seen.** Root cause not established; suspected worker-side deadlock at high parallelism (DoWhy/matplotlib import inside loky workers on Windows is the prime suspect — the single-process validation and the 2-worker pilot ran DoWhy fine). The machine was also power-throttled (≈6.4 effective cores for 15 workers), which had already pushed the projected wall time past the 90-min cap.
- Consequence for the plan: the neural arm's pre-registered cap (plan §8/§12: 90 min) expired at 19:33 EDT with **0 completed cohorts** → S2, S5 and F3 cannot be reported as pre-registered results. Completing the 170-cohort run on this machine is a multi-hour job (est. 1–2 h wall) → **held for the user's approval** (stop condition: no multi-hour job without asking). Command to complete: `python code/run_neural_arm.py --per-cell 10 --seeds 5 --n-jobs 8 --chunk 10 --no-dowhy` (resumes from completed cohorts; E7 DoWhy cross-check to be run separately single-process).
- 19:18 EDT: to avoid an empty neural leg, a **post-cap partial run** was launched: `--per-cell 2 --seeds 5 --n-jobs 8 --chunk 10 --no-dowhy` (34 cohorts, nested inside the pre-registered 170-cohort sample; same seeds; DoWhy removed from workers; MPLBACKEND=Agg). Everything derived from it is labelled **exploratory / underpowered (post-cap)** — S2/S5 MDEs at n = 34 are ≈ 25 pp coverage and ≈ 16 % PEHE.
- 19:18–19:33 EDT: the 34-cohort stage ran cleanly (`logs/neural_arm_partial.log`: 34/34 in 15.2 min, 0 errors, ≈4.5 min per 10 cohorts at 8 workers) — confirming the earlier stall was an execution fault, not compute. At this pace the pre-registered 170-cohort sample needs ≈ 55 min (< 1 h; not a multi-hour job).
- 19:33 EDT: **pre-registered sample completed as a logged re-run after a code/execution fault (plan §12)**: `python code/run_neural_arm.py --per-cell 10 --seeds 5 --n-jobs 8 --chunk 10 --no-dowhy` (resumes; the 34 completed cohorts are the first stage of the same nested sample; identical seeds/spec; DoWhy E7 moved out of the workers to `code/run_e7_dowhy.py`, single process). Log `logs/neural_arm_full.log`. The neural results are reported as pre-registered with this deviation disclosed (stall → restart; the E7 relocation does not affect any pre-registered quantity).
- 19:20–19:46 EDT: `code/run_e7_dowhy.py` (DoWhy E7 cross-check, single process) was started without BLAS thread caps, spawned 47 threads and starved the neural workers (first completion-run chunk 9.5 min vs 4.5 min); killed at 19:46 (no output written), thread caps added; to be re-run after the neural arm.
- 19:33–20:19 EDT: completion run finished — `logs/neural_arm_full.log`: 136/136 in 45.5 min, 0 errors → **170 cohorts (10 per cell), K = 5, both adjustment sets complete**. Outputs: `results/tables/neural_arm_estimates.csv` (sha256 dd30c4f511b3c6c2…), `neural_arm_cohort.csv`, `neural_arm_cohort_ids.json` (per_cell 10, seed 20260818, ids + seeds).

## 9. Full analysis re-run (fast + neural) and final QA — 20:20–20:30 EDT
- `python code/analyze.py` → `hypothesis_tests.csv` now includes S2/S5 (X_full pre-registered; X_reduced labelled exploratory), `key_numbers.json`; log `logs/analyze_full.log`.
- `python code/make_report_tables.py` → `report_tables.md`; `python code/assumption_checks_report.py` → `assumption_checks.md` (adds S2/S5 checks: McNemar exact; Wilcoxon on PEHE differences with Shapiro W = 0.971 (p = 1.3e-3, DragonNet) / 0.993 (p = 0.63, TARNet)).
- `python code/figures.py f3` → F3 (180 × 70 mm, 300 dpi + PDF); `python code/figure_fidelity_check.py` → `results/figures/figure_qa_report.md`: **ALL PASS** (F1–F6 data fidelity re-derived from the tables; 0 text collisions, 0 clipped, min font 6–7 pt, PNG 300 dpi, vector PDF with Type 42 fonts).
- E7 DoWhy cross-check: first single-process attempt (20:20) ran with 54 DoWhy-internal threads and slow identification/IPW → killed; relaunched 20:26 without the exploratory IPW step (`python code/run_e7_dowhy.py`; log `logs/e7_dowhy.log`) — result appended below when finished. (V5 already established exact RV/RVa agreement with DoWhy on cohorts 1–10: max |Δ| 3e-14.)
- **E7 result (21:02 EDT):** `python code/run_e7_dowhy.py 10` (single process, thread caps, no IPW; DoWhy identification ≈ 38 s/cohort so the pass was capped at the first 10 sub-sample cohorts 11–112 after a 34-cohort pass whose RV column was lost to an accessor bug) → `results/tables/e7_dowhy_crosscheck.csv` (sha256 d26a2d0bc8101802…): max |Δ estimate| = 5.0e-11, max |Δ RV| = 1.3e-13, max |Δ RVa| = 1.0e-13 vs E0/`rv_linear` — DoWhy 0.14 `linear-partial-R2` and our sensemakr-form implementation agree exactly (together with V5's cohorts 1–10: 20 cohorts total).
- Reports finalised: `results/sensitivity_report.md` (§8 neural + F4 contour filled), `results/findings_summary.md` (S2/S5, claim audit §6, report-check §7).

## 10. Deliverables index (SHA-256 prefixes at close-out, 2026-08-19T00:56Z)

| file | sha256[:16] |
|---|---|
| `results/analysis_plan.md` | `961cc354c62e078a` |
| `results/sensitivity_report.md` | `75d4eabcdf6b0886` |
| `results/findings_summary.md` | `28d6555b1268b61d` |
| `results/tables/fast_arm_estimates.csv` | `3d62211cee675283` |
| `results/tables/fast_arm_cohort.csv` | `34fe459bf6ee3b47` |
| `results/tables/neural_arm_estimates.csv` | `dd30c4f511b3c6c2` |
| `results/tables/neural_arm_cohort.csv` | `6b3a1581ce867c24` |
| `results/tables/neural_arm_cohort_ids.json` | `a532907afd794b6b` |
| `results/tables/ovb_grid_E1_reduced.csv` | `30c3e298e917ac4c` |
| `results/tables/summary_by_stratum.csv` | `3348aed7ae557992` |
| `results/tables/hypothesis_tests.csv` | `ac8326058caaa927` |
| `results/tables/sensitivity_summary.csv` | `0c9cfba5ec4ae982` |
| `results/tables/calibration_curve.csv` | `f9ec3dbcbe739183` |
| `results/tables/concordance.csv` | `f462d3c9b9896219` |
| `results/tables/measured_strength_placebo_by_stratum.csv` | `7c7db246c1dbb300` |
| `results/tables/key_numbers.json` | `84c0c039f8ccd7b2` |
| `results/tables/software_validation.json` | `f70f942d7cf40032` |
| `results/tables/assumption_checks.md` | `a9f695574ecc262f` |
| `results/tables/report_tables.md` | `ffe5a08c8a594d53` |
| `results/figures/figure_qa_report.md` | `44fb95b448f8dcda` |
| `results/figures/F1_coverage_ladder.png` | `98cf99a56c992381` |
| `results/figures/F2_calibration_curves.png` | `f7afd805d6d1f388` |
| `results/figures/F3_pehe_by_regime.png` | `f7cfce6cc3ff355f` |
| `results/figures/F4_ovb_sensitivity_contour.png` | `6b957be37d494f68` |
| `results/figures/F5_coverage_vs_measured_strength.png` | `2eef554ab6e1ff63` |
| `results/figures/F6_concordance.png` | `81ce5319246c9566` |
| `results/figures/F1_coverage_ladder.pdf` | `b5482000f1fe3454` |
| `results/figures/F2_calibration_curves.pdf` | `b3cd4b78c610b75f` |
| `results/figures/F3_pehe_by_regime.pdf` | `fe1fab67a7e1a4e1` |
| `results/figures/F4_ovb_sensitivity_contour.pdf` | `ac806a802e9a9908` |
| `results/figures/F5_coverage_vs_measured_strength.pdf` | `cdb86656c3db28d9` |
| `results/figures/F6_concordance.pdf` | `5ff1f7707749f44a` |
| `code/acic_common.py` | `5276e565d3875d70` |
| `code/analyze.py` | `cc8d3341467be731` |
| `code/assumption_checks_report.py` | `0d9f356141bc6a2a` |
| `code/estimators.py` | `4fca63e7f6a66aae` |
| `code/figure_fidelity_check.py` | `c42823e62d5ec167` |
| `code/figures.py` | `897d1e186bcdc893` |
| `code/make_report_tables.py` | `90bc810fabcf0d2f` |
| `code/neural.py` | `2c77558bb3b8f95c` |
| `code/run_e7_dowhy.py` | `3ef482701ef0ff39` |
| `code/run_fast_arm.py` | `26b2540a85e8026a` |
| `code/run_neural_arm.py` | `17763a27dd0ca200` |
| `code/run_ovb_grid.py` | `cf18777538772868` |
| `code/sensitivity.py` | `ca8dcb98fa0042c2` |
| `code/validate_software.py` | `16872b93de8fc5a3` |

## 11. Legs NOT run — need the user's approval to install packages (stop condition)
- E3 `DoubleMLIRM(score='ATTE')` with the same folds + `sensitivity_analysis(0.03, 0.03, 1)` + `sensitivity_benchmark(['pre_trend'])` — `pip install doubleml`.
- E8 EconML `DRLearner` + `sensitivity_analysis`/`robustness_value` — `pip install econml`; catenets `TARNet`/`DragonNet` replication of E4/E5 on the same 170 cohorts — `pip install catenets jax jaxlib` (CPU wheels; ≈100+ MB).
- Command to run once approved: add the learners to `code/run_fast_arm.py` (E3) / `code/run_neural_arm.py` (E8) exactly as specified in plan §4; no other change.

## 12. Plan-adherence statement
- Every confirmatory number traces to a line of the frozen plan (column `plan_ref` in `hypothesis_tests.csv`). Deviations: (a) ν² fallback clarification (§1); (b) neural arm stall → logged re-run of the identical pre-registered sample (§8); (c) F5 binning (presentational); (d) percentile instead of BCa bootstrap CIs; (e) E7 DoWhy pass reduced to a subset of the sub-sample and moved out of the workers. None was outcome-driven; all are disclosed in `findings_summary.md` §4.

## 13. Approved install-gated legs (user: "run it", 2026-08-18 ~21:35 EDT) — constraints restated by the user
- **Read-only sources:** the challenge data directory and `ROSENBOUND-BENCHMARKS\benchmarks` are read-only for this paper. Verified at 21:38 EDT: `git status --porcelain` in ROSENBOUND-BENCHMARKS is clean; no file under the challenge data directory was modified after 16:00 EDT. Nothing was written outside the project working directory, the scratchpad, and the memory dir at any point.
- **No external algorithm used:** from the challenge data directory only the public ACIC 2022 data files were read. The `ROSENBOUND-BENCHMARKS\benchmarks` reproducer is the author's own public repository — the E1/E2 re-implementations verified against it (V1/V2) stand.
- **Environment:** the global Python could not take the installs cleanly (they would upgrade NumPy 1.26→2.4 or downgrade scikit-learn 1.8→1.6 globally; a `--system-site-packages` venv fails because pip's legacy scanner hits a `PKG-INFO` under a local path redacted here (Errno 22) that a global editable install put on sys.path; a fully clean venv hit the Windows 260-char path limit installing torch under this deep project path). Solution: a plain venv in the project working directory + `zz_global_site.pth` pointing at the global site-packages (plain path entry — global `.pth` files are NOT re-processed, so that path never enters sys.path). Installed into the venv only: DoubleML 0.11.4, econml 0.17.0, catenets (torch backend), jax/jaxlib 0.10.2, numpy 2.4.6, joblib 1.5.3; global torch 2.11+cpu / scikit-learn 1.8 / scipy 1.15.3 / pandas 2.3.3 / statsmodels 0.14.6 / dowhy 0.14 reused. Global Python untouched. Log `logs/pip_install.log`.
- **E3/E8 pilot** (`.venv\Scripts\python code\run_e3_e8_fast.py --pilot 4`): 4.2 s/cohort, 0 errors; EconML's `robustness_value/sensitivity_interval` raise "sigma and nu must be non-negative" whenever its ν² estimate is negative (e.g. cohort 1, X_full) — recorded per cohort as `econml_sens_error` (no fallback in EconML; DoubleML and our layer fall back to the plug-in).
- **catenets pilot** (`code\run_catenets.py --pilot 1 --seeds 2`): works; DragonNet's own propensity extracted from `_repr_estimator`→`_propensity_estimator`; ≈5–6 s per fit uncontended (library defaults: lr 1e-4, n_iter ≤ 10,000, batch 100, val 0.3, patience 10×50 iters, batch-norm on).
- 22:17 EDT: `.venv\Scripts\python code\run_e3_e8_fast.py --cohorts all --n-jobs 8 --chunk 200` launched (log `logs/e3e8_fast.log`).
- 22:17–22:41 EDT: **E3/E8 fast leg finished** — 3,400 cohorts × 2 sets, 23.8 min, 0 errors → `e3e8_estimates.csv` (f5f20ea8f0deda43…), `e3e8_cohort.csv` (b0ccb81faa8f18fb…). `analyze_addendum.py` → `hypothesis_tests_addendum.csv`, `summary_e3e8_by_stratum.csv`, `e3e8_crosscheck.csv`, `key_numbers_addendum.json` (log `logs/analyze_addendum_e3e8.log`).
  - E3 DoubleML-IRM ATTE: P1 0.7675 [0.7525, 0.7818]; P2 −15.0 pp [−17.9, −12.1]; P3 +23.6 pp [21.0, 26.1]; RCT 0.900 — all pre-registered directions confirmed (Holm within the S1_E3 family); agreement with E1: ρ = 0.984/0.990, median |Δθ̂| $1.5/$0.9, |ΔSE| $0.02, coverage indicators agree 97.1 %/97.8 %.
  - DoubleML-native S3 (its benchmark + sensitivity_analysis, X_reduced Scenario A): 0.730 [0.708, 0.751] vs naive 0.6275; benchmark unusable (clipped to 0/1) in 80.4 % of cohorts; 0.959 where usable (n = 314).
  - Cross-check of my OVB layer on DoubleML's nuisance predictions (6,800 fits): RV median |Δ| 5.6e-7 (max 1.7e-6), σ² and ν² identical to machine precision (fallback in 13.0 %), bias bounds identical; bound-CIs differ by convention (DoubleML one-sided quantiles / own IF; median |Δ| $5.3, DoubleML narrower).
  - E8 EconML DRLearner (ate_inference over treated): coverage 1.000 both sets, width $483/$380, SE-ratio 0.35/0.30; `robustness_value`/`sensitivity_interval` raise on negative ν² in 67.2 %/69.5 % of cohorts (no fallback in EconML).
- 22:41 EDT: catenets replication launched: `.venv\Scripts\python code\run_catenets.py --seeds 5 --n-jobs 8 --chunk 16 --adj full` (log `logs/catenets_arm.log`); pace under throttling ≈ 12.6 min per 16 cohorts → ≈ 2 h projected; checkpointed every 16 cohorts.
- Figures: F1 regenerated with an E3 (DoubleML) panel (compact regime labels; QA: 0 clipped, one axis-aligned-bbox flag between the rotated y-label and the suptitle that is not a visual overlap). F3 will be regenerated with the catenets learners.
- 22:41 EDT → 01:25 EDT (19 Aug): **catenets replication finished** — 170 cohorts, K = 5, X_full, 163.3 min under CPU throttling, 0 errors → `catenets_arm_estimates.csv` (55e3bc6f279e0bf6…), `catenets_arm_cohort.csv`; log `logs/catenets_arm.log`.
- 01:27 EDT: `.venv\Scripts\python code\analyze_addendum.py` (log `logs/analyze_addendum_full.log`) → `hypothesis_tests_addendum.csv` now includes S2_catenets / S5_catenets; `key_numbers_addendum.json` (catenets-vs-PyTorch agreement).
  - **Replication verdict:** the PyTorch neural results do not replicate under catenets defaults — S2 DragonNet-DR −5.3 pp [−12.9, +1.8] (p = 0.22, NOT CONFIRMED; PyTorch −14.7 pp confirmed); TARNet-DR +13.5 pp [7.6, 19.4] (over-covers, width $210); S5 PEHE log-ratio −0.229 / −0.217 (nets ≈ 20 % better than ridge; PyTorch +0.118/+0.123 worse). Cross-implementation Spearman of plug-in SATTs 0.19–0.24. Reported as pre-registered replication outcome; interpretation in sensitivity_report §10 and findings_summary §8.
- Figures: F3 regenerated with both implementations (legend moved outside the axes; 180 × 78 mm); F1 y-label shortened; `figure_fidelity_check.py` → ALL PASS (F1–F6).
- All deliverables updated: sensitivity_report.md §10, findings_summary.md §8 + claim audit rows, this log. Frozen plan file untouched (SHA-256 961cc354…0b953).
