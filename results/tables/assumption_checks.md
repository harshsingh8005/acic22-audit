# Assumption checks for the pre-registered tests

Unit of replication for every test = cohort (independent ACIC 2022 datasets, generated independently by the organisers). 
No cohort contributes twice to any between-cohort comparison; within-cohort comparisons (adjustment set, estimator) use paired tests.


## E1_AIPW

- **P1 one-proportion score test vs 0.95** — n = 3200, successes = 2494; expected counts under H0: n·0.95 = 3040, n·0.05 = 160 (both ≫ 5 → normal approximation valid; Wilson CI reported). Independence: cohorts are independent datasets.
- **P2 two-proportion z (pooled) + Newcombe CI** — n_Strong = 1600 (successes 1127), n_Weak = 1600 (successes 1367); all cell counts ≥ 5 → asymptotics valid; independent groups (different cohorts).
- **P3 Welch t on d = cover_full − cover_reduced (values −1/0/+1)** — n_A = 1600, n_B = 1600; distribution of d in A: {-1: 0.009, 0: 0.731, 1: 0.259}, in B: {-1: 0.03, 0: 0.943, 1: 0.028}; SD_A = 0.454, SD_B = 0.240; Levene F = 284.62, p = 3.06e-61 (variances unequal → Welch, as pre-specified). Normality: d is discrete, so Shapiro–Wilk is not meaningful; with n = 1,600 per group the sampling distribution of the mean is normal by the CLT (bounded variable). Robustness re-computations (labelled, not replacements): permutation test of the mean difference (10,000 label permutations, seed 20260818) p = 0.0000 (observed diff = 0.2525); Mann–Whitney U p = 3.48e-78. Same conclusion as the pre-specified Welch t.

## E2_DRATT

- **P1 one-proportion score test vs 0.95** — n = 3200, successes = 2472; expected counts under H0: n·0.95 = 3040, n·0.05 = 160 (both ≫ 5 → normal approximation valid; Wilson CI reported). Independence: cohorts are independent datasets.
- **P2 two-proportion z (pooled) + Newcombe CI** — n_Strong = 1600 (successes 1156), n_Weak = 1600 (successes 1316); all cell counts ≥ 5 → asymptotics valid; independent groups (different cohorts).
- **P3 Welch t on d = cover_full − cover_reduced (values −1/0/+1)** — n_A = 1600, n_B = 1600; distribution of d in A: {-1: 0.019, 0: 0.755, 1: 0.226}, in B: {-1: 0.029, 0: 0.943, 1: 0.028}; SD_A = 0.450, SD_B = 0.239; Levene F = 237.16, p = 1.12e-51 (variances unequal → Welch, as pre-specified). Normality: d is discrete, so Shapiro–Wilk is not meaningful; with n = 1,600 per group the sampling distribution of the mean is normal by the CLT (bounded variable). Robustness re-computations (labelled, not replacements): permutation test of the mean difference (10,000 label permutations, seed 20260818) p = 0.0000 (observed diff = 0.2094); Mann–Whitney U p = 8.30e-58. Same conclusion as the pre-specified Welch t.

## E0_OLS

- **P1 one-proportion score test vs 0.95** — n = 3200, successes = 2186; expected counts under H0: n·0.95 = 3040, n·0.05 = 160 (both ≫ 5 → normal approximation valid; Wilson CI reported). Independence: cohorts are independent datasets.
- **P2 two-proportion z (pooled) + Newcombe CI** — n_Strong = 1600 (successes 937), n_Weak = 1600 (successes 1249); all cell counts ≥ 5 → asymptotics valid; independent groups (different cohorts).
- **P3 Welch t on d = cover_full − cover_reduced (values −1/0/+1)** — n_A = 1600, n_B = 1600; distribution of d in A: {-1: 0.006, 0: 0.763, 1: 0.231}, in B: {-1: 0.021, 0: 0.961, 1: 0.018}; SD_A = 0.431, SD_B = 0.198; Levene F = 285.33, p = 2.20e-61 (variances unequal → Welch, as pre-specified). Normality: d is discrete, so Shapiro–Wilk is not meaningful; with n = 1,600 per group the sampling distribution of the mean is normal by the CLT (bounded variable). Robustness re-computations (labelled, not replacements): permutation test of the mean difference (10,000 label permutations, seed 20260818) p = 0.0000 (observed diff = 0.2287); Mann–Whitney U p = 4.21e-75. Same conclusion as the pre-specified Welch t.

## S3 (OVB bound calibration) — paired binary outcomes within cohort (OVB-adjusted CI covers / naive CI covers) → exact McNemar (no distributional assumption; discordant pairs reported in hypothesis_tests.csv). Coverage proportion CI: Wilson (n = 1,600).
- Measured strengths used: cf_y* clipped at 0 in 405 of 1600 Scenario-A cohorts (E1); cf_d* clipped at 0 in 443. ν² plug-in fallback used in the short model in 6.5% of Scenario-A cohorts.

## S4 (RCT gate) — one-proportion Wilson CI (n = 200; expected counts 190/10 under 0.95 → adequate); trend test = logistic regression of cover on ordinal score (RCT=0, Weak=1, Strong=2), Wald z on the slope (n = 3,400; no separation).

## S6 (Spearman concordance) — rank-based, no distributional assumption; heavy ties at the floor (Γ* = 1 in ~66% of cohorts) handled by average ranks; bootstrap CIs (B = 500 resamples of cohorts).

## S2/S5 (E5_DragonNet, X_full) — n cohorts = 170. S2: exact McNemar on paired coverage indicators (no distributional assumption). S5: Wilcoxon signed-rank on paired PEHE differences (pre-specified non-parametric): median diff = 6.98, mean = 7.19, skew = 0.50; Shapiro–Wilk on differences W = 0.971, p = 1.30e-03 (non-normal; the rank test was pre-specified regardless and is not affected; symmetry of differences is approximate — median log-ratio with bootstrap CI reported as the effect size).

## S2/S5 (E4_TARNet, X_full) — n cohorts = 170. S2: exact McNemar on paired coverage indicators (no distributional assumption). S5: Wilcoxon signed-rank on paired PEHE differences (pre-specified non-parametric): median diff = 7.01, mean = 6.88, skew = 0.09; Shapiro–Wilk on differences W = 0.993, p = 6.34e-01 (consistent with normality; the rank test was pre-specified regardless and is not affected; symmetry of differences is approximate — median log-ratio with bootstrap CI reported as the effect size).