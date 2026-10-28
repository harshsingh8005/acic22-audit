### T1 — Overall performance by estimator and adjustment set (all 3,400 cohorts)

| estimator | adj | n | coverage | cov_lo | cov_hi | bias_mean | rmse | width_mean | se_ratio | nu2_fallback_rate | n_trimmed_mean |
|---|---|---|---|---|---|---|---|---|---|---|---|
| E0_OLS | full | 3400 | 0.693 | 0.677 | 0.708 | 16.903 | 25.794 | 52.551 | 1.454 |  | 0.000 |
| E0_OLS | reduced | 3400 | 0.588 | 0.572 | 0.605 | 21.929 | 29.509 | 51.931 | 1.491 |  | 0.000 |
| E1_AIPW | full | 3400 | 0.786 | 0.772 | 0.800 | 30.183 | 49.777 | 123.085 | 1.261 | 0.181 | 0.000 |
| E1_AIPW | reduced | 3400 | 0.671 | 0.655 | 0.686 | 27.233 | 39.790 | 83.470 | 1.363 | 0.079 | 0.000 |
| E2_DRATT | full | 3400 | 0.780 | 0.766 | 0.794 | 19.100 | 30.977 | 75.953 | 1.259 | 0.000 | 126.779 |
| E2_DRATT | reduced | 3400 | 0.684 | 0.668 | 0.699 | 24.124 | 34.050 | 71.706 | 1.314 | 0.000 | 109.708 |

### T2 — Coverage ladder by estimator × adjustment × strength × source

| estimator | adj | strength | source | n | coverage | cov_lo | cov_hi | bias_mean | rmse | width_mean | se_ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| E0_OLS | full | RCT | RCT | 200 | 0.850 | 0.794 | 0.893 | 0.839 | 14.549 | 43.153 | 1.323 |
| E0_OLS | full | Strong | A | 800 | 0.613 | 0.578 | 0.646 | 21.901 | 29.597 | 56.761 | 1.376 |
| E0_OLS | full | Strong | B | 800 | 0.559 | 0.524 | 0.593 | 26.008 | 33.680 | 59.129 | 1.420 |
| E0_OLS | full | Weak | A | 800 | 0.853 | 0.826 | 0.875 | 8.649 | 16.553 | 47.280 | 1.171 |
| E0_OLS | full | Weak | B | 800 | 0.709 | 0.676 | 0.739 | 15.069 | 22.147 | 49.383 | 1.289 |
| E0_OLS | reduced | RCT | RCT | 200 | 0.850 | 0.794 | 0.893 | 0.717 | 14.656 | 43.658 | 1.318 |
| E0_OLS | reduced | Strong | A | 800 | 0.362 | 0.330 | 0.396 | 33.496 | 38.238 | 54.125 | 1.337 |
| E0_OLS | reduced | Strong | B | 800 | 0.571 | 0.537 | 0.605 | 26.246 | 34.098 | 59.791 | 1.428 |
| E0_OLS | reduced | Weak | A | 800 | 0.651 | 0.618 | 0.683 | 17.989 | 22.701 | 45.953 | 1.182 |
| E0_OLS | reduced | Weak | B | 800 | 0.703 | 0.670 | 0.733 | 15.288 | 22.513 | 49.922 | 1.299 |
| E1_AIPW | full | RCT | RCT | 200 | 0.900 | 0.851 | 0.934 | 1.702 | 15.568 | 52.687 | 1.154 |
| E1_AIPW | full | Strong | A | 800 | 0.816 | 0.788 | 0.842 | 47.919 | 71.244 | 186.708 | 1.108 |
| E1_AIPW | full | Strong | B | 800 | 0.593 | 0.558 | 0.626 | 36.128 | 51.039 | 106.000 | 1.334 |
| E1_AIPW | full | Weak | A | 800 | 0.935 | 0.916 | 0.950 | 23.042 | 40.890 | 133.026 | 0.996 |
| E1_AIPW | full | Weak | B | 800 | 0.774 | 0.743 | 0.801 | 20.762 | 33.420 | 84.203 | 1.220 |
| E1_AIPW | reduced | RCT | RCT | 200 | 0.910 | 0.862 | 0.942 | 1.009 | 15.449 | 52.088 | 1.163 |
| E1_AIPW | reduced | Strong | A | 800 | 0.502 | 0.468 | 0.537 | 41.327 | 50.900 | 96.373 | 1.209 |
| E1_AIPW | reduced | Strong | B | 800 | 0.603 | 0.568 | 0.636 | 35.081 | 49.012 | 101.927 | 1.317 |
| E1_AIPW | reduced | Weak | A | 800 | 0.749 | 0.718 | 0.778 | 19.290 | 25.186 | 62.119 | 1.022 |
| E1_AIPW | reduced | Weak | B | 800 | 0.769 | 0.738 | 0.797 | 19.791 | 32.275 | 81.308 | 1.230 |
| E2_DRATT | full | RCT | RCT | 200 | 0.900 | 0.851 | 0.934 | 1.739 | 15.579 | 52.480 | 1.159 |
| E2_DRATT | full | Strong | A | 800 | 0.756 | 0.725 | 0.785 | 23.288 | 34.942 | 85.530 | 1.195 |
| E2_DRATT | full | Strong | B | 800 | 0.689 | 0.656 | 0.720 | 27.444 | 39.397 | 84.745 | 1.308 |
| E2_DRATT | full | Weak | A | 800 | 0.897 | 0.875 | 0.917 | 11.075 | 21.446 | 70.894 | 1.016 |
| E2_DRATT | full | Weak | B | 800 | 0.748 | 0.716 | 0.776 | 18.932 | 28.008 | 68.509 | 1.182 |
| E2_DRATT | reduced | RCT | RCT | 200 | 0.910 | 0.862 | 0.942 | 0.943 | 15.478 | 51.996 | 1.168 |
| E2_DRATT | reduced | Strong | A | 800 | 0.504 | 0.469 | 0.538 | 36.903 | 43.796 | 78.026 | 1.186 |
| E2_DRATT | reduced | Strong | B | 800 | 0.689 | 0.656 | 0.720 | 27.217 | 39.394 | 84.817 | 1.317 |
| E2_DRATT | reduced | Weak | A | 800 | 0.735 | 0.703 | 0.764 | 19.536 | 24.954 | 60.681 | 1.004 |
| E2_DRATT | reduced | Weak | B | 800 | 0.751 | 0.720 | 0.780 | 18.634 | 27.836 | 68.227 | 1.189 |

### T3 — Pre-registered hypothesis tests (Holm within family)

| family | id | test | n | effect | ci_lo | ci_hi | statistic | p_raw | p_holm | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| PRIMARY | P1 | E1_AIPW X_full 95% coverage < 0.95 (score test vs 0.95; Wilson CI) | 3200 | 0.779 | 0.765 | 0.793 | -44.286 | 0.00e+00 | 0.00e+00 | CONFIRMED |
| PRIMARY | P2 | E1_AIPW X_full coverage(Strong) − coverage(Weak) < 0 (two-proportion z; Newcombe CI) | 3200 | -0.150 | -0.178 | -0.122 | -10.231 | 1.43e-24 | 1.43e-24 | CONFIRMED |
| PRIMARY | P3 | E1_AIPW mean(cover_full−cover_reduced | A) − mean(· | B) > 0 (Welch t) | 3400 | 0.253 | 0.227 | 0.278 | 19.660 | 5.70e-80 | 1.14e-79 | CONFIRMED |
| S1_E2 | P1 | E2_DRATT X_full 95% coverage < 0.95 (score test vs 0.95; Wilson CI) | 3200 | 0.772 | 0.758 | 0.787 | -46.071 | 0.00e+00 | 0.00e+00 | CONFIRMED |
| S1_E2 | P2 | E2_DRATT X_full coverage(Strong) − coverage(Weak) < 0 (two-proportion z; Newcombe CI) | 3200 | -0.100 | -0.129 | -0.071 | -6.747 | 1.51e-11 | 1.51e-11 | CONFIRMED |
| S1_E2 | P3 | E2_DRATT mean(cover_full−cover_reduced | A) − mean(· | B) > 0 (Welch t) | 3400 | 0.209 | 0.184 | 0.234 | 16.457 | 8.79e-58 | 1.76e-57 | CONFIRMED |
| S1_E0 | P1 | E0_OLS X_full 95% coverage < 0.95 (score test vs 0.95; Wilson CI) | 3200 | 0.683 | 0.667 | 0.699 | -69.269 | 0.00e+00 | 0.00e+00 | CONFIRMED |
| S1_E0 | P2 | E0_OLS X_full coverage(Strong) − coverage(Weak) < 0 (two-proportion z; Newcombe CI) | 3200 | -0.195 | -0.226 | -0.163 | -11.855 | 2.04e-32 | 2.04e-32 | CONFIRMED |
| S1_E0 | P3 | E0_OLS mean(cover_full−cover_reduced | A) − mean(· | B) > 0 (Welch t) | 3400 | 0.229 | 0.205 | 0.252 | 19.270 | 1.07e-76 | 2.13e-76 | CONFIRMED |
| S3 | S3_E1 | E1_AIPW X_reduced Scenario A: OVB-adjusted CI at measured (cf_y*,cf_d*,rho=1) covers truth; Wilson lower ≥ 0.93; McNemar vs naive CI | 1600 | 0.807 | 0.787 | 0.825 | 0.000 | 1.01e-87 | 2.01e-87 | NOT MET |
| S3 | S3_E2 | E2_DRATT X_reduced Scenario A: OVB-adjusted CI at measured (cf_y*,cf_d*,rho=1) covers truth; Wilson lower ≥ 0.93; McNemar vs naive CI | 1600 | 0.778 | 0.757 | 0.798 | 0.000 | 6.91e-77 | 6.91e-77 | NOT MET |
| S4 | S4_E1_AIPW | E1_AIPW X_full RCT coverage (Wilson CI; upper<0.95 ⇒ non-confounding under-coverage); trend z = logistic slope RCT→Weak→Strong | 200 | 0.900 | 0.851 | 0.934 | -10.560 | 4.57e-26 | 9.14e-26 | RCT under-covers (non-confounding under-coverage present) |
| S4 | S4_E2_DRATT | E2_DRATT X_full RCT coverage (Wilson CI; upper<0.95 ⇒ non-confounding under-coverage); trend z = logistic slope RCT→Weak→Strong | 200 | 0.900 | 0.851 | 0.934 | -7.929 | 2.21e-15 | 2.21e-15 | RCT under-covers (non-confounding under-coverage present) |
| S4 | S4_E0_OLS | E0_OLS X_full RCT coverage (Wilson CI; upper<0.95 ⇒ non-confounding under-coverage); trend z = logistic slope RCT→Weak→Strong | 200 | 0.850 | 0.794 | 0.893 | -12.364 | 4.12e-35 | 1.23e-34 | RCT under-covers (non-confounding under-coverage present) |
| S2 | S2_E5_DragonNet_full | E5_DragonNet neural-DR interval coverage < E1 AIPW coverage, same cohorts, X_full (exact McNemar; paired bootstrap CI) | 170 | -0.147 | -0.206 | -0.088 | 2.000 | 1.62e-06 | 3.25e-06 | CONFIRMED |
| S5 | S5_E5_DragonNet_full | E5_DragonNet PEHE-on-treated vs ridge T-learner, X_full (paired Wilcoxon; median log-ratio + bootstrap CI); no direction pre-declared | 170 | 0.123 | 0.097 | 0.146 | 1928.000 | 9.72e-17 | 9.72e-17 | reported (no direction pre-declared) |
| S2 | S2_E4_TARNet_full | E4_TARNet neural-DR interval coverage < E1 AIPW coverage, same cohorts, X_full (exact McNemar; paired bootstrap CI) | 170 | -0.035 | -0.076 | 0.000 | 3.000 | 0.146 | 0.146 | NOT CONFIRMED |
| S5 | S5_E4_TARNet_full | E4_TARNet PEHE-on-treated vs ridge T-learner, X_full (paired Wilcoxon; median log-ratio + bootstrap CI); no direction pre-declared | 170 | 0.118 | 0.089 | 0.145 | 1444.000 | 1.29e-19 | 2.58e-19 | reported (no direction pre-declared) |
| S2_reduced(exploratory) | S2_E5_DragonNet_reduced | E5_DragonNet neural-DR interval coverage < E1 AIPW coverage, same cohorts, X_reduced (exact McNemar; paired bootstrap CI) | 170 | -0.141 | -0.200 | -0.088 | 2.000 | 3.03e-06 | 6.06e-06 | CONFIRMED |
| S5_reduced(exploratory) | S5_E5_DragonNet_reduced | E5_DragonNet PEHE-on-treated vs ridge T-learner, X_reduced (paired Wilcoxon; median log-ratio + bootstrap CI); no direction pre-declared | 170 | 0.096 | 0.065 | 0.111 | 3158.000 | 1.61e-10 | 1.61e-10 | reported (no direction pre-declared) |
| S2_reduced(exploratory) | S2_E4_TARNet_reduced | E4_TARNet neural-DR interval coverage < E1 AIPW coverage, same cohorts, X_reduced (exact McNemar; paired bootstrap CI) | 170 | -0.088 | -0.141 | -0.035 | 4.000 | 0.003 | 0.003 | CONFIRMED |
| S5_reduced(exploratory) | S5_E4_TARNet_reduced | E4_TARNet PEHE-on-treated vs ridge T-learner, X_reduced (paired Wilcoxon; median log-ratio + bootstrap CI); no direction pre-declared | 170 | 0.086 | 0.066 | 0.127 | 2385.000 | 3.03e-14 | 6.06e-14 | reported (no direction pre-declared) |

### T4 — Sensitivity diagnostics by stratum (medians; fractions above floor)

| estimator | adj | strength | source | n | median_rv | median_rva | median_gamma_star | frac_gamma_gt1 | median_e_pt | median_e_ci | frac_e_ci_gt1 | naive_coverage | ovb03_coverage | nu2_fallback_rate | median_n_pairs | median_proxy_gamma |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E0_OLS | full | RCT | RCT | 200 | 0.061 | 0.000 |  | 0.000 | 1.414 | 1.000 | 0.360 | 0.850 |  |  |  | 1.000 |
| E0_OLS | full | Strong | A | 800 | 0.079 | 0.000 |  | 0.000 | 1.609 | 1.000 | 0.463 | 0.613 |  |  |  | 1.065 |
| E0_OLS | full | Strong | B | 800 | 0.080 | 0.000 |  | 0.000 | 1.648 | 1.000 | 0.470 | 0.559 |  |  |  | 1.125 |
| E0_OLS | full | Weak | A | 800 | 0.066 | 0.000 |  | 0.000 | 1.473 | 1.000 | 0.341 | 0.853 |  |  |  | 1.000 |
| E0_OLS | full | Weak | B | 800 | 0.074 | 0.000 |  | 0.000 | 1.531 | 1.000 | 0.421 | 0.709 |  |  |  | 1.017 |
| E0_OLS | reduced | RCT | RCT | 200 | 0.063 | 0.000 |  | 0.000 | 1.428 | 1.000 | 0.350 | 0.850 |  |  |  | 1.000 |
| E0_OLS | reduced | Strong | A | 800 | 0.112 | 0.029 |  | 0.000 | 1.774 | 1.276 | 0.598 | 0.362 |  |  |  | 1.543 |
| E0_OLS | reduced | Strong | B | 800 | 0.080 | 0.000 |  | 0.000 | 1.650 | 1.000 | 0.460 | 0.571 |  |  |  | 1.117 |
| E0_OLS | reduced | Weak | A | 800 | 0.075 | 0.000 |  | 0.000 | 1.507 | 1.000 | 0.449 | 0.651 |  |  |  | 1.055 |
| E0_OLS | reduced | Weak | B | 800 | 0.075 | 0.000 |  | 0.000 | 1.537 | 1.000 | 0.414 | 0.703 |  |  |  | 1.016 |
| E1_AIPW | full | RCT | RCT | 200 | 0.060 | 0.000 | 1.000 | 0.380 | 1.419 | 1.000 | 0.320 | 0.900 | 0.970 | 0.000 | 196.000 | 1.000 |
| E1_AIPW | full | Strong | A | 800 | 0.063 | 0.000 | 1.000 | 0.335 | 2.042 | 1.000 | 0.237 | 0.816 | 0.956 | 0.246 | 98.000 | 1.000 |
| E1_AIPW | full | Strong | B | 800 | 0.040 | 0.000 | 1.000 | 0.345 | 1.844 | 1.000 | 0.420 | 0.593 | 0.870 | 0.129 | 79.000 | 1.024 |
| E1_AIPW | full | Weak | A | 800 | 0.067 | 0.000 | 1.000 | 0.345 | 1.606 | 1.000 | 0.165 | 0.935 | 0.991 | 0.292 | 162.000 | 1.000 |
| E1_AIPW | full | Weak | B | 800 | 0.060 | 0.000 | 1.000 | 0.335 | 1.622 | 1.000 | 0.320 | 0.774 | 0.907 | 0.100 | 136.000 | 1.000 |
| E1_AIPW | reduced | RCT | RCT | 200 | 0.063 | 0.000 | 1.000 | 0.375 | 1.432 | 1.000 | 0.305 | 0.910 | 0.980 | 0.000 | 196.000 | 1.000 |
| E1_AIPW | reduced | Strong | A | 800 | 0.068 | 0.000 | 1.021 | 0.529 | 1.912 | 1.000 | 0.491 | 0.502 | 0.741 | 0.111 | 113.000 | 1.169 |
| E1_AIPW | reduced | Strong | B | 800 | 0.040 | 0.000 | 1.000 | 0.345 | 1.822 | 1.000 | 0.430 | 0.603 | 0.884 | 0.116 | 79.500 | 1.017 |
| E1_AIPW | reduced | Weak | A | 800 | 0.077 | 0.000 | 1.000 | 0.470 | 1.535 | 1.000 | 0.350 | 0.749 | 0.881 | 0.019 | 178.000 | 1.000 |
| E1_AIPW | reduced | Weak | B | 800 | 0.057 | 0.000 | 1.000 | 0.339 | 1.615 | 1.000 | 0.333 | 0.769 | 0.907 | 0.087 | 136.000 | 1.000 |
| E2_DRATT | full | RCT | RCT | 200 | 0.061 | 0.000 | 1.000 | 0.375 | 1.422 | 1.000 | 0.330 | 0.900 | 0.970 | 0.000 | 196.000 | 1.000 |
| E2_DRATT | full | Strong | A | 800 | 0.083 | 0.000 | 1.000 | 0.346 | 1.667 | 1.000 | 0.324 | 0.756 | 0.870 | 0.000 | 85.000 | 1.000 |
| E2_DRATT | full | Strong | B | 800 | 0.083 | 0.000 | 1.000 | 0.338 | 1.696 | 1.000 | 0.379 | 0.689 | 0.810 | 0.000 | 68.000 | 1.000 |
| E2_DRATT | full | Weak | A | 800 | 0.068 | 0.000 | 1.000 | 0.340 | 1.495 | 1.000 | 0.236 | 0.897 | 0.971 | 0.000 | 156.000 | 1.000 |
| E2_DRATT | full | Weak | B | 800 | 0.073 | 0.000 | 1.000 | 0.341 | 1.591 | 1.000 | 0.371 | 0.748 | 0.879 | 0.000 | 125.000 | 1.000 |
| E2_DRATT | reduced | RCT | RCT | 200 | 0.063 | 0.000 | 1.000 | 0.370 | 1.432 | 1.000 | 0.310 | 0.910 | 0.980 | 0.000 | 196.000 | 1.000 |
| E2_DRATT | reduced | Strong | A | 800 | 0.110 | 0.000 | 1.000 | 0.495 | 1.851 | 1.027 | 0.502 | 0.504 | 0.691 | 0.000 | 101.000 | 1.194 |
| E2_DRATT | reduced | Strong | B | 800 | 0.080 | 0.000 | 1.000 | 0.352 | 1.703 | 1.000 | 0.381 | 0.689 | 0.810 | 0.000 | 68.000 | 1.000 |
| E2_DRATT | reduced | Weak | A | 800 | 0.078 | 0.000 | 1.000 | 0.470 | 1.539 | 1.000 | 0.355 | 0.735 | 0.874 | 0.000 | 175.000 | 1.000 |
| E2_DRATT | reduced | Weak | B | 800 | 0.074 | 0.000 | 1.000 | 0.352 | 1.598 | 1.000 | 0.365 | 0.751 | 0.882 | 0.000 | 125.000 | 1.000 |

### T5 — Measured strength of the withheld confounder (E1) and pre-period placebo effect, by stratum

| strength | source | n | E1_cf_y_star_median | E1_cf_d_star_median | E1_rho_hat_star_median | E1_rho_true_median | E1_B_star_median | placebo_att_mean | placebo_att_sd | placebo_se_mean | E1_full_minus_red_att_mean | placebo_att_ci_lo | placebo_att_ci_hi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RCT | RCT | 200 | 0.0220 | 0.0000 | -0.3353 | -0.1181 | 0.0000 | -1.0873 | 16.9961 | 17.2781 | 0.6938 | -3.4429 | 1.2682 |
| Strong | A | 800 | 0.0125 | 0.3042 | 0.0994 | 0.5359 | 17.1659 | 140.2553 | 35.4219 | 31.1424 | 6.5921 | 137.8007 | 142.7099 |
| Strong | B | 800 | 0.0204 | 0.0000 | 0.0168 | 1.9501 | 0.0000 | 1.5550 | 48.1610 | 37.0035 | 1.0466 | -1.7823 | 4.8924 |
| Weak | A | 800 | 0.0169 | 0.2334 | 0.2050 | 0.4563 | 8.8803 | 94.8904 | 19.8981 | 20.1809 | 3.7526 | 93.5116 | 96.2693 |
| Weak | B | 800 | 0.0209 | 0.0000 | -0.0414 | 1.3480 | 0.0000 | 1.4359 | 30.1630 | 27.9555 | 0.9709 | -0.6543 | 3.5261 |

### T6 — Cross-family Spearman concordance (X_full)

| estimator | pair | spearman | ci_lo | ci_hi | n |
|---|---|---|---|---|---|
| E1_AIPW | rv~gamma_star | 0.472 | 0.447 | 0.498 | 3400 |
| E1_AIPW | rv~e_pt | 0.810 | 0.797 | 0.821 | 3400 |
| E1_AIPW | gamma_star~e_pt | 0.466 | 0.439 | 0.494 | 3400 |
| E2_DRATT | rv~gamma_star | 0.634 | 0.612 | 0.658 | 3400 |
| E2_DRATT | rv~e_pt | 0.975 | 0.973 | 0.977 | 3400 |
| E2_DRATT | gamma_star~e_pt | 0.616 | 0.593 | 0.638 | 3400 |
| E0_OLS | rv~e_pt | 0.990 | 0.989 | 0.991 | 3400 |

### T7 — Test detail columns

| id | cov_strong | cov_weak | mean_d_A | mean_d_B | cov_full_A | cov_red_A | cov_full_B | cov_red_B | naive_coverage_same_cohorts | n_ovb_only | n_naive_only | median_cf_y_star | median_cf_d_star | median_rho_true | trend_slope_logit | cov_RCT | cov_neural_dr | cov_e1 | width_neural_dr | width_e1 | rmse_neural_plugin | rmse_neural_dr | rmse_e1 | cov_seed_ensemble_interval(exploratory) | n_discordant | median_pehe_neural | median_pehe_ridge |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P1 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| P2 | 0.704 | 0.854 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| P3 |  |  | 0.250 | -0.003 | 0.876 | 0.626 | 0.683 | 0.686 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| P1 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| P2 | 0.723 | 0.823 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| P3 |  |  | 0.207 | -0.002 | 0.827 | 0.619 | 0.718 | 0.720 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| P1 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| P2 | 0.586 | 0.781 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| P3 |  |  | 0.226 | -0.003 | 0.733 | 0.507 | 0.634 | 0.637 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| S3_E1 |  |  |  |  |  |  |  |  | 0.626 | 290 | 0 | 0.015 | 0.270 | 0.498 |  |  |  |  |  |  |  |  |  |  |  |  |  |
| S3_E2 |  |  |  |  |  |  |  |  | 0.619 | 254 | 0 | 0.046 | 0.063 | 0.877 |  |  |  |  |  |  |  |  |  |  |  |  |  |
| S4_E1_AIPW | 0.704 | 0.854 |  |  |  |  |  |  |  |  |  |  |  |  | -0.831 | 0.900 |  |  |  |  |  |  |  |  |  |  |  |
| S4_E2_DRATT | 0.723 | 0.823 |  |  |  |  |  |  |  |  |  |  |  |  | -0.590 | 0.900 |  |  |  |  |  |  |  |  |  |  |  |
| S4_E0_OLS | 0.586 | 0.781 |  |  |  |  |  |  |  |  |  |  |  |  | -0.850 | 0.850 |  |  |  |  |  |  |  |  |  |  |  |
| S2_E5_DragonNet_full |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 0.676 | 0.824 | 61.153 | 118.982 | 31.085 | 34.802 | 48.047 | 0.618 | 29 |  |  |
| S5_E5_DragonNet_full |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 60.770 | 54.388 |
| S2_E4_TARNet_full |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 0.788 | 0.824 | 89.538 | 118.982 | 30.401 | 39.708 | 48.047 | 0.559 | 12 |  |  |
| S5_E4_TARNet_full |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 61.009 | 54.388 |
| S2_E5_DragonNet_reduced |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 0.576 | 0.718 | 58.858 | 82.393 | 35.481 | 37.800 | 38.361 | 0.482 | 28 |  |  |
| S5_E5_DragonNet_reduced |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 61.081 | 55.274 |
| S2_E4_TARNet_reduced |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 0.629 | 0.718 | 69.754 | 82.393 | 35.690 | 37.936 | 38.361 | 0.353 | 23 |  |  |
| S5_E4_TARNet_reduced |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 61.470 | 55.274 |