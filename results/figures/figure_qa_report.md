# Figure QA report (alterlab-figure-qa)

- **F1 data fidelity**: PASS — 30 plotted coverage points re-derived from fast_arm_estimates.csv equal summary_by_stratum.csv (mismatches: 0); all 5 strata × 2 adjustment sets × 3 estimators present
- **F2 data fidelity**: PASS — 90 calibration points recomputed (mismatches: 0); RCT n=200, Weak/Strong n=1600 per curve
- **F4 data fidelity**: PASS — grid table has 3400 cohorts; diagonal coverage A@0.03 = 0.811, B@0.03 = 0.896 (match sensitivity_report §5); measured points = 3,400 cohort rows of fast_arm_cohort.csv
- **F5 data fidelity**: PASS — 3,200 confounded cohorts; zero-strength bin n = 1790 (figure shows 1790); positive quintiles n = 282 each
- **F6 data fidelity**: PASS — 3 Spearman coefficients recomputed from fast_arm_estimates.csv match concordance.csv (n = 3400 each); floor mass Γ*=1: 0.658
- **F3 data fidelity**: PASS — PEHE rows per learner: {'E1_AIPW': 0, 'E4_TARNet': 170, 'E5_DragonNet': 170, 'E6_ridgeT': 170}; cohorts = 170

## Export / legibility (from *_provenance.json 'qa')

- **F1_coverage_ladder**: PASS — text collisions 0, clipped 0, min font 6.0 pt, PNG [2125, 732] px @ 300 dpi, size [180.0, 62.0] mm; PDF (vector, Type 42 fonts) alongside
- **F2_calibration_curves**: PASS — text collisions 0, clipped 0, min font 7.0 pt, PNG [2125, 732] px @ 300 dpi, size [180.0, 62.0] mm; PDF (vector, Type 42 fonts) alongside
- **F3_pehe_by_regime**: PASS — text collisions 0, clipped 0, min font 6.5 pt, PNG [2125, 921] px @ 300 dpi, size [180.0, 78.0] mm; PDF (vector, Type 42 fonts) alongside
- **F4_ovb_sensitivity_contour**: PASS — text collisions 0, clipped 0, min font 6.0 pt, PNG [2125, 850] px @ 300 dpi, size [180.0, 72.0] mm; PDF (vector, Type 42 fonts) alongside
- **F5_coverage_vs_measured_strength**: PASS — text collisions 0, clipped 0, min font 6.0 pt, PNG [2125, 767] px @ 300 dpi, size [180.0, 65.0] mm; PDF (vector, Type 42 fonts) alongside
- **F6_concordance**: PASS — text collisions 0, clipped 0, min font 7.0 pt, PNG [2125, 708] px @ 300 dpi, size [180.0, 60.0] mm; PDF (vector, Type 42 fonts) alongside

**Overall: ALL PASS**