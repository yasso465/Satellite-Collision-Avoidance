# Reproducing the study

## Environment

- MATLAB R2026b was used during development.
- Aerospace Toolbox is required for the TLE/SGP4 workflow.
- Check the individual MATLAB script headers for their expected input filenames.

## Required inputs (not yet in this repository)

- `debris_matlab.tle`: filtered, MATLAB-readable 3-line TLE catalog (10,866 objects).
- `debris_leo.tle`: original 3-line TLE catalog (16,883 objects).
- `debris_filtered.tle`: retained original records, required for the excluded-object validation script.
- `robust09_high_precision_validation.csv`: previous candidate's local-encounter validation, used by the early-burn search.
- `early_burn_full_catalog_results.csv`: output from nominal validation, required for the paired Monte Carlo script.

## Intended script sequence

1. `search_early_burn_times.m` — exploratory screening; requires the previous candidate's validation CSV.
2. `validate_early_burn_full_catalog.m` — nominal full **retained-catalog** validation of C3, C63, C123.
3. `monte_carlo_early_burn_10000.m` — paired local Monte Carlo for C3/C63.
4. `validate_C63_excluded.m` — separate validation of objects excluded during filtering; review the resulting log before making full-original-catalog claims.

These scripts are **not yet published in this repository**. The sequence documents the existing workflow; it is not an executable quick-start until code and input files are uploaded and checked.

## Scientific caution

The 10.5 km threshold is a conservative screening metric used in this experiment, not a collision radius. The Monte Carlo uncertainty model is hypothetical and does not represent catalog-derived covariance. Numerical results have not been independently reproduced from the GitHub repository yet.
