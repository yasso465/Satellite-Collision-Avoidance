# Data provenance and publication status

Source archive supplied by the project author: `adaptive_conjunctions.zip` (65 entries, including intermediate experiments).

## Core reproducibility files

| Source file | Target directory | Role |
|---|---|---|
| `debris_leo.tle` | `data/` | Original 16,883-object 3-line TLE catalog |
| `debris_filtered.tle` | `data/` | Original-format retained-object catalog |
| `debris_matlab.tle` | `data/` | MATLAB-compatible 10,866-object TLE catalog |
| `early_burn_full_catalog_results.csv` | `results/` | Detailed retained-catalog screening |
| `early_burn_full_catalog_summary.csv` | `results/` | Summary of C3, C63, C123 |
| `early_burn_mc_10000_trials.csv` | `results/` | Trial-by-trial paired Monte Carlo |
| `early_burn_mc_10000_summary.csv` | `results/` | Paired Monte Carlo summary |
| `C63_excluded_validation.csv` | `results/` | C63-specific excluded-catalog audit |
| `robust09_high_precision_validation.csv` | `results/` | Input for early-burn search |
| `C63_visual_debris_10.tle` | `data/` | Small visualization subset |

## Scientific cautions

- The original catalog and the MATLAB-compatible catalog are not interchangeable.
- The 10,866 retained objects were screened for C63, and 6,014 of 6,017 excluded objects were propagated successfully in a separate C63 audit. Three excluded objects could not be added to the SGP4 scenario.
- Monte Carlo perturbations use assumed error scales; the observed 10.5 km threshold crossing rate is not a collision probability.
- Do not use the older `excluded_two_maneuvers_validation.csv` as proof of C63 validation; use `C63_excluded_validation.csv`.

Large data files should be reviewed for provenance, permitted redistribution, and reproducibility before public publication. The project README must not claim they are available on GitHub until they have actually been committed.
