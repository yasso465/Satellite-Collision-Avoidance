# Satellite Collision Avoidance — Two-Impulse Maneuver Optimization

MATLAB numerical study for the MathWorks Satellite Collision Avoidance challenge (#225).

> **Status: research prototype.** MATLAB source files, TLE catalogs and detailed result CSVs have been uploaded; a clean-checkout MATLAB reproduction test and final animations are still pending. The 10.5 km distance is an **experimental screening threshold**, not a physical collision criterion.

## Objective

Find low-total-ΔV impulsive maneuvers that avoid close approaches to catalogued orbital debris and return a spacecraft to approximately its initial orbital elements over a three-hour simulation.

- Epoch: **2026-10-12 00:00 UTC**
- Duration: **10,800 s**
- Initial spacecraft orbit: circular, **500 km** altitude, **51.6°** inclination
- Spacecraft model: two-body Keplerian propagation with impulsive tangential burns
- Debris model: SGP4 propagation from three-line TLE records
- Return-orbit tolerances tested: |Δa| ≤ **5 km**, |Δe| ≤ **0.001**

## Selected candidate: C63

| Quantity | Value |
|---|---:|
| First burn time | 2,600 s (00:43:20 UTC) |
| First tangential ΔV | −0.40 m/s |
| Second burn time | 8,000 s (02:13:20 UTC) |
| Second tangential ΔV | +0.20 m/s |
| Total ΔV | **0.60 m/s** |
| Final semi-major-axis difference | **0.361384 km** |
| Final eccentricity difference | **0.000057205** |
| Closest nominal approach in 10,866 retained objects (locally refined to 0.1 s) | **11.735 km** |
| Object at closest nominal approach | NORAD **44804** |
| Local Monte Carlo trials | **10,000** |
| Trials with minimum below 10.5 km | **0 observed** |
| Smallest separation observed in local Monte Carlo | **10.608 km** |

Compared with the earlier 0.90 m/s reference maneuver, C63 uses **33.3% less ΔV**. This is not a proof of global optimality.

## Published result files

- [Nominal retained-catalog validation summary](results/early_burn_full_catalog_summary.csv)
- [Paired 10,000-trial Monte Carlo summary](results/early_burn_mc_10000_summary.csv)
- [C63 excluded-object validation summary](docs/C63_EXCLUDED_VALIDATION.md)

The full row-level CSVs and TLE inputs are now included in [`results/`](results/) and [`data/`](data/).

## Method

1. Filter the original 16,883-entry TLE catalog to identify candidate close approaches.
2. Propagate the nominal spacecraft trajectory and SGP4 debris trajectories.
3. Screen two-burn maneuvers for minimum separation, total ΔV, and return-orbit tolerances.
4. Validate promising candidates against the **10,866 retained objects**, using 1-second global sampling and 0.1-second refinement around close encounters.
5. Compare candidates using paired Monte Carlo perturbations in **29 local encounter windows**.

### Uncertainty assumptions

The local Monte Carlo assumed independent Gaussian perturbations: debris position σ=250 m per inertial axis, debris velocity σ=0.05 m/s per axis, burn timing σ=1 s, tangential burn magnitude σ=0.01 m/s, and radial/normal burn components σ=0.005 m/s per axis.

These are **assumed uncertainties**, not covariance estimates inferred from the TLEs. Consequently, the observed fraction of runs below the experimental 10.5 km threshold is **not** a physical collision probability.

## Validation scope and limitations

- The retained-catalog validation covered **10,866 objects**, with no SGP4 errors reported in that run.
- The original catalog contains **16,883 objects**, including **6,017 initially excluded**. A separate **C63-specific** excluded-object CSV was supplied and checked: **6,014 of 6,017** excluded objects propagated successfully, with a minimum 1-second sampled separation of **60.872769 km** (NORAD 47945); three objects (58274, 58285, 60491) could not be added to the scenario. See [excluded validation summary](docs/C63_EXCLUDED_VALIDATION.md). Thus 16,880 of 16,883 catalog entries were numerically checked across the retained and excluded sets; three remain unvalidated.
- Three SGP4 failures were previously identified for NORAD **58274, 58285, 60491**. Propagation failure does not independently establish atmospheric reentry.
- Monte Carlo perturbations were assessed in selected **local windows**, not as an exhaustive perturbed full-catalog scan.
- Numerical sampling at 1 s and 0.1 s does not establish a rigorous continuous-time lower bound.
- The project does **not** currently include a trained AI model.

## Quick start: standalone spacecraft demo

Run `src/demo_C63_two_body.m` in MATLAB. This standalone example plots the baseline orbit and the C63 two-burn trajectory, and prints final orbital-element differences. It does **not** propagate debris or validate closest approaches. The validation scripts and TLE inputs are present in this repository, but their clean-checkout execution has not yet been verified.

## Reproduction

**Software:** MATLAB R2026b, with Aerospace Toolbox for SGP4/TLE workflows. Additional dependencies will be documented alongside the final source files.

The MATLAB scripts, required TLE files, and generated result CSVs are now available. See [running instructions](docs/RUNNING.md) and [data manifest](docs/DATA_MANIFEST.md). **Important:** the scripts have not yet been re-executed from a clean checkout, and some may expect input files in the MATLAB current folder rather than in `data/` or `results/`. A tested one-command reproduction workflow is still pending.

## Repository layout

```text
Satellite-Collision-Avoidance/
├── README.md
├── LICENSE
├── src/          # MATLAB simulations, search, validation
├── data/         # TLE input and example subset
├── results/      # Numerical CSV summaries
├── figures/      # Plots
└── animations/   # 3D scenario and close-approach videos
```

## License

MIT (see LICENSE when added).
