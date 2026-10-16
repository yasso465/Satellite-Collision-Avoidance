# C63 — validation of objects excluded during filtering

Source: `C63_excluded_validation.csv` supplied with the research data archive.

| Metric | Result |
|---|---:|
| Original TLE catalog | 16,883 objects |
| Retained and separately validated | 10,866 objects |
| Excluded objects checked with C63 | 6,017 |
| Excluded objects successfully propagated | 6,014 |
| Excluded objects with propagation errors | 3 |
| Smallest 1 s sampled separation among successful excluded objects | 60.872769 km |
| NORAD ID of nearest excluded object | 47945 |
| Excluded encounters refined at 0.1 s | 0 |

Three excluded objects (NORAD **58274, 58285, 60491**) failed with `Unable to add satellite to the satelliteScenario.`. Their trajectories were **not validated**, and this failure alone does not establish reentry.

This check applies specifically to **candidate C63**, not merely an earlier reference maneuver. It is separate from the retained-catalog validation. The minimum is based on 1-second sampling; no continuous-time separation guarantee is implied.

For full row-level traceability, consult the source CSV in the original data archive until that larger CSV is uploaded to the repository.
