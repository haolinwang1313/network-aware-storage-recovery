# Reproducibility

## Model and evidence scope

The constant-condition formulation uses DC power flow, charge-only recovery, and event-feasible storage states to determine an exact continuous minimum recovery duration. The chronological formulation evaluates charging availability within a fixed 18-hour interval before the next event and includes generator ramping and bidirectional storage. Its charging-window variable does not move the event onset.

Weather-driven cases use observed temperatures from three NOAA stations over September 1-10, 2022. The mapping to the IEEE RTS-24 system and electrical temperature-response coefficients are controlled assumptions. The 54 cases are parameter combinations within this weather interval. Recovery decisions are common across specified future temperature scenarios, with the event trajectory revealed at event onset.

The saved chronological sweep contains 50 cases with a zero-unserved-energy recovery solution and four without one in the available interval. At the same minimum charging window, the myopic baseline has positive unserved energy in 12 of the 50 recoverable cases. These comparisons depend on the stated model, information, and intertemporal constraints.

The saved AC records found no complete AC realization of the tested RTS-24 trajectories. A solver's failure to find a realization does not prove physical infeasibility. The independent analytical two-bus AC example addresses the location dependence of stored energy and is a separate construction.

## Scientific question to output mapping

| Question or check | Saved output | Generator |
|---|---|---|
| Constant-condition parameter sweep | `results/summary.csv` | `src/run_experiments.py` |
| Two-bus equal-energy comparison | `results/equal_soc_counterexample.json` | `src/run_experiments.py` |
| 54 chronological parameter combinations | `results_v2/weather_summary.json`, `results_v2/weather_summary.csv` | `src/run_validation_v2.py` |
| Complete hourly trajectories | `results_v2/all_weather_trajectories.json.gz` | `src/run_validation_v2.py` |
| All charging windows for September 5-6 | `results_v2/all_windows_day5.json`, `results_v2/all_windows_day5.csv` | `src/run_validation_v2.py` |
| Storage LP/MILP boundary comparison | `results_v2/lp_milp_boundary.json` | `src/run_validation_v2.py` |
| Parameter sensitivities | `results_v2/sensitivities.json`, `results_v2/sensitivities.csv` | `src/run_validation_v2.py` |
| Forecast perturbation comparison | `results_v2/forecast_envelope.json`, `results_v2/forecast_envelope.csv` | `src/run_validation_v2.py` |
| AC trajectory replay and static back-check | `results_v2/ac_replay.json`, `results_v2/legacy_ac_check.json` | `src/run_validation_v2.py --ac` |
| Additional AC trajectory study | `results_v2/ac_followup.json` | `tools/ac_followup.py` |
| Independent analytical AC construction | `results_v2/ac_analytic_example.json` | `tools/ac_analytic_example.py` |
| Execution environment and input identity | `results/run_metadata.json`, `results_v2/run_metadata.json` | Corresponding experiment driver |

## Figure mapping

The following companion plots are produced by `src/plot_results.py`. Each plot has an editable SVG and a PNG rendering in `results_v2/figures/`.

| Plot stem | Numerical input | Interpretation |
|---|---|---|
| `observed_temperature` | Prepared `data/public/weather_hourly.csv` | Observed station temperature traces |
| `charging_window` | `results_v2/all_windows_day5.json` | Unserved energy versus charging availability before a fixed event onset |
| `forecast_envelope` | `results_v2/forecast_envelope.json` | Nominal and specified-envelope designs under future temperature perturbations |

The saved two-bus energy-location comparison is available in `results/equal_soc_counterexample.json`. A new run of `src/run_experiments.py` also produces `summary.json` and `two_bus_detail.json` in its output directory; these additional files are not included among the distributed results.

## Reproduction sequence

1. Install the pinned dependencies in the Python environment described in [usage](usage.md).
2. Retrieve the source inputs and prepare the reduced RTS-24 case using the two documented preparation scripts.
3. Run the model tests and analytical AC checks.
4. Run full experiments into new `generated/` directories.
5. Compare result statuses, energy quantities, and charging-window boundaries with the distributed outputs.
6. Generate companion plots from the saved or newly generated numerical files.

Input retrieval requires network access. Once preparation is complete, model execution uses the local input files. Saved numerical outputs and figures are directly inspectable without the input preparation step.

Record the software commit, Python and solver versions, preparation output, and any local parameter changes when reporting a new run. Run metadata includes timestamps and platform details, so byte-for-byte equality of newly generated metadata is not an expected test. Numerical comparisons should account for solver tolerances; retain unsuccessful and no-advantage cases when summarizing a sweep.

## Local packaging checks

The 29 model tests and 3 analytic AC tests passed with existing input fixtures in a separate test copy, using NumPy 1.26.4 and SciPy 1.13.0. The original experiment environment is recorded in `results_v2/run_metadata.json`; it was not recreated for these checks. Windows additionally needs `tzdata` for the named Pacific time zone. The acquired reduced case helper reproduces every original physical field. Network acquisition was not exercised in this check.
