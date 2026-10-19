# Minimum Charging Time Windows for Battery Energy Storage Between Heatwave Demand Peaks

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/downloads/)

Research software and saved numerical results for determining how much charging time battery energy storage systems need between successive heatwave demand peaks. For given storage locations and ratings, the models determine whether inter-peak operation can prepare stored-energy states that meet demand in every specified next-peak scenario.

**Authors:** Haolin Wang, Zhi Wu, Wei Gu, and Shuai Yao.

## Models and scope

- **Unchanged operating conditions:** an exact minimum-duration linear program for charge-only recharging under convex DC network constraints.
- **Hourly operating conditions:** joint generation and storage scheduling determines the minimum charging window before a fixed next-peak onset, including generator ramping and bidirectional storage operation.
- **Saved comparisons:** two-bus examples and 54 RTS-24 parameter combinations. In the representative case, the proposed schedule requires 9 h, compared with 10 h for a one-step baseline using current-hour information only. Within the available 18 h, 50 cases admit zero-unserved-energy schedules; the baseline has positive unserved energy at the same minimum window in 12 of those 50 cases.

The cases use observed temperatures from September 1-10, 2022, with assumed station-to-bus mapping and electrical temperature responses. They are parameter combinations from one weather interval, not 54 independent heatwaves. Next-peak scenarios share one inter-peak schedule and are revealed at peak onset. The charging window does not shift that onset.

The network-scale scheduling results use DC power flow. Saved RTS-24 AC checks found no complete AC realization of the tested trajectories; unsuccessful searches do not prove AC infeasibility. The independent analytical two-bus AC construction is a separate example.

## Repository structure

```text
.
|-- README.md             # Overview, usage, and result mappings
|-- CITATION.cff          # Software citation metadata
|-- LICENSE              # MIT license
|-- requirements.txt     # Python dependencies
|-- src/                 # Models, experiment drivers, and plotting
|-- tests/               # Model checks
|-- tools/               # Input preparation and additional AC studies
|-- data/                # Source catalog and data documentation
|-- results/             # Saved constant-condition results
`-- results_v2/          # Saved chronological results, AC checks, and plots
```

## Quick start

Run commands from the repository root. Saved results can be inspected without installing dependencies or downloading inputs. To execute the models, create an environment:

```bash
uv venv --python 3.13
uv pip install -r requirements.txt
```

Activate it with `source .venv/bin/activate` on Linux/macOS or `.\.venv\Scripts\Activate.ps1` in PowerShell. Then retrieve inputs and run the checks:

```bash
python tools/fetch_public_inputs.py
python tools/prepare_reduced_case.py
python -m unittest discover -s tests -v
python -m unittest discover -s tools -p test_ac_analytic.py -v
```

To reproduce the experiment families and generate companion plots:

```bash
python src/run_experiments.py --out generated/results
python src/run_validation_v2.py --out generated/results_v2 --ac
python src/plot_results.py --results generated/results_v2
```

Use a new, empty output directory for each run; the drivers do not remove stale outputs. `--ac` adds AC replay and the static-case back-check, while `--quick` runs only a development subset. Full chronological runs can be substantially more expensive than reading saved results. Plotting requires the prepared weather input. The separate AC follow-up tools use fixed result paths and should be run only in a separate working copy.

## Results and manuscript mapping

| Manuscript item or check | Saved numerical source |
|---|---|
| Fig. 1: two-bus equal-energy comparison | [equal_soc_counterexample.json](results/equal_soc_counterexample.json) |
| Fig. 2: charging-window comparison | [all_windows_day5.csv](results_v2/all_windows_day5.csv), [JSON](results_v2/all_windows_day5.json) |
| Fig. 3: temperature-offset comparison | [forecast_envelope.csv](results_v2/forecast_envelope.csv), [JSON](results_v2/forecast_envelope.json) |
| Table I and the 54-case summary | [weather_summary.csv](results_v2/weather_summary.csv), [JSON](results_v2/weather_summary.json) |
| Constant-condition sweep | [summary.csv](results/summary.csv) |
| Ramping and temperature-response sensitivities | [sensitivities.csv](results_v2/sensitivities.csv) |
| LP/MILP boundary comparison | [lp_milp_boundary.json](results_v2/lp_milp_boundary.json) |
| Complete chronological trajectories | [all_weather_trajectories.json.gz](results_v2/all_weather_trajectories.json.gz) |

Table I selects the original line-rating cases with 200-MW chargers. The [companion plots](results_v2/figures/) show observed temperatures, charging windows, and temperature offsets; they do not reproduce the manuscript's figure numbering or final layout. The plotter generates these companion plots, not the manuscript's two-bus schematic. A new constant-condition run also writes `summary.json` and `two_bus_detail.json`, which are not distributed reference files.

Additional AC records are [ac_replay.json](results_v2/ac_replay.json), [legacy_ac_check.json](results_v2/legacy_ac_check.json), [ac_followup.json](results_v2/ac_followup.json), and [ac_analytic_example.json](results_v2/ac_analytic_example.json). Retain unsuccessful cases when interpreting these records.

## Data and reproducibility

Input preparation retrieves NOAA ISD-Lite observations for Sacramento, Fresno, and San Francisco and the MATPOWER RTS-24 case. It requires network access and prepares 240 hourly observations per station. See the [input documentation](data/README.md), [source catalog](data/catalog.yaml), and [data dictionary](data/data_dictionary.md) for provenance, checksums, and field definitions. Downloaded inputs are not bundled in this repository.

The original experiments used Python 3.13, NumPy 2.3.5, and SciPy 1.17.0; their environments and input identities are recorded in [constant-condition metadata](results/run_metadata.json) and [chronological metadata](results_v2/run_metadata.json). These are historical run records, including the original source commit identifiers and scope assessments.

Separate packaging checks passed 29 model tests and 3 analytical AC tests using existing input fixtures with NumPy 1.26.4 and SciPy 1.13.0. Those checks did not recreate the original experiment environment or exercise network acquisition. Compare numerical quantities using solver tolerances; retrieval times and platform metadata need not match byte for byte. Passing the independent two-bus AC tests does not establish AC feasibility of the RTS-24 schedules.

## License

The project code is released under the [MIT License](LICENSE). Retrieved third-party inputs remain subject to their source terms; the software license does not grant redistribution rights for those inputs.
