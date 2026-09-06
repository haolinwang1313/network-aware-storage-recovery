# A Network-Constrained Recovery Criterion for Energy Storage During Heatwaves

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/downloads/)

This repository is the research software companion to a study of storage recovery between successive heatwave stress episodes. It contains network-constrained recovery models, deterministic experiment drivers, tests, saved numerical results, and editable plots. The central question is which storage states can supply a specified next event and how soon the system can reach those states.

## Overview

- **Recovery criterion:** event-feasible storage states account for both energy availability and network constraints.
- **Constant recovery conditions:** an exact continuous-time linear program determines the minimum recovery duration under charge-only operation.
- **Time-varying conditions:** chronological optimization determines the minimum charging window before a fixed next-event onset, including generator ramping and bidirectional storage operation.
- **Saved evidence:** synthetic two-bus and RTS-24 tests, 54 weather-driven parameter combinations, charging-window comparisons, sensitivities, forecast perturbations, and AC replay records.

The weather-driven cases combine observed station temperatures with an assumed mapping to the RTS-24 system and assumed electrical temperature responses. They describe controlled model tests. The recovery results use DC power flow; the saved RTS-24 AC checks found no complete AC realization of the tested trajectories. See [reproducibility notes](docs/reproducibility.md) for interpretation and result-file mappings.

## Repository Structure

```text
.
|-- README.md                         # Project overview and entry points
|-- CITATION.cff                      # Software citation metadata
|-- LICENSE                           # MIT license
|-- requirements.txt                  # Pinned Python dependencies
|-- src/
|   |-- recovery.py                   # Storage and network recovery model
|   |-- minimum_time.py               # Exact constant-condition duration LP
|   |-- chronology.py                 # Fixed-clock chronological optimization
|   |-- cases.py                      # Synthetic two-bus and RTS-24 cases
|   |-- weather_cases.py              # Weather-to-grid scenario construction
|   |-- ac_validation.py              # AC trajectory checks
|   |-- run_experiments.py            # Constant-condition experiment driver
|   |-- run_validation_v2.py          # Weather-driven validation driver
|   `-- plot_results.py               # SVG and PNG result plots
|-- tests/                            # Model and solver checks
|-- tools/
|   |-- fetch_public_inputs.py        # Retrieve public weather and RTS inputs
|   |-- prepare_reduced_case.py       # Prepare the reduced RTS-24 input
|   |-- ac_followup.py                # Additional AC trajectory checks
|   |-- ac_analytic_example.py        # Independent two-bus AC construction
|   `-- test_ac_analytic.py           # Analytic AC example checks
|-- data/
|   |-- README.md                     # Input sources and preparation instructions
|   |-- catalog.yaml                  # Source and reference checksum catalog
|   `-- data_dictionary.md            # Input and result fields
|-- results/                          # Saved constant-condition results
|-- results_v2/                       # Saved chronological and AC results
|   `-- figures/                      # Editable SVG and rendered PNG plots
`-- docs/
    |-- usage.md                      # Installation and execution commands
    `-- reproducibility.md            # Result mappings and model scope
```

## Dependencies & Installation

The original experiments used Python 3.13. Install the pinned dependencies in an isolated environment; the model uses open-source solvers. Run commands from the repository root.

```bash
uv venv --python 3.13
uv pip install -r requirements.txt
```

Use `.venv/bin/python` on Linux/macOS or `.venv\Scripts\python.exe` on Windows for the `python` commands below, or activate the environment first. See [usage](docs/usage.md) for activation instructions.

## Usage

Saved JSON, CSV, and figure files can be inspected directly. To run tests and experiments, first retrieve and prepare the inputs:

```bash
python tools/fetch_public_inputs.py
python tools/prepare_reduced_case.py
python -m unittest discover -s tests -v
python -m unittest discover -s tools -v
```

Run the two experiment families into separate generated-output directories:

See [experiment and plotting commands](docs/usage.md#run-experiments). The commands write to fresh directories under `generated/` so the distributed reference results remain available.

The full chronological experiment includes optimization and AC replay and can be substantially more expensive than reading the saved results. [Usage](docs/usage.md) also explains how to regenerate plots from saved numerical outputs.

## Data Availability

The input preparation script retrieves NOAA ISD-Lite observations and the MATPOWER RTS-24 case. The weather interval is September 1-10, 2022, with three stations and 240 hourly observations per station. Network access is required for input retrieval. Run `prepare_reduced_case.py` after retrieval to create the reduced case used by the constant-condition experiments. Consult [data documentation](data/README.md) for source URLs, attribution, and preparation details.

The saved numerical outputs and figures are available in `results/` and `results_v2/`. The complete chronological trajectories are stored as `results_v2/all_weather_trajectories.json.gz`.

## License

The project code is released under the [MIT License](LICENSE). Retrieved third-party inputs remain subject to their respective source terms; see [data documentation](data/README.md).
