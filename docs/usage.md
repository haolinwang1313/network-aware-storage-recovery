# Usage

Run all commands from the repository root. Saved CSV, JSON, SVG, and PNG files can be read without installing the optimization environment or downloading the input data.

## Environment

The original experiments used Python 3.13. Create an isolated environment with the pinned dependencies:

```bash
uv venv --python 3.13
uv pip install -r requirements.txt
```

Activate it in PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Or in a POSIX shell:

```bash
source .venv/bin/activate
```

The remaining commands use this environment's `python` executable.

## Prepare inputs

```bash
python tools/fetch_public_inputs.py
python tools/prepare_reduced_case.py
```

The first command requires network access and retrieves public weather and full RTS-24 inputs. The second prepares `data/rts24.json` for the reduced DC case. Review [input documentation](../data/README.md) for source provenance and retrieval requirements. Full model tests and experiment execution require these inputs.

## Run tests

```bash
python -m unittest discover -s tests -v
python -m unittest discover -s tools -p test_ac_analytic.py -v
```

The second command checks the independent two-bus AC example. Passing that example does not establish AC realizability of the RTS-24 trajectories.

## Run experiments

```bash
python src/run_experiments.py --out generated/results
python src/run_validation_v2.py --out generated/results_v2 --ac
python src/plot_results.py --results generated/results_v2
```

`run_experiments.py` writes the constant-condition sweep and analytical example. `run_validation_v2.py` writes the 54-case chronological sweep, charging-window boundary checks, sensitivity and forecast results, trajectories, and run metadata. `--ac` adds AC replay and the static-case AC back-check. The separate `ac_followup.json` study is described below.

Use an empty output directory for each run. The drivers write incremental results and do not remove stale files left by earlier executions. Their `--quick` option runs a development subset, so use the full commands above when comparing with the saved full-study results.

## Regenerate plots from saved results

After preparing the weather input, copy the saved JSON result inputs into a new plotting directory. This command works in either shell:

```bash
python -c "from pathlib import Path; import shutil; p=Path('generated/plots'); p.mkdir(parents=True, exist_ok=True); [shutil.copy2(Path('results_v2')/name,p/name) for name in ('all_windows_day5.json','forecast_envelope.json')]"
python src/plot_results.py --results generated/plots
```

The plotter reads those JSON files and `data/public/weather_hourly.csv`, then writes `observed_temperature`, `charging_window`, and `forecast_envelope` in SVG and PNG formats to `generated/plots/figures/`, along with a figure manifest. These are the companion result plots.

## Additional AC studies

The saved additional AC study is `results_v2/ac_followup.json`; the independent construction is `results_v2/ac_analytic_example.json`. `tools/ac_followup.py` and `tools/ac_analytic_example.py` use fixed result paths. To rerun them, first make a separate working copy of the repository with prepared input data, then execute them there. This keeps the distributed reference results available for comparison.

## Inspect outputs

See [reproducibility](reproducibility.md) for the mapping between scientific questions and files. Compare statuses and numerical quantities with solver tolerances rather than expecting identical run timestamps or platform metadata. An unsuccessful AC replay is a recorded scientific outcome, not evidence of AC infeasibility.
