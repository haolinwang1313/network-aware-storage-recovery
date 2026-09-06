# Data dictionary

| File or table | Fields and units |
|---|---|
| `data/public/weather_hourly.csv` | `utc`, `local`: ISO timestamps; city `_c` columns: temperature (degrees Celsius); `_filled`: interpolation flags |
| `data/public/rts24_full.json` | MATPOWER version 2 bus, branch, generator and cost matrices; `baseMVA` and upstream provenance |
| `data/rts24.json` | Active load and bus-aggregated generation (MW); branches: 1-based endpoints, reactance (p.u.), rateA interpreted as DC MW limit, transformer tap |
| `results_v2/weather_summary.csv` | Date pair, line multiplier, charge rating, status, minimum available window (h), energy not served (MWh) |
| `results_v2/all_windows_day5.csv` | Proposed and baseline ENS for every available window from 0 to 18 h |
| `results_v2/forecast_envelope.csv` | Temperature offsets, design windows, and ENS |
| `results_v2/sensitivities.csv` | Controlled parameter variations and outcomes |

Status and numerical tolerances are defined by the model implementation. All 54 combinations are retained, including four without feasible recovery in the available interval. Time-window length is charging availability before a fixed event start, not a shifted event time.
