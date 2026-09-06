# Input data

The study uses the IEEE RTS-24 network and 240 hourly temperatures for Sacramento, Fresno, and San Francisco from September 1–10, 2022 (Pacific daylight time). The electrical mapping and temperature-response coefficients are synthetic modeling assumptions.

## Acquisition

Run from the repository root, after installing dependencies:

```bash
python tools/fetch_public_inputs.py
python tools/prepare_reduced_case.py
```

The first command contacts NOAA NCEI and MATPOWER and writes `data/public/`. It uses the current upstream MATPOWER file; compare the reference SHA-256 and Git blob before interpreting a rerun as an exact reproduction. The second command writes `data/rts24.json` from the downloaded full network.

`catalog.yaml` records reference file sizes, hashes, source URLs, and acquisition commands. `acquisition_path` gives the local destination; `path: TODO` marks input acquisition pending. Redistribution permission is pending verification; consult upstream terms before redistributing downloaded inputs. The software MIT license does not assign a license to upstream data.

The source network is MATPOWER `case24_ieee_rts.m`, Git blob `862c195c98e319e858829ebc0dea0f928694b71d`, credited to the IEEE Reliability Test System Task Force, RTS-96 report, and MATPOWER case contributors. NOAA station IDs are `724830-23232`, `723890-93193`, and `724940-23234`.

Reference weather CSV uses UTC and local ISO timestamps, temperature in degrees Celsius, and flags for filled hours. Downloader-generated metadata has a new retrieval time. The numerical tests verify the weather CSV checksum; a mismatch must be investigated rather than updating the expected value automatically.
