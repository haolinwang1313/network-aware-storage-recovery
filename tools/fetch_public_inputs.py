"""Fetch public NOAA observations and the exact upstream MATPOWER case.

Run explicitly; never download during numerical tests. All source URLs and
SHA-256 values are saved. No empirical grid-response coefficients are inferred.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import io
import json
import re
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.request import Request, urlopen
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'public'
STATIONS = {'sacramento': '724830-23232', 'fresno': '723890-93193',
            'san_francisco': '724940-23234'}
TZ = ZoneInfo('America/Los_Angeles')
START = datetime(2022, 9, 1, tzinfo=TZ)
HOURS = 240


def get(url: str) -> bytes:
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={'User-Agent': 'PES-Letter-Research/1.0'}), timeout=60) as r:
                return r.read()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError('Unreachable')


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    history_url = 'https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv'
    history_raw = get(history_url)
    history = list(csv.DictReader(io.StringIO(history_raw.decode())))
    metadata = {'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                'timezone': 'America/Los_Angeles', 'start_local': START.isoformat(),
                'hours': HOURS, 'history_url': history_url,
                'history_sha256': hashlib.sha256(history_raw).hexdigest(),
                'interpretation': 'Observed station temperatures; assignment to RTS buses and grid sensitivities are synthetic.',
                'stations': {}}
    series = {}
    timestamps = [START.astimezone(timezone.utc) + timedelta(hours=i) for i in range(HOURS)]
    for key, station in STATIONS.items():
        matching = [v for v in history if v['USAF'] + '-' + v['WBAN'] == station]
        if not matching:
            raise RuntimeError(f'Cannot verify station metadata: {station}')
        url = f'https://www.ncei.noaa.gov/pub/data/noaa/isd-lite/2022/{station}-2022.gz'
        raw = get(url)
        lines = gzip.decompress(raw).decode('ascii').splitlines()
        records = {}
        kept = []
        for line in lines:
            fields = line.split()
            stamp = datetime(*map(int, fields[:4]), tzinfo=timezone.utc)
            if stamp in timestamps:
                if stamp in records:
                    raise RuntimeError(f'Duplicate hourly record: {station} {stamp}')
                value = int(fields[4])
                records[stamp] = None if value == -9999 else value / 10
                kept.append(line)
        (OUT / f'{key}_isd_lite_subset.txt').write_text('\n'.join(kept) + '\n')
        values = [records.get(t) for t in timestamps]
        raw_values = values.copy()
        missing = [i for i, v in enumerate(values) if v is None]
        for i in missing:
            left = i - 1
            while left >= 0 and raw_values[left] is None:
                left -= 1
            right = i + 1
            while right < HOURS and raw_values[right] is None:
                right += 1
            if left < 0 or right == HOURS or right - left - 1 > 2:
                raise RuntimeError(f'Gap exceeds declared two-hour interpolation limit at {station}:{i}')
            values[i] = raw_values[left] + (raw_values[right] - raw_values[left]) * (i-left) / (right-left)
        if any(not (-50 <= v <= 60) for v in values):
            raise RuntimeError('Invalid temperature or unit conversion.')
        series[key] = (values, set(missing))
        metadata['stations'][key] = {'station_id': station, 'metadata': matching,
            'url': url, 'download_sha256': hashlib.sha256(raw).hexdigest(),
            'raw_hourly_records': len(records), 'filled_hours': missing,
            'min_c': min(values), 'max_c': max(values)}
    with (OUT / 'weather_hourly.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        keys = list(STATIONS)
        writer.writerow(['utc','local'] + [k+'_c' for k in keys] + [k+'_filled' for k in keys])
        for i, stamp in enumerate(timestamps):
            writer.writerow([stamp.isoformat(),stamp.astimezone(TZ).isoformat()] +
                [series[k][0][i] for k in keys] + [int(i in series[k][1]) for k in keys])
    metadata['processed_sha256'] = hashlib.sha256((OUT/'weather_hourly.csv').read_bytes()).hexdigest()
    save(OUT / 'weather_metadata.json', metadata)
    url = 'https://raw.githubusercontent.com/MATPOWER/matpower/master/data/case24_ieee_rts.m'
    raw = get(url)
    blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
    if blob != '862c195c98e319e858829ebc0dea0f928694b71d':
        raise RuntimeError(f'Upstream case changed: {blob}; review required, not silently accepted.')
    text = raw.decode()
    case = {'version': '2', 'baseMVA': 100.0, 'source_url': url, 'source_git_blob': blob,
            'source_sha256': hashlib.sha256(raw).hexdigest()}
    for name in ['bus','gen','branch','gencost']:
        block = re.search(r'mpc\.' + name + r'\s*=\s*\[(.*?)\];', text, re.S).group(1)
        clean = '\n'.join(line.split('%')[0] for line in block.splitlines())
        case[name] = [list(map(float, row.split())) for row in clean.split(';') if row.strip()]
    save(OUT / 'rts24_full.json', case)
    print(json.dumps(metadata, indent=2), flush=True)
    print('FULL_CASE', {k: len(case[k]) for k in ['bus','gen','branch','gencost']}, flush=True)


if __name__ == '__main__':
    main()
