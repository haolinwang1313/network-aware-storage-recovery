"""Build the DC input table from a locally acquired MATPOWER version 2 case."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def reduce_case(full):
    buses = full['bus']
    generation = {int(row[0]): 0.0 for row in buses}
    for row in full['gen']:
        if row[7] > 0:
            generation[int(row[0])] += row[8]
    return {
        'name': 'IEEE RTS-24 DC adaptation',
        'base_mva': full['baseMVA'],
        'source': full.get('source_url', 'MATPOWER case24_ieee_rts.m'),
        'source_git_blob': full.get('source_git_blob', 'unverified'),
        'load_mw': [row[2] for row in buses],
        'generation_mw': [generation[int(row[0])] for row in buses],
        'branch_columns': ['from_bus_1based', 'to_bus_1based', 'reactance_pu',
                           'rateA_mva_used_as_mw', 'tap'],
        'branches': [[int(row[0]), int(row[1]), row[3], row[5], row[8] or 1.0]
                     for row in full['branch'] if row[10] > 0],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'data/public/rts24_full.json')
    parser.add_argument('--output', type=Path, default=ROOT / 'data/rts24.json')
    args = parser.parse_args()
    result = reduce_case(json.loads(args.input.read_text(encoding='utf-8')))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(args.output)


if __name__ == '__main__':
    main()
