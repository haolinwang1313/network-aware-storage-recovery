"""Declared follow-up AC checks after the unbuffered schedules failed.

Retains every tested point. Dates 3, 6 and 7 extend, rather than replace, the
initial dates 4 and 5. Margins are explicit planning assumptions, not calibrated
AC feasibility guarantees. The numerical source modules remain unchanged.
"""
from __future__ import annotations
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from cases import ROOT
from run_validation_v2 import ac_case
from run_experiments import plain, save


def main():
    out=ROOT/'results_v2';out.mkdir(exist_ok=True)
    rows=[];full=[]
    for day in [3,6,7]:
        for loss,line in [(0.,1.),(.03,.98),(.05,.95)]:
            settings=dict(day=day,line_scale=1.,charger_mw=200.,loss_margin=loss,line_margin=line)
            r=ac_case(settings);full.append(plain(r))
            segments=[r['prefix']]+r['events'] if 'prefix' in r else []
            points=[p for s in segments for p in s['rows'] if p['witness']]
            row=dict(**settings,status=r['status'],window_h=r.get('window_h'),
                intervals=sum(s['hours'] for s in segments),witnessed=sum(s['witness_hours'] for s in segments),
                max_balance_mva=max((p['audit']['balance_mva'] for p in points),default=None))
            rows.append(row);save(out/'ac_followup.json',full)
            with (out/'ac_followup_summary.csv').open('w',newline='') as stream:
                w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            print(json.dumps(row),flush=True)
    save(out/'ac_followup_metadata.json',dict(completed_utc=datetime.now(timezone.utc).isoformat(),
        rationale='Follow-up after initial AC non-convergence; no failures removed.',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted((ROOT/'src').glob('*.py'))},designs=len(rows)))


if __name__=='__main__':main()
