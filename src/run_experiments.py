"""Run the complete declared grid of synthetic experiments; retain null/fail cases."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import platform
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy
from cases import ROOT, two_bus, rts24
from minimum_time import minimum_time
from recovery import (TOL, recourse, first_ready, recover_monolithic,
                      full_recovery_time, time_to_go, _master)


def plain(value):
    if isinstance(value,np.ndarray): return value.tolist()
    if isinstance(value,np.generic): return value.item()
    if isinstance(value,dict): return {k:plain(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [plain(v) for v in value]
    return value


def save(path, value):
    path.write_text(json.dumps(plain(value),indent=2,allow_nan=False)+'\n',encoding='utf-8')


def summarize(factory, **settings):
    start=time.perf_counter()
    grid,st,gap,events,initial=factory(**settings)
    result=dict(settings=settings,initial_energy_mwh=initial.tolist(),
                event_names=[v.name for v in events])
    # This zero-charge feasibility problem confirms that current service can
    # already be restored even when the future-event feasible set is unreachable.
    lp,_,cost=_master(grid,st,gap.prefix(1),initial,[])
    result['current_service_feasible']=bool(lp.solve(cost).success)
    result['full_state_worst_ens_mwh']=max(recourse(grid,st,v,st.emax)['ens'] for v in events)
    opt=first_ready(grid,st,gap,events,initial)
    result['ready_h']=opt['ready_h']
    result['status']=opt.get('reason','event_ready')
    if opt['ready_h'] is None:
        result['baseline_ready_h']=None
        result['runtime_s']=time.perf_counter()-start
        return result,None
    sol=opt['solution']; h=opt['steps']
    direct=recover_monolithic(grid,st,gap.prefix(h),events,initial)
    previous=recover_monolithic(grid,st,gap.prefix(h-1),events,initial) if h else None
    if not direct['feasible'] or (previous and previous['feasible']):
        raise AssertionError('Direct chronology failed first-ready boundary check.')
    if abs(direct['charge_mwh']-sol['charge_mwh'])>1e-5:
        raise AssertionError('Cut and monolithic optimum mismatch.')
    continuous=minimum_time(grid,st,gap,events,initial)
    separated=minimum_time(grid,st,gap,events,initial,method='cuts')
    if not continuous['feasible'] or not separated['feasible']:
        raise AssertionError('Continuous-time feasibility mismatch.')
    if abs(continuous['minimum_h']-separated['minimum_h'])>1e-5:
        raise AssertionError('Continuous-time formulations disagree.')
    if np.ceil(continuous['minimum_h']/gap.dt-1e-7)*gap.dt!=opt['ready_h']:
        raise AssertionError('Continuous-time and discrete horizon disagree.')
    baseline=time_to_go(grid,st,gap,initial)
    baseline_ens=[max(recourse(grid,st,v,e)['ens'] for v in events) for e in baseline['energies']]
    bh=next((j for j,q in enumerate(baseline_ens) if q<=TOL),None)
    result.update(continuous_ready_h=continuous['minimum_h'],
                  continuous_separated_h=separated['minimum_h'],
                  baseline_ready_h=None if bh is None else bh*gap.dt,
                  baseline_ens_at_ready_mwh=baseline_ens[h],
                  planned_energy_mwh=sol['energy'],
                  baseline_energy_at_ready_mwh=baseline['energies'][h],
                  charge_input_mwh=sol['charge_mwh'],
                  monolithic_charge_input_mwh=direct['charge_mwh'],
                  full_soc_recovery_h=full_recovery_time(grid,st,gap,initial,st.emax),
                  cut_count=sol['cuts'],cut_iterations=sol['iterations'],
                  max_qp_optimality_gap=max(baseline['optimality_gaps'],default=0),
                  max_physical_residual=max([*sol['audit'].values(),
                    *[val for r in sol['event_checks'] for val in r['audit'].values()]]))
    detail=dict(result=result,grid=vars(grid),storage=vars(st),gap=vars(gap),
                events=[vars(v) for v in events],recovery_solution=sol,
                continuous=continuous,continuous_separated=separated,
                direct_recovery_solution=direct,baseline=baseline,
                baseline_worst_ens_mwh=baseline_ens,scanned_horizons=opt['scanned'])
    result['runtime_s']=time.perf_counter()-start
    return result,detail


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'results')
    parser.add_argument('--quick',action='store_true',help='Run toy plus two RTS cases, not full sweep.')
    args=parser.parse_args(); out=args.out; out.mkdir(parents=True,exist_ok=True)
    started=datetime.now(timezone.utc).isoformat()
    rows=[]
    toy,detail=summarize(two_bus)
    toy['case']='two_bus'; rows.append(plain(toy)); save(out/'two_bus_detail.json',detail)
    grid,st,gap,events,e0=two_bus()
    states=[np.array([6.,0.]),np.array([0.,6.]),np.array([8.,0.])]
    save(out/'equal_soc_counterexample.json',[
        dict(energy_mwh=e,total_mwh=e.sum(),ens_mwh=recourse(grid,st,events[0],e)['ens']) for e in states])
    print('toy',toy['ready_h'],toy['baseline_ready_h'],flush=True)
    pairs=[(1.,1.),(1.,.8)] if args.quick else [
        (severity,line) for severity in [.6,.8,1.,1.2] for line in [.55,.65,.8,1.,1.25]]
    for severity,line in pairs:
        row,detail=summarize(rts24,severity=severity,line_scale=line)
        row['case']='rts24'; rows.append(plain(row))
        save(out/'summary.json',rows)  # incremental, including all null/failure cases
        print(f'severity={severity} line={line} ready={row["ready_h"]} '
              f'baseline={row["baseline_ready_h"]} status={row["status"]}',flush=True)
        if (severity,line) in [(1.,1.),(1.,.8)]:
            save(out/f'rts24_line_{line:.2f}_detail.json',detail)
        if (severity,line)==(1.,.8) and detail:
            planned=np.asarray(row['planned_energy_mwh'])
            sensitivity=[]
            for observed in [1.,1.05,1.1]:
                g,s,_,ev,_=rts24(severity=observed,line_scale=.8)
                sensitivity.append(dict(actual_severity=observed,
                    worst_ens_mwh=max(recourse(g,s,v,planned)['ens'] for v in ev)))
            save(out/'forecast_miss_sensitivity.json',sensitivity)
    fields=['case','severity','line_scale','status','ready_h','baseline_ready_h',
            'baseline_ens_at_ready_mwh','full_state_worst_ens_mwh','full_soc_recovery_h','continuous_ready_h',
            'charge_input_mwh','monolithic_charge_input_mwh','max_physical_residual']
    with (out/'summary.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader()
        for row in rows:
            flattened={**row,**row['settings']}
            writer.writerow({k:flattened.get(k) for k in fields})
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for folder in ['src','tests','data'] for p in sorted((ROOT/folder).glob('*'))
            if p.is_file() and p.suffix in ['.py','.json']}
    save(out/'run_metadata.json',dict(started_utc=started,
        completed_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
        mode='quick' if args.quick else 'complete 4x5 grid plus analytical case',
        evidence_status='synthetic DC proof-of-concept; not publication-ready validation',
        input_sha256=hashes))
    feasible=[r for r in rows if r['case']=='rts24' and r['ready_h'] is not None]
    earlier=[r for r in feasible if r['baseline_ready_h'] is None or r['ready_h']<r['baseline_ready_h']]
    print('RTS feasible/earlier/no-advantage/full-state-fail:',len(feasible),len(earlier),
          len(feasible)-len(earlier),len(rows)-1-len(feasible),flush=True)


if __name__=='__main__': main()
