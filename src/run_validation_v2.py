"""Bounded completion suite: observed chronology, sensitivities, MILP and AC.

Run: python src/run_validation_v2.py --out results_v2 --ac
A failed hypothesis/AC search is a recorded scientific outcome, not a hidden
exception. Broken invariants and solver failures DO raise and fail the workflow.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import time
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from cases import ROOT, rts24
from recovery import recover_cuts
from chronology import EPS, Limits, optimise, minimum_window, balanced_recovery
from weather_cases import episode_case
from run_experiments import plain, save


def prepare(settings,limits=None):
    c=episode_case(**settings);g,s=c['grid'],c['storage'];limits=limits or Limits()
    before=optimise(g,s,[c['first']],c['initial'],limits)
    if not before['feasible']:raise RuntimeError('Unrestricted shedding predecessor was infeasible.')
    p=before['events'][0];e=p['energy'][-1].copy();pg=p['g'][-1].copy()
    return c,before,e,pg


def one_case(settings,limits=None,retain=False):
    start=time.perf_counter();limits=limits or Limits()
    c,before,e0,p0=prepare(settings,limits);g,s=c['grid'],c['storage']
    solution=minimum_window(g,s,c['gap'],c['events'],e0,limits,p0)
    full=optimise(g,s,c['events'],s.emax,limits)
    row=dict(**c['settings'],**vars(limits),first_ens_mwh=before['ens'],
             initial_recovery_mwh=float(e0.sum()),window_h=solution['window_h'],
             full_soc_event_ens_mwh=full['ens'],
             status='ready' if solution['window_h'] is not None else solution['reason'])
    sol=solution['solution'];baseline=None;baseline_event=None
    if sol['feasible']:
        row.update(planned_worst_ens_mwh=sol['ens'],planned_terminal_mwh=float(sol['energy'].sum()),
            max_simultaneous_mw=max(a['simultaneous_mw'] for a in sol['audits']),
            max_physical_residual=max(v for a in sol['audits'] for k,v in a.items() if k!='simultaneous_mw'))
        if row['max_simultaneous_mw']>EPS:
            raise AssertionError('Relaxed schedule must be verified physical, not silently accepted.')
    if solution['window_h'] is not None:
        window=int(round(solution['window_h']/c['gap'].dt))
        baseline=balanced_recovery(g,s,c['gap'],e0,window,replace(limits,integer_storage=False),p0)
        row['baseline_status']='feasible' if baseline['feasible'] else baseline['reason']
        if baseline['feasible']:
            baseline_event=optimise(g,s,c['events'],baseline['energy'],limits,g0=baseline['generator_state'])
            row.update(baseline_ens_mwh=baseline_event['ens'],baseline_terminal_mwh=float(baseline['energy'].sum()),
                       baseline_qp_gap=max(baseline['qp_gaps']))
        row['preceding_window_ens_mwh']=solution.get('preceding_ens')
        row['bisection_checked_windows']=solution.get('checked')
    row['runtime_s']=time.perf_counter()-start
    detail=dict(row=row,settings={k:str(v) if isinstance(v,Path) else v for k,v in settings.items()},limits=vars(limits),initial_before_first=c['initial'],
        initial_recovery=e0,initial_recovery_generation=p0,first=before,minimum_window=solution,
        baseline=baseline,baseline_event=baseline_event,
        traces={'first':vars(c['first']),'gap':vars(c['gap']),'events':[vars(v) for v in c['events']]})
    return plain(row),detail


def write_csv(path,rows):
    fields=list(dict.fromkeys(k for r in rows for k,v in r.items() if not isinstance(v,(list,dict))))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)


def ac_case(settings,limits=Limits()):
    from ac_validation import replay
    c,before,e0,p0=prepare(settings,limits);g,s=c['grid'],c['storage']
    m=minimum_window(g,s,c['gap'],c['events'],e0,limits,p0)
    if m['window_h'] is None:
        return dict(settings=settings,status='no_dc_ready_schedule',window_h=None,
                    reason=m['reason'],planned_ens=m['solution'].get('ens'))
    sol=m['solution']
    prefix=replay(g,s,[('prior_episode',c['actual_first'],before['events'][0]),
                       ('recovery',c['actual_gap'],sol['gap'])],ramp_fraction=limits.ramp_fraction,
                       minimum_fraction=limits.minimum_fraction)
    checks=[]
    for label,event,path in zip(['observed','plus_1C','plus_2C'],c['actual_events'],sol['events']):
        checks.append(replay(g,s,[(label,event,path)],prefix['terminal_pg'],
                             limits.ramp_fraction,limits.minimum_fraction))
    all_ok=prefix['all_hours_witnessed'] and all(r['all_hours_witnessed'] for r in checks)
    return dict(settings=settings,status='ac_witnessed' if all_ok else 'no_complete_ac_witness',
                window_h=m['window_h'],prefix=prefix,events=checks,
                fixed_storage_schedule=True,generation_redispatch_allowed=True,
                no_reactive_support_from_storage=True,
                interpretation='AC witness for these trajectories, not AC optimality or general guarantee.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=ROOT/'results_v2')
    p.add_argument('--ac',action='store_true')
    p.add_argument('--quick',action='store_true')
    p.add_argument('--weather-path',type=Path,default=None,help='Explicit development subset, never a silent fallback.')
    args=p.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
    started=datetime.now(timezone.utc).isoformat();rows=[];details=[]
    base={} if args.weather_path is None else {'weather_path':args.weather_path}
    days=[5] if args.quick else list(range(1,10))
    configs=[dict(day=d,line_scale=l,charger_mw=r,**base) for d in days
             for l in [.9,.95,1.] for r in [100.,200.]]
    # Grid includes null and non-restorable regimes. Pilot development used day 5;
    # remaining dates test transfer within the same observed heatwave, not independent climate events.
    for settings in configs:
        row,detail=one_case(settings);rows.append(row);details.append(plain(detail))
        print('WEATHER',json.dumps(row),flush=True)
        save(out/'weather_summary.json',rows);write_csv(out/'weather_summary.csv',rows)
    with gzip.open(out/'all_weather_trajectories.json.gz','wt',encoding='utf-8') as f:
        json.dump(details,f,allow_nan=False,separators=(',',':'))
    checks=[]
    representative=dict(day=5,line_scale=1.,charger_mw=200.,**base)
    c,before,e0,p0=prepare(representative);g,s=c['grid'],c['storage']
    optimal=minimum_window(g,s,c['gap'],c['events'],e0,g0=p0)
    assert optimal['window_h'] is not None
    w=int(optimal['window_h'])
    for ww in [w-1,w]:
        continuous=optimise(g,s,c['events'],e0,Limits(),c['gap'],ww,p0)
        integer=optimise(g,s,c['events'],e0,Limits(integer_storage=True),c['gap'],ww,p0)
        delta=abs(integer['primary_ens']-continuous['primary_ens'])
        if delta>1e-4:raise AssertionError('Storage relaxation mismatch; reconsider conclusions.')
        checks.append(dict(window_h=ww,lp_ens=continuous['primary_ens'],milp_ens=integer['primary_ens'],difference=delta))
    save(out/'lp_milp_boundary.json',checks)
    baseline_windows=[]
    for ww in range(c['gap'].h+1):
        b=balanced_recovery(g,s,c['gap'],e0,ww,g0=p0)
        row={'window_h':ww,'baseline_feasible':b['feasible']}
        if b['feasible']:
            q=optimise(g,s,c['events'],b['energy'],g0=b['generator_state'])
            row.update(baseline_ens_mwh=q['ens'],baseline_total_mwh=float(b['energy'].sum()))
        # Independent exhaustive boundary check at every unchanged clock horizon.
        q=optimise(g,s,c['events'],e0,gap=c['gap'],window=ww,g0=p0)
        row.update(proposed_feasible=q['feasible'],proposed_ens_mwh=q.get('ens'))
        baseline_windows.append(row)
    first=next(r['window_h'] for r in baseline_windows if r['proposed_feasible'] and r['proposed_ens_mwh']<=EPS)
    if first!=w:raise AssertionError('Bisection and exhaustive chronological boundary differ.')
    save(out/'all_windows_day5.json',baseline_windows);write_csv(out/'all_windows_day5.csv',baseline_windows)
    sens=[]
    for key,values in [('eta',[.9,.98]),('mapping_shift',[1,2]),('demand_slope',[.009,.015]),
                       ('loss_margin',[.02,.04])]:
        for value in values:
            row,_=one_case({**representative,key:value});row.update(sensitivity=key,value=value);sens.append(row)
    for ramp in [.1,1.]:
        row,_=one_case(representative,Limits(ramp_fraction=ramp));row.update(sensitivity='ramp_fraction',value=ramp);sens.append(row)
    null,_=one_case({**representative,'demand_slope':0.,'generation_slope':0.,'rating_slope':0.})
    null.update(sensitivity='no_temperature_effect',value=0.);sens.append(null)
    save(out/'sensitivities.json',sens);write_csv(out/'sensitivities.csv',sens)
    # Forecast-error envelope: nominal/robust share the SAME initial recovery
    # state and fixed actual 18 h interval. Errors affect only the future episode.
    forecast=[]
    for design,label in [(c['events'][:1],'nominal'),(c['events'],'0_to_2C_envelope')]:
        m=minimum_window(g,s,c['gap'],design,e0,g0=p0)
        if m['window_h'] is None:continue
        q=m['solution']
        for delta in [-1.,0.,.5,1.,1.5,2.,2.5,3.]:
            tr=c['trace_builder'](c['next_start'],6,delta)
            evaluation=optimise(g,s,[tr],q['energy'],g0=q['generator_state'])
            forecast.append(dict(design=label,window_h=m['window_h'],temperature_error_c=delta,
                                 ens_mwh=evaluation['ens']))
    save(out/'forecast_envelope.json',forecast);write_csv(out/'forecast_envelope.csv',forecast)
    ac=[]
    if args.ac:
        # Declared margin grid, including the unbuffered model. No failed row is removed.
        for day in ([5] if args.quick else [4,5]):
            for margin,linemargin in [(0.,1.),(.02,.98),(.04,.95)]:
                result=ac_case(dict(day=day,line_scale=1.,charger_mw=200.,
                                   loss_margin=margin,line_margin=linemargin,**base))
                ac.append(plain(result));save(out/'ac_replay.json',ac)
                print('AC',day,margin,linemargin,result['status'],flush=True)
        # Back-check previous headline static example, without repairing its schedule.
        from ac_validation import replay,legacy_path
        g1,s1,gap1,ev1,e1=rts24(line_scale=.8)
        r=recover_cuts(g1,s1,gap1.prefix(8),ev1,e1)
        prefix=replay(g1,s1,[('legacy_recovery',gap1.prefix(8),legacy_path(r['path'],'charge'))],ramp_fraction=1.)
        events=[replay(g1,s1,[(v.name,v,legacy_path(q['path'],'discharge'))],
                       prefix['terminal_pg'],ramp_fraction=1.) for v,q in zip(ev1,r['event_checks'])]
        save(out/'legacy_ac_check.json',dict(prefix=prefix,events=events))
    try:sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    except (subprocess.CalledProcessError,FileNotFoundError):sha='local-uncommitted-execution'
    hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
            for folder in ['src','tests','data','tools'] for f in sorted((ROOT/folder).rglob('*'))
            if f.is_file() and '__pycache__' not in str(f) and f.suffix in ['.py','.json','.csv','.txt']}
    counts=dict(total=len(rows),ready=sum(r['status']=='ready' for r in rows),
                baseline_positive_ens_at_minimum_window=sum(r.get('baseline_ens_mwh',0)>EPS for r in rows),
                no_ready_schedule=sum(r['status']!='ready' for r in rows))
    metadata=dict(started_utc=started,completed_utc=datetime.now(timezone.utc).isoformat(),
        source_commit=sha,python=platform.python_version(),
        dependencies={k:importlib.metadata.version(k) for k in ['numpy','scipy']+(['PYPOWER'] if args.ac else [])},
        input_sha256=hashes,mode='development' if args.quick or args.weather_path else 'full',
        experiment_counts=counts,ac_executed=args.ac,
        interpretation='Observed weather, synthetic electrical mapping/response; finite-scenario two-stage information.')
    save(out/'run_metadata.json',metadata)
    print('COMPLETION',json.dumps(counts),flush=True)


if __name__=='__main__':main()
