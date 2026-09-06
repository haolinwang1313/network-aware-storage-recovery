"""AC witness search with fixed storage schedules and independent circuit audit.

Uses PYPOWER AC OPF only as a feasibility-witness finder. Thermal generators
may redispatch; their original reactive bounds and bus voltage limits remain.
Storage injects no reactive power. Subsequent successful hours enforce per-unit
ramps from the preceding AC witness. A failed solve is NOT proof of infeasibility.
"""
from __future__ import annotations
import copy
import json
import warnings
from pathlib import Path
import numpy as np
from cases import ROOT
from recovery import Grid, Storage, Trace


def load_full_case(path: Path | None=None):
    p=ROOT/'data/public/rts24_full.json' if path is None else path
    raw=json.loads(p.read_text())
    return {k:np.asarray(raw[k],float) if k in ['bus','gen','branch','gencost'] else raw[k]
            for k in ['version','baseMVA','bus','gen','branch','gencost']}


def circuit_audit(case: dict, solution: dict) -> dict[str,float]:
    """Build Ybus and branch currents directly; do not reuse solver matrices."""
    base=case['baseMVA'];bus=case['bus'];branch=case['branch'];n=len(bus)
    ybus=np.zeros((n,n),complex); yseries=[]
    for l,line in enumerate(branch):
        a,b=int(line[0])-1,int(line[1])-1
        y=complex(1)/(line[2]+1j*line[3])*line[10]
        bc=1j*line[4]*line[10]/2
        tap=(line[8] if line[8] else 1.)*np.exp(1j*np.deg2rad(line[9]))
        yff=(y+bc)/(tap*np.conj(tap)); yft=-y/np.conj(tap)
        ytf=-y/tap; ytt=y+bc
        ybus[a,a]+=yff;ybus[a,b]+=yft;ybus[b,a]+=ytf;ybus[b,b]+=ytt
        yseries.append((a,b,yff,yft,ytf,ytt))
    ybus[np.arange(n),np.arange(n)]+=(bus[:,4]+1j*bus[:,5])/base
    voltage=solution['bus'][:,7]*np.exp(1j*np.deg2rad(solution['bus'][:,8]))
    snet=voltage*np.conj(ybus@voltage)*base
    specified=-bus[:,2]-1j*bus[:,3]
    for unit in solution['gen']:
        specified[int(unit[0])-1]+=unit[1]+1j*unit[2]
    flow=[]
    for a,b,yff,yft,ytf,ytt in yseries:
        sf=voltage[a]*np.conj(yff*voltage[a]+yft*voltage[b])*base
        st=voltage[b]*np.conj(ytf*voltage[a]+ytt*voltage[b])*base
        flow.append(max(abs(sf),abs(st)))
    pos=lambda x:float(np.maximum(x,0).max(initial=0))
    sg=solution['gen']; cg=case['gen']; vm=abs(voltage)
    out=dict(balance_mva=float(np.max(np.abs(snet-specified))),
             voltage_pu=max(pos(vm-bus[:,11]),pos(bus[:,12]-vm)),
             thermal_mva=pos(np.asarray(flow)-branch[:,5]),
             active_mw=max(pos(sg[:,1]-cg[:,8]),pos(cg[:,9]-sg[:,1])),
             reactive_mvar=max(pos(sg[:,2]-cg[:,3]),pos(cg[:,4]-sg[:,2])),
             min_voltage_pu=float(vm.min()),max_voltage_pu=float(vm.max()),
             max_line_loading=float(np.max(np.asarray(flow)/branch[:,5])))
    return out


def check_hour(grid: Grid, st: Storage, trace: Trace, schedule: dict, t: int,
               previous: np.ndarray | None=None,ramp_fraction=.30,minimum_fraction=0.):
    from pypower.api import runopf, ppoption
    original=load_full_case();case=copy.deepcopy(original)
    gen=case['gen'];bus=case['bus'];branch=case['branch']
    served=np.maximum(0,trace.demand[t]-schedule['shed'][t])
    original_pd=original['bus'][:,2]
    ratio=np.divide(served,original_pd,out=np.zeros_like(served),where=original_pd>0)
    bus[:,2]=served;bus[:,3]=original['bus'][:,3]*ratio
    for j,b in enumerate(st.bus):bus[b,2]+=schedule['charge'][t,j]-schedule['discharge'][t,j]
    for j in range(len(gen)):
        b=int(gen[j,0])-1
        weight=original['gen'][j,8]/grid.generation[b] if grid.generation[b]>0 else 0.
        gen[j,8]=trace.capacity[t,b]*weight
        gen[j,9]=minimum_fraction*gen[j,8]
        if previous is not None:
            ramp=ramp_fraction*original['gen'][j,8]*trace.dt
            gen[j,8]=min(gen[j,8],previous[j]+ramp)
            gen[j,9]=max(gen[j,9],previous[j]-ramp)
        if gen[j,9]>gen[j,8]+1e-8:
            return dict(witness=False,reason='ramp_bounds_inconsistent'),None
        gen[j,1]=np.clip(schedule['g'][t,b]*weight,gen[j,9],gen[j,8])
    branch[:,5]=trace.rating[t]
    branch[:,6:8]=trace.rating[t][:,None]
    # Convex generation cost is only a numerical witness-search objective.
    # Economic costs from this feasibility check are not used in comparisons.
    case['gencost']=np.tile([2.,0.,0.,3.,.01,1.,0.],(len(gen),1))
    opt=ppoption(VERBOSE=0,OUT_ALL=0,OPF_ALG=560,PDIPM_MAX_IT=150,
                 OPF_VIOLATION=1e-7,PDIPM_FEASTOL=1e-8)
    records=[]
    for attempt in range(2):
        trial=copy.deepcopy(case)
        if attempt:
            trial['bus'][:,7]=1.;trial['bus'][:,8]=0.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                solution=runopf(trial,opt)
            if not bool(solution['success']):
                records.append('optimizer_not_converged');continue
            a=circuit_audit(case,solution)
            valid=(a['balance_mva']<=1e-4 and a['voltage_pu']<=1e-6 and
                   max(a['thermal_mva'],a['active_mw'],a['reactive_mvar'])<=1e-4)
            if valid:
                return dict(witness=True,attempts=attempt+1,audit=a,
                            max_generator_change_from_dc_mw=float(np.max(np.abs(
                                np.bincount(solution['gen'][:,0].astype(int)-1,
                                            weights=solution['gen'][:,1],minlength=grid.n)-schedule['g'][t])))),solution['gen'][:,1].copy()
            records.append({'audit_failure':a})
        except Exception as exc:
            records.append(type(exc).__name__+': '+str(exc))
    return dict(witness=False,reason='no_ac_witness_found',attempts=records),None


def replay(grid,st,segments,previous=None,ramp_fraction=.30,minimum_fraction=0.):
    rows=[];chain=True;pg=previous
    for label,tr,path in segments:
        for t in range(tr.h):
            outcome,next_pg=check_hour(grid,st,tr,path,t,pg,ramp_fraction,minimum_fraction)
            outcome.update(segment=label,step=t,linked_to_preceding_ac_witness=pg is not None)
            rows.append(outcome)
            if not outcome['witness']:chain=False
            pg=next_pg
    return dict(hours=len(rows),witness_hours=sum(r['witness'] for r in rows),
                all_hours_witnessed=chain,rows=rows,terminal_pg=pg)


def legacy_path(p,mode):
    z=np.zeros_like(p['power'])
    return dict(g=p['g'],charge=p['power'] if mode=='charge' else z,
                discharge=p['power'] if mode=='discharge' else z,
                shed=p['shed'],energy=p['energy'])
