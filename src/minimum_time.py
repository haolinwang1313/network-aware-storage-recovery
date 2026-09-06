"""Exact continuous-time recovery LP for STATIONARY recovery conditions.

The time-scaling/perspective transformation is standard. Its applicability here
rests on convex instantaneous DC constraints, charge-only operation, and no
intertemporal generation constraints or storage self-discharge. A constant
average dispatch realizes every feasible integrated solution. Nonstationary
traces must use the multi-period solver in recovery.py instead.
"""
from __future__ import annotations
import numpy as np
from recovery import (LP, TOL, Grid, Storage, Trace, _check_inputs,
                      add_stage, decode, audit, recourse)


def scaled_master(grid: Grid, st: Storage, gap: Trace, initial: np.ndarray, cuts: list):
    _check_inputs(grid,st,gap,initial)
    if gap.h < 1:
        raise ValueError('Need at least one stationary recovery sample.')
    for values in (gap.demand,gap.capacity,gap.rating):
        if not np.allclose(values,values[0],rtol=0,atol=1e-10):
            raise ValueError('Continuous-time LP requires stationary recovery conditions.')
    lp=LP()
    tau=lp.var((1,))[0]
    energy=lp.var((st.k,),hi=st.emax)
    charge=lp.var((st.k,))
    generation=lp.var((grid.n,))
    theta=lp.var((grid.n,),lo=-np.inf)
    flow=lp.var((grid.m,),lo=-np.inf)
    lp.equal({int(theta[0]):1})
    for j in range(st.k):
        lp.equal({int(energy[j]):1,int(charge[j]):-st.eta_c[j]},initial[j])
        lp.upper({int(charge[j]):1,int(tau):-st.pc[j]},0)
    for i in range(grid.n):
        lp.upper({int(generation[i]):1,int(tau):-gap.capacity[0,i]},0)
        row={int(generation[i]):1,int(tau):-gap.demand[0,i]}
        for j in np.where(st.bus==i)[0]: row[int(charge[j])]=-1
        for l,(a,b,*_) in enumerate(grid.branches):
            if int(a)==i: row[int(flow[l])]=row.get(int(flow[l]),0)-1
            if int(b)==i: row[int(flow[l])]=row.get(int(flow[l]),0)+1
        lp.equal(row)
    for l,(a,b,x,_,tap) in enumerate(grid.branches):
        susceptance=grid.base_mva/(x*tap)
        lp.equal({int(flow[l]):1,int(theta[int(a)]):-susceptance,
                  int(theta[int(b)]):susceptance})
        lp.upper({int(flow[l]):1,int(tau):-gap.rating[0,l]},0)
        lp.upper({int(flow[l]):-1,int(tau):-gap.rating[0,l]},0)
    for a,b in cuts:
        lp.upper({int(energy[j]):-float(a[j]) for j in range(st.k)},-b)
    index=dict(tau=int(tau),energy=energy,charge=charge,g=generation,theta=theta,flow=flow)
    return lp,index


def audit_scaled(grid,st,gap,initial,index,x):
    tau=float(x[index['tau']]); energy=x[index['energy']]
    if tau<TOL:
        if np.max(np.abs(energy-initial))>TOL:
            raise AssertionError('Energy changed at zero recovery time.')
        return dict(zero_time_energy_residual=float(np.max(np.abs(energy-initial))))
    trace=Trace('constant-average-recovery',gap.demand[:1],gap.capacity[:1],gap.rating[:1],tau)
    path=dict(g=x[index['g']][None,:]/tau,theta=x[index['theta']][None,:]/tau,
              flow=x[index['flow']][None,:]/tau,power=x[index['charge']][None,:]/tau,
              energy=np.vstack([initial,energy]),shed=np.zeros((1,grid.n)))
    return audit(grid,st,trace,path,'charge',initial)


def minimum_time(grid: Grid, st: Storage, gap: Trace, events: list[Trace],
                 initial: np.ndarray, method='monolithic') -> dict:
    if method not in ('monolithic','cuts') or not events:
        raise ValueError('Choose a method and at least one event.')
    cuts=[]
    for iteration in range(200):
        lp,index=scaled_master(grid,st,gap,initial,cuts)
        event_indices=[]
        if method=='monolithic':
            event_indices=[add_stage(lp,grid,st,v,None,'discharge',link=index['energy'])
                           for v in events]
        objective=np.zeros(len(lp.lo)); objective[index['tau']]=1
        sol=lp.solve(objective)
        if not sol.success:
            return dict(feasible=False)
        energy=sol.x[index['energy']]
        checks=[recourse(grid,st,v,energy) for v in events]
        if max(r['ens'] for r in checks)<=TOL:
            return dict(feasible=True,minimum_h=float(sol.fun),energy_mwh=energy,
                        charge_input_mwh=float(sol.x[index['charge']].sum()),
                        average_charge_mw=sol.x[index['charge']]/sol.fun if sol.fun>TOL else np.zeros(st.k),
                        iterations=iteration+1,cuts=len(cuts),
                        audit=audit_scaled(grid,st,gap,initial,index,sol.x),
                        max_event_audit=max(v for r in checks for v in r['audit'].values()),
                        scenario_ens_mwh=[r['ens'] for r in checks])
        if method=='monolithic':
            raise AssertionError('No-shedding monolithic result failed independent recourse.')
        for r in checks:
            if r['ens']>TOL:
                a=-r['gradient']; cuts.append((a,float(r['ens']+a@energy)))
    raise RuntimeError('Continuous-time cut iteration limit exceeded.')
