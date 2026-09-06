"""Chronological operating check; mature LP/MILP, not a new algorithm.

The second stress episode retains its actual clock time. `window` permits
charging only during the last W hours of the INTER-EPISODE interval; no weather
samples or hours are removed. Storage may discharge during recovery. Ramping
makes the restored state (energy, generator output), not energy alone.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp, minimize
from recovery import LP, Grid, Storage, Trace, _check_inputs

EPS = 1e-5

@dataclass
class Limits:
    ramp_fraction: float = .30  # per hour times installed aggregate Pmax
    minimum_fraction: float = .0  # aggregate flexibility approximation
    integer_storage: bool = False


def solve_lp(lp: LP, c: np.ndarray, integer: list[int]):
    if not integer:
        return lp.solve(c)
    integrality = np.zeros(len(c), int)
    integrality[integer] = 1
    constraints = [LinearConstraint(lp.matrix(lp.eq), np.asarray(lp.rhs), np.asarray(lp.rhs))]
    if lp.ub:
        constraints.append(LinearConstraint(lp.matrix(lp.ub), -np.inf, np.asarray(lp.urhs)))
    out = milp(c, integrality=integrality, bounds=Bounds(lp.lo,lp.hi),
               constraints=constraints, options={'mip_rel_gap':1e-8,'time_limit':120.})
    if out.status not in (0,2):
        raise RuntimeError(f'MILP did not establish optimum/infeasibility: {out.status} {out.message}')
    return out


def block(lp: LP, grid: Grid, st: Storage, tr: Trace, limits: Limits,
          initial: np.ndarray | None = None, link: dict | None = None,
          g0: np.ndarray | None = None, allow_shed: bool = True,
          charge_mask: np.ndarray | None = None) -> dict[str,Any]:
    _check_inputs(grid,st,tr,initial)
    if (initial is None)==(link is None):
        raise ValueError('Exactly one energy initialisation is required.')
    if tr.h < 1:
        raise ValueError('A chronological block must contain at least one interval.')
    if limits.ramp_fraction <= 0 or not 0 <= limits.minimum_fraction <= 1:
        raise ValueError('Invalid generator operating limits.')
    if g0 is not None:
        g0=np.asarray(g0,float)
        if g0.shape!=(grid.n,) or not np.all(np.isfinite(g0)) or np.any(g0 < -EPS) or np.any(g0 > grid.generation+EPS):
            raise ValueError('Invalid preceding generator state.')
    h,n,k=tr.h,grid.n,st.k
    mask=np.ones(h) if charge_mask is None else np.asarray(charge_mask,float)
    if mask.shape!=(h,) or np.any((mask<0)|(mask>1)):
        raise ValueError('Invalid charging availability.')
    g=lp.var((h,n),lo=limits.minimum_fraction*tr.capacity,hi=tr.capacity)
    theta=lp.var((h,n),lo=-np.inf)
    flow=lp.var((h,grid.m),lo=-tr.rating,hi=tr.rating)
    ch=lp.var((h,k),hi=mask[:,None]*st.pc)
    dis=lp.var((h,k),hi=st.pd)
    e=lp.var((h+1,k),hi=st.emax)
    shed=lp.var((h,n),hi=tr.demand if allow_shed else 0)
    integers=[]
    if limits.integer_storage:
        z=lp.var((h,k),hi=1)
        integers=z.ravel().tolist()
        for t in range(h):
            for j in range(k):
                lp.upper({int(ch[t,j]):1,int(z[t,j]):-st.pc[j]},0)
                lp.upper({int(dis[t,j]):1,int(z[t,j]):st.pd[j]},st.pd[j])
    for j in range(k):
        row={int(e[0,j]):1.}
        if link is not None: row[int(link['energy'][-1,j])]=-1.
        lp.equal(row,0 if link is not None else initial[j])
    for t in range(h):
        lp.equal({int(theta[t,0]):1.})
        for l,(a,b,x,_,tap) in enumerate(grid.branches):
            bval=grid.base_mva/(x*tap)
            lp.equal({int(flow[t,l]):1,int(theta[t,int(a)]):-bval,int(theta[t,int(b)]):bval})
        for i in range(n):
            row={int(g[t,i]):1,int(shed[t,i]):1}
            for j in np.where(st.bus==i)[0]:
                row[int(ch[t,j])]=-1; row[int(dis[t,j])]=1
            for l,(a,b,*_) in enumerate(grid.branches):
                if int(a)==i: row[int(flow[t,l])]=row.get(int(flow[t,l]),0)-1
                if int(b)==i: row[int(flow[t,l])]=row.get(int(flow[t,l]),0)+1
            lp.equal(row,tr.demand[t,i])
            prev = int(g[t-1,i]) if t else (int(link['g'][-1,i]) if link is not None else None)
            ramp=limits.ramp_fraction*grid.generation[i]*tr.dt
            if prev is not None:
                lp.upper({int(g[t,i]):1,prev:-1},ramp)
                lp.upper({int(g[t,i]):-1,prev:1},ramp)
            elif g0 is not None:
                lp.upper({int(g[t,i]):1},g0[i]+ramp)
                lp.upper({int(g[t,i]):-1},ramp-g0[i])
        for j in range(k):
            lp.equal({int(e[t+1,j]):1,int(e[t,j]):-1,
                      int(ch[t,j]):-st.eta_c[j]*tr.dt,int(dis[t,j]):tr.dt/st.eta_d[j]})
    return dict(g=g,theta=theta,flow=flow,charge=ch,discharge=dis,energy=e,shed=shed,
                integers=integers,mask=mask)


def decoded(idx: dict, x: np.ndarray) -> dict[str,np.ndarray]:
    return {key:x[idx[key]] for key in ['g','theta','flow','charge','discharge','energy','shed']}


def physical_audit(grid: Grid, st: Storage, tr: Trace, p: dict, initial: np.ndarray,
                   limits: Limits, g0: np.ndarray | None=None, mask: np.ndarray | None=None) -> dict:
    g,f,theta,ch,dis,e,shed=[p[key] for key in ['g','flow','theta','charge','discharge','energy','shed']]
    maximum=lambda x:float(np.max(np.abs(x),initial=0))
    positive=lambda x:float(np.max(np.maximum(x,0),initial=0))
    balance=g+shed-tr.demand
    for j,b in enumerate(st.bus): balance[:,b]+=dis[:,j]-ch[:,j]
    dc=np.zeros_like(f)
    for l,(a,b,x,_,tap) in enumerate(grid.branches):
        a,b=int(a),int(b); balance[:,a]-=f[:,l]; balance[:,b]+=f[:,l]
        dc[:,l]=f[:,l]-grid.base_mva/(x*tap)*(theta[:,a]-theta[:,b])
    prev=np.vstack([g0,g[:-1]]) if g0 is not None else g[:-1]
    current=g if g0 is not None else g[1:]
    availability=np.ones(tr.h) if mask is None else mask
    bounds=[positive(-ch),positive(ch-availability[:,None]*st.pc),positive(-dis),positive(dis-st.pd),
            positive(-e),positive(e-st.emax),positive(-shed),positive(shed-tr.demand),
            positive(limits.minimum_fraction*tr.capacity-g),positive(g-tr.capacity),positive(np.abs(f)-tr.rating)]
    out=dict(balance_mw=maximum(balance),dc_mw=maximum(dc),initial_mwh=maximum(e[0]-initial),
             dynamics_mwh=maximum(e[1:]-e[:-1]-tr.dt*(ch*st.eta_c-dis/st.eta_d)),
             bounds=max(bounds),ramp_mw=positive(np.abs(current-prev)-limits.ramp_fraction*grid.generation*tr.dt),
             simultaneous_mw=float(np.minimum(ch,dis).max(initial=0)))
    if max(v for k,v in out.items() if k!='simultaneous_mw')>EPS:
        raise AssertionError(f'Chronological physical audit failed: {out}')
    if limits.integer_storage and out['simultaneous_mw']>EPS:
        raise AssertionError('Integer storage had simultaneous charging/discharging.')
    return out


def optimise(grid: Grid,st: Storage,events: list[Trace],initial: np.ndarray,
             limits: Limits | None=None,gap: Trace | None=None,window: int | None=None,
             g0: np.ndarray | None=None,fixed_gap_path: dict | None=None) -> dict:
    """Shared recovery followed by scenario-specific recourse; all event paths
    are revealed at event onset. This is not hourly causal stochastic recourse.
    Objective is worst-case ENS, then minimum throughput on its optimal face.
    """
    limits=Limits() if limits is None else limits
    if not events: raise ValueError('At least one event is required.')
    if gap is not None and (not isinstance(window,(int,np.integer)) or not 0<=window<=gap.h):
        raise ValueError('Charging window must lie within the unchanged gap.')
    lp=LP(); gapidx=None; integers=[]
    if gap is not None:
        mask=(np.arange(gap.h)>=gap.h-window).astype(float)
        gapidx=block(lp,grid,st,gap,limits,initial=initial,g0=g0,allow_shed=False,charge_mask=mask)
        integers+=gapidx['integers']
        if fixed_gap_path is not None:
            for key in ['charge','discharge','g']:
                for idx,value in zip(gapidx[key].ravel(),fixed_gap_path[key].ravel()):
                    lp.equal({int(idx):1},value)
    eventidx=[]
    for ev in events:
        idx=block(lp,grid,st,ev,limits,initial=initial if gap is None else None,
                  link=gapidx,g0=g0 if gap is None else None,allow_shed=True)
        eventidx.append(idx); integers+=idx['integers']
    z=int(lp.var((1,))[0])
    for ev,idx in zip(events,eventidx):
        row={int(j):ev.dt for j in idx['shed'].ravel()}; row[z]=-1.
        lp.upper(row,0.)
    objective=np.zeros(len(lp.lo)); objective[z]=1
    primary=solve_lp(lp,objective,integers)
    if not primary.success: return {'feasible':False,'reason':'recovery_or_ramp_infeasible'}
    best=float(primary.fun)
    lp.upper({z:1.},best+1e-7)
    secondary=np.zeros(len(lp.lo))
    for ev,idx in zip(([gap] if gap is not None else [])+events,
                      ([gapidx] if gapidx is not None else [])+eventidx):
        secondary[idx['charge'].ravel()]=ev.dt
        secondary[idx['discharge'].ravel()]=ev.dt
    result=solve_lp(lp,secondary,integers)
    if not result.success: raise RuntimeError('Lexicographic operating solve failed.')
    gap_path=decoded(gapidx,result.x) if gap is not None else None
    ready=gap_path['energy'][-1] if gap_path is not None else initial
    pg=gap_path['g'][-1] if gap_path is not None else g0
    paths=[decoded(idx,result.x) for idx in eventidx]
    audits=[physical_audit(grid,st,ev,p,ready,limits,pg) for ev,p in zip(events,paths)]
    if gap is not None: audits.append(physical_audit(grid,st,gap,gap_path,initial,limits,g0,gapidx['mask']))
    ens=[float(p['shed'].sum()*ev.dt) for p,ev in zip(paths,events)]
    return dict(feasible=True,ready=max(ens)<=EPS,ens=max(ens),primary_ens=best,
                scenario_ens=ens,energy=ready,generator_state=pg,gap=gap_path,events=paths,
                audits=audits,throughput=float(result.fun),integer_storage=limits.integer_storage)


def minimum_window(grid,st,gap,events,initial,limits=None,g0=None) -> dict:
    """Minimum suffix charging availability, not a moved event start.
    Opening more charging hours nests the feasible sets, allowing bisection.
    All 18 gap hours, energy states, and ramp transitions remain in every solve.
    """
    cache={}
    def at(w):
        if w not in cache: cache[w]=optimise(grid,st,events,initial,limits,gap,w,g0)
        return cache[w]['feasible'] and cache[w]['ready']
    if not at(gap.h):
        return dict(window_h=None,reason='not_ready_with_full_gap',solution=cache[gap.h])
    if at(0):
        return dict(window_h=0.,solution=cache[0],checked=sorted(cache))
    lo,hi=0,gap.h
    while hi-lo>1:
        mid=(hi+lo)//2
        if at(mid): hi=mid
        else: lo=mid
    return dict(window_h=hi*gap.dt,window_steps=hi,solution=cache[hi],
                preceding_ens=cache[lo].get('ens'),checked=sorted(cache))


def balanced_recovery(grid,st,gap,initial,window,limits=None,g0=None) -> dict:
    """Network/ramp-aware myopic time-to-go QP adaptation. Same actual gap;
    no second-event information. Convex LP oracle checks QP suboptimality.
    The continuous model is accepted only after no-simultaneous audit.
    """
    limits=Limits() if limits is None else limits
    if limits.integer_storage: raise ValueError('Baseline QP is continuous; audit simultaneity instead.')
    current=initial.copy(); pg=None if g0 is None else g0.copy()
    delta=.5*float(np.min(st.pd/st.emax)); paths=[]; audits=[]; gaps=[]
    for t in range(gap.h):
        tr=gap.at(t); mask=np.array([float(t>=gap.h-window)])
        lp=LP(); idx=block(lp,grid,st,tr,limits,initial=current,g0=pg,allow_shed=False,charge_mask=mask)
        end=idx['energy'][-1]
        seed_c=np.zeros(len(lp.lo)); seed_c[end]=-1
        seed=lp.solve(seed_c)
        if not seed.success: return dict(feasible=False,reason='myopic_current_service_infeasible',failed_step=t)
        scale=max(float(st.emax.sum()),1.)
        def fun(x):
            e=x[end]
            return float((-e.sum()+delta*np.sum(e*e/st.pd)/2)/scale)
        def jac(x):
            g=np.zeros(len(x));g[end]=(-1+delta*x[end]/st.pd)/scale
            return g
        vertices=[seed.x]; weights=np.ones(1); x=seed.x.copy()
        for iteration in range(150):
            oracle=lp.solve(jac(x)); bound=max(0.,float(jac(x)@x-oracle.fun)*scale)
            if bound<=1e-4: break
            vertices.append(oracle.x); vv=np.asarray(vertices); weights=np.r_[weights,0.]
            out=minimize(lambda w:fun(w@vv),weights,jac=lambda w:vv@jac(w@vv),method='SLSQP',
                         bounds=[(0.,1.)]*len(weights),
                         constraints=[{'type':'eq','fun':lambda w:w.sum()-1.,'jac':lambda w:np.ones(len(w))}],
                         options={'maxiter':500,'ftol':1e-14})
            if not out.success: raise RuntimeError(f'QP correction failed: {out.message}')
            keep=out.x>1e-10;weights=out.x[keep];weights/=weights.sum()
            vertices=[v for v,use in zip(vertices,keep) if use];x=weights@np.asarray(vertices)
        else: raise RuntimeError('QP oracle bound did not converge.')
        p=decoded(idx,x)
        audit=physical_audit(grid,st,tr,p,current,limits,pg,mask)
        if audit['simultaneous_mw']>EPS:
            raise AssertionError('Baseline used nonphysical simultaneous storage operation.')
        paths.append(p);audits.append(audit);gaps.append(bound)
        current=p['energy'][-1];pg=p['g'][-1]
    combined={key:np.concatenate([p[key] for p in paths],axis=0) for key in paths[0] if key!='energy'}
    combined['energy']=np.vstack([initial]+[p['energy'][-1] for p in paths])
    return dict(feasible=True,energy=current,generator_state=pg,path=combined,audits=audits,qp_gaps=gaps)
