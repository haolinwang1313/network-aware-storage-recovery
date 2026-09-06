"""Network-aware recovery to a forecast-event feasible set.

Research prototype: lossless DC grid; continuous curtailable generation; no
commitment/ramping; charge-only recovery and discharge-only stress windows.
Every scenario shares ONE recovery trajectory. The complete stress trace is
revealed at its start; this is NOT an intra-event nonanticipative policy.
Units: MW, MWh, hours, radians. No proprietary solver or network access required.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import numpy as np
from numpy.typing import NDArray
from scipy.optimize import linprog, minimize
from scipy.sparse import coo_matrix

Array = NDArray[np.float64]
TOL = 1e-6


@dataclass
class Grid:
    name: str
    base_mva: float
    load: Array
    generation: Array
    branches: Array  # from/to ZERO-based, x, rating, tap

    def __post_init__(self) -> None:
        self.load = np.asarray(self.load, float)
        self.generation = np.asarray(self.generation, float)
        self.branches = np.asarray(self.branches, float)
        if self.load.ndim != 1 or self.generation.shape != self.load.shape:
            raise ValueError('Grid load and generation must have equal 1-D shape.')
        if self.base_mva <= 0 or not np.isfinite(self.base_mva):
            raise ValueError('base_mva must be finite and positive.')
        if self.branches.ndim != 2 or self.branches.shape[1] != 5:
            raise ValueError('Expected branches [from,to,x,rating,tap].')
        if not all(np.all(np.isfinite(a)) for a in (self.load, self.generation, self.branches)):
            raise ValueError('Non-finite grid data.')
        if min(self.load.min(), self.generation.min()) < 0:
            raise ValueError('Negative load or generation.')
        if np.any(self.branches[:, 2:] <= 0):
            raise ValueError('Reactance, rating and tap must be positive.')
        bus = self.branches[:, :2]
        if np.any(bus != bus.astype(int)) or np.any(bus < 0) or np.any(bus >= len(self.load)):
            raise ValueError('Invalid branch bus index.')

    @property
    def n(self) -> int:
        return len(self.load)

    @property
    def m(self) -> int:
        return len(self.branches)


@dataclass
class Storage:
    bus: NDArray[np.int64]
    emax: Array
    pc: Array
    pd: Array
    eta_c: Array
    eta_d: Array

    def __post_init__(self) -> None:
        raw = np.asarray(self.bus)
        if raw.ndim != 1 or np.any(raw != raw.astype(int)):
            raise ValueError('Storage buses must be integers.')
        self.bus = raw.astype(int)
        for key in ('emax', 'pc', 'pd', 'eta_c', 'eta_d'):
            a = np.asarray(getattr(self, key), float)
            if a.shape != self.bus.shape or not np.all(np.isfinite(a)) or np.any(a <= 0):
                raise ValueError(f'Invalid storage parameter {key}.')
            setattr(self, key, a)
        if np.any(self.eta_c > 1) or np.any(self.eta_d > 1):
            raise ValueError('Efficiency cannot exceed one.')

    @property
    def k(self) -> int:
        return len(self.bus)


@dataclass
class Trace:
    name: str
    demand: Array
    capacity: Array
    rating: Array
    dt: float = 1.0

    def __post_init__(self) -> None:
        for key in ('demand', 'capacity', 'rating'):
            a = np.asarray(getattr(self, key), float)
            if a.ndim != 2 or not np.all(np.isfinite(a)) or np.any(a < 0):
                raise ValueError(f'Invalid trace {key}.')
            setattr(self, key, a)
        if self.demand.shape != self.capacity.shape or self.rating.shape[0] != self.h:
            raise ValueError('Trace dimensions differ.')
        if self.dt <= 0 or not np.isfinite(self.dt):
            raise ValueError('dt must be finite and positive.')

    @property
    def h(self) -> int:
        return len(self.demand)

    def prefix(self, h: int) -> Trace:
        if not isinstance(h, int) or h < 0 or h > self.h:
            raise ValueError('Invalid prefix length.')
        return Trace(self.name, self.demand[:h], self.capacity[:h], self.rating[:h], self.dt)

    def at(self, t: int) -> Trace:
        return Trace(self.name, self.demand[t:t+1], self.capacity[t:t+1], self.rating[t:t+1], self.dt)


class LP:
    def __init__(self) -> None:
        self.lo: list[float] = []
        self.hi: list[float] = []
        self.eq: list[dict[int, float]] = []
        self.rhs: list[float] = []
        self.ub: list[dict[int, float]] = []
        self.urhs: list[float] = []

    def var(self, shape: tuple[int, ...], lo: Any = 0, hi: Any = np.inf) -> NDArray[np.int64]:
        start = len(self.lo)
        count = int(np.prod(shape))
        self.lo.extend(np.broadcast_to(lo, shape).astype(float).ravel())
        self.hi.extend(np.broadcast_to(hi, shape).astype(float).ravel())
        return np.arange(start, start + count).reshape(shape)

    def equal(self, row: dict[int, float], rhs: float = 0.0) -> int:
        self.eq.append(row)
        self.rhs.append(float(rhs))
        return len(self.eq) - 1

    def upper(self, row: dict[int, float], rhs: float) -> None:
        self.ub.append(row)
        self.urhs.append(float(rhs))

    def matrix(self, rows: list[dict[int, float]]):
        rr, cc, vv = [], [], []
        for i, row in enumerate(rows):
            for j, val in row.items():
                rr.append(i); cc.append(j); vv.append(val)
        return coo_matrix((vv, (rr, cc)), shape=(len(rows), len(self.lo))).tocsr()

    def solve(self, cost: Array | None = None):
        c = np.zeros(len(self.lo)) if cost is None else cost
        result = linprog(c, A_ub=self.matrix(self.ub) if self.ub else None,
                         b_ub=np.asarray(self.urhs) if self.ub else None,
                         A_eq=self.matrix(self.eq), b_eq=np.asarray(self.rhs),
                         bounds=list(zip(self.lo, self.hi)), method='highs',
                         options={'primal_feasibility_tolerance': 1e-8,
                                  'dual_feasibility_tolerance': 1e-8})
        if result.status not in (0, 2):
            raise RuntimeError(f'LP failed, status={result.status}: {result.message}')
        return result


def _check_inputs(grid: Grid, st: Storage, trace: Trace, initial: Array | None) -> None:
    if trace.demand.shape[1] != grid.n or trace.rating.shape[1] != grid.m:
        raise ValueError('Grid and trace dimensions differ.')
    if np.any(st.bus < 0) or np.any(st.bus >= grid.n):
        raise ValueError('Storage bus outside grid.')
    if initial is not None:
        initial = np.asarray(initial, float)
        if initial.shape != (st.k,) or not np.all(np.isfinite(initial)):
            raise ValueError('Invalid initial energy vector.')
        if np.any(initial < -TOL) or np.any(initial > st.emax + TOL):
            raise ValueError('Initial energy outside capacity bounds.')


def add_stage(lp: LP, grid: Grid, st: Storage, trace: Trace, initial: Array | None,
              mode: str, allow_shed: bool = False, link: NDArray | None = None) -> dict:
    _check_inputs(grid, st, trace, initial)
    if mode not in ('charge', 'discharge') or (initial is None) == (link is None):
        raise ValueError('Choose a mode and exactly one initial state specification.')
    h, n, k = trace.h, grid.n, st.k
    g = lp.var((h, n), hi=trace.capacity)
    theta = lp.var((h, n), lo=-np.inf)
    f = lp.var((h, grid.m), lo=-trace.rating, hi=trace.rating)
    p = lp.var((h, k), hi=st.pc if mode == 'charge' else st.pd)
    e = lp.var((h+1, k), hi=st.emax)
    shed = lp.var((h, n), hi=trace.demand if allow_shed else 0)
    init_rows = []
    for j in range(k):
        row = {int(e[0,j]): 1.0}
        if link is not None:
            row[int(link[j])] = -1.0
        init_rows.append(lp.equal(row, 0 if link is not None else float(initial[j])))
    for t in range(h):
        lp.equal({int(theta[t,0]): 1.0})
        for line, (a,b,x,_,tap) in enumerate(grid.branches):
            a,b = int(a),int(b)
            susceptance = grid.base_mva / (x*tap)
            lp.equal({int(f[t,line]): 1.0, int(theta[t,a]): -susceptance,
                      int(theta[t,b]): susceptance})
        for i in range(n):
            row = {int(g[t,i]): 1.0, int(shed[t,i]): 1.0}
            for j in np.where(st.bus == i)[0]:
                row[int(p[t,j])] = -1.0 if mode == 'charge' else 1.0
            for line, (a,b,*_) in enumerate(grid.branches):
                if int(a) == i:
                    row[int(f[t,line])] = row.get(int(f[t,line]),0) - 1
                if int(b) == i:
                    row[int(f[t,line])] = row.get(int(f[t,line]),0) + 1
            lp.equal(row, trace.demand[t,i])
        for j in range(k):
            coefficient = -st.eta_c[j]*trace.dt if mode == 'charge' else trace.dt/st.eta_d[j]
            lp.equal({int(e[t+1,j]):1.0, int(e[t,j]):-1.0, int(p[t,j]):coefficient})
    return dict(g=g, theta=theta, flow=f, power=p, energy=e, shed=shed,
                init_rows=np.asarray(init_rows), mode=mode)


def decode(index: dict, x: Array) -> dict:
    return {key: x[index[key]] for key in ('g','theta','flow','power','energy','shed')}


def audit(grid: Grid, st: Storage, trace: Trace, path: dict, mode: str,
          initial: Array) -> dict[str, float]:
    """Recompute physical equations directly, without using the LP matrix."""
    g,th,f,p,e,shed = (path[x] for x in ('g','theta','flow','power','energy','shed'))
    balance = g + shed - trace.demand
    for j,bus in enumerate(st.bus):
        balance[:,bus] += (-1 if mode == 'charge' else 1) * p[:,j]
    dc = np.zeros_like(f)
    for l,(a,b,x,_,tap) in enumerate(grid.branches):
        a,b = int(a),int(b)
        balance[:,a] -= f[:,l]; balance[:,b] += f[:,l]
        dc[:,l] = f[:,l] - grid.base_mva/(x*tap)*(th[:,a]-th[:,b])
    inc = p*st.eta_c*trace.dt if mode == 'charge' else -p/st.eta_d*trace.dt
    maxabs = lambda a: float(np.max(np.abs(a), initial=0))
    pos = lambda a: float(np.maximum(a,0).max(initial=0))
    limits = [pos(-g),pos(g-trace.capacity),pos(np.abs(f)-trace.rating),
              pos(-p),pos(p-(st.pc if mode=='charge' else st.pd)),
              pos(-e),pos(e-st.emax),pos(-shed),pos(shed-trace.demand)]
    result = dict(balance_mw=maxabs(balance),dc_mw=maxabs(dc),
                  dynamics_mwh=maxabs(e[1:]-e[:-1]-inc),
                  initial_mwh=maxabs(e[0]-initial),bounds=max(limits),
                  reference_angle=maxabs(th[:,0]))
    if max(result.values()) > 5e-6:
        raise AssertionError(f'Independent physical audit failed: {result}')
    return result


def recourse(grid: Grid, st: Storage, trace: Trace, energy: Array) -> dict:
    lp = LP()
    idx = add_stage(lp,grid,st,trace,np.asarray(energy),'discharge',True)
    cost = np.zeros(len(lp.lo)); cost[idx['shed'].ravel()] = trace.dt
    sol = lp.solve(cost)
    if not sol.success:
        raise RuntimeError('Recourse should be feasible through load shedding.')
    path = decode(idx,sol.x)
    gradient = np.asarray(sol.eqlin.marginals[idx['init_rows']])
    return dict(ens=float(sol.fun),gradient=gradient,path=path,
                audit=audit(grid,st,trace,path,'discharge',np.asarray(energy)))


def _master(grid: Grid, st: Storage, gap: Trace, initial: Array, cuts: list,
            target: Array | None = None, aggregate: bool = False):
    lp = LP(); idx = add_stage(lp,grid,st,gap,initial,'charge')
    end = idx['energy'][-1]
    for a,b in cuts:
        lp.upper({int(end[j]):-float(a[j]) for j in range(st.k)},-b)
    if target is not None:
        for j in range(st.k):
            lp.upper({int(end[j]):-1},-float(target[j]))
    cost = np.zeros(len(lp.lo))
    cost[idx['power'].ravel()] = gap.dt
    if aggregate:
        cost[:] = 0; cost[end] = -1
    return lp,idx,cost


def recover_cuts(grid: Grid, st: Storage, gap: Trace, events: list[Trace], initial: Array,
                 max_iterations: int = 200) -> dict:
    """Standard LP/Benders feasibility separation, NOT a new optimizer."""
    if not events:
        raise ValueError('At least one event is required.')
    cuts: list[tuple[Array,float]] = []
    history = []
    for iteration in range(max_iterations):
        lp,idx,cost = _master(grid,st,gap,initial,cuts)
        sol = lp.solve(cost)
        if not sol.success:
            return dict(feasible=False,iterations=iteration+1,cuts=len(cuts),history=history)
        path = decode(idx,sol.x); energy = path['energy'][-1]
        checks = [recourse(grid,st,event,energy) for event in events]
        worst = max(r['ens'] for r in checks)
        history.append(dict(energy=energy.tolist(),worst_ens=worst))
        if worst <= TOL:
            return dict(feasible=True,energy=energy,charge_mwh=float(sol.fun),
                        iterations=iteration+1,cuts=len(cuts),history=history,
                        path=path,event_checks=checks,
                        cut_coefficients=[{'a':a.tolist(),'b':b} for a,b in cuts],
                        audit=audit(grid,st,gap,path,'charge',initial))
        for r in checks:
            if r['ens'] > TOL:
                a = -r['gradient']; b = float(r['ens'] + a@energy)
                cuts.append((a,b))
    raise RuntimeError('Cut iteration limit reached; not treated as success or infeasibility.')


def recover_monolithic(grid: Grid, st: Storage, gap: Trace, events: list[Trace],
                       initial: Array) -> dict:
    if not events:
        raise ValueError('At least one event is required.')
    lp = LP(); recovery = add_stage(lp,grid,st,gap,initial,'charge')
    end = recovery['energy'][-1]
    idxs = [add_stage(lp,grid,st,event,None,'discharge',link=end) for event in events]
    cost = np.zeros(len(lp.lo)); cost[recovery['power'].ravel()] = gap.dt
    sol = lp.solve(cost)
    if not sol.success:
        return dict(feasible=False)
    path = decode(recovery,sol.x)
    checks = [audit(grid,st,event,decode(idx,sol.x),'discharge',path['energy'][-1])
              for event,idx in zip(events,idxs)]
    return dict(feasible=True,energy=path['energy'][-1],charge_mwh=float(sol.fun),
                audit=audit(grid,st,gap,path,'charge',initial),event_audits=checks,path=path)


def first_ready(grid: Grid, st: Storage, gap: Trace, events: list[Trace],
                initial: Array, method: str = 'cuts') -> dict:
    if method not in ('cuts','monolithic'):
        raise ValueError('Unknown method.')
    fn = recover_cuts if method=='cuts' else recover_monolithic
    if not events:
        raise ValueError('At least one event is required.')
    full_ens = max(recourse(grid,st,event,st.emax)['ens'] for event in events)
    if full_ens > TOL:
        return dict(ready_h=None,steps=None,scanned=[],reason='full_state_not_event_feasible',
                    full_state_worst_ens=full_ens)
    scanned = []
    for h in range(gap.h+1):
        result = fn(grid,st,gap.prefix(h),events,initial)
        scanned.append(dict(h=h,feasible=result['feasible']))
        if result['feasible']:
            return dict(ready_h=h*gap.dt,steps=h,solution=result,scanned=scanned)
    return dict(ready_h=None,steps=None,scanned=scanned,reason='not_reached_within_search_horizon')


def full_recovery_time(grid: Grid, st: Storage, gap: Trace, initial: Array,
                       target: Array) -> float | None:
    target = np.asarray(target,float)
    _check_inputs(grid,st,gap,target)
    for h in range(gap.h+1):
        lp,idx,cost = _master(grid,st,gap.prefix(h),initial,[],target)
        sol = lp.solve(cost)
        if sol.success:
            audit(grid,st,gap.prefix(h),decode(idx,sol.x),'charge',initial)
            return h*gap.dt
    return None


def time_to_go(grid: Grid, st: Storage, gap: Trace, initial: Array) -> dict:
    """Network-constrained one-step QP value-function adaptation of Evans &
    Tindemans (2020), Sec. III, Eqs. (5)-(7). Uses this prototype's charge-only
    phase and efficiencies; not a full reproduction of their risk-sizing study.
    Value = sum(e - delta*e^2/(2*pd)); delta keeps all marginal values positive.
    """
    delta = 0.5*float(np.min(st.pd/st.emax))
    current = np.asarray(initial,float).copy()
    energies = [current.copy()]; audits = []; objective_gaps = []
    for t in range(gap.h):
        trace = gap.at(t); lp,idx,_ = _master(grid,st,trace,current,[])
        end = idx['energy'][-1]
        seed_cost = np.zeros(len(lp.lo)); seed_cost[end] = -1
        seed = lp.solve(seed_cost)
        if not seed.success:
            raise RuntimeError('Current load cannot be served without discharging storage.')
        scale = max(float(np.sum(st.emax)),1.0)
        def fun(x):
            e = x[end]
            return float((-e.sum()+delta*np.sum(e*e/st.pd)/2)/scale)
        def jac(x):
            grad = np.zeros(len(x)); grad[end] = (-1+delta*x[end]/st.pd)/scale
            return grad
        # Fully-corrective conditional-gradient solution of this convex QP.
        # Each vertex comes from a feasible LP; correction optimizes only
        # convex-combination weights. Avoids redundant DC equalities in QP
        # solvers. The LP oracle supplies a genuine first-order error bound.
        vertices = [seed.x.copy()]
        weights = np.array([1.0])
        x = seed.x.copy()
        for inner in range(100):
            oracle = lp.solve(jac(x))
            gap_bound = max(0.0,float(jac(x)@x-oracle.fun)*scale)
            if gap_bound <= 5e-5:
                break
            vertices.append(oracle.x.copy())
            vv = np.asarray(vertices)
            weights = np.r_[weights,0.0]
            correction = minimize(lambda w: fun(w@vv),weights,
                                  jac=lambda w: vv@jac(w@vv),method='SLSQP',
                                  bounds=[(0.,1.)]*len(weights),
                                  constraints=[{'type':'eq','fun':lambda w:w.sum()-1,
                                                'jac':lambda w:np.ones(len(w))}],
                                  options={'maxiter':500,'ftol':1e-14})
            if not correction.success:
                raise RuntimeError(f'QP correction failed: {correction.message}')
            weights = correction.x
            keep = weights > 1e-10
            weights = weights[keep]; weights /= weights.sum()
            vertices = [v for v,take in zip(vertices,keep) if take]
            x = weights@np.asarray(vertices)
        else:
            raise RuntimeError(f'QP optimality bound failed: {gap_bound}')
        path = decode(idx,x)
        audits.append(audit(grid,st,trace,path,'charge',current))
        objective_gaps.append(gap_bound)
        current = path['energy'][-1].copy(); energies.append(current)
    return dict(energies=np.asarray(energies),audits=audits,optimality_gaps=objective_gaps)
