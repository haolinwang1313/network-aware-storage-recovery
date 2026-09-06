"""Deterministic, explicitly synthetic engineering stress tests."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from recovery import Grid, Storage, Trace

ROOT = Path(__file__).resolve().parents[1]


def constant_trace(name, grid, hours, demand, capacity, ratings, dt=1.0):
    return Trace(name,np.tile(demand,(hours,1)),np.tile(capacity,(hours,1)),
                 np.tile(ratings,(hours,1)),dt)


def two_bus(eta=1.0, event_rating=2.0, gap_hours=14):
    grid = Grid('Two-bus analytical example',100,np.array([2.,2.]),np.array([6.,0.]),
                np.array([[0,1,.1,4.,1.]]))
    st = Storage(np.array([0,1]),np.array([12.,12.]),np.array([2.,2.]),
                 np.array([3.,3.]),np.full(2,eta),np.full(2,eta))
    gap = constant_trace('recovery',grid,gap_hours,grid.load,grid.generation,[4.])
    stress = constant_trace('heat-stress',grid,2,[2.,5.],[5.,0.],[event_rating])
    return grid,st,gap,[stress],np.zeros(2)


def rts24(severity=1.0, line_scale=1.0, eta=.95, event_hours=4, gap_hours=24):
    data = json.loads((ROOT/'data/rts24.json').read_text())
    branches = np.asarray(data['branches'],float); branches[:,:2] -= 1
    grid = Grid(data['name'],float(data['base_mva']),np.asarray(data['load_mw']),
                np.asarray(data['generation_mw']),branches)
    buses = np.array([3,8,19,23])-1
    st = Storage(buses,np.full(4,1200.),np.full(4,300.),np.full(4,300.),
                 np.full(4,eta),np.full(4,eta))
    gap = constant_trace('limited-headroom-recovery',grid,gap_hours,grid.load,
                         .9*grid.generation,.95*line_scale*grid.branches[:,3])
    # Temperatures are controlled test inputs, NOT an observed heatwave. The
    # slopes are sensitivity assumptions, not empirically calibrated laws.
    events=[]
    for label,region,extra in [('south',np.arange(10),8.),('central',np.arange(12,20),8.),
                               ('widespread',np.arange(24),2.)]:
        temperature=np.full(grid.n,25.+10.*severity)
        temperature[region]+=extra*severity
        excess=np.maximum(temperature-25.,0)
        demand=grid.load*(1.+.012*excess)
        capacity=.9*grid.generation*(1.-.005*excess)
        line_excess=np.array([(excess[int(a)]+excess[int(b)])/2 for a,b,*_ in grid.branches])
        rating=line_scale*grid.branches[:,3]*(1.-.01*line_excess)
        events.append(constant_trace(label,grid,event_hours,demand,capacity,rating))
    initial=.1*st.emax
    return grid,st,gap,events,initial
