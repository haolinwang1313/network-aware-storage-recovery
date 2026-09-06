"""Observed-weather-driven *synthetic* RTS-24 experiments.

Temperatures are observed, but geographic assignment, temperature-to-electric
coefficients, charging availability and equipment modifications are assumptions.
The six-hour stress episodes are afternoon operating windows WITHIN a heatwave,
not six-hour meteorological heatwaves. Local clocks never shift in optimisation.
"""
from __future__ import annotations
import csv
from datetime import datetime, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
from cases import ROOT, rts24
from recovery import Trace

TZ=ZoneInfo('America/Los_Angeles')
SHAPE=np.array([.76,.73,.71,.70,.70,.72,.76,.81,.85,.89,.92,.95,.97,.99,1.,1.,1.,.99,.97,.94,.91,.87,.83,.79])


def load_weather(path: Path | None=None):
    path=ROOT/'data/public/weather_hourly.csv' if path is None else path
    rows=list(csv.DictReader(path.open()))
    by_time={datetime.fromisoformat(r['local']):np.array([float(r[k+'_c']) for k in
             ['sacramento','fresno','san_francisco']]) for r in rows}
    if len(rows)!=len(by_time): raise ValueError('Duplicate weather timestamps.')
    return by_time


def episode_case(day=5,line_scale=.8,charger_mw=100.,eta=.95,loss_margin=0.,line_margin=1.,
                 demand_slope=.012,generation_slope=.005,rating_slope=.01,
                 initial_fraction=.9,weather_path=None,mapping_shift=0):
    if not 1<=day<=9 or not 0<charger_mw<=300 or not 0<=initial_fraction<=1:
        raise ValueError('Invalid experiment inputs.')
    if loss_margin<0 or not 0<line_margin<=1:
        raise ValueError('Margins must not increase planning flexibility.')
    grid,st,_,_,_=rts24(line_scale=line_scale,eta=eta)
    st.pc[:]=charger_mw
    # Electrical areas, NOT coordinates of the standard test system.
    groups=np.r_[np.zeros(10,int),np.ones(6,int),np.full(8,2,int)]
    groups=(groups+mapping_shift)%3
    weather=load_weather(weather_path)
    def trace(start,hours,delta=0.,planning=True):
        timestamps=[start+timedelta(hours=i) for i in range(hours)]
        missing=[t.isoformat() for t in timestamps if t not in weather]
        if missing: raise ValueError('Missing weather hours: '+str(missing[:3]))
        temp=np.array([weather[t][groups]+delta for t in timestamps])
        excess=np.maximum(temp-25.,0)
        shape=np.array([SHAPE[t.hour] for t in timestamps])
        d=grid.load*shape[:,None]*(1+demand_slope*excess)
        c=.9*grid.generation*(1-generation_slope*excess)
        branch_excess=np.array([(excess[:,int(a)]+excess[:,int(b)])/2 for a,b,*_ in grid.branches]).T
        r=line_scale*grid.branches[:,3]*(1-rating_slope*branch_excess)
        if planning: d=d*(1+loss_margin); r=r*line_margin
        return Trace(start.isoformat()+f'_delta{delta}',d,c,r)
    first_start=datetime(2022,9,day,14,tzinfo=TZ)
    gap_start=first_start+timedelta(hours=6)
    next_start=first_start+timedelta(days=1)
    first=trace(first_start,6); gap=trace(gap_start,18)
    events=[trace(next_start,6,delta) for delta in [0.,1.,2.]]
    actual_events=[trace(next_start,6,delta,False) for delta in [0.,1.,2.]]
    return dict(grid=grid,storage=st,first=first,gap=gap,events=events,
                initial=st.emax*initial_fraction,actual_first=trace(first_start,6,planning=False),
                actual_gap=trace(gap_start,18,planning=False),actual_events=actual_events,
                trace_builder=trace,next_start=next_start,
                settings=dict(day=day,line_scale=line_scale,charger_mw=charger_mw,eta=eta,
                              loss_margin=loss_margin,line_margin=line_margin,demand_slope=demand_slope,
                              generation_slope=generation_slope,rating_slope=rating_slope,
                              initial_fraction=initial_fraction,mapping_shift=mapping_shift))
