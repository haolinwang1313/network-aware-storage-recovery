"""Reproducible editable SVG and 300 dpi PNG figures; data remain in CSV.

Contracts: (1) observed hourly temperature, no electrical-effect claim;
(2) fixed-clock charging-window comparison, zero-based ENS axis;
(3) forecast envelope, both design window lengths shown in legend.
Default matplotlib colours; line styles distinguish series in monochrome.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from cases import ROOT


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=ROOT/'results_v2')
    p.add_argument('--weather',type=Path,default=ROOT/'data/public/weather_hourly.csv')
    args=p.parse_args();out=args.results/'figures';out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    manifest=[]
    def save(fig,name,source):
        fig.tight_layout()
        for suffix in ['svg','png']:
            path=out/(name+'.'+suffix);fig.savefig(path,dpi=300,bbox_inches='tight')
            manifest.append(dict(path=str(path.relative_to(args.results)),source=source,
                                 sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        plt.close(fig)
    rows=list(csv.DictReader(args.weather.open()))
    stamps=[datetime.fromisoformat(r['local']) for r in rows]
    fig,ax=plt.subplots(figsize=(8.2,3.5))
    for name,style in zip(['sacramento','fresno','san_francisco'],['-','--',':']):
        ax.plot(stamps,[float(r[name+'_c']) for r in rows],linestyle=style,label=name.replace('_',' ').title())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d',tz=ZoneInfo('America/Los_Angeles')))
    ax.set(ylabel='Air temperature (°C)',xlabel='Local date (America/Los_Angeles)',
           title='NOAA station temperatures, September 2022')
    ax.legend(ncol=3,frameon=False,loc='upper left');ax.grid(axis='y',alpha=.2)
    save(fig,'observed_temperature','data/public/weather_hourly.csv')
    rows=json.loads((args.results/'all_windows_day5.json').read_text())
    fig,ax=plt.subplots(figsize=(7.2,3.7))
    for key,label,style,marker in [('proposed_ens_mwh','Event-conditioned recovery','-','o'),
                                  ('baseline_ens_mwh','Myopic network-aware balancing','--','s')]:
        y=[np.nan if r.get(key) is None else (0 if r[key]<=1e-5 else r[key]) for r in rows]
        ax.plot([r['window_h'] for r in rows],y,linestyle=style,marker=marker,markersize=3.5,label=label)
    ax.set(xlabel='Available charging window before fixed event onset (h)',ylabel='Worst-case ENS (MWh)',
           title='RTS-24, September 5–6 chronology; three future temperature scenarios',xlim=(0,18),ylim=(0,None))
    ax.set_xticks(range(0,19,3));ax.legend(frameon=False);ax.grid(alpha=.2)
    save(fig,'charging_window','results_v2/all_windows_day5.csv')
    rows=json.loads((args.results/'forecast_envelope.json').read_text())
    fig,ax=plt.subplots(figsize=(7.2,3.7))
    for name,label,style in [('nominal','Nominal design','--'),('0_to_2C_envelope','0–2 °C envelope','-')]:
        part=[r for r in rows if r['design']==name]
        ax.plot([r['temperature_error_c'] for r in part],[0 if r['ens_mwh']<=1e-5 else r['ens_mwh'] for r in part],
                linestyle=style,marker='o',markersize=4,label=f"{label}: {part[0]['window_h']:g} h charging window")
    ax.set(xlabel='Future-episode temperature perturbation (°C)',ylabel='ENS (MWh)',ylim=(0,None),
           title='Specified forecast envelope and out-of-envelope stress tests')
    ax.legend(frameon=False);ax.grid(alpha=.2)
    save(fig,'forecast_envelope','results_v2/forecast_envelope.csv')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':main()
