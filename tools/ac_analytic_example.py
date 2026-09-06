"""Constructive nonzero-flow AC counterpart of the spatial-energy example.

No nonlinear solver is used. This is an explicitly modified two-bus example,
not a claim that the original DC toy (3 MW inverter) is AC feasible. Both loads
retain a fixed 0.3 Mvar receiver demand; the storage inverter supplies zero Q.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ac_validation import circuit_audit


def operating_point(receiver_net_p,receiver_q,source_local_p,limit,pmax):
    base=100.; r=.01; x=.1; p=receiver_net_p/base; q=receiver_q/base
    # Exact receiving-end squared-voltage equation for a two-bus pi-free line.
    coefficient=2*(r*p+x*q)-1.
    roots=np.roots([1.,coefficient,(r*r+x*x)*(p*p+q*q)])
    if not np.all(np.isreal(roots)):raise AssertionError('No high-voltage solution.')
    vr=float(np.sqrt(np.max(roots.real)))
    current=(p-1j*q)/vr
    vs=vr+(r+1j*x)*current
    rotation=np.exp(-1j*np.angle(vs))
    voltages=np.array([vs,complex(vr)])*rotation
    sending=voltages[0]*np.conj((voltages[0]-voltages[1])/(r+1j*x))*base
    bus=np.zeros((2,13));bus[:,0]=[1,2];bus[:,1]=[3,1]
    bus[:,2]=[source_local_p,receiver_net_p];bus[:,3]=[0.,receiver_q]
    bus[:,7]=abs(voltages);bus[:,8]=np.rad2deg(np.angle(voltages));bus[:,11]=1.05;bus[:,12]=.95
    gen=np.zeros((1,21));gen[0,0]=1;gen[0,1]=source_local_p+sending.real;gen[0,2]=sending.imag
    gen[0,3]=5.;gen[0,4]=-5.;gen[0,7]=1.;gen[0,8]=pmax
    branch=np.array([[1,2,r,x,0.,limit,limit,limit,0.,0.,1.,-360.,360.]])
    case={'baseMVA':base,'bus':bus,'gen':gen,'branch':branch}
    audit=circuit_audit(case,{'bus':bus,'gen':gen})
    if audit['balance_mva']>1e-8 or max(audit[k] for k in ['voltage_pu','thermal_mva','active_mw','reactive_mvar'])>1e-8:
        raise AssertionError(audit)
    return {'audit':audit,'bus_voltage_complex':[[float(v.real),float(v.imag)] for v in voltages],
            'source_pg_mw':float(gen[0,1]),'source_qg_mvar':float(gen[0,2]),
            'line_sending_p_mw':float(sending.real),'line_sending_q_mvar':float(sending.imag)}


def main():
    event=operating_point(1.9,.3,2.,2.,5.)
    recovery=operating_point(3.9,.3,2.,4.,6.)
    duration=2.;discharge=3.1;energy=duration*discharge;charge=1.9
    # At zero receiver storage, serving 5 MW would require receiving >2 MW
    # through a 2 MVA-rated line. This is impossible regardless of sender energy.
    out={'type':'constructive AC example; not the original 3-MW DC toy',
         'storage_discharge_rating_mw':4.,'storage_charge_rating_mw':2.,'efficiency':1.,
         'event_duration_h':duration,'receiver_event_load_mw':5.,'receiver_q_mvar':.3,
         'receiver_discharge_mw':discharge,'required_receiver_energy_mwh':energy,
         'receiver_only_state_mwh':[0.,energy],'sender_only_equal_energy_state_mwh':[energy,0.],
         'sender_only_unserved_energy_lower_bound_mwh':(5.-2.)*duration,
         'constructed_recovery_h':energy/charge,
         'interpretation':'Constructed feasible duration, NOT the AC minimum recovery time.',
         'event':event,'recovery':recovery}
    path=Path(__file__).resolve().parents[1]/'results_v2/ac_analytic_example.json'
    path.parent.mkdir(exist_ok=True);path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
