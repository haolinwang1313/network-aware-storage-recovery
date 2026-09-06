"""Physical and data invariants for the completed chronological experiments."""
import csv
import hashlib
import json
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from cases import ROOT,two_bus,rts24
from recovery import Trace
from chronology import Limits,optimise,minimum_window,balanced_recovery,physical_audit
from ac_validation import circuit_audit,load_full_case
from weather_cases import episode_case


class ChronologyTests(unittest.TestCase):
    def test_suffix_window_does_not_delete_hours(self):
        g,s,gap,events,e0=two_bus()
        m=minimum_window(g,s,gap,events,e0,Limits(ramp_fraction=1.))
        self.assertEqual(m['window_h'],3.)
        path=m['solution']['gap']
        self.assertEqual(len(path['g']),gap.h)
        self.assertLess(np.max(path['charge'][:-3]),1e-6)
        self.assertLess(m['solution']['ens'],1e-5)
        self.assertGreater(m['preceding_ens'],1.)

    def test_generator_ramp_can_lengthen_the_window(self):
        g,s,gap,events,e0=two_bus()
        self.assertEqual(minimum_window(g,s,gap,events,e0,Limits(ramp_fraction=.3))['window_h'],4.)

    def test_integer_storage_matches_continuous_toy(self):
        g,s,gap,events,e0=two_bus()
        for w in [2,3,6]:
            a=optimise(g,s,events,e0,Limits(),gap,w)
            b=optimise(g,s,events,e0,Limits(integer_storage=True),gap,w)
            self.assertAlmostEqual(a['primary_ens'],b['primary_ens'],places=5)
            self.assertLess(max(x['simultaneous_mw'] for x in b['audits']),1e-6)

    def test_ramping_requires_generation_state(self):
        g,s,gap,events,_=two_bus()
        e=np.array([0.,6.]);limits=Limits(ramp_fraction=.05)
        low=optimise(g,s,events,e,limits,g0=np.array([0.,0.]))
        high=optimise(g,s,events,e,limits,g0=np.array([4.,0.]))
        self.assertGreater(low['ens'],high['ens']+1.)
        self.assertLess(high['ens'],1e-5)

    def test_discharge_during_recovery_is_allowed(self):
        g,s,_,_,_=two_bus()
        gap=Trace('two-phase',np.array([[0,3],[0,0]]),np.array([[2,0],[2,0]]),np.full((2,1),2.))
        event=Trace('next',np.array([[0,2]]),np.zeros((1,2)),np.full((1,1),2.))
        out=optimise(g,s,[event],np.array([0.,2.]),Limits(ramp_fraction=1),gap,2)
        self.assertTrue(out['ready'])
        self.assertGreater(out['gap']['discharge'][0].sum(),.99)
        self.assertGreater(out['gap']['charge'][1].sum(),.99)

    def test_nonstationary_windows_are_nested(self):
        g,s,gap,events,e0=two_bus()
        gap.capacity[:,0]=np.linspace(5,7,gap.h)
        values=[optimise(g,s,events,e0,Limits(),gap,w)['ens'] for w in range(gap.h+1)]
        self.assertTrue(np.all(np.diff(values)<=1e-5))

    def test_physical_audit_rejects_corrupted_energy(self):
        g,s,gap,events,e0=two_bus()
        out=optimise(g,s,events,e0,Limits(),gap,3)
        out['gap']['energy'][2,0]+=1
        with self.assertRaises(AssertionError):
            physical_audit(g,s,gap,out['gap'],e0,Limits())

    def test_invalid_window_is_rejected(self):
        g,s,gap,events,e0=two_bus()
        for w in [-1,gap.h+1]:
            with self.assertRaises(ValueError):optimise(g,s,events,e0,gap=gap,window=w)

    def test_ac_audit_uses_complex_network_equations(self):
        bus=np.zeros((2,13));bus[:,0]=[1,2];bus[:,7]=1.;bus[:,11]=1.05;bus[:,12]=.95
        gen=np.zeros((1,21));gen[0,0]=1;gen[0,3]=100;gen[0,4]=-100;gen[0,8]=100
        branch=np.array([[1,2,.01,.1,0,100,100,100,0,0,1,-360,360]])
        case={'baseMVA':100.,'bus':bus,'gen':gen,'branch':branch}
        solution={'bus':bus.copy(),'gen':gen.copy()}
        self.assertLess(circuit_audit(case,solution)['balance_mva'],1e-10)
        solution['bus'][1,7]=.9
        audit=circuit_audit(case,solution)
        self.assertGreater(audit['balance_mva'],1.)
        self.assertGreater(audit['voltage_pu'],.04)


class PublicDataTests(unittest.TestCase):
    def test_public_weather_is_complete_and_frozen(self):
        path=ROOT/'data/public/weather_hourly.csv'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         'ed5afc7dfca4bb515c20cf7673cb50cb2a190cc0bc3ad91daae3e661ea03c230')
        rows=list(csv.DictReader(path.open()))
        self.assertEqual(len(rows),240)
        self.assertEqual(len({r['utc'] for r in rows}),240)
        for key in ['sacramento','fresno','san_francisco']:
            self.assertEqual(sum(int(r[key+'_filled']) for r in rows),0)

    def test_full_source_matches_dc_subset(self):
        full=load_full_case();g,*_=rts24()
        np.testing.assert_allclose(full['bus'][:,2],g.load)
        grouped=np.bincount(full['gen'][:,0].astype(int)-1,weights=full['gen'][:,8],minlength=24)
        np.testing.assert_allclose(grouped,g.generation)
        np.testing.assert_allclose(full['branch'][:,3],g.branches[:,2])
        self.assertEqual(full['gen'].shape[0],33)

    def test_events_have_fixed_clock_and_eighteen_hour_gap(self):
        c=episode_case(day=5)
        self.assertEqual(c['first'].h,6)
        self.assertEqual(c['gap'].h,18)
        self.assertEqual(c['next_start'].isoformat(),'2022-09-06T14:00:00-07:00')
        self.assertFalse(np.allclose(c['gap'].demand,c['gap'].demand[0]))

    def test_unmodified_rts_ac_solver_and_audit(self):
        from pypower.api import runopf,ppoption
        full=load_full_case()
        result=runopf(full,ppoption(VERBOSE=0,OUT_ALL=0))
        self.assertTrue(result['success'])
        a=circuit_audit(full,result)
        self.assertLess(a['balance_mva'],1e-3)
        self.assertLess(a['thermal_mva'],1e-3)
        self.assertLess(a['voltage_pu'],1e-5)


if __name__=='__main__':unittest.main(verbosity=2)
