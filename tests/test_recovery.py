"""Scientific invariants, analytic counterexamples and formulation checks."""
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from cases import two_bus, rts24
from minimum_time import minimum_time
from recovery import (TOL, recourse, recover_cuts, recover_monolithic, first_ready,
                      full_recovery_time, time_to_go, Storage, Trace, audit)


class RecoveryTests(unittest.TestCase):
    def test_01_equal_total_energy_is_not_sufficient(self):
        g,s,gap,events,e0=two_bus()
        q=lambda e:recourse(g,s,events[0],np.asarray(e))['ens']
        self.assertAlmostEqual(q([6,0]),6)
        self.assertAlmostEqual(q([0,6]),0)
        self.assertAlmostEqual(q([8,0]),6)

    def test_02_analytical_time_and_cut(self):
        g,s,gap,events,e0=two_bus()
        out=first_ready(g,s,gap,events,e0)
        self.assertEqual(out['ready_h'],3)
        self.assertEqual(out['solution']['cuts'],1)
        np.testing.assert_allclose(out['solution']['energy'],[0,6],atol=TOL)
        cut=out['solution']['cut_coefficients'][0]
        np.testing.assert_allclose(cut['a'],[0,1],atol=TOL)
        self.assertAlmostEqual(cut['b'],6)
        self.assertEqual(full_recovery_time(g,s,gap,e0,s.emax),12)

    def test_03_efficiency_and_hourly_rounding(self):
        for eta in [.8,.9,.95,1.]:
            with self.subTest(eta=eta):
                g,s,gap,events,e0=two_bus(eta=eta)
                expected=int(np.ceil(6/(2*eta*eta)-1e-10))
                self.assertEqual(first_ready(g,s,gap,events,e0)['ready_h'],expected)

    def test_04_uncongested_aggregate_null_case(self):
        g,s,gap,events,e0=two_bus(event_rating=100.)
        rng=np.random.default_rng(20260905)
        for _ in range(40):
            energy=rng.uniform(0,6,2)
            q=recourse(g,s,events[0],energy)['ens']
            self.assertAlmostEqual(q,max(0,4-energy.sum()),places=7)

    def test_05_monolithic_equals_cuts_all_toy_horizons(self):
        g,s,gap,events,e0=two_bus()
        for h in range(7):
            a=recover_cuts(g,s,gap.prefix(h),events,e0)
            b=recover_monolithic(g,s,gap.prefix(h),events,e0)
            self.assertEqual(a['feasible'],b['feasible'])
            if a['feasible']:
                self.assertAlmostEqual(a['charge_mwh'],b['charge_mwh'],places=6)

    def test_06_dual_support_inequalities(self):
        g,s,gap,events,e0=rts24(line_scale=.8)
        rng=np.random.default_rng(42)
        for event in events:
            for _ in range(8):
                anchor=rng.uniform(0,1,s.k)*s.emax
                r=recourse(g,s,event,anchor)
                self.assertTrue(np.all(r['gradient']<=1e-7))
                other=rng.uniform(0,1,s.k)*s.emax
                q=recourse(g,s,event,other)['ens']
                support=r['ens']+r['gradient']@(other-anchor)
                self.assertLessEqual(support,q+1e-5)

    def test_07_rts24_reproduction_and_nonanticipativity(self):
        g,s,gap,events,e0=rts24(line_scale=.8)
        a=first_ready(g,s,gap,events,e0)
        b=first_ready(g,s,gap,events,e0,method='monolithic')
        self.assertEqual(a['ready_h'],b['ready_h'])
        self.assertEqual(a['ready_h'],8)
        self.assertAlmostEqual(a['solution']['charge_mwh'],b['solution']['charge_mwh'],places=5)
        common=a['solution']['energy']
        for r in a['solution']['event_checks']:
            np.testing.assert_allclose(r['path']['energy'][0],common,atol=TOL)
        clairvoyant=[first_ready(g,s,gap,[event],e0)['ready_h'] for event in events]
        self.assertGreaterEqual(a['ready_h'],max(clairvoyant))

    def test_08_network_qp_baseline_is_feasible_and_changes_decision(self):
        g,s,gap,events,e0=two_bus()
        q=time_to_go(g,s,gap.prefix(6),e0)
        np.testing.assert_allclose(q['energies'][3],[3,3],atol=1e-5)
        self.assertAlmostEqual(recourse(g,s,events[0],q['energies'][3])['ens'],3,places=5)
        self.assertLess(recourse(g,s,events[0],q['energies'][6])['ens'],TOL)
        self.assertLess(max(q['optimality_gaps']),5e-5)

    def test_09_full_soc_can_still_fail_with_installed_resources(self):
        g,s,gap,events,e0=rts24(line_scale=.55)
        result=first_ready(g,s,gap,events,e0)
        self.assertIsNone(result['ready_h'])
        self.assertEqual(result['reason'],'full_state_not_event_feasible')
        self.assertGreater(result['full_state_worst_ens'],0)

    def test_10_more_initial_energy_cannot_worsen_stress_ens(self):
        g,s,gap,events,e0=rts24(line_scale=.8)
        rng=np.random.default_rng(7)
        for event in events:
            for _ in range(5):
                e=rng.random(s.k)*s.emax
                higher=e+rng.random(s.k)*(s.emax-e)
                self.assertLessEqual(recourse(g,s,event,higher)['ens'],
                                     recourse(g,s,event,e)['ens']+TOL)

    def test_11_invalid_inputs_fail_loudly(self):
        g,s,gap,events,e0=two_bus()
        with self.assertRaises(ValueError):
            recourse(g,s,events[0],np.array([-1,0]))
        with self.assertRaises(ValueError):
            recourse(g,s,events[0],np.array([13,0]))
        with self.assertRaises(ValueError):
            Trace('invalid',np.ones((2,2)),np.ones((1,2)),np.ones((2,1)))
        with self.assertRaises(ValueError):
            Storage(np.array([0]),np.ones(1),np.ones(1),np.ones(1),np.array([1.1]),np.ones(1))
        with self.assertRaises(ValueError):
            first_ready(g,s,gap,[],e0)

    def test_12_independent_audit_catches_corruption(self):
        g,s,gap,events,e0=two_bus()
        r=recourse(g,s,events[0],np.array([0.,6.]))
        r['path']['flow'][0,0]+=1
        with self.assertRaises(AssertionError):
            audit(g,s,events[0],r['path'],'discharge',np.array([0.,6.]))

    def test_13_rts_source_dimensions(self):
        g,*_=rts24()
        self.assertEqual((g.n,g.m),(24,38))
        self.assertAlmostEqual(g.load.sum(),2850.)
        self.assertAlmostEqual(g.generation.sum(),3405.)
        self.assertEqual(np.count_nonzero(g.branches[:,4]!=1),5)

    def test_14_ready_initial_state_needs_zero_recovery(self):
        g,s,gap,events,e0=two_bus()
        self.assertEqual(first_ready(g,s,gap,events,np.array([0.,6.]))['ready_h'],0)

    def test_15_time_scaling_exactness(self):
        for factory,settings in [(two_bus,{}),(two_bus,{'eta':.9}),
                                 (rts24,{}),(rts24,{'line_scale':.8})]:
            g,s,gap,events,e0=factory(**settings)
            a=minimum_time(g,s,gap,events,e0)
            b=minimum_time(g,s,gap,events,e0,method='cuts')
            discrete=first_ready(g,s,gap,events,e0)
            self.assertTrue(a['feasible'] and b['feasible'])
            self.assertAlmostEqual(a['minimum_h'],b['minimum_h'],places=5)
            self.assertEqual(np.ceil(a['minimum_h']-1e-7),discrete['ready_h'])
        g,s,gap,events,e0=two_bus(eta=.9)
        self.assertAlmostEqual(minimum_time(g,s,gap,events,e0)['minimum_h'],6/(2*.9*.9))

    def test_16_time_scaling_rejects_nonstationary_gap(self):
        g,s,gap,events,e0=two_bus()
        gap.demand[1,0]+=.1
        with self.assertRaises(ValueError):
            minimum_time(g,s,gap,events,e0)


if __name__=='__main__':
    unittest.main(verbosity=2)
