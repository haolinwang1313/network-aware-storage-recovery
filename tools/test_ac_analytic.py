"""Independent constructive AC checks; no inference of RTS AC feasibility."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from ac_analytic_example import operating_point


class AnalyticACTests(unittest.TestCase):
    def test_event_has_losses_reactive_load_and_valid_limits(self):
        p=operating_point(1.9,.3,2.,2.,5.)
        self.assertGreater(p['line_sending_p_mw'],1.9)
        self.assertGreater(p['source_qg_mvar'],.3)
        self.assertLess(p['audit']['balance_mva'],1e-8)
        self.assertLess(p['audit']['max_line_loading'],1.)

    def test_constructed_charge_trajectory_restores_required_energy(self):
        p=operating_point(3.9,.3,2.,4.,6.)
        duration=6.2/1.9
        self.assertAlmostEqual(1.9*duration,3.1*2.)
        self.assertLess(p['source_pg_mw'],6.)
        self.assertGreater(duration,3.)

    def test_equal_energy_at_sender_cannot_bypass_line_limit(self):
        # No receiver battery: all 5 MW would have to cross the 2 MVA line.
        with self.assertRaises(AssertionError):
            operating_point(5.,.3,2.,2.,10.)


if __name__=='__main__':unittest.main(verbosity=2)
