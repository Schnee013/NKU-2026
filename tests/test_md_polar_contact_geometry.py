"""Known geometry and time boundaries for the shared polar-contact extension."""
import sys
from pathlib import Path
import unittest
import numpy as np

from src.modeling.md_analysis.polar_contacts import short_linear_geometry, fill_short_gaps, window_specs, summarize_mask


class PolarContactGeometryTest(unittest.TestCase):
    limits = dict(donor_acceptor_max_A=3.2, H_acceptor_max_A=2.4, DHA_min_deg=150.)

    def test_known_linear_and_bent_geometries(self):
        for acceptor, wanted_angle, expected in [([2.7,0,0],180.,True), ([1,2,0],90.,False), ([3.3,0,0],180.,False)]:
            x=np.array([[0.,0,0],[1.,0,0],acceptor])
            da,ha,angle,valid=short_linear_geometry(x,[[0,1,2]],self.limits)
            self.assertAlmostEqual(angle[0],wanted_angle)
            self.assertEqual(bool(valid[0]),expected)
        x=np.array([[0.,0,0],[.7,0,0],[3.2,0,0]])
        self.assertFalse(short_linear_geometry(x,[[0,1,2]],self.limits)[3][0])

    def test_invalid_wrapped_covalent_bond_is_not_a_hbond(self):
        with self.assertRaisesRegex(ValueError,'covalent'):
            short_linear_geometry(np.array([[0.,0,0],[10.,0,0],[11.,0,0]]),[[0,1,2]],self.limits)

    def test_gap_tolerance_never_fills_endpoints_or_long_absences(self):
        mask=np.array([0,1,0,0,1,0,0,0,1,0],dtype=bool)
        result=fill_short_gaps(mask,.01,.02)
        self.assertEqual(result.tolist(),[False,True,True,True,True,False,False,False,True,False])
        self.assertEqual(mask.tolist(),[False,True,False,False,True,False,False,False,True,False])

    def test_full_window_and_block_denominators_and_strict_spans(self):
        screen={'primary_window_ns':[50,100],'secondary_windows_ns':[[0,100]]+[[i,i+10] for i in range(0,100,10)]}
        times=np.arange(10001)*.01;mask=np.ones(10001,dtype=bool);mask[5001]=False
        rows=summarize_mask(mask,times,screen,.05)
        lookup={r['window']:r for r in rows}
        self.assertEqual(lookup['primary']['n_frames'],5001)
        self.assertEqual(lookup['full']['n_frames'],10001)
        self.assertEqual(sum(r['n_frames'] for r in rows if r['window'].startswith('block_')),10001)
        self.assertAlmostEqual(lookup['primary']['longest_uninterrupted_ns'],49.98)
        self.assertAlmostEqual(lookup['primary']['longest_gap_tolerant_ns'],50.)
        self.assertAlmostEqual(lookup['primary']['occupancy'],5000/5001)


if __name__=='__main__':
    unittest.main()
