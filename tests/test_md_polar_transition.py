import importlib.util
from pathlib import Path
import unittest
import numpy as np

from src.modeling.md_analysis import polar_transition as module


class EndpointMotionTests(unittest.TestCase):
    def test_exact_additivity_with_acceptor_exchange(self):
        d0, d1 = np.array([0.,0,0]), np.array([3.,1,0])
        a0, a1 = np.array([[4.,0,0],[7,0,0]]), np.array([[8.,0,0],[5,0,0]])
        r = module.endpoint_decomposition(d0,a0,d1,a1)
        self.assertAlmostEqual(r['donor_motion']+r['acceptor_motion'],r['delta'])
        self.assertEqual(r,module.endpoint_decomposition(d0,a0[::-1],d1,a1[::-1]))

    def test_rigid_coordinate_invariance_and_fixed_endpoint(self):
        d0,d1=np.array([0.,0,0]),np.array([1.,0,0]);a=np.array([[4.,0,0],[5.,0,0]])
        r=module.endpoint_decomposition(d0,a,d1,a)
        self.assertAlmostEqual(r['acceptor_motion'],0.)
        self.assertAlmostEqual(r['donor_motion'],-1.)
        q=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        shifted=module.endpoint_decomposition(d0@q+8,a@q+8,d1@q+8,a@q+8)
        for k in r:self.assertAlmostEqual(r[k],shifted[k])

    def test_fit_recovers_rigid_structure_without_reflection(self):
        x=np.array([[0,0,0],[1,0,0],[0,2,0],[0,0,3]],float)
        q=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        np.testing.assert_allclose(module.fit(x@q+5,np.arange(4),x),x,atol=1e-12)


if __name__=='__main__':unittest.main()
