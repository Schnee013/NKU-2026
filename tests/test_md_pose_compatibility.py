import sys
from pathlib import Path
import unittest
import json
import os
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
from src.modeling.md_analysis.pose_compatibility import energy_gradient, mapped_names, stereo_valid


class CompatibilityGeometryTests(unittest.TestCase):
    def test_steric_gradient_matches_finite_difference(self):
        class ZeroField:
            def CalcEnergy(self,x):return 0.
            def CalcGrad(self,x):return np.zeros(len(x))
        x=np.array([[0.,0.,0.],[2.,1.,0.]])
        env=np.array([[1.,0.,0.],[3.,0.,1.]])
        free=np.array([1]);moving=np.array([1]);limits=np.array([[3.,3.]])
        energy,gradient=energy_gradient(x,ZeroField(),free,moving,env,limits,100.)
        numerical=[]
        for k in range(3):
            a=x.copy();b=x.copy();a[1,k]+=1e-6;b[1,k]-=1e-6
            numerical.append((energy_gradient(a,ZeroField(),free,moving,env,limits,100.)[0]-energy_gradient(b,ZeroField(),free,moving,env,limits,100.)[0])/2e-6)
        np.testing.assert_allclose(gradient,numerical,rtol=1e-6,atol=1e-6)

    def test_real_chemical_maps_and_stereo(self):
        from rdkit import Chem
        base=ROOT/'studies/receptor-engineering/molecular-dynamics/configs'
        contracts=json.loads((base/'ligand_geometry_contract.json').read_text())
        for name in ['ascr3','ascr18']:
            c=contracts['v368r_'+name];mol=Chem.SDMolSupplier(str(base/(name+'.sdf')),removeHs=False)[0]
            names=mapped_names(mol,c);x=mol.GetConformer().GetPositions()
            self.assertTrue(stereo_valid(x,names,c))
            mirror=x.copy();mirror[:,0]*=-1;self.assertFalse(stereo_valid(mirror,names,c))
            Chem.AssignStereochemistryFrom3D(mol,confId=0,replaceExistingTags=True)
            for s in c['stereocenters']:self.assertEqual(mol.GetAtomWithIdx(names.index(s['center'])).GetProp('_CIPCode'),s['source_cip'])


if __name__=='__main__':unittest.main()
