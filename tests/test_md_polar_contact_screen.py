"""Scientific counting contracts for all-residue polar screening."""
import sys
from pathlib import Path
import unittest
import numpy as np

from src.modeling.md_analysis.polar_contact_screen import aggregate_flags,distinct_residues
from src.modeling.md_analysis.analyze_ligand_contact_patterns import frequent_itemsets,runs


class PolarScreenTest(unittest.TestCase):
    def test_two_oxygens_and_multiple_hydrogens_count_as_one_contact(self):
        flags=aggregate_flags(np.array([True,True,False,True]),[0,0,1,2],4)
        self.assertEqual(flags.tolist(),[True,False,True,False])

    def test_salt_and_hbond_to_the_same_residue_are_not_two_anchors(self):
        defs=[{'resid':345},{'resid':345},{'resid':181}]
        self.assertFalse(distinct_residues((0,1),defs))
        self.assertTrue(distinct_residues((0,2),defs))

    def test_frequent_singles_can_have_no_persistent_joint_pattern(self):
        m=np.array([[1,0],[1,0],[1,1],[0,1],[0,1]],dtype=bool)
        results={i for i,_ in frequent_itemsets(m,3)}
        self.assertEqual(results,{(0,),(1,)})

    def test_no_order_cap_and_no_duration_from_missing_frames(self):
        m=np.ones((5,4),dtype=bool)
        self.assertEqual(len(list(frequent_itemsets(m,4))),15)
        self.assertEqual(runs(np.array([1,1,0,1,1],dtype=bool),np.arange(5)*.01),
                         [(0.,.01,.01),(.03,.04,.010000000000000002)])


if __name__=='__main__':unittest.main()
