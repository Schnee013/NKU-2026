"""Sequence identity, scoring definition and failure-boundary checks."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from src.modeling.gprotein.design import STUDY, designs, read_fasta
from src.modeling.gprotein.parse_af3_results import compute_interface_pae, find_best_model

class GproteinTests(unittest.TestCase):
    def test_all_grafts_equal_retained_sequences(self):
        rows=list(designs())
        self.assertEqual(len(rows),11)
        for name,seq,donor,n,end in rows:
            with self.subTest(name=name):
                self.assertEqual(seq,read_fasta(STUDY/'chimeras'/(name+'.fasta')))
                self.assertEqual(len(seq),472)
                self.assertEqual(end+n,472)
    def test_all_prediction_jobs_use_retained_chimeras(self):
        seqs={name:seq for name,seq,*_ in designs()}
        jobs=[]
        with (STUDY/'prediction-jobs.csv').open() as f:
            records=list(csv.DictReader(f))
        for record in records:
            jobs.append(json.loads((STUDY/record['input_path']).read_text()))
        self.assertEqual(len(jobs),42)
        self.assertEqual(len({j['name'] for j in jobs}),42)
        for j in jobs:
            matching=[name for name in seqs if name in j['name']]
            self.assertEqual(len(matching),1,j['name'])
            self.assertEqual(j['dialect'],'alphafold3')
            self.assertTrue(all(type(seed) is int for seed in j['modelSeeds']))
            self.assertIn(seqs[matching[0]],[x['protein']['sequence'] for x in j['sequences']],j['name'])
    def test_prediction_plan_has_unique_score_identities(self):
        with (STUDY/'prediction-jobs.csv').open() as f:records=list(csv.DictReader(f))
        self.assertEqual(len(records),42)
        self.assertEqual(len({r['score_id'] for r in records}),42)
        # test_af3.py exercises the parser with synthetic outputs for all 42 jobs
        # in both native and numbered layouts.
    def test_cross_chain_all_directed_token_pairs(self):
        # Independent arithmetic: (2+4+6+8)/4 = 5; within-chain zeroes excluded.
        self.assertEqual(compute_interface_pae({'token_chain_ids':['A','A','B'],'pae':[[0,0,2],[0,0,4],[6,8,0]]}),5)
    def test_single_chain_has_no_interface(self):
        self.assertEqual(compute_interface_pae({'token_chain_ids':['A'],'pae':[[0]]}),float('inf'))
    def test_model_selection_uses_ranking_and_first_tie(self):
        with tempfile.TemporaryDirectory() as d:
            for idx,score in enumerate([0.1,0.6,0.8,0.8,0.4]):
                Path(d,f'fold_summary_confidences_{idx}.json').write_text(json.dumps({'ranking_score':score}))
            self.assertEqual(find_best_model(Path(d)).model_index,2)
    def test_missing_predictions_fail_without_output(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'output'
            p=subprocess.run([sys.executable,'-m','src.modeling.gprotein.parse_af3_results','--downloads',str(Path(d)/'missing'),'--output',str(out)],capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0)
            self.assertIn('does not exist',p.stderr)
            self.assertFalse(out.exists())
    def test_design_output_conflict_preserves_previous(self):
        with tempfile.TemporaryDirectory() as d:
            marker=Path(d)/'existing';marker.write_text('keep')
            p=subprocess.run([sys.executable,'-m','src.modeling.gprotein.design','--output',d],capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0)
            self.assertEqual(marker.read_text(),'keep')
if __name__=='__main__':
    unittest.main()
