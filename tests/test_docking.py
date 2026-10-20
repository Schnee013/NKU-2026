"""Synthetic score/coordinate/contact joins and invalid docking inputs."""
import csv,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from src.modeling.docking.summarize import reconstruct, read_csv, verify_scientific_bodies, read_pose_models
from src.modeling.docking.output_safety import protect_output
from src.modeling.docking.run_vina import parse_config
from examples.synthetic_data import write_docking_dataset
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'studies/receptor-engineering/receptor-comparison'

class DockingReferenceTests(unittest.TestCase):
    def test_synthetic_paired_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            data=write_docking_dataset(Path(tmp)/'data')
            rows,stats=reconstruct(data)
            self.assertEqual(stats,dict(docking_runs=2,receptor_models=1,score_rows=2))
            self.assertEqual(rows,[dict(model_rel='toy',ascr3_best_kcal_mol=-6.,ascr18_best_kcal_mol=-7.,ascr18_minus_ascr3_kcal_mol=-1.)])
    def test_synthetic_coordinate_and_contact_materials(self):
        with tempfile.TemporaryDirectory() as tmp:
            data=write_docking_dataset(Path(tmp)/'data')
            self.assertEqual(verify_scientific_bodies(data),dict(docking_runs=2,coordinate_poses=2,contact_features=2))
            table=next(data.glob('reference/prolif/**/interaction_frequency.csv'))
            table.write_text(table.read_text().replace('1.0,2','0.5,2'))
            with self.assertRaisesRegex(ValueError,'Contact summary mismatch'):verify_scientific_bodies(data)
    def test_contact_residue_identity_matches_retained_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            data=write_docking_dataset(Path(tmp)/'data')
            fasta=next(data.glob('inputs/receptors/fasta/*.fasta'))
            fasta.write_text('>synthetic changed sequence\nG\n')
            with self.assertRaisesRegex(ValueError,'IFP atom/residue identity mismatch'):verify_scientific_bodies(data)
    def test_constructs_preserve_model_set_and_original_numbering(self):
        rows=read_csv(DATA/'inputs/constructs.csv');groups={}
        for r in rows:groups.setdefault((r['receptor'],r['construct']),[]).append(r)
        self.assertEqual(len(groups),7)
        for key,items in groups.items():self.assertEqual({int(x['model_index']) for x in items},set(range(5)))
        for r in rows:
            if r['receptor']=='gpr2_twf154_wt' and r['construct']=='deltaN_aligned':self.assertEqual((r['residues'],r['present_original_residue_ranges']),('370','22-163;179-406'))
            if r['receptor']=='gpr3_twf154_wt' and r['construct']=='deltaN_aligned':self.assertEqual((r['residues'],r['present_original_residue_ranges']),('380','8-149;164-401'))
    def test_reference_and_existing_output_rejected(self):
        with self.assertRaises(ValueError):protect_output(DATA/'reference/overwrite')
        with tempfile.TemporaryDirectory() as d:
            (Path(d)/'partial').write_text('incomplete')
            with self.assertRaises(FileExistsError):protect_output(d,empty_directory=True)
    def write_score_fixture(self,root,runs=None,scores=None):
        if runs is None:
            runs=[{'result_id':'model/'+lig,'model_rel':'model','ligand':lig,'logged_scores':1}
                  for lig in ['ascr3_deprotonated','ascr18_deprotonated']]
        if scores is None:
            scores=[{'result_id':r['result_id'],'pose':1,'affinity_kcal_mol':-6-i,'rmsd_lb_A':0,'rmsd_ub_A':0,
                     'pose_affinity_kcal_mol':'','pose_rmsd_lb_A':'','pose_rmsd_ub_A':'','pose_present':'false'}
                    for i,r in enumerate(runs)]
        (root/'reference').mkdir(exist_ok=True)
        for name,rows in [('runs.csv',runs),('pose-scores.csv',scores)]:
            with (root/'reference'/name).open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        return runs,scores
    def test_score_summary_does_not_require_coordinate_materials(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.write_score_fixture(root)
            rows,stats=reconstruct(root)
            self.assertEqual(rows[0]['ascr18_minus_ascr3_kcal_mol'],-1)
            self.assertEqual(stats,{'docking_runs':2,'receptor_models':1,'score_rows':2})
            with self.assertRaises(FileNotFoundError):verify_scientific_bodies(root)
    def test_duplicate_run_and_pair_identities_rejected(self):
        for duplicate_result in [True,False]:
            with self.subTest(duplicate_result=duplicate_result),tempfile.TemporaryDirectory() as d:
                root=Path(d);runs,scores=self.write_score_fixture(root)
                extra=dict(runs[0])
                if not duplicate_result:extra['result_id']='different-result-id'
                self.write_score_fixture(root,runs+[extra],scores)
                with self.assertRaisesRegex(ValueError,'Duplicate docking run identity'):reconstruct(root)
    def test_duplicate_pose_identity_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);runs,scores=self.write_score_fixture(root)
            runs[0]['logged_scores']=2
            self.write_score_fixture(root,runs,scores+[dict(scores[0])])
            with self.assertRaisesRegex(ValueError,'Incomplete/duplicate score keys'):reconstruct(root)
    def test_nonfinite_scores_rejected(self):
        for value in ['NaN','Infinity','-Infinity']:
            with self.subTest(value=value),tempfile.TemporaryDirectory() as d:
                root=Path(d);runs,scores=self.write_score_fixture(root);scores[0]['affinity_kcal_mol']=value
                self.write_score_fixture(root,runs,scores)
                with self.assertRaisesRegex(ValueError,'Nonfinite'):reconstruct(root)
    def test_coordinate_model_requires_end_and_atoms(self):
        for text in ['MODEL 1\nREMARK VINA RESULT: -6 0 0\n', 'MODEL 1\nREMARK VINA RESULT: -6 0 0\nENDMDL\n']:
            with self.subTest(text=text),tempfile.TemporaryDirectory() as d:
                path=Path(d)/'poses.pdbqt';path.write_text(text)
                with self.assertRaises(ValueError):read_pose_models(path)
    def test_duplicate_config_field_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'config';p.write_text('seed = 1919\nseed = 0\n')
            with self.assertRaises(ValueError):parse_config(p)
if __name__=='__main__':unittest.main()
