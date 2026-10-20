"""Mutation-screen scientific joins, numerical comparisons, and output protection."""
import csv
import importlib.util
from decimal import Decimal
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace
from src.modeling.mutant_screen.summarize import reconstruct, validate_output, verify_scientific_bodies, read_csv, POSE_COLUMNS, read_contacts, load_runs, load_pose_scores, select_minimum
from src.modeling.mutant_screen.prepare import EXPECTED, prepare

ROOT=Path(__file__).resolve().parents[1]
STUDY=ROOT/'studies/receptor-engineering/mutant-screen'


class MutantStudyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from examples.synthetic_data import write_mutant_dataset
        cls.temporary=tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.study=write_mutant_dataset(Path(cls.temporary.name)/'mutants')
        cls.result=reconstruct(cls.study)
        cls.records={r['record_id']:r for r in cls.result['records']}

    def test_independent_scores_parent_contrasts_and_saturation(self):
        self.assertEqual(self.result['counts'],dict(receptors=81,ligand_runs=162,score_rows=162,
            single_structure_records=77,unique_single_designs=76,source_WT_single_batches=3,
            saturation_designs=76,batch08_combinations=1,batch08_parents=2))
        for site,stats in self.result['saturation_statistics'].items():
            self.assertEqual((stats['n'],stats['D_positive'],stats['D_zero']),(19,19,0))
            self.assertEqual(Decimal(stats['D_median']),Decimal(1))
        combo=self.records['08/synthetic/N84A_V368A/af0/seed1919']
        self.assertEqual((combo['deltaS3'],combo['deltaS18'],combo['D']),('-2','1','3'))
        self.assertEqual(combo['epsilon']['05'],dict(S3='0',S18='1',D='1'))
        for parent in combo['parents'].values():
            self.assertEqual((parent['score_change3'],parent['score_change18'],parent['D_change']),('-1','1','2'))
        wt=self.records['05/synthetic/baseline/af0/seed1919']
        self.assertEqual(wt['D'],'0')
        self.assertEqual(len([r for r in self.records.values() if r['variant']=='W338L']),2)

    def test_complete_synthetic_pose_and_contact_joins(self):
        self.assertEqual(verify_scientific_bodies(self.study,load_runs(self.study)),dict(ligand_runs=162,coordinate_poses=162))

    def test_parent_geometry_has_known_two_angstrom_translation(self):
        from src.modeling.mutant_screen.geometry import reconstruct as geometry
        rows=geometry(self.study)
        self.assertEqual(len(rows),4)
        for row in rows:
            self.assertEqual(row['atom_count'],3)
            self.assertAlmostEqual(row['heavy_atom_rmsd_no_fit_A'],2.)
            self.assertAlmostEqual(row['centroid_shift_A'],2.)

    def test_fresh_output_required(self):
        with tempfile.TemporaryDirectory() as d:
            for dest in (self.study,self.study/'new',self.study.parent,Path(d)):
                with self.subTest(dest=dest),self.assertRaises(ValueError):validate_output(self.study,dest)
            validate_output(self.study,Path(d)/'fresh')

    def test_missing_input_is_an_error(self):
        from examples.synthetic_data import write_mutant_dataset
        with tempfile.TemporaryDirectory() as d:
            study=write_mutant_dataset(Path(d)/'mutants')
            next(study.glob('results/contacts/**/contact_residues.csv')).unlink()
            with self.assertRaises(FileNotFoundError):reconstruct(study)


class ScoreInputTests(unittest.TestCase):
    def fixture(self,study):
        run={'run_id':'example/ascr3','pose_count':2,'selected_mode':1,
             'selected_pose':{'tied_minimum_modes':[1,2]}}
        rows=[dict(zip(POSE_COLUMNS,['example/ascr3',str(mode),'-6.000','0','0','-6.000','0','0'])) for mode in [1,2]]
        (study/'reference').mkdir(exist_ok=True)
        self.write_rows(study,rows)
        return run,rows

    def write_rows(self,study,rows):
        with (study/'reference/pose-scores.csv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=POSE_COLUMNS);writer.writeheader();writer.writerows(rows)

    def test_every_numeric_field_including_last_rejects_nonfinite(self):
        for field in POSE_COLUMNS[2:]:
            for value in ['NaN','Infinity','-Infinity']:
                with self.subTest(field=field,value=value),tempfile.TemporaryDirectory() as d:
                    study=Path(d);run,rows=self.fixture(study);rows[0][field]=value;self.write_rows(study,rows)
                    with self.assertRaisesRegex(ValueError,'nonfinite'):load_pose_scores(study,[run])

    def test_complete_mode_set_and_actual_tie_set(self):
        with tempfile.TemporaryDirectory() as d:
            study=Path(d);run,rows=self.fixture(study)
            self.assertEqual(select_minimum(run,rows)['mode'],'1')
            run['selected_pose']['tied_minimum_modes']=[1]
            with self.assertRaisesRegex(ValueError,'tie set'):select_minimum(run,rows)
            run['selected_pose']['tied_minimum_modes']=[1,2];run['selected_mode']=2
            with self.assertRaisesRegex(ValueError,'selection'):select_minimum(run,rows)
            rows[1]['mode']='3';self.write_rows(study,rows)
            with self.assertRaisesRegex(ValueError,'mode set'):load_pose_scores(study,[run])


class ContactInputTests(unittest.TestCase):
    def test_synthetic_fingerprints_keep_their_pose_denominators(self):
        from examples.synthetic_data import write_mutant_dataset
        with tempfile.TemporaryDirectory() as d:
            study=write_mutant_dataset(Path(d)/'mutants')
            for run in load_runs(study):
                self.assertEqual(len(read_contacts(study/run['contacts_path'],run)),1)
            run=load_runs(study)[0];run['selected_pose']['analyzed_pose_count']=2
            with self.assertRaisesRegex(ValueError,'denominator'):read_contacts(study/run['contacts_path'],run)

    def test_zero_contacts_and_missing_values_remain_distinct(self):
        run={'contact_parameters':{'interactions':['VdWContact']},'selected_pose':{'analyzed_pose_count':1,'pose_index_in_ifp':0}}
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            (folder/'contact_residues.csv').write_text('residue,original_number,interaction_types\n')
            (folder/'ifp_presence.csv').write_text('ligand,UNL1\nprotein,ALA22.A\ninteraction,VdWContact\nFrame,\n0,False\n')
            self.assertEqual(read_contacts(folder,run),{})
            (folder/'ifp_presence.csv').write_text('ligand,UNL1\nprotein,ALA22.A\ninteraction,VdWContact\nFrame,\n0,\n')
            with self.assertRaisesRegex(ValueError,'Missing fingerprint value'):
                read_contacts(folder,run)


class PreparationTests(unittest.TestCase):
    def prepare_receptor(self, change_charge=False, missing_pdbqt=False, wrong_sequence=False, route='saved-construct'):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);study=root/'study';saved=root/'user-inputs'
            saved.mkdir(parents=True)
            (study/'inputs').mkdir(parents=True)
            receptor=''.join(f'ATOM  {i:5d}  CA  ALA A{n:4d}    {1.0:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00  0.00    -0.100 C\n'
                             for i,n in enumerate(sorted(EXPECTED),1))
            (saved/'input.pdb').write_text(receptor.replace('ALA A  22','GLY A  22') if wrong_sequence else receptor)
            (saved/'receptor.pdbqt').write_text(receptor)
            if missing_pdbqt:(saved/'receptor.pdbqt').unlink()
            (study/'inputs/sequence.fasta').write_text('>example\n'+'A'*406+'\n')
            (study/'inputs/receptors.json').write_text(json.dumps({'construct_definition':{
                'id':'test-construct','chain':'A','prediction_residue_count':406,'docking_residue_count':370,
                'retained_original_ranges':[[22,163],[179,406]],'deleted_original_ranges':[[1,21],[164,178]],
                'renumbered':False},'receptors':[{
                'construct_id':'test-construct','record_id':'example','batch':'05','alias':'example','variant':'baseline',
                'prediction':{'fasta_paths':['inputs/sequence.fasta']},
                'preparation_route':'full-length AF3, deletion and TM alignment'}]}))
            output=root/'prepared'
            def run_command(*args,**kwargs):
                # Formatting outside ATOM records does not alter the prepared molecule.
                generated='REMARK prepared receptor\n'+receptor
                if change_charge:generated=generated.replace('-0.100','-0.200')
                (output/'receptor.pdbqt').write_text(generated)
                (output/'receptor.json').write_text('{}')
            polymer=SimpleNamespace(get_valid_monomers=lambda:{f'A:{n}' for n in EXPECTED})
            meeko=SimpleNamespace(Polymer=SimpleNamespace(from_json=lambda text:polymer))
            with patch.dict('sys.modules',{'meeko':meeko}),patch('src.modeling.mutant_screen.prepare.subprocess.run',side_effect=run_command):
                try:
                    return prepare(study,'example',None,output,route=route,
                                   construct_pdb=saved/'input.pdb',construct_pdbqt=saved/'receptor.pdbqt')
                finally:
                    if missing_pdbqt or wrong_sequence or route=='full-prediction':
                        self.assertFalse(output.exists())

    def test_preparation_accepts_identical_chemistry_with_different_headers(self):
        result=self.prepare_receptor()
        self.assertEqual(result['prepared_coordinate_validation']['max_heavy_atom_delta_A'],0.0)
        self.assertEqual(result['charge_model'],'gasteiger')

    def test_preparation_rejects_changed_heavy_atom_charge(self):
        with self.assertRaisesRegex(ValueError,'charge or AutoDock type'):
            self.prepare_receptor(change_charge=True)

    def test_missing_user_input_fails_before_creating_output(self):
        with self.assertRaisesRegex(ValueError,'user-supplied construct PDBQT'):
            self.prepare_receptor(missing_pdbqt=True)

    def test_preparation_rejects_construct_for_wrong_record(self):
        with self.assertRaisesRegex(ValueError,'Construct sequence does not match selected record'):
            self.prepare_receptor(wrong_sequence=True)

    def test_prediction_route_requires_user_reference(self):
        with self.assertRaisesRegex(ValueError,'common alignment reference'):
            self.prepare_receptor(route='full-prediction')


if __name__=='__main__':unittest.main()
