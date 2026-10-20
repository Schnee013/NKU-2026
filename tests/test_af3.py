"""Native AF3 identity, ranking precision and structural confidence checks."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from src.modeling.af3 import best_prediction, ca_plddts, predictions, single_protein_sequence
from src.modeling.gprotein.parse_af3_results import STUDY, compute_interface_pae
from src.modeling.mutant_screen.prepare import prediction_sequence


def cif(values):
    # Deliberately reordered columns; C-alpha and all-atom averages differ.
    return ('data_test\nloop_\n_atom_site.B_iso_or_equiv\n_atom_site.label_atom_id\n'
            '_atom_site.group_PDB\n' + ''.join(f'{value} CA ATOM\n10 CB ATOM\n' for value in values) + '#\n')


def native_job(folder, request, scores, values=(90, 70)):
    folder.mkdir()
    (folder/'test_data.json').write_text(json.dumps(request))
    with (folder/'test_ranking_scores.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['seed','sample','ranking_score'])
        for seed,sample,score in scores:
            writer.writerow([seed,sample,score])
            prefix=f'test_seed-{seed}_sample-{sample}'
            sample_dir=folder/f'seed-{seed}_sample-{sample}';sample_dir.mkdir()
            (sample_dir/f'{prefix}_model.cif').write_text(cif(values))
            (sample_dir/f'{prefix}_summary_confidences.json').write_text(json.dumps({
                'ranking_score':round(score,2),'ptm':0.6,'iptm':0.5,'fraction_disordered':0.1,
                'has_clash':0.0,'chain_pair_pae_min':[[0,2],[6,0]]}))
            (sample_dir/f'{prefix}_confidences.json').write_text(json.dumps({
                'token_chain_ids':['A','A','B'],'pae':[[0,0,2],[0,0,4],[6,8,0]]}))


class NativeAF3Tests(unittest.TestCase):
    def request(self):
        return {'name':'test','dialect':'alphafold3','version':4,'modelSeeds':[7,11],
                'sequences':[{'protein':{'id':'A','sequence':'ACD','unpairedMsa':'>q\nACD',
                                          'pairedMsa':'','templates':[]}}]}

    def test_native_and_processed_input_reaches_construct_preparation(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'request.json').write_text(json.dumps(self.request()))
            identity={'prediction':{'job_request_path':'request.json'}}
            self.assertEqual(prediction_sequence(root,identity),'ACD')
            (root/'sequence.fasta').write_text('>q\nACE\n')
            identity['prediction']['fasta_paths']=['sequence.fasta']
            with self.assertRaisesRegex(ValueError,'inconsistent'):
                prediction_sequence(root,identity)

    def test_single_chain_copy_and_older_input(self):
        self.assertEqual(single_protein_sequence([{'sequences':[{'proteinChain':{'sequence':'ACD','count':1}}]}]),'ACD')
        request=self.request();request['sequences'][0]['protein']['id']=['A','B']
        with self.assertRaisesRegex(ValueError,'one explicit'):
            single_protein_sequence(request)

    def test_multiple_seeds_full_precision_and_first_exact_tie(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)/'job'
            native_job(folder,self.request(),[(7,0,0.8211),(11,0,0.8249),(7,1,0.8249)])
            chosen=best_prediction(predictions(folder))
            self.assertEqual((chosen.seed,chosen.sample),(11,0))
            self.assertEqual(chosen.ranking_score,0.8249)
            self.assertEqual(chosen.summary()['ranking_score'],0.82)
            self.assertEqual(chosen.identity(),{'best_model':'','best_seed':11,'best_sample':0})

    def test_incomplete_or_ambiguous_native_outputs_fail(self):
        for fault in ['missing_csv','duplicate_row','missing_cif','extra_sample','duplicate_summary']:
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as d:
                folder=Path(d)/'job'
                request=self.request();request['modelSeeds']=[7]
                native_job(folder,request,[(7,0,0.82)])
                if fault=='missing_csv':(folder/'test_ranking_scores.csv').unlink()
                elif fault=='duplicate_row':
                    with (folder/'test_ranking_scores.csv').open('a') as f:f.write('7,0,0.82\n')
                elif fault=='missing_cif':next(folder.rglob('*.cif')).unlink()
                elif fault=='extra_sample':(folder/'seed-7_sample-1').mkdir()
                else:
                    original=next(folder.rglob('*summary_confidences.json'))
                    (original.parent/'extra_summary_confidences.json').write_bytes(original.read_bytes())
                with self.assertRaises(ValueError):predictions(folder)

    def test_ca_values_follow_column_names(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'model.cif';p.write_text(cif([90,70]))
            values=ca_plddts(p)
            self.assertEqual(values,[90,70])
            self.assertEqual(sum(values)/len(values),80)

    def test_pae_uses_token_chains_and_rejects_wrong_shape(self):
        full={'token_chain_ids':['A','A','B'],'atom_chain_ids':['A']*6+['B']*5,
              'pae':[[0,0,2],[0,0,4],[6,8,0]]}
        self.assertEqual(compute_interface_pae(full),5)
        full['pae']=[[0,2],[6,0]]
        with self.assertRaisesRegex(ValueError,'dimensions'):compute_interface_pae(full)

    def test_all_42_numbered_exports_through_cli(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);downloads=root/'inputs';downloads.mkdir()
            with (STUDY/'prediction-jobs.csv').open() as f:records=list(csv.DictReader(f))
            for record in records:
                folder=downloads/record['score_id'].lower().replace('-','_');folder.mkdir()
                for index,score in enumerate([0.1,0.6,0.8,0.8,0.4]):
                    (folder/f'fold_summary_confidences_{index}.json').write_text(json.dumps({
                        'ranking_score':score,'ptm':0.6,'iptm':0.5,'fraction_disordered':0.1,
                        'has_clash':0.0,'chain_pair_pae_min':[[0,2],[6,0]]}))
                (folder/'fold_model_2.cif').write_text(cif([80]*472))
                (folder/'fold_full_data_2.json').write_text(json.dumps({
                    'token_chain_ids':['A','A','B'],'pae':[[0,0,2],[0,0,4],[6,8,0]]}))
            out=root/'scores'
            result=subprocess.run([sys.executable,'-m','src.modeling.gprotein.parse_af3_results',
                                   '--downloads',str(downloads),'--output',str(out)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            rows=[]
            for filename in ['monomer_scores.csv','complex_scores.csv']:
                with (out/filename).open() as f:rows.extend(csv.DictReader(f))
            self.assertEqual(len(rows),42)
            self.assertTrue(all(r['best_model']=='2' and r['best_seed']=='' and r['best_sample']=='' for r in rows))
            self.assertTrue(all(r['ranking_score']=='0.8' for r in rows))

    def test_all_42_native_jobs_through_cli_and_input_identity_checks(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);downloads=root/'inputs';downloads.mkdir()
            with (STUDY/'prediction-jobs.csv').open() as f:records=list(csv.DictReader(f))
            for index,record in enumerate(records):
                request=json.loads((STUDY/record['input_path']).read_text());seed=request['modelSeeds'][0]
                native_job(downloads/f'output_{index}',request,[(seed,0,0.8211),(seed,1,0.8249)],
                           values=[80]*467+[60]*5)
            command=[sys.executable,'-m','src.modeling.gprotein.parse_af3_results','--downloads',str(downloads)]
            out=root/'scores'
            result=subprocess.run(command+['--output',str(out)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            with (out/'monomer_scores.csv').open() as f:monomers=list(csv.DictReader(f))
            with (out/'complex_scores.csv').open() as f:complexes=list(csv.DictReader(f))
            self.assertEqual((len(monomers),len(complexes)),(11,31))
            self.assertTrue(all(r['best_sample']=='1' and r['best_model']=='' for r in monomers+complexes))
            self.assertTrue(all(r['plddt_global']=='79.79' for r in monomers))
            self.assertTrue(all(r['plddt_cterminal']=='60.0' for r in monomers if r['graft_len']=='5'))
            self.assertTrue(all(r['ipae_mean']=='5.0' for r in complexes))
            self.assertIn('Ste2_Gpa1-WT_complex',{r['name'] for r in complexes})
            data_path=downloads/'output_0/test_data.json';request=json.loads(data_path.read_text())
            request['sequences'][0]['protein']['sequence']='WRONG'
            data_path.write_text(json.dumps(request))
            out=root/'invalid'
            result=subprocess.run(command+['--output',str(out)],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('sequences differ',result.stderr)
            self.assertFalse(out.exists())


if __name__=='__main__':
    unittest.main()
