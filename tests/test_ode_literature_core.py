"""Check reaction mass balances, output algebra and portable CLI."""
import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
from scipy.integrate import solve_ivp

from src.modeling.ode.literature_core import LiteratureCoreModel, PROJECT_PATH, ROOT, SOURCE_PATH, RECONSTRUCTED


def full_rhs(t, y, model, ligand):
    """Full mass balances, sharing rates/input/output with the reduced model."""
    H,A,D,B,S,Q,F,Fp,Z,Zp,C = y
    p,o = model.rates,model.output
    act = model.drive(ligand)*H
    hyd = (p['k4']+p['k41']*Zp)*A
    re = p['k5']*D*B
    ass = p['k6']*S*B
    dis = (p['k7']+p['k71']*Zp)*Q
    fon = p['k8']*Q*F
    foff = (p['k9']+p['k91']*Z)*Fp
    zon = p['k10']*Fp*Z
    zoff = p['k11']*Zp
    m = max(Fp/model.pools[3],0.)
    dc = (o['b']+o['a']*m**o['n']/(o['K']**o['n']+m**o['n'])-o['delta']*C)/3600
    return np.array([-act+re,act-hyd,hyd-re,act-re-ass+dis,-ass+dis,
                     ass-dis,-fon+foff,fon-foff,-zon+zoff,zon-zoff,dc])


class LiteratureCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = LiteratureCoreModel()
        cls.times = np.unique(np.r_[0,np.geomspace(.001,900,80),np.linspace(900,28800,80)])

    def test_reduction_matches_all_full_states(self):
        for ligand in (0,1e-6,.1,1,10):
            with self.subTest(ligand=ligand):
                reduced=self.model.solve(ligand,self.times)
                full=solve_ivp(lambda t,y:full_rhs(t,y,self.model,ligand),
                               (0,self.times[-1]),self.model.reconstruct(self.model.initial),
                               t_eval=self.times,method='BDF',rtol=1e-9,atol=1e-12)
                self.assertTrue(full.success)
                np.testing.assert_allclose(self.model.reconstruct(reduced.y),full.y,atol=1e-7,rtol=0)

    def test_physical_domain_and_conservation(self):
        y=self.model.reconstruct(self.model.solve(1,self.times).y)
        self.assertGreaterEqual(float(y.min()),-1e-9)
        pools=np.array([y[0]+y[1]+y[2],y[0]+y[3]+y[5],y[4]+y[5],y[6]+y[7],y[8]+y[9]])
        np.testing.assert_allclose(pools,np.repeat(self.model.pools[:,None],len(self.times),axis=1),atol=1e-12,rtol=0)

    def test_zero_input_color_matches_hour_unit_analytical_solution(self):
        o=self.model.output
        sol=self.model.solve(0,self.times)
        expected=o['b']/o['delta']*(-np.expm1(-o['delta']*self.times/3600))
        np.testing.assert_allclose(sol.y[5],expected,atol=1e-8,rtol=0)
        np.testing.assert_allclose(sol.y[:5],
            np.repeat(self.model.initial[:5,None],len(self.times),axis=1),
            atol=1e-12,rtol=0)

    def test_output_is_exact_merge_of_components(self):
        t={'basal_activity':.01,'efficacy':.85,
           'half_activation':.25,'hill_coefficient':2.}
        g={'copy_number':1,'leak_rate':.005,'gain_rate':.4,'half_activation':.25}
        c={'production_rate_per_hour':1.,'loss_rate_per_hour':.12}
        m=np.linspace(0,1,1001)
        tf=t['basal_activity']+t['efficacy']*m**t['hill_coefficient']/(t['half_activation']**t['hill_coefficient']+m**t['hill_coefficient'])
        direct=c['production_rate_per_hour']*(g['leak_rate']+g['copy_number']*g['gain_rate']*tf/(g['half_activation']+tf))
        o=self.model.output
        merged=o['b']+o['a']*m**o['n']/(o['K']**o['n']+m**o['n'])
        np.testing.assert_allclose(merged,direct,atol=1e-14,rtol=0)
        self.assertEqual(o['delta'],c['loss_rate_per_hour'])

    def test_changed_endpoints_and_zero_loss(self):
        config=json.loads(PROJECT_PATH.read_text())
        replacements={'input':{'u_max_per_second':.002,'K_L_normalized':.3},
                      'output':{'b':.01,'a':.15,'K':.4,'n':1.5,'delta':0.}}
        for group,values in replacements.items():
            for key,value in values.items():
                config[group][key]['value']=value
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'endpoints.json';path.write_text(json.dumps(config))
            model=LiteratureCoreModel(path)
            reduced=model.solve(.2,self.times)
            full=solve_ivp(lambda t,y:full_rhs(t,y,model,.2),(0,self.times[-1]),
                           model.reconstruct(model.initial),t_eval=self.times,
                           method='BDF',rtol=1e-9,atol=1e-12)
            self.assertTrue(full.success)
            np.testing.assert_allclose(model.reconstruct(reduced.y),full.y,atol=1e-7,rtol=0)
            np.testing.assert_allclose(model.solve(0,self.times).y[5],.01*self.times/3600,atol=1e-8,rtol=0)

    def test_invalid_inputs_and_times_rejected(self):
        for ligand in (-1,np.nan,np.inf):
            with self.assertRaises(ValueError): self.model.solve(ligand,[0,1])
        for times in ([],[0],[1,2],[0,0],[0,np.nan],[[0,1],[2,3]]):
            with self.assertRaises(ValueError): self.model.solve(1,times)

    def test_invalid_parameters_and_source_rejected(self):
        base=json.loads(PROJECT_PATH.read_text())
        source=json.loads(SOURCE_PATH.read_text())
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'config.json'
            for value in (0,-1,float('nan')):
                config=copy.deepcopy(base);config['input']['K_L_normalized']['value']=value
                path.write_text(json.dumps(config))
                with self.assertRaises(ValueError): LiteratureCoreModel(path)
            config=copy.deepcopy(base);config['output']['b']['source']=''
            path.write_text(json.dumps(config))
            with self.assertRaises(ValueError): LiteratureCoreModel(path)
            source['time_unit']='h';path.write_text(json.dumps(source))
            with self.assertRaises(ValueError): LiteratureCoreModel(source_path=path)

    def test_cli_works_outside_repo_and_records_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'run'
            result=subprocess.run([sys.executable,str(ROOT/'src/modeling/ode/run_literature_core.py'),
                                   '--hours','.1','--ligand','0','.1','--output',str(dest)],
                                  cwd=tmp,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            record=json.loads((dest/'run.json').read_text())
            self.assertEqual(record['hours'],.1)
            self.assertEqual(record['ligand_normalized'],[0.,.1])
            self.assertEqual(record['resolved_endpoint_parameters'],self.model.parameters)
            self.assertEqual(record['fixed_source_rates'],self.model.rates)
            self.assertEqual(record['fixed_source_pools_uM'],self.model.pools.tolist())
            with (dest/'timecourses.csv').open() as f: rows=list(csv.DictReader(f))
            self.assertTrue(rows)
            self.assertEqual({row['ligand_normalized'] for row in rows},{'0.0','0.1'})

    def test_repeat_run_preserves_existing_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'run'
            command=[sys.executable,str(ROOT/'src/modeling/ode/run_literature_core.py'),
                     '--output',str(dest)]
            first=subprocess.run(command,cwd=tmp,capture_output=True,text=True)
            self.assertEqual(first.returncode,0,first.stderr)
            plot=subprocess.run([sys.executable,str(ROOT/'src/modeling/ode/plot_literature_core.py'),
                                 str(dest/'timecourses.csv'),'--output',str(dest/'timecourses.png')],cwd=tmp,capture_output=True,text=True)
            self.assertEqual(plot.returncode,0,plot.stderr)
            before={p.name:p.read_bytes() for p in dest.iterdir()}
            repeat=subprocess.run(command+['--hours','.1','--ligand','0'],
                                  cwd=tmp,capture_output=True,text=True)
            self.assertNotEqual(repeat.returncode,0)
            self.assertIn('choose a new directory',repeat.stderr)
            self.assertEqual(before,{p.name:p.read_bytes() for p in dest.iterdir()})

    def test_existing_empty_directory_or_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            for kind in ('directory','file'):
                dest=Path(tmp)/kind
                if kind=='directory': dest.mkdir()
                else: dest.write_text('preserve me')
                result=subprocess.run([sys.executable,str(ROOT/'src/modeling/ode/run_literature_core.py'),
                                       '--output',str(dest)],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('already exists',result.stderr)
                if kind=='directory': self.assertEqual(list(dest.iterdir()),[])
                else: self.assertEqual(dest.read_text(),'preserve me')

    def test_runner_protects_inputs_and_references(self):
        protected=ROOT/'studies/signal-model/pathway'
        for dest in (protected/'reference'/'should-not-be-created',
                     protected/'configs'/'should-not-be-created'):
            result=subprocess.run([sys.executable,str(ROOT/'src/modeling/ode/run_literature_core.py'),
                                   '--output',str(dest)],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('outside input and reference',result.stderr)
            self.assertFalse(dest.exists())

    def test_missing_input_does_not_create_a_successful_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'run'
            result=subprocess.run([sys.executable,str(ROOT/'src/modeling/ode/run_literature_core.py'),
                                   '--parameters',str(Path(tmp)/'missing.json'),
                                   '--output',str(dest)],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse(dest.exists())

    def test_failed_solver_leaves_no_result_or_success_metadata(self):
        from src.modeling.ode import run_literature_core
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'run'
            with mock.patch.object(sys,'argv',['run_literature_core','--output',str(dest)]), \
                 mock.patch.object(LiteratureCoreModel,'solve',side_effect=RuntimeError('interrupted solve')):
                with self.assertRaisesRegex(RuntimeError,'interrupted solve'):
                    run_literature_core.main()
            self.assertTrue(dest.is_dir())
            self.assertFalse((dest/'timecourses.csv').exists())
            self.assertFalse((dest/'run.json').exists())

    def test_plot_rejects_study_and_existing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv_path=Path(tmp)/'curves.csv'
            csv_path.write_text('ligand_normalized,time_s,color_au\n0,0,0\n')
            existing=Path(tmp)/'existing.png';existing.write_bytes(b'preserve')
            protected=ROOT/'studies/signal-model/pathway/should-not-be-created.png'
            for dest in (existing,protected):
                result=subprocess.run([sys.executable,str(ROOT/'src/modeling/ode/plot_literature_core.py'),str(csv_path),'--output',str(dest)],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                if dest==existing:self.assertEqual(dest.read_bytes(),b'preserve')
                else:self.assertFalse(dest.exists())

    def test_generated_zero_input_curves_match_analytic_solution(self):
        for name,config in [('default',PROJECT_PATH),('alternate',PROJECT_PATH.with_name('alternate_project.json'))]:
            with self.subTest(name=name),tempfile.TemporaryDirectory() as tmp:
                output=Path(tmp)/'curves'
                result=subprocess.run([sys.executable,str(ROOT/'src/modeling/ode/run_literature_core.py'),
                    '--parameters',str(config),'--ligand','0','--output',str(output)],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                data=np.genfromtxt(output/'timecourses.csv',delimiter=',',names=True)
                model=LiteratureCoreModel(config)
                parameters=json.loads(config.read_text())['output']
            b=parameters['b']['value'];delta=parameters['delta']['value']
            expected=b*data['time_s']/3600 if delta==0 else b/delta*(-np.expm1(-delta*data['time_s']/3600))
            self.assertEqual(len(data),400)
            np.testing.assert_allclose(data['color_au'],expected,atol=1e-8,rtol=0)
            np.testing.assert_allclose(data['MAPK_fraction'],0,atol=1e-14,rtol=0)

    def test_source_parameter_correspondence(self):
        # Supady Model B (PDF p49) and Delft Tables 4.1–4.2 define these
        # identifiers, reactions and rates independently of the configuration.
        expected = (
            ('K1', 'k1', 'receptor activation', 0.071, 'uM^-1 s^-1'),
            ('K2', 'k2', 'receptor inactivation', 2.457, 's^-1'),
            ('K3', 'k12', 'active receptor loss', 1.199, 's^-1'),
            ('K4', 'k3', 'G-protein activation', 0.724, 'uM^-1 s^-1'),
            ('K5', 'k4', 'GTP hydrolysis', 1.259, 's^-1'),
            ('K6', 'k41', 'active Sst2-assisted hydrolysis', 1.124, 'uM^-1 s^-1'),
            ('K7', 'k5', 'G-protein reassociation', 0.002, 'uM^-1 s^-1'),
            ('K8', 'k6', 'scaffold assembly', 0.021, 'uM^-1 s^-1'),
            ('K9', 'k7', 'scaffold disassembly', 0.004, 's^-1'),
            ('K10', 'k71', 'active Sst2-assisted scaffold disassembly', 0.261, 'uM^-1 s^-1'),
            ('K11', 'k8', 'Fus3 activation', 0.818, 'uM^-1 s^-1'),
            ('K12', 'k91', 'Sst2-dependent Fus3 inactivation', 0.017, 'uM^-1 s^-1'),
            ('K13', 'k10', 'Sst2-state conversion', 0.342, 'uM^-1 s^-1'),
            ('K14', 'k11', 'Sst2-state recovery', 0.643, 's^-1'),
        )
        source = self.model.source
        parameters = source['parameters']
        self.assertEqual(set(parameters), {row[1] for row in expected} | {'k9'})
        self.assertEqual(source['evidence_status_default'], 'source_fitted')
        for delft_id, name, reaction, value, unit in expected:
            parameter = parameters[name]
            self.assertEqual(parameter['delft_id'], delft_id)
            self.assertEqual(parameter['reaction'], reaction)
            self.assertEqual(parameter['value'], value)
            self.assertEqual(parameter['unit'], unit)
            self.assertEqual(parameter['source_pointer'],
                             f'Supady 2009 PDF p49, Model B; Delft Table4.2 {delft_id}')
            self.assertEqual(parameter.get('evidence_status', source['evidence_status_default']),
                             'source_fitted')
        self.assertIn('Sst2 and active Fus3', parameters['k91']['unit_note'])
        self.assertEqual(parameters['k9']['value'], 0.)
        self.assertEqual(parameters['k9']['unit'], 's^-1')
        self.assertEqual(parameters['k9']['reaction'], 'basal Fus3 inactivation')
        self.assertEqual(parameters['k9']['evidence_status'], 'source_fitted_zero')
        self.assertEqual(parameters['k9']['source_pointer'], 'Supady 2009 PDF p49, Model B k9')
        self.assertNotIn('delft_id', parameters['k9'])
        self.assertEqual(self.model.rates['k9'], 0.)

if __name__=='__main__':
    unittest.main()
