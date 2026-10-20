"""Input qualification, complete comparison sets, and protected outputs."""
import csv
from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from src.modeling.dual_channel.calibration import fit, predict, eco_scenarios
from src.modeling.dual_channel.replay import PROTOCOL_ID, read_pairs, run, serialize_predictions
from src.modeling.economic.eco_model_pkg.config import load_model

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'studies/signal-model/dual-channel'
ECO = ROOT / 'studies/economic-model/configs/eco_model_example_model.json'


class ContractsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from src.modeling.dual_channel.generate_inputs import generate
        cls.temporary=tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.inputs=Path(cls.temporary.name)/'generated'
        generate(FIXTURE/'inputs',cls.inputs,ROOT/'studies/signal-model/pathway/configs/literature_core_source.json')

    def setUp(self):
        rows = read_pairs(self.inputs / 'calibration.csv')
        self.model = fit([[r['x3'], r['x18']] for r in rows], [r['target'] for r in rows],
                         sample_ids=[r['sample_id'] for r in rows], protocol_id=PROTOCOL_ID,
                         evidence_status='simulated', reference_basis='spike_count')

    def test_unknown_qc_and_saturation_are_rejected(self):
        prediction = predict(self.model, [.02, .003], protocol_id=PROTOCOL_ID)
        self.assertIn('qc_status_unknown', prediction['reasons'])
        self.assertIn('saturation_status_unknown', prediction['reasons'])
        self.assertFalse(prediction['within_calibration_range'])
        prediction = predict(self.model, [.02, .003], protocol_id=PROTOCOL_ID,
                             saturated='false', qc_pass='pass')
        self.assertIn('qc_status_unknown', prediction['reasons'])

    def test_regression_matches_independent_normal_equations(self):
        rows=read_pairs(self.inputs/'calibration.csv')
        x=np.array([[1.,r['x3'],r['x18']] for r in rows]);y=np.array([r['target'] for r in rows])
        beta=np.linalg.solve(x.T@x,x.T@y)
        variance=float(np.sum((y-x@beta)**2)/(len(y)-3))
        covariance=variance*np.linalg.inv(x.T@x)
        np.testing.assert_allclose(self.model.coefficients,beta,rtol=1e-7,atol=1e-9)
        np.testing.assert_allclose(self.model.covariance,covariance,rtol=1e-7,atol=1e-9)
        self.assertAlmostEqual(self.model.residual_variance,variance,places=7)
        self.assertEqual(self.model.df,len(y)-3)
        self.assertEqual(list(self.model.density_range),[0.,120.])
        self.assertEqual(self.model.evidence_status,'simulated')
        equations=self.model.hull_equations
        self.assertTrue(np.all(x[:,1:]@equations[:,:2].T+equations[:,2]<=1e-10))

    def test_declared_reversed_channel_order_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'channel order'):
            predict(self.model, [.003, .02], protocol_id=PROTOCOL_ID,
                    saturated=False, qc_pass=True, channel_order=('x18', 'x3'))

    def test_unit_mismatch_and_nonnegative_density_range_are_rejected(self):
        row = read_pairs(self.inputs / 'calibration.csv')[4]
        signal = [row['x3'], row['x18']]
        prediction = predict(self.model, signal, protocol_id=PROTOCOL_ID, saturated=False, qc_pass=True)
        self.assertTrue(prediction['within_calibration_range'])
        prediction['density_unit'] = 'J2/g soil'
        self.assertIn('unit_mismatch', eco_scenarios(prediction, load_model(ECO))['reasons'])
        limited = replace(self.model, density_range=(0., 1.))
        prediction = predict(limited, signal, protocol_id=PROTOCOL_ID, saturated=False, qc_pass=True)
        self.assertGreater(prediction['density'], 1.)
        self.assertIn('outside_density_calibration', prediction['reasons'])

    def test_missing_columns_and_missing_file_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            malformed = Path(temporary) / 'bad.csv'
            malformed.write_text('sample_id,target,background,C3,x18,data_origin\na,0,0,0,0,synthetic_ODE\n')
            with self.assertRaisesRegex(ValueError, 'missing paired columns'):
                read_pairs(malformed)
            with self.assertRaises(FileNotFoundError):
                read_pairs(Path(temporary) / 'absent.csv')

    def test_zero_density_is_explicit_and_missing_is_rejected(self):
        rows = read_pairs(self.inputs / 'calibration.csv')
        self.assertEqual(rows[0]['target'], 0.)
        self.assertEqual(rows[0]['background'], 0.)
        with tempfile.TemporaryDirectory() as temporary:
            missing = Path(temporary) / 'missing.csv'
            missing.write_text('sample_id,target,background,x3,x18,data_origin\na,,0,0,0,synthetic_ODE\n')
            with self.assertRaises(ValueError):
                read_pairs(missing)

    def test_endpoint_units_and_matched_blank(self):
        with (self.inputs / 'ode-endpoints.csv').open() as handle:
            endpoints = list(csv.DictReader(handle))
        paired = {r['sample_id']: r for r in read_pairs(self.inputs / 'calibration.csv')}
        self.assertEqual(len(endpoints), 18)
        for row in endpoints:
            pigment = float(row['ode_pigment_uM'])
            corrected = float(row['ode_A480_blank_corrected'])
            np.testing.assert_allclose(float(row['ode_A480']), .048 * pigment, rtol=1e-10, atol=1e-12)
            np.testing.assert_allclose(corrected, .048 * (pigment - .5), rtol=1e-7, atol=1e-9)
            self.assertEqual(paired[row['sample_id']]['x' + row['channel']], corrected)

    def test_full_replay_metrics_counts_and_integration(self):
        with tempfile.TemporaryDirectory() as temporary:
            output=Path(temporary)/'replay'
            summary=run(self.inputs,output,ECO,PROTOCOL_ID)
            batch=json.loads((output/'predictions.json').read_text())
            self.assertEqual(batch['metadata'],dict(protocol_id=PROTOCOL_ID,evidence_status='simulated',
                reference_basis='spike_count',calibration_density_range=[0.,120.],density_unit='J2/100 cm^3 soil'))
            calibration=read_pairs(self.inputs/'calibration.csv')
            evaluation=read_pairs(self.inputs/'evaluation.csv')
            self.assertEqual([r['sample_id'] for r in batch['samples']],[r['sample_id'] for r in evaluation])
            for name,columns in [('channel3',['x3']),('channel18',['x18']),('dual',['x3','x18'])]:
                x=np.array([[1.,*[r[k] for k in columns]] for r in calibration])
                y=np.array([r['target'] for r in calibration])
                beta=np.linalg.solve(x.T@x,x.T@y)
                expected=np.array([[1.,*[r[k] for k in columns]] for r in evaluation])@beta
                np.testing.assert_allclose([r['estimates'][name] for r in batch['samples']],expected,rtol=1e-7,atol=1e-7)
                for ratio,metrics in summary['scenario_metrics'].items():
                    indices=[i for i,r in enumerate(evaluation) if r['mi_ascr3_to_ascr18']==ratio]
                    residual=expected[indices]-np.array([evaluation[i]['target'] for i in indices])
                    self.assertEqual(metrics[name]['n'],12)
                    self.assertAlmostEqual(metrics[name]['rmse'],float(np.sqrt(np.mean(residual**2))),places=6)
                    self.assertAlmostEqual(metrics[name]['mae'],float(np.mean(abs(residual))),places=6)
            rows=batch['samples']
            self.assertEqual(summary['counts'],dict(evaluation=36,valid_point=30,eco_conditional=30))
            self.assertEqual(sum(r['estimates']['dual']<0 for r in rows),6)
            for row in rows:
                if row['estimates']['dual']<0:
                    self.assertIn('negative_density',row['prediction']['reasons'])
                    self.assertEqual(row['eco']['status'],'blocked')
                if row['eco']['status']=='conditional_demo_only':self.assertTrue(row['eco']['scenarios'])
            before={p.name:p.read_bytes() for p in output.iterdir()}
            with self.assertRaises(FileExistsError):run(self.inputs,output,ECO,PROTOCOL_ID)
            self.assertEqual(before,{p.name:p.read_bytes() for p in output.iterdir()})

    def test_batch_serialization_preserves_complete_api_outputs(self):
        row = read_pairs(self.inputs / 'calibration.csv')[4]
        prediction = predict(self.model, [row['x3'], row['x18']], protocol_id=PROTOCOL_ID,
                             saturated=False, qc_pass=True)
        eco = eco_scenarios(prediction, load_model(ECO))
        self.assertEqual(eco['status'], 'conditional_demo_only')
        self.assertEqual(prediction['protocol_id'], PROTOCOL_ID)
        self.assertEqual(prediction['calibration_density_range'], [0., 120.])
        for key, expected in (('evidence_status', 'simulated'), ('reference_basis', 'spike_count')):
            self.assertEqual(prediction[key], expected)
            self.assertEqual(eco[key], expected)
        self.assertEqual(prediction['density_unit'], 'J2/100 cm^3 soil')
        self.assertTrue(all(item['density_unit'] == prediction['density_unit'] for item in eco['scenarios']))
        results = [dict(sample_id=row['sample_id'], prediction=prediction, eco=eco)]
        before = deepcopy(results)
        batch = serialize_predictions(results)
        self.assertEqual(results, before)
        self.assertEqual(batch['samples'][0]['prediction']['prediction_interval_95'],
                         prediction['prediction_interval_95'])
        self.assertEqual(batch['samples'][0]['eco']['physical_interval'], eco['physical_interval'])
        self.assertEqual(batch['samples'][0]['eco']['interval_truncated_at_zero'],
                         eco['interval_truncated_at_zero'])
        self.assertFalse(set(batch['metadata']) & set(batch['samples'][0]['prediction']))
        self.assertTrue(all('density_unit' not in item for item in batch['samples'][0]['eco']['scenarios']))
        for key, replacement in (('protocol_id', 'other'), ('evidence_status', 'measured_unvalidated'),
                                 ('reference_basis', 'independent_field_count'),
                                 ('calibration_density_range', [0., 80.]), ('density_unit', 'J2/g soil')):
            mixed = deepcopy(results) + deepcopy(results)
            mixed[1]['prediction'][key] = replacement
            with self.assertRaisesRegex(ValueError, 'inconsistent prediction metadata'):
                serialize_predictions(mixed)
        mixed = deepcopy(results)
        mixed[0]['eco']['scenarios'][0]['density_unit'] = 'J2/g soil'
        with self.assertRaisesRegex(ValueError, 'inconsistent economic density unit'):
            serialize_predictions(mixed)
        mixed = deepcopy(results)
        mixed[0]['eco']['reference_basis'] = 'independent_field_count'
        with self.assertRaisesRegex(ValueError, 'inconsistent economic metadata'):
            serialize_predictions(mixed)

    def test_output_protection_and_failure_does_not_publish_partial_results(self):
        with self.assertRaises(ValueError):
            run(self.inputs, FIXTURE / 'reference/new-output', ECO, PROTOCOL_ID)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'interrupted'
            with patch('src.modeling.dual_channel.replay.plot', side_effect=RuntimeError('injected interruption')):
                with self.assertRaisesRegex(RuntimeError, 'injected interruption'):
                    run(self.inputs, output, ECO, PROTOCOL_ID)
            self.assertFalse(output.exists())
            self.assertEqual(list(Path(temporary).iterdir()), [])


if __name__ == '__main__':
    unittest.main()
