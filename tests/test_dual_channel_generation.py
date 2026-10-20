"""Check paired-channel input generation and regression from model configurations."""
import csv
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np

from src.modeling.dual_channel.generate_inputs import generate
from src.modeling.dual_channel.replay import run

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT/'studies/signal-model/dual-channel'


class GenerationTests(unittest.TestCase):
    def test_parameters_generate_disjoint_synthetic_inputs_and_regression(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inputs=root/'parameters';inputs.mkdir()
            for name in ('scenario.json','generation-design.json','channel3.json','channel18.json'):
                shutil.copy2(FIXTURES/'inputs'/name,inputs/name)
            output=root/'generated'
            record=generate(inputs,output,ROOT/'studies/signal-model/pathway/configs/literature_core_source.json')
            self.assertEqual({n:v['rows'] for n,v in record['outputs'].items()},
                             {'ode-endpoints.csv':18,'calibration.csv':9,'evaluation.csv':36})
            tables={}
            for name in record['outputs']:
                with (output/name).open() as handle:tables[name]=list(csv.DictReader(handle))
                self.assertTrue(all(r['data_origin']=='synthetic_ODE' for r in tables[name]))
            calibration=tables['calibration.csv'];evaluation=tables['evaluation.csv']
            self.assertFalse({r['sample_id'] for r in calibration}&{r['sample_id'] for r in evaluation})
            self.assertFalse({(float(r['target']),float(r['background'])) for r in calibration}&
                             {(float(r['target']),float(r['background'])) for r in evaluation})
            blank=calibration[0]
            np.testing.assert_allclose([float(blank['x3']),float(blank['x18'])],[0,0],atol=1e-12)
            before={p.name:p.read_bytes() for p in output.iterdir()}
            with self.assertRaises(FileExistsError):generate(inputs,output,ROOT/'studies/signal-model/pathway/configs/literature_core_source.json')
            self.assertEqual(before,{p.name:p.read_bytes() for p in output.iterdir()})
            summary=run(output,root/'replay',ROOT/'studies/economic-model/configs/eco_model_example_model.json',
                        'SIM_V4_MI1TO10_CE10TO1_Q1_8H_A480')
            self.assertEqual(summary['evidence_status'],'simulated')
            self.assertEqual(summary['counts']['evaluation'],36)
            for ratio in summary['scenario_metrics']:
                for model in ('channel3','channel18','dual'):
                    self.assertTrue(np.isfinite(summary['scenario_metrics'][ratio][model]['rmse']))
