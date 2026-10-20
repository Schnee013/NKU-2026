"""Input/output boundaries and representative-coordinate serialization."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from src.modeling.md_analysis.polar_io import load_config
from src.modeling.md_analysis.analyze_local_mechanism import write_representative_pdb


class PolarInputOutputTests(unittest.TestCase):
    def test_output_cannot_cover_input_or_published_material(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            case = dict(id='one', topology='data/input/solution.tpr',
                        reference='data/input/reference.pdb', trajectory='data/input/input.xtc',
                        output='runs/one')
            cfg = dict(cases=[case], output_root='runs')
            path = root/'config.json'
            with patch('src.modeling.md_analysis.polar_io.REPO_ROOT', root):
                path.write_text(json.dumps(cfg))
                self.assertEqual(load_config(path)['cases'][0]['output'], str(root/'runs/one'))
                for bad in ('data/input', 'studies/results', 'data/input/nested'):
                    case['output'] = bad
                    path.write_text(json.dumps(cfg))
                    with self.assertRaises(ValueError):
                        load_config(path)
                case['output'] = 'runs/one'
                cfg['cases'].append(dict(case))
                path.write_text(json.dumps(cfg))
                with self.assertRaises(ValueError):
                    load_config(path)

    def test_snapshot_keeps_identity_time_precision_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            template = root/'template.pdb'
            template.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n')
            coordinates = np.array([[[1.23456, -2.34567, 3.45678]]])
            output = root/'snapshot.pdb'
            write_representative_pdb(output, template, coordinates, [25.])
            lines = output.read_text().splitlines()
            atom = next(line for line in lines if line.startswith('ATOM'))
            self.assertEqual(atom[:30], template.read_text().splitlines()[0][:30])
            np.testing.assert_allclose([float(atom[n:n+8]) for n in (30,38,46)],
                                       [1.235,-2.346,3.457], atol=1e-12)
            self.assertIn('REMARK   2 TIME_NS 25', lines)
            with self.assertRaises(FileExistsError):
                write_representative_pdb(output, template, coordinates, [25.])
            with self.assertRaises(ValueError):
                write_representative_pdb(root/'invalid.pdb', template,
                                         np.full((1,1,3), np.nan), [0.])
            self.assertFalse((root/'invalid.pdb').exists())


if __name__ == '__main__':
    unittest.main()
