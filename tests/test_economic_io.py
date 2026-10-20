"""Public economic references, independent accounting and output protection."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import io
from contextlib import redirect_stderr

from src.modeling.economic.eco_model_pkg import load_model, evaluate
from src.modeling.economic.run_mvp import clean_rows
from src.modeling.economic.economic_model import evaluate_config

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "studies/economic-model"


class EconomicPublicIOTests(unittest.TestCase):
    def run_cli(self, args, cwd):
        env = dict(os.environ, PYTHONPATH=str(ROOT), PYTHONDONTWRITEBYTECODE="1")
        return subprocess.run([sys.executable, "-m", *args], cwd=cwd, env=env,
                              capture_output=True, text=True)

    def test_independent_comparison_formulas(self):
        model = load_model(STUDY / "configs/eco_model_example_model.json")
        # Independent source equation and explicit cash-flow algebra. Tiers are
        # checked at boundaries, including cost discontinuities.
        for density in [0.0, 6.0, 6.8, 10.0, 79.999999, 80.0, 100.0, 869.6]:
            with self.subTest(density=density):
                result = evaluate(density, model)
                relative = 1 if density <= 6 else .34 + .66 * .95 ** (density / 6 - 1)
                revenue = 113500 * .690044053
                loss = revenue * (1 - relative)
                cost, effect = (0, 0) if density < 10 else (220, .586805) if density < 80 else (2700, .691234)
                self.assertAlmostEqual(result.program_net_benefit_vs_no_control, loss * effect - cost - 5, places=8)
                self.assertAlmostEqual(result.detection_net_benefit_vs_uniform, loss * (effect - .586805) - cost - 5 + 220, places=8)

    def test_framework_cli_writes_complete_finite_scenario_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'new/results.csv'
            result=self.run_cli(['src.modeling.economic.run_mvp','--out',str(output)],tmp)
            self.assertEqual(result.returncode,0,result.stderr)
            with output.open() as handle:rows=list(csv.DictReader(handle))
            self.assertEqual(len(rows),21)
            keys={(r['scenario_id'],r['intervention_effect'],r['early_program_cost_per_area']) for r in rows}
            self.assertEqual(len(keys),21)
            for row in rows:
                for value in row.values():
                    try:number=float(value)
                    except ValueError:continue
                    self.assertTrue(math.isfinite(number))

    def test_cli_outside_checkout_and_nonempty_output_rejected(self):
        with tempfile.TemporaryDirectory(prefix="economics external ") as tmp:
            out = Path(tmp) / 'report with spaces'
            args = ["src.modeling.economic.eco_model", "100", "--report-dir", str(out)]
            first = self.run_cli(args, tmp)
            self.assertEqual(first.returncode, 0, first.stderr)
            before = {p.name: p.read_bytes() for p in out.iterdir()}
            second = self.run_cli(args, tmp)
            self.assertNotEqual(second.returncode, 0)
            self.assertIn("absent or empty", second.stderr)
            self.assertEqual(before, {p.name: p.read_bytes() for p in out.iterdir()})
            report = json.loads((out / 'diagnostic_report.json').read_text())
            self.assertAlmostEqual(report['diagnosis']['detection_net_benefit_vs_uniform'], 496.2593745327813)

    def test_reference_output_path_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.run_cli(["src.modeling.economic.eco_model", "100", "--report-dir", str(STUDY / 'reference')], tmp)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn('outside the packaged', run.stderr)

    def test_failed_input_has_no_success_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'failed-output'
            for extra in [["-1"], ["100", "--config", str(Path(tmp) / 'missing.json')]]:
                run = self.run_cli(["src.modeling.economic.eco_model", *extra, "--report-dir", str(out)], tmp)
                self.assertNotEqual(run.returncode, 0)
                self.assertFalse(out.exists())

    def test_report_render_failure_leaves_no_partial_final_report(self):
        from src.modeling.economic.eco_model_pkg.cli import main
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'report'
            def fail_after_partial(report, directory):
                (Path(directory) / 'diagnostic_report.txt').write_text('partial')
                raise OSError('injected renderer failure')
            with patch('sys.argv', ['economic', '100', '--report-dir', str(out)]), patch(
                'src.modeling.economic.eco_model_pkg.cli.write_diagnostic_report', side_effect=fail_after_partial
            ), redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main()
            self.assertEqual(list(out.iterdir()), [])
            self.assertEqual(list(Path(tmp).glob('.economic-report-*')), [])

    def test_framework_cli_refuses_nonempty_parent(self):
        with tempfile.TemporaryDirectory() as tmp:
            sentinel = Path(tmp) / 'unrelated.txt'
            sentinel.write_text('keep')
            out = Path(tmp) / 'framework.csv'
            run = self.run_cli(["src.modeling.economic.run_mvp", "--out", str(out)], tmp)
            self.assertNotEqual(run.returncode, 0)
            self.assertFalse(out.exists())
            self.assertEqual(sentinel.read_text(), 'keep')


if __name__ == '__main__':
    unittest.main()
