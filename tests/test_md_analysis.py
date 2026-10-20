"""Known-geometry checks for each supplied basic MD driver and local primitives."""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import MDAnalysis as mda
from examples.md_geometry import create_inputs, run_drivers
from src.modeling.md_analysis.common import equal_weight_kabsch_coordinates, ensure_dir
from src.modeling.md_analysis.independent_geometry import independent_fit
from src.modeling.md_analysis.verify_tables import load_case_summary
from src.modeling.md_analysis.analyze_local_mechanism import first_run, stretches, summary, validate_coordinates
from src.modeling.md_analysis.analyze_hbond import summarize_events, aggregate_pair_occupancy, scan_direction
from src.modeling.md_analysis.compare_ligand_paths import (
    arm_inputs, parse_named_values, summarize_contact_occupancy, validate_arm_mappings,
)
from src.modeling.md_analysis.common import analysis_context

ROOT=Path(__file__).resolve().parents[1]




class MDGeometryTests(unittest.TestCase):
    def test_rigid_fit_and_reflection(self):
        ref=np.array([[0.,0,0],[1,0,0],[0,2,0],[0,0,3]])
        R=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        mobile=ref@R+np.array([3,5,7])
        fit,rmsd=equal_weight_kabsch_coordinates(mobile,ref,mobile)
        np.testing.assert_allclose(fit,ref,atol=1e-10)
        np.testing.assert_allclose(independent_fit(mobile,ref,mobile),ref,atol=1e-10)
        self.assertLess(rmsd,1e-7)
        reflected=ref*np.array([-1,1,1])
        _,value=equal_weight_kabsch_coordinates(reflected,ref,reflected)
        self.assertGreater(value,.1)
        with self.assertRaises(ValueError):equal_weight_kabsch_coordinates(ref[:2],ref,ref)

    def test_local_event_and_selection_rules(self):
        times=np.arange(21)*.1
        self.assertEqual(first_run(np.ones(21,dtype=bool),times),0.)
        self.assertIsNone(first_run(np.arange(21)%2==0,times))
        self.assertEqual(first_run(np.arange(21)>=5,times),.5)
        self.assertEqual(stretches({1,2,3,5,6,8},3),[1,2,3])
        self.assertEqual(summary(np.arange(21),times,[1,2])['late_mean'],15.)

    def test_coordinate_shape_and_finite_values(self):
        valid = np.zeros((2, 3))
        np.testing.assert_array_equal(validate_coordinates(valid, 2), valid)
        for invalid in (np.ones((3, 3)), np.full((2, 3), np.nan)):
            with self.assertRaises(ValueError):
                validate_coordinates(invalid, 2)

    def test_hbond_denominators_and_replicate_aggregation(self):
        events = [dict(rep='rep1', direction='forward', pair='donor-acceptor',
                       distance_A=2.8, angle_deg=160., time_ns=50.)]
        _, table = summarize_events(events, {('rep1', 'forward'): 4})
        self.assertEqual(table.iloc[0].occupancy, .25)
        one = aggregate_pair_occupancy(table, ['rep1'])
        two = aggregate_pair_occupancy(table, ['rep1', 'rep2'])
        self.assertEqual(one.iloc[0].occupancy, .25)
        self.assertEqual(two.iloc[0].occupancy, .125)
        self.assertEqual(two.iloc[0].configured_reps, 2)

    def test_contact_window_keeps_early_and_excludes_later_frames(self):
        times = np.array([0., 49.5, 50., 100., 100.5])
        residue = (84, 'ASN')
        rows = summarize_contact_occupancy(
            'WT', times, [{residue}, set(), {residue}, set(), {residue}],
            50., 100., 4.,
        )
        self.assertEqual(rows[0]['early_occupancy'], .5)
        self.assertEqual(rows[0]['late_occupancy'], .5)
        self.assertEqual(rows[0]['early_frames'], 2)
        self.assertEqual(rows[0]['late_frames'], 2)
        self.assertTrue(rows[0]['initial_contact'])

    def test_arm_system_overrides_global_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / 'input'
            processed = payload / 'processed/rep1'
            processed.mkdir(parents=True)
            (processed / 'protein_ref.pdb').touch()
            (processed / 'prod_protein.xtc').touch()
            config = root / 'systems.json'
            config.write_text(json.dumps({'default': 'WT', 'systems': {'WT': {}, 'N84V': {}}}))
            systems = parse_named_values(['wt=WT', 'mutant=N84V'], '--arm-system')
            with patch.dict(os.environ, {'MD_SYSTEM': 'WT', 'MD_OUTPUT_ROOT': str(root / 'output')}):
                context, _, _ = arm_inputs(config, payload, 1, systems['mutant'])
                self.assertEqual(context.system, 'N84V')
                self.assertEqual(analysis_context(config_path=config, payload_base=payload).system, 'WT')
            with self.assertRaises(ValueError):
                validate_arm_mappings({'wt': config}, {'wt': payload}, systems, 'wt')

    def test_hbond_empty_candidates_keep_sampled_denominator(self):
        universe = mda.Universe.empty(1, trajectory=True)
        universe.add_TopologyAttr('names', ['C'])
        universe.add_TopologyAttr('resids', [1])
        universe.add_TopologyAttr('resnames', ['ALA'])
        universe.add_TopologyAttr('segids', ['A'])
        universe.add_TopologyAttr('elements', ['C'])
        universe.load_new(np.zeros((3, 1, 3)), order='fac', dt=500.)
        events, count = scan_direction(
            universe, universe.atoms, universe.atoms, 'rep1', 'forward', 1, (0., 1.),
        )
        self.assertEqual(events, [])
        self.assertEqual(count, 3)

    def test_combined_summary_requires_complete_unique_scientific_keys(self):
        header = 'case,site,neighbor,late_mean\n'
        rows = ['WT,84,171,0\n', 'WT,368,70,2.3\n',
                'N84V,84,171,4.2\n', 'N84V,368,70,0\n']
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'pairs.csv'
            def read(body):
                path.write_text(header + ''.join(body))
                return load_case_summary(
                    path, ['WT', 'N84V'], ['site', 'neighbor'], [(84, 171), (368, 70)],
                )
            table = read(rows)
            self.assertEqual(len(table), 4)
            self.assertEqual(table.iloc[0].late_mean, 0.)
            for invalid in (
                rows + [rows[0]],  # duplicate scientific key
                rows[:2],  # missing case
                rows[:-1],  # missing pair within a present case
                rows[:3] + ['N84V,368,77,0\n'],  # unconfigured pair
                rows[:3] + ['N84V,368,,0\n'],  # missing identity
            ):
                with self.subTest(rows=invalid), self.assertRaises(ValueError):
                    read(invalid)

    def test_output_protection(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'result';ensure_dir(p);(p/'old').write_text('old')
            with self.assertRaises(FileExistsError):ensure_dir(p)
        with self.assertRaises(ValueError):ensure_dir(ROOT/'studies/forbidden')

    def test_seven_basic_drivers_on_known_geometry(self):
        with tempfile.TemporaryDirectory(prefix='md geometry ') as tmp:
            p=Path(tmp)
            config=create_inputs(p/'input')
            run_drivers(config,p/'output',p/'mpl')
            out=p/'output/synthetic/analysis'
            q=pd.read_csv(out/'selection_rmsd/selection_rmsd.csv')
            np.testing.assert_allclose(q.selection_rmsd_A,[0,1,2],atol=1e-4)
            q=pd.read_csv(out/'convergence/rmsd_summary.csv');self.assertLess(abs(q.backbone_rmsd_A).max(),1e-6)
            q=pd.read_csv(out/'rmsf/rmsf.csv');self.assertLess(abs(q.rmsf_A).max(),1e-6)
            q=pd.read_csv(out/'contacts/contact_persistence.csv');np.testing.assert_allclose(q.retained_fraction,[1,1,0])
            hb=pd.read_csv(out/'hbond/hbond_summary.csv')
            self.assertTrue(((hb.events==1) & np.isclose(hb.occupancy,1/3)).any())
            com=pd.read_csv(out/'com_distance/com.csv');self.assertTrue(np.all(np.diff(com.distance_A)>0))
            bs=pd.read_csv(out/'bsa/bsa.csv');self.assertTrue((bs.bsa_A2>=-1e-3).all());self.assertGreater(bs.bsa_A2.iloc[0],bs.bsa_A2.iloc[-1])
            # Refuse an existing result without overwriting it.
            before={str(f.relative_to(out)):f.read_bytes() for f in out.rglob('*') if f.is_file()}
            with self.assertRaises(subprocess.CalledProcessError):run_drivers(config,p/'output',p/'mpl')
            self.assertEqual(before,{str(f.relative_to(out)):f.read_bytes() for f in out.rglob('*') if f.is_file()})


if __name__=='__main__':unittest.main()
