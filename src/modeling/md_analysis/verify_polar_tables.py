"""Numerical and identity checks of saved polar-contact materials."""
import json
import numpy as np
import pandas as pd


def verify_polar_tables(study):
    count = 0

    def check(ok, message):
        nonlocal count
        if not bool(ok):
            raise ValueError(message)
        count += 1

    def close(a, b, message, atol=1e-10):
        check(np.allclose(a, b, atol=atol, rtol=0), message)

    cfg = json.loads((study/'configs/polar-anchors.json').read_text())
    summary = pd.read_csv(study/'results/polar-summary/all_windows_polar_summary.csv')
    keys = ['case_id', 'feature', 'metric', 'window']
    check(not summary.duplicated(keys).any(), 'Duplicate polar summary identity')
    check(set(summary.case_id) == {c['id'] for c in cfg['cases']}, 'Incomplete polar cohort')
    primary = pd.read_csv(study/'results/polar-summary/primary_polar_summary.csv')
    pd.testing.assert_frame_equal(primary.reset_index(drop=True),
                                  summary[summary.window == 'primary'].reset_index(drop=True))
    for case in cfg['cases']:
        cid = case['id']
        root = study/'results/polar'/cid/'polar'
        series = pd.read_csv(root/'polar_timeseries.csv.gz')
        check(len(series) == 10001, cid+': polar frame count')
        times = series.time_ns.to_numpy()
        close(times, np.arange(10001)*.01, cid+': polar time grid', 1e-7)
        definitions = json.loads((root/'polar_atom_definitions.json').read_text())
        expected = {item['feature']: item['status'] for item in definitions}
        check(len(summary[summary.case_id == cid]) == len(expected)*4*12,
              cid+': complete targeted metrics and windows')
        for row in summary[summary.case_id == cid].itertuples():
            check(row.status == expected[row.feature], cid+': applicability differs')
            if row.status != 'APPLICABLE':
                check(pd.isna(row.occupancy), cid+': inapplicable contact must remain missing')
                continue
            mask = series[row.feature+'__'+row.metric].to_numpy(dtype=bool)
            window = (times >= row.start_ns) & ((times <= row.end_ns) if row.end_inclusive else (times < row.end_ns))
            check(int(window.sum()) == row.n_frames, cid+': polar denominator')
            check(int(mask[window].sum()) == row.present_frames, cid+': polar present frames')
            close(mask[window].mean(), row.occupancy, cid+': polar occupancy')
            if row.metric.startswith('salt_'):
                cutoff = float(row.metric.removeprefix('salt_').removesuffix('A'))
                check(np.array_equal(mask, series[row.feature+'__NO_min_A'].to_numpy() <= cutoff), cid+': salt threshold')
        for item in json.loads((root/'structures.json').read_text()):
            pdb = root/item['filename']
            atoms = [line for line in pdb.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
            check(len(atoms) == item['atoms'], cid+': representative atom count')
            check(abs(item['index']*.01-item['time_ns']) < 1e-9, cid+': representative time')

    screen = json.loads((study/'configs/polar-screen.json').read_text())
    definition = json.loads((study/'configs/polar-screen-definition.json').read_text())
    full = pd.read_csv(study/'results/polar-screen/comparison/single_feature_summary.csv')
    check(not full.duplicated(['case_id','feature','mode','window']).any(), 'Duplicate all-residue summary keys')
    check(set(full.case_id) == {c['id'] for c in screen['cases']}, 'Incomplete all-residue cohort')
    for case in screen['cases']:
        cid = case['id'];root = study/'results/polar-screen'/cid
        features = pd.read_csv(root/'features.csv')
        with np.load(root/'polar_features.npz', allow_pickle=False) as arrays:
            times = arrays['times_ns']
            matrices = {mode: arrays[mode] for mode in ('primary', 'tight', 'relaxed')}
            close(times, np.arange(10001)*.01, cid+': all-residue grid', 1e-7)
            check(list(arrays['features']) == list(features.feature), cid+': feature identities')
            feature_index = {name: i for i, name in enumerate(arrays['features'])}
            subset = full[full.case_id == cid]
            applicable = subset[subset.status == 'APPLICABLE']
            check(set(applicable.feature) == set(feature_index), cid+': full feature coverage')
            check(len(applicable) == len(feature_index)*16, cid+': all modes and windows')
            for row in subset.itertuples():
                if row.status != 'APPLICABLE':
                    check(row.feature not in feature_index and pd.isna(row.occupancy), cid+': chemical absence')
                    continue
                window = (times >= row.start_ns) & ((times <= row.end_ns) if row.end_inclusive else (times < row.end_ns))
                values = matrices[row.mode][window, feature_index[row.feature]]
                check(len(values) == row.n_frames and int(values.sum()) == row.present_frames, cid+': screen denominator')
                close(values.mean(), row.occupancy, cid+': screen occupancy')
            for threshold in ('primary', 'tight', 'relaxed'):
                flags = matrices[threshold]
                check(flags.dtype == bool and flags.shape == (10001, len(features)), cid+': feature matrix')
                salt = arrays['salt_feature_indices']
                cutoff = definition['salt_cutoffs_A'][threshold]
                distances = arrays['salt_min_NO_A']
                mismatch = flags[:, salt] != (distances <= cutoff)
                check(np.all(np.abs(distances[mismatch]-cutoff) <= np.spacing(np.float32(cutoff))), cid+': salt flags')

    patterns = pd.read_csv(study/'results/polar-screen/comparison/all_frequent_patterns.csv')
    pattern_keys = ['target_case', 'mode', 'window', 'features']
    check(not patterns.duplicated(pattern_keys).any(), 'Duplicate polar pattern identities')
    passing = pd.read_csv(study/'results/polar-screen-summary/qualifying_patterns.csv')
    pd.testing.assert_frame_equal(passing.reset_index(drop=True),
                                  patterns[patterns.qualifies].reset_index(drop=True))
    coverage = json.loads((study/'results/polar-screen-summary/screen_coverage.json').read_text())
    for window in coverage['windows']:
        subset = patterns[(patterns.target_case == window['target_case']) &
                          (patterns['mode'] == window['mode']) & (patterns.window == window['window'])]
        check(len(subset) == window['target_frequent_patterns'], 'Complete failed/passing pattern count')
        check(int(subset.qualifies.sum()) == window['qualifying_patterns'], 'Qualifying pattern count')

    transitions = json.loads((study/'configs/polar-transition.json').read_text())
    for case in transitions['cases']:
        cid = case['id'];root = study/'results/polar-transition'/cid
        table = pd.read_csv(root/'timeline.csv')
        check(len(table) == 10001, cid+': transition rows')
        with np.load(root/'geometry.npz', allow_pickle=False) as arrays:
            close(arrays['times_ns'], np.arange(10001)*.01, cid+': transition grid', 1e-7)
            for name in ('core_xyz_A', 'local_xyz_A', 'neighbor_distances_A'):
                check(np.isfinite(arrays[name]).all(), cid+': finite '+name)
    selected = pd.read_csv(study/'results/polar-transition-summary/compatibility_selected_models.csv')
    candidates = pd.read_csv(study/'results/polar-transition-summary/compatibility_candidates.csv')
    check(len(selected) == 10 and set(selected.ligand) == {'ascr3', 'ascr18'}, 'Static model comparison coverage')
    check(not selected.duplicated(['ligand', 'state']).any(), 'Duplicate static model identities')
    check(not candidates.duplicated(['ligand', 'state', 'conformer']).any(), 'Duplicate static candidate identities')
    for row in selected.itertuples():
        subset = candidates[(candidates.ligand == row.ligand) & (candidates.state == row.state)]
        check(len(subset) == row.refined_candidates, 'Static candidate denominator')
        check(int(subset.geometry_pass.sum()) == row.accepted_candidates, 'Static accepted candidates')
        check((study/'results/polar-transition/structures'/row.filename).is_file(), 'Missing constructed structure')
    return count
