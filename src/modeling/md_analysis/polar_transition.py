"""Config-driven follow-up of a polar-contact transition in existing trajectories.

Fits, endpoint-motion attribution and one-intermediate proximity checks are
descriptive geometry, not evidence of causality or a new dynamics calculation.
"""
from pathlib import Path
import json
import hashlib
import warnings
import csv
import numpy as np
from .common import ensure_dir
from .polar_io import load_config


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def table(path, rows):
    rows = list(rows)
    with open(path, 'w') as handle:
        if rows:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)


def fit(x, mobile, fixed):
    center, origin = x[mobile].mean(0), fixed.mean(0)
    a, _, b = np.linalg.svd((x[mobile] - center).T @ (fixed - origin))
    rot = a @ np.diag([1, 1, np.linalg.det(a @ b)]) @ b
    return (x - center) @ rot + origin


def min_distance(donor, acceptors):
    return np.linalg.norm(np.asarray(acceptors) - donor, axis=-1).min(axis=-1)


def endpoint_decomposition(d0, a0, d1, a1):
    """Symmetric two-endpoint attribution, exact even if nearest O changes."""
    before, after = min_distance(d0, a0), min_distance(d1, a1)
    ligand = .5 * (min_distance(d1, a0) - before + after - min_distance(d0, a1))
    return dict(before=float(before), after=float(after), delta=float(after-before),
                donor_motion=float(ligand), acceptor_motion=float(after-before-ligand))


def run_transition_case(cfg, case):
    import MDAnalysis as mda
    from MDAnalysis.lib.formats.libmdaxdr import XTCFile
    from MDAnalysis.lib.distances import distance_array, calc_dihedrals, calc_angles, calc_bonds
    from .polar_contacts import short_linear_geometry

    root = Path(cfg['output_root']); out = ensure_dir(root / case['id'])
    if (out / 'completion.json').exists():
        raise RuntimeError('Completed output exists: ' + str(out))
    identities = {}
    for key in ['topology', 'reference', 'trajectory']:
        digest = sha(case[key]); assert digest == case['expected_' + key + '_sha256']
        identities[key] = dict(path=case[key], sha256=digest, bytes=Path(case[key]).stat().st_size)
    u, ref = mda.Universe(case['topology']), mda.Universe(case['reference'])
    assert [(a.name, int(a.resid), a.resname[:4]) for a in u.atoms] == [(a.name, int(a.resid), a.resname) for a in ref.atoms]
    assert len(u.select_atoms('protein').residues) == 406
    for residue, name in case['expected_sites'].items():
        assert u.select_atoms('protein and resid ' + str(residue)).residues.resnames.tolist() == [name]
    assert abs(u.select_atoms('not protein').charges.sum() + 1) < .001
    tcfg = cfg['transition']; site, mutation = tcfg['acceptor_resid'], tcfg['mutation_resid']
    donor = u.select_atoms('not protein and name ' + tcfg['ligand_donor']).indices
    acceptor = u.select_atoms(f'protein and resid {site} and name ' + ' '.join(tcfg['acceptor_names'])).indices
    assert len(donor) == 1 and len(acceptor) == 2
    hs = [a.index for a in u.atoms[int(donor[0])].bonded_atoms if a.mass < 2]
    assert len(hs) == 1
    triplets = np.array([[int(donor[0]), hs[0], int(a)] for a in acceptor])
    def ids(query): return u.select_atoms(query).indices
    ligand = ids('not protein and not name H*')
    ring = ids('not protein and name ' + ' '.join(tcfg['ring_names']))
    carbox = ids('not protein and name O5 O6')
    d70 = ids('protein and resid 70 and name OD1 OD2')
    r345 = ids('protein and resid 345 and name NE NH1 NH2')
    r189 = ids('protein and resid 189 and name NE NH1 NH2')
    mutant = ids(f'protein and resid {mutation} and name NE NH1 NH2')
    dside = ids(f'protein and resid {site} and not backbone and not name H*')
    mside = ids(f'protein and resid {mutation} and not backbone and not name H*')
    tracked = np.unique(np.r_[ligand, hs, acceptor, dside, mside,
                             ids(f'protein and resid {site} {mutation} and backbone')])
    remap = {int(a): i for i, a in enumerate(tracked)}
    fixed = mda.Universe(cfg['alignment_reference'])
    assert sha(cfg['alignment_reference']) == cfg['alignment_reference_sha256']
    queries = dict(core='protein and backbone and resid ' + ' '.join(map(str, cfg['alignment_core_resids'])),
                   local='protein and backbone and resid ' + ' '.join(map(str, tcfg['local_fit_resids'])))
    fits = {}
    for name, query in queries.items():
        a, b = u.select_atoms(query), fixed.select_atoms(query)
        assert list(zip(a.resids, a.names)) == list(zip(b.resids, b.names))
        fits[name] = (a.indices, b.positions.astype(float))
    chi_definitions = {f'acceptor_chi1': ids(f'protein and resid {site} and name N CA CB CG')}
    chi_definitions['acceptor_chi1'] = np.array([ids(f'protein and resid {site} and name {n}')[0] for n in ['N','CA','CB','CG']])
    chi_definitions['acceptor_chi2'] = np.array([ids(f'protein and resid {site} and name {n}')[0] for n in ['CA','CB','CG','OD1']])
    if len(mutant):
        for k, names in enumerate([['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','NE'],['CG','CD','NE','CZ']], 1):
            chi_definitions[f'mutation_chi{k}'] = np.array([ids(f'protein and resid {mutation} and name {n}')[0] for n in names])
    protein_heavy = u.select_atoms('protein and not name H*')
    starts = np.r_[0, np.flatnonzero(np.diff(protein_heavy.resids)) + 1]
    residue_numbers = protein_heavy.resids[starts]
    prior = np.load(Path(case['prior_polar_output']) / 'polar_features.npz')
    prior_idx = prior['features'].tolist().index(tcfg['feature'])
    prior_flags = prior['primary'][:, prior_idx].copy()
    prior.close()
    rows, coords, neighbor_rows, neighbors = [], {name: [] for name in fits}, [], []
    independent_errors = []; chosen_xyz = {}
    with XTCFile(case['trajectory'], 'r') as handle:
        for i, fr in enumerate(handle):
            x = (fr.x * 10).astype(float); time = float(fr.time / 1000)
            da, ha, angles, valid = short_linear_geometry(x, triplets, tcfg['hbond'])
            best = int(np.argmin(np.where(valid, da, da + 100)))
            row = dict(frame=i, time_ns=time, hbond=int(valid.any()), min_DA_A=float(da.min()),
                       best_DA_A=float(da[best]), best_HA_A=float(ha[best]), best_DHA_deg=float(angles[best]),
                       distance_pass=int(((da <= 3.2) & (ha <= 2.4)).any()),
                       angle_pass=int((angles >= 150).any()),
                       arg345_carbox_A=float(distance_array(x[r345], x[carbox]).min()),
                       arg189_carbox_A=float(distance_array(x[r189], x[carbox]).min()),
                       mutation_acceptor_heavy_A=float(distance_array(x[mside], x[dside]).min()))
            if len(mutant):
                row.update(arg368_asp70_NO_A=float(distance_array(x[mutant], x[d70]).min()),
                           arg368_asp89_NO_A=float(distance_array(x[mutant], x[acceptor]).min()),
                           arg368_ligand_min_A=float(distance_array(x[mutant], x[ligand]).min()))
            for name, indices in chi_definitions.items():
                row[name + '_deg'] = float(np.degrees(calc_dihedrals(*x[indices])))
            assert row['hbond'] == int(prior_flags[i]), ('prior flag differs', i)
            rows.append(row)
            for name, (mobile, reference) in fits.items():
                coords[name].append(fit(x, mobile, reference)[tracked])
            if i % tcfg['neighbor_stride'] == 0:
                ds = distance_array(x[dside], x[protein_heavy.indices]).min(0)
                ms = distance_array(x[mside], x[protein_heavy.indices]).min(0)
                pair = np.array([np.minimum.reduceat(ds, starts), np.minimum.reduceat(ms, starts)])
                neighbors.append(pair); neighbor_rows.append(i)
            if i % 101 == 0:
                independent_errors.append(max(float(np.max(np.abs(calc_bonds(x[triplets[:,0]], x[triplets[:,2]]) - da))),
                                              float(np.max(np.abs(np.degrees(calc_angles(x[triplets[:,0]], x[triplets[:,1]], x[triplets[:,2]])) - angles)))))
    assert len(rows) == cfg['expected_frames']
    times = np.array([r['time_ns'] for r in rows]); assert np.allclose(times, np.arange(len(rows)) * .01, atol=1e-5)
    assert max(independent_errors) < .001
    coords = {name: np.array(value, dtype=np.float32) for name, value in coords.items()}
    neighbors = np.array(neighbors); neighbor_frames = np.array(neighbor_rows)
    table(out / 'timeline.csv', rows)
    np.savez_compressed(out / 'geometry.npz', times_ns=times, atom_indices=tracked,
                        atom_names=u.atoms[tracked].names.astype(str), atom_resids=u.atoms[tracked].resids,
                        neighbor_frames=neighbor_frames, neighbor_resids=residue_numbers,
                        neighbor_distances_A=neighbors, **{f'{name}_xyz_A': value for name,value in coords.items()})
    phases = []
    for name, lo, hi in tcfg['phases']:
        mask = (times >= lo) & (times < hi)
        p = dict(case_id=case['id'], phase=name, start_ns=lo, end_ns=hi, frames=int(mask.sum()))
        for key in rows[0]:
            if key not in ['frame','time_ns'] and not key.endswith('_deg'):
                values = np.array([r[key] for r in rows])[mask]
                p[key + ('_fraction' if key in ['hbond','distance_pass','angle_pass'] else '_median')] = float(np.mean(values) if key in ['hbond','distance_pass','angle_pass'] else np.median(values))
        for label, key in [('arg345_salt','arg345_carbox_A'),('arg368_Asp70_salt','arg368_asp70_NO_A')]:
            if key in rows[0]:p[label + '_fraction'] = float(np.mean([r[key] <= 4 for r,m in zip(rows,mask) if m]))
        phases.append(p)
    table(out / 'phase_summary.csv', phases)
    paths = []
    for name, lo, hi in tcfg['phases']:
        mask = (times[neighbor_frames] >= lo) & (times[neighbor_frames] < hi)
        for j, r in enumerate(residue_numbers):
            if abs(r-site) <= 2 or abs(r-mutation) <= 2:continue
            a, b = neighbors[mask,0,j] <= tcfg['proximity_A'], neighbors[mask,1,j] <= tcfg['proximity_A']
            if a.any() or b.any():
                paths.append(dict(case_id=case['id'],phase=name,intermediate=int(r),
                                  acceptor_neighbor_fraction=float(a.mean()),mutation_neighbor_fraction=float(b.mean()),
                                  simultaneous_fraction=float((a&b).mean())))
    table(out / 'one_intermediate_proximity.csv', paths)
    if case['id'] == tcfg['target_case']:
        structures = ensure_dir(root / 'structures')
        hb = np.array([r['hbond'] for r in rows], bool); chosen = {}
        for name, lo, hi in tcfg['phases']:
            if name not in tcfg['export_phases']:continue
            ix = np.flatnonzero((times >= lo) & (times < hi))
            keep = hb[ix] if hb[ix].mean() >= .5 else ~hb[ix]
            ix = ix[keep]
            atoms = [remap[int(a)] for a in np.r_[ring,acceptor]]
            cloud = coords['local'][ix][:,atoms]
            pick = ix[np.argmin(np.mean(np.sum((cloud-cloud.mean(0))**2,axis=2),axis=1))]
            chosen[name] = int(pick)
        for target_time in tcfg['compatibility_times_ns']:
            ix = np.flatnonzero(hb & (np.abs(times-target_time) <= tcfg['compatibility_tolerance_ns']))
            assert len(ix), ('no eligible real frame',target_time)
            chosen[f'compat_{target_time:g}ns'] = int(ix[np.argmin(np.abs(times[ix]-target_time))])
        records = []
        with XTCFile(case['trajectory'], 'r') as handle:
            for i, fr in enumerate(handle):
                keys = [name for name, frame in chosen.items() if frame == i]
                if not keys:continue
                x = fit((fr.x*10).astype(float), *fits['core']); u.atoms.positions = x
                for name in keys:
                    path = structures / (name+'.pdb')
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore');u.atoms.write(str(path))
                    reread=mda.Universe(path)
                    assert len(reread.atoms)==len(u.atoms) and np.max(np.abs(reread.atoms.positions-x))<.0006
                    records.append(dict(state=name, frame=i, time_ns=times[i], filename=path.name,
                                        sha256=sha(path), measurements=rows[i]))
        dump(structures/'representatives.json',dict(case_id=case['id'], records=records,
             selection='Phase: actual majority-state frame nearest local-fit mean sugar/Asp oxygen coordinates; compatibility: nearest eligible hydrogen-bonded frame to fixed times.'))
        decomposition = []
        for label, before, after in tcfg['transitions']:
            bounds = {p[0]:p[1:] for p in tcfg['phases']}
            ia = np.flatnonzero((times>=bounds[before][0])&(times<bounds[before][1]))[::10]
            ib = np.flatnonzero((times>=bounds[after][0])&(times<bounds[after][1]))[::10]
            for alignment, cloud in coords.items():
                stats=[]
                for a in ia:
                    for b in ib:
                        stats.append(endpoint_decomposition(cloud[a,remap[int(donor[0])]],cloud[a,[remap[int(v)] for v in acceptor]],
                                                            cloud[b,remap[int(donor[0])]],cloud[b,[remap[int(v)] for v in acceptor]]))
                decomposition.append(dict(transition=label,alignment=alignment,pairs=len(stats),
                    **{k+'_mean_A':float(np.mean([r[k] for r in stats])) for k in stats[0]},
                    sugar_center_shift_median_A=float(np.median([np.linalg.norm(cloud[b,[remap[int(v)] for v in ring]].mean(0)-cloud[a,[remap[int(v)] for v in ring]].mean(0)) for a in ia for b in ib]))))
        table(out/'motion_decomposition.csv',decomposition)
        absent = (times>=20)&(times<48)&~hb
        dump(out/'mid_state_interruptions.json',dict(absent_frames=int(absent.sum()),
            distance_failed=int(sum(r['distance_pass']==0 for r,m in zip(rows,absent) if m)),
            distance_pass_but_angle_failed=int(sum(r['distance_pass']==1 for r,m in zip(rows,absent) if m)),
            window_ns=[20,48], note='Distance tests include DA and HA on the same triplet; absent with distance pass requires an angle failure for that triplet.'))
    dump(out/'completion.json',dict(status='complete',frames=len(rows),input_identity=identities,
        prior_hbond_flags_exact_match=True,independent_geometry_max_error=max(independent_errors),
        config_sha256=cfg['_config_sha256'],module_sha256=sha(__file__)))
    print('TRANSITION_COMPLETE',case['id'],flush=True)


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--case')
    args=parser.parse_args();cfg=load_config(args.config)
    if args.case and args.case not in {c['id'] for c in cfg['cases']}:parser.error('Unknown case ID')
    for case in cfg['cases']:
        if not args.case or args.case==case['id']:run_transition_case(cfg,case)

if __name__=='__main__':main()
