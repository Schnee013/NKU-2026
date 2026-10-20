"""Polar-contact helpers used by the shared ligand contact-pattern driver.

All cutoffs come from a hash-bound screen configuration. Missing chemical
donors are non-applicable, not zero-occupancy observations.
"""
from pathlib import Path
import json
import numpy as np
from .common import ensure_dir


def short_linear_geometry(x, triplets, criteria):
    triplets = np.asarray(triplets, dtype=int)
    d, h, a = (x[triplets[:, i]] for i in range(3))
    dh, ah = d - h, a - h
    da = np.linalg.norm(d - a, axis=1)
    ha = np.linalg.norm(ah, axis=1)
    nh = np.linalg.norm(dh, axis=1)
    if np.any(nh < .7) or np.any(nh > 1.3) or np.any(ha == 0):
        raise ValueError('Invalid covalent N-H geometry in processed coordinates')
    angle = np.degrees(np.arccos(np.clip((dh * ah).sum(axis=1) / (nh * ha), -1, 1)))
    valid = ((da <= criteria['donor_acceptor_max_A']) &
             (ha <= criteria['H_acceptor_max_A']) & (angle >= criteria['DHA_min_deg']))
    return da, ha, angle, valid


def fill_short_gaps(mask, dt_ns, max_gap_ns):
    """Fill only bounded false runs of at most max_gap_ns sampled absence."""
    result = np.asarray(mask, dtype=bool).copy()
    edges = np.diff(np.r_[False, ~result, False].astype(int))
    for start, stop in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
        if start > 0 and stop < len(result) and (stop - start) * dt_ns <= max_gap_ns + 1e-9:
            result[start:stop] = True
    return result


def window_specs(screen):
    lo, hi = screen['primary_window_ns']
    yield 'primary', lo, hi, True
    for a, b in screen['secondary_windows_ns']:
        is_full = a == 0 and b == hi
        yield ('full' if is_full else f'block_{a:g}_{b:g}'), a, b, is_full or b == hi


def summarize_mask(mask, times, screen, max_gap_ns):
    from .analyze_ligand_contact_patterns import runs
    rows = []
    for label, lo, hi, inclusive in window_specs(screen):
        selected = (times >= lo) & (times <= hi if inclusive else times < hi)
        t, v = times[selected], mask[selected]
        intervals = runs(v, t)
        tolerant = runs(fill_short_gaps(v, times[1] - times[0], max_gap_ns), t)
        absent = np.flatnonzero(~v)
        rows.append(dict(window=label, start_ns=lo, end_ns=hi, end_inclusive=inclusive,
            n_frames=len(v), present_frames=int(v.sum()), occupancy=float(v.mean()),
            interval_count=len(intervals), longest_uninterrupted_ns=max((r[2] for r in intervals), default=0.),
            median_uninterrupted_ns=float(np.median([r[2] for r in intervals])) if intervals else 0.,
            present_at_window_start=bool(v[0]), present_at_window_end=bool(v[-1]),
            initial_uninterrupted_ns=float(t[absent[0]-1]-t[0]) if v[0] and len(absent) else float(t[-1]-t[0]) if v[0] else 0.,
            gap_tolerance_ns=max_gap_ns,
            longest_gap_tolerant_ns=max((r[2] for r in tolerant), default=0.)))
    return rows


def run_polar_case(cfg, case):
    import MDAnalysis as mda
    from MDAnalysis.lib.formats.libmdaxdr import XTCFile
    from MDAnalysis.lib.distances import calc_bonds, calc_angles
    from .analyze_local_mechanism import heavy, ligand_groups, hbond_candidates
    from .analyze_ligand_contact_patterns import sha, dump, table, runs

    screen_path = Path(cfg['polar_screen'])
    assert sha(screen_path) == cfg['polar_screen_sha256'], 'Frozen polar screen changed'
    screen = json.loads(screen_path.read_text())
    out = ensure_dir(Path(case['output']))
    if (out / 'polar_completion.json').exists():
        raise RuntimeError(f'Refuse to overwrite completed polar analysis: {out}')
    identity = {}
    for key in ['topology', 'reference', 'trajectory']:
        digest = sha(case[key])
        assert digest == case['expected_' + key + '_sha256'], (case['id'], key, 'input SHA changed')
        identity[key] = dict(path=case[key], sha256=digest, bytes=Path(case[key]).stat().st_size)
    assert identity['trajectory']['bytes'] == case['expected_trajectory_bytes']
    u, ref = mda.Universe(case['topology']), mda.Universe(case['reference'])
    assert len(u.atoms) == len(ref.atoms)
    assert all((a.name, int(a.resid), a.resname[:4]) == (b.name, int(b.resid), b.resname)
               for a, b in zip(u.atoms, ref.atoms)), 'TPR/reference order mismatch'
    protein, ligand = u.select_atoms('protein'), u.select_atoms('not protein')
    assert len(ligand.residues) == 1 and len(protein.residues) == 406
    assert protein.residues.resids.tolist() == list(range(1, 407))
    for resid, name in case['expected_sites'].items():
        assert protein.select_atoms('resid ' + str(resid)).residues.resnames.tolist() == [name]
    assert np.isclose(ligand.charges.sum(), -1., atol=1e-3), 'Ligand formal charge mismatch'
    lh = heavy(ligand)
    groups = ligand_groups(u, lh)
    carbox = lh[[groups[a.index] == 'carboxylate' and a.name.startswith('O') for a in lh]]
    assert len(carbox) == 2
    if case['ligand'] == cfg['target_ligand']:
        assert set(carbox.names) == set(screen['carboxylate_atoms']), 'Target carboxylate atom mapping mismatch'
    table(out / 'ligand_atom_mapping.csv', [dict(index=int(a.index), resid=int(a.resid), name=a.name,
        resname=a.resname, charge=float(a.charge), group=groups.get(a.index, 'hydrogen'),
        bonded_indices=' '.join(str(b.index) for b in a.bonded_atoms)) for a in ligand])
    table(out / 'sequence.csv', [dict(resid=int(r.resid), resname=r.resname) for r in protein.residues])

    definitions, metadata = {}, []
    targets = list(screen['targets']) + list(cfg.get('context_anchors', []))
    specifications = [(f'R{r}_ligand', r, carbox, 'carboxylate',
                       'predefined' if r in screen['targets'] else 'context') for r in targets]
    asp = protein.select_atoms(f"resid {cfg['competitor_resid']}")
    assert asp.residues.resnames.tolist() == ['ASP'] and np.isclose(asp.charges.sum(), -1, atol=1e-3)
    ao = asp.select_atoms('name OD1 OD2'); assert len(ao) == 2
    specifications.append((f"R{cfg['competing_anchor']}_D{cfg['competitor_resid']}", cfg['competing_anchor'], ao, 'Asp carboxylate', 'competitor'))
    for key, resid, acceptors, role, purpose in specifications:
        residue = protein.select_atoms('resid ' + str(resid))
        info = dict(feature=key, resid=resid, residue_name=residue.residues.resnames[0],
                    acceptor_role=role, analysis_role=purpose)
        if residue.residues.resnames.tolist() != ['ARG']:
            info.update(status='NOT_APPLICABLE', reason='Corresponding residue has no arginine guanidinium donor; no zero occupancy assigned.')
            metadata.append(info)
            continue
        assert np.isclose(residue.charges.sum(), 1., atol=1e-3)
        donors = residue.select_atoms('name ' + ' '.join(screen['arginine_donor_atoms']))
        assert set(donors.names) == set(screen['arginine_donor_atoms'])
        dids, aids = set(donors.indices), set(acceptors.indices)
        triplets = np.array([(d, h, a) for d, h, a, direction in hbond_candidates(u, residue, acceptors)
                             if direction == 'protein_to_ligand' and d in dids and a in aids], dtype=int)
        assert triplets.shape == (10, 3), (key, 'unexpected arginine donor/H/acceptor mapping')
        definitions[key] = (donors.indices, acceptors.indices, triplets)
        info.update(status='APPLICABLE', donor_residue_charge=float(residue.charges.sum()),
            donor_indices=donors.indices.tolist(), donor_names=donors.names.tolist(),
            acceptor_indices=acceptors.indices.tolist(), acceptor_names=acceptors.names.tolist(),
            triplets=triplets.tolist())
        metadata.append(info)
    dump(out / 'polar_atom_definitions.json', metadata)

    assert sha(cfg['alignment_reference']) == cfg['alignment_reference_sha256']
    alignment = mda.Universe(cfg['alignment_reference'])
    query = 'protein and backbone and resid ' + ' '.join(map(str, cfg['alignment_core_resids']))
    mobile, fixed = u.select_atoms(query), alignment.select_atoms(query)
    assert list(zip(mobile.resids, mobile.names)) == list(zip(fixed.resids, fixed.names))
    fixed_xyz = fixed.positions.astype(float); fixed_center = fixed_xyz.mean(axis=0)
    def align(x, ids):
        center = x[mobile.indices].mean(axis=0)
        left, _, right = np.linalg.svd((x[mobile.indices]-center).T @ (fixed_xyz-fixed_center))
        rotation = left @ np.diag([1., 1., np.linalg.det(left @ right)]) @ right
        return (x[ids] - center) @ rotation + fixed_center

    n = cfg['expected_frames']; times = np.empty(n)
    ligand_xyz = np.empty((n, len(lh), 3), dtype=np.float32)
    values = {k: np.empty((n, 7), dtype=float) for k in definitions}
    audit, checks = [], dict(max_heavy_bond_A=0., max_adjacent_CA_A=0., max_relative_ligand_step_A=0.)
    hi = set(heavy(u.atoms).indices)
    bonds = np.array([b.indices for b in u.bonds if all(i in hi for i in b.indices)])
    ca = protein.select_atoms('name CA').indices; previous = None
    audit_frames = set(np.linspace(0, n-1, 21, dtype=int)); count = 0
    with XTCFile(case['trajectory'], 'r') as trajectory:
        for i, frame in enumerate(trajectory):
            assert i < n
            x = (frame.x * 10).astype(float)
            assert x.shape == (len(u.atoms), 3) and np.isfinite(x).all()
            times[i] = frame.time / 1000.
            ligand_xyz[i] = align(x, lh.indices)
            checks['max_heavy_bond_A'] = max(checks['max_heavy_bond_A'], float(np.linalg.norm(x[bonds[:,0]]-x[bonds[:,1]], axis=1).max()))
            checks['max_adjacent_CA_A'] = max(checks['max_adjacent_CA_A'], float(np.linalg.norm(np.diff(x[ca], axis=0),axis=1).max()))
            relative = x[lh.indices].mean(axis=0) - x[ca].mean(axis=0)
            if previous is not None:
                checks['max_relative_ligand_step_A'] = max(checks['max_relative_ligand_step_A'], float(np.linalg.norm(relative-previous)))
            previous = relative
            for key, (donors, acceptors, triplets) in definitions.items():
                ds = np.sqrt(((x[donors, None] - x[acceptors][None]) ** 2).sum(axis=2))
                da, ha, angle, valid = short_linear_geometry(x, triplets, screen['short_linear_hbond'])
                eligible = np.flatnonzero(valid) if valid.any() else np.arange(len(triplets))
                best = int(eligible[np.argmin(da[eligible])])
                values[key][i] = [ds.min(), int(valid.any()), da[best], ha[best], angle[best], best, int(np.argmin(ds))]
                if i in audit_frames:
                    d, h, a = triplets.T
                    other_da = calc_bonds(x[d], x[a]); other_ha = calc_bonds(x[h], x[a])
                    other_angle = np.degrees(calc_angles(x[d], x[h], x[a]))
                    other_valid = ((other_da <= screen['short_linear_hbond']['donor_acceptor_max_A']) &
                        (other_ha <= screen['short_linear_hbond']['H_acceptor_max_A']) &
                        (other_angle >= screen['short_linear_hbond']['DHA_min_deg']))
                    other_no = min(float(calc_bonds(x[di][None], x[ai][None])[0]) for di in donors for ai in acceptors)
                    error = max(float(np.max(np.abs(da-other_da))), float(np.max(np.abs(ha-other_ha))), abs(float(ds.min())-other_no))
                    angle_error = float(np.max(np.abs(angle-other_angle)))
                    assert error < 1e-4 and angle_error < .001 and np.array_equal(valid, other_valid)
                    audit.append(dict(frame=i,time_ns=float(times[i]),feature=key,max_distance_error_A=error,
                        max_angle_error_deg=angle_error, all_triplet_flags_match=True))
            count += 1
    assert count == n and np.allclose(times, np.arange(n)*cfg['expected_step_ps']/1000, atol=1e-6)
    assert checks['max_heavy_bond_A'] < 3 and checks['max_adjacent_CA_A'] < 4.5 and checks['max_relative_ligand_step_A'] < 10
    table(out / 'independent_geometry_samples.csv', audit)
    masks, series = {}, {'time_ns': times}
    for key, data in values.items():
        for j, label in enumerate(['NO_min_A','short_linear_hbond','best_DA_A','best_HA_A','best_DHA_deg','best_triplet','nearest_pair']):
            series[key+'__'+label] = data[:,j]
        masks[key+'__short_linear_hbond'] = data[:,1].astype(bool)
        for cut in screen['salt_bridge_sensitivity_A']:
            masks[key+f'__salt_{cut:g}A'] = data[:,0] <= cut
    for key, mask in masks.items(): series[key] = mask.astype(int)
    import pandas as pd
    pd.DataFrame(series).to_csv(out/'polar_timeseries.csv.gz',index=False)
    summaries, episodes = [], []
    for key, mask in masks.items():
        feature, metric = key.split('__')
        summaries.extend(dict(case_id=case['id'], feature=feature, metric=metric, status='APPLICABLE', **r)
                         for r in summarize_mask(mask,times,screen,cfg['gap_tolerance_ns']))
        for a,b,duration in runs(mask,times):
            episodes.append(dict(feature=feature,metric=metric,start_ns=a,end_ns=b,duration_ns=duration,
                left_censored=a==times[0],right_censored=b==times[-1]))
    for meta in metadata:
        if meta['status']=='NOT_APPLICABLE':
            for metric in ['short_linear_hbond'] + [f'salt_{c:g}A' for c in screen['salt_bridge_sensitivity_A']]:
                for label,lo,hi,inclusive in window_specs(screen):
                    summaries.append(dict(case_id=case['id'],feature=meta['feature'],metric=metric,window=label,
                        start_ns=lo,end_ns=hi,end_inclusive=inclusive,status='NOT_APPLICABLE',reason=meta['reason']))
    table(out/'polar_summary.csv',summaries)
    pd.DataFrame(episodes,columns=['feature','metric','start_ns','end_ns','duration_ns','left_censored','right_censored']).to_csv(out/'uninterrupted_episodes.csv',index=False)
    competition=[]
    ak=f"R{cfg['competing_anchor']}_ligand__salt_{screen['salt_bridge_max_NO_A']:g}A"
    bk=f"R{cfg['competing_anchor']}_D{cfg['competitor_resid']}__salt_{screen['salt_bridge_max_NO_A']:g}A"
    if ak in masks and bk in masks:
        a,b=masks[ak],masks[bk]
        for label,lo,hi,inclusive in window_specs(screen):
            selected=(times>=lo)&(times<=hi if inclusive else times<hi)
            for state,mask in {'ligand_only':a&~b,'Asp_only':b&~a,'both':a&b,'neither':~a&~b}.items():
                competition.append(dict(case_id=case['id'],window=label,state=state,n_frames=int(selected.sum()),
                    state_frames=int(mask[selected].sum()),fraction=float(mask[selected].mean())))
    pd.DataFrame(competition,columns=['case_id','window','state','n_frames','state_frames','fraction']).to_csv(out/'competition_states.csv',index=False)

    # Choose actual frames from the most frequent late contact state, then the
    # nearest core-aligned ligand coordinates to that state's mean.
    state_keys=[k for k in [ak,bk,'R189_ligand__salt_4A'] if k in masks]
    state=np.column_stack([masks[k] for k in state_keys])
    lo,hi=screen['primary_window_ns']; late=(times>=lo)&(times<=hi)
    patterns,counts=np.unique(state[late],axis=0,return_counts=True)
    modal=patterns[np.argmax(counts)]; eligible=late&np.all(state==modal,axis=1)
    center=ligand_xyz[eligible].mean(axis=0)
    candidates=np.flatnonzero(eligible)
    representative=int(candidates[np.argmin(((ligand_xyz[candidates]-center)**2).mean(axis=(1,2)))])
    choices=[dict(label='start',index=0),dict(label='late_modal',index=representative),dict(label='end',index=n-1)]
    u.add_TopologyAttr('chainIDs',['A']*len(u.atoms));ligand.chainIDs='L';ligand.residues.resnames='LIG'
    u.add_TopologyAttr('elements',[a.name[0] for a in u.atoms])
    structures=[]
    with XTCFile(case['trajectory'],'r') as trajectory:
        for i,frame in enumerate(trajectory):
            selected=[c for c in choices if c['index']==i]
            if not selected:continue
            x=align((frame.x*10).astype(float),np.arange(len(u.atoms)));u.atoms.positions=x
            for choice in selected:
                filename=case['id']+'_'+choice['label']+'.pdb'
                with mda.Writer(str(out/filename),n_atoms=len(u.atoms),bonds='all') as writer:writer.write(u.atoms)
                distances=[]
                for key,(donors,acceptors,triplets) in definitions.items():
                    arr=values[key][i];ai,bi=np.unravel_index(int(arr[6]),(len(donors),len(acceptors)))
                    d,h,a=map(int,triplets[int(arr[5])])
                    distances.append(dict(feature=key,NO_min_A=float(arr[0]),nearest_donor_index=int(donors[ai]),
                        nearest_acceptor_index=int(acceptors[bi]),short_linear_hbond=bool(arr[1]),donor_index=d,hydrogen_index=h,
                        acceptor_index=a,DA_A=float(arr[2]),HA_A=float(arr[3]),DHA_deg=float(arr[4])))
                structures.append(dict(**choice,time_ns=float(times[i]),filename=filename,atoms=len(u.atoms),
                    distances=distances,selection='Fixed endpoint' if choice['label']!='late_modal' else 'Actual frame nearest mean core-aligned ligand coordinates within most frequent primary-window contact state',
                    modal_state=dict(zip(state_keys,map(bool,modal))),eligible_late_frames=int(eligible.sum())))
    dump(out/'structures.json',structures)
    dump(out/'input_identity.json',identity)
    dump(out/'polar_validation.json',dict(status='PASS',frames=count,first_ns=float(times[0]),last_ns=float(times[-1]),
        atoms=len(u.atoms),protein_residues=len(protein.residues),ligand_atoms=len(ligand),
        ligand_carboxylate_oxygen_names=carbox.names.tolist(),independent_geometry_samples=len(audit),**checks))
    dump(out/'polar_completion.json',dict(status='TERMINAL_NORMAL',config_sha256=sha(cfg['_path']),
        screen_sha256=sha(screen_path),source_sha256=sha(__file__),input_identity=identity,
        files={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
