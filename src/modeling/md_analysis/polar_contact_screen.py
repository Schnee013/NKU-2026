"""Whole-receptor polar screening for the shared ligand-contact driver.

Uses topology bonds, chemically eligible N/O hydrogen bonds, charged cation
groups, and the existing exact frequent-itemset implementation. No residue or
initial-pocket preselection is made. A geometric hydrogen bond is not an energy.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from .common import ensure_dir


def aggregate_flags(valid, feature_ids, n_features):
    """Count a residue/group once even when several H/O atom pairs qualify."""
    flags = np.zeros(n_features, dtype=bool)
    flags[np.asarray(feature_ids, dtype=int)[valid]] = True
    return flags


def distinct_residues(itemset, definitions):
    ids = [definitions[i]['resid'] for i in itemset]
    return len(ids) == len(set(ids))


def definitions_for(u, screen):
    from .analyze_local_mechanism import heavy, ligand_groups, hbond_candidates
    protein, ligand = u.select_atoms('protein'), u.select_atoms('not protein')
    lh = heavy(ligand)
    groups = ligand_groups(u, lh)
    lid = set(ligand.indices)
    triplets, feature_ids, features, lookup, mapping = [], [], [], {}, []

    def feature(kind, residue, group, region, direction, **extra):
        key = f'{kind}:{residue.resid}:{region}:{direction}:{group}'
        if key not in lookup:
            lookup[key] = len(features)
            features.append(dict(feature=key, kind=kind, resid=int(residue.resid),
                resname=residue.resname, ligand_group=group, region=region,
                direction=direction, **extra))
        return lookup[key]

    for d, h, a, direction in hbond_candidates(u, protein, ligand):
        # Strong N/O geometry is the requested screen; sulfur bonds need a
        # separate definition and are not silently judged with N/O cutoffs.
        if u.atoms[d].name[0] not in 'NO' or u.atoms[a].name[0] not in 'NO':
            continue
        la, pa = (a, d) if direction == 'protein_to_ligand' else (d, a)
        assert la in lid and pa not in lid
        atom = u.atoms[pa]
        region = 'backbone' if atom.name in ('N', 'O', 'OXT', 'OT1', 'OT2') else 'sidechain'
        fi = feature('HB', atom.residue, groups[la], region, direction)
        triplets.append((d, h, a)); feature_ids.append(fi)
        mapping.append(dict(triplet=len(triplets)-1, feature=features[fi]['feature'],
            donor_index=int(d), hydrogen_index=int(h), acceptor_index=int(a),
            donor=f'{u.atoms[d].resname}{u.atoms[d].resid}:{u.atoms[d].name}',
            hydrogen=u.atoms[h].name,
            acceptor=f'{u.atoms[a].resname}{u.atoms[a].resid}:{u.atoms[a].name}'))

    carbox = lh[[groups[a.index] == 'carboxylate' and a.name.startswith('O') for a in lh]]
    assert len(carbox) == 2 and all(not any(b.name.startswith('H') for b in a.bonded_atoms) for a in carbox)
    assert np.isclose(ligand.charges.sum(), -1, atol=1e-3)
    salt_pairs, salt_feature_ids = [], []
    for residue in protein.residues:
        names = {'ARG': ['NE','NH1','NH2'], 'LYS': ['NZ'], 'HSP': ['ND1','NE2'],
                 'HIP': ['ND1','NE2']}.get(residue.resname)
        if names:
            donors = residue.atoms.select_atoms('name ' + ' '.join(names))
            assert set(donors.names) == set(names)
            # A terminal residue can also carry its terminal charge.
            assert float(residue.atoms.charges.sum()) > .5
            fi = feature('SB', residue, 'carboxylate', 'sidechain', 'cation_to_ligand',
                         donor_indices=donors.indices.tolist(), acceptor_indices=carbox.indices.tolist())
            for d in donors.indices:
                for a in carbox.indices:
                    salt_pairs.append((d,a)); salt_feature_ids.append(fi)
        # Include a protonated N terminus when the topology actually contains it.
        if int(residue.resid) == int(protein.residues[0].resid):
            nitrogens = residue.atoms.select_atoms('name N')
            if len(nitrogens) == 1 and sum(a.name.startswith('H') for a in nitrogens[0].bonded_atoms) >= 3:
                fi = feature('SB', residue, 'carboxylate', 'backbone', 'cation_to_ligand',
                    donor_indices=nitrogens.indices.tolist(), acceptor_indices=carbox.indices.tolist())
                for a in carbox.indices:
                    salt_pairs.append((nitrogens[0].index,a)); salt_feature_ids.append(fi)
    return dict(features=features, triplets=np.asarray(triplets,dtype=int),
        hbond_feature_ids=np.asarray(feature_ids,dtype=int), triplet_mapping=mapping,
        salt_pairs=np.asarray(salt_pairs,dtype=int), salt_feature_ids=np.asarray(salt_feature_ids,dtype=int),
        ligand_groups=groups, carboxylate_indices=carbox.indices.tolist())


def scan_all_polar(cfg, case):
    import MDAnalysis as mda
    from MDAnalysis.lib.formats.libmdaxdr import XTCFile
    from MDAnalysis.lib.distances import calc_bonds, calc_angles
    from .analyze_local_mechanism import heavy
    from .analyze_ligand_contact_patterns import sha, dump, table
    from .polar_contacts import short_linear_geometry

    screen_path = Path(cfg['screen_path'])
    assert sha(screen_path) == cfg['screen_sha256']
    screen = json.loads(screen_path.read_text())
    out = ensure_dir(Path(case['output']))
    if (out/'all_polar_completion.json').exists():
        raise RuntimeError(f'Refuse to overwrite completed analysis: {out}')
    identity = {}
    for key in ('topology','reference','trajectory'):
        digest = sha(case[key])
        assert digest == case['expected_'+key+'_sha256'], (case['id'],key,'SHA mismatch')
        identity[key] = dict(path=case[key],sha256=digest,bytes=Path(case[key]).stat().st_size)
    assert identity['trajectory']['bytes'] == case['expected_trajectory_bytes']
    u, ref = mda.Universe(case['topology']), mda.Universe(case['reference'])
    assert len(u.atoms) == len(ref.atoms)
    assert all((a.name,int(a.resid),a.resname[:4]) == (b.name,int(b.resid),b.resname)
               for a,b in zip(u.atoms,ref.atoms)), 'TPR/reference atom-order mismatch'
    protein, ligand = u.select_atoms('protein'),u.select_atoms('not protein')
    assert len(ligand.residues) == 1 and protein.residues.resids.tolist() == list(range(1,407))
    for resid,name in case['expected_sites'].items():
        assert protein.select_atoms('resid '+str(resid)).residues.resnames.tolist() == [name]
    defs = definitions_for(u,screen)
    features, triplets = defs['features'],defs['triplets']
    hfi, sfi, pairs = defs['hbond_feature_ids'],defs['salt_feature_ids'],defs['salt_pairs']
    nf, n = len(features),cfg['expected_frames']
    matrices = {name:np.zeros((n,nf),dtype=bool) for name in ('primary','tight','relaxed')}
    times = np.empty(n)
    definitions = {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in defs.items() if k!='ligand_groups'}
    dump(out/'atom_definitions.json',definitions)
    table(out/'features.csv',features)
    table(out/'ligand_atom_mapping.csv',[dict(index=int(a.index),name=a.name,resname=a.resname,
        group=defs['ligand_groups'].get(a.index,'hydrogen'),charge=float(a.charge),
        bonds=' '.join(str(b.index) for b in a.bonded_atoms)) for a in ligand])
    table(out/'sequence.csv',[dict(resid=int(r.resid),resname=r.resname) for r in protein.residues])
    heavy_ids=set(heavy(u.atoms).indices)
    bonds=np.array([b.indices for b in u.bonds if all(i in heavy_ids for i in b.indices)])
    ca=protein.select_atoms('name CA').indices; lh=heavy(ligand).indices
    checks=dict(max_heavy_bond_A=0.,max_adjacent_CA_A=0.,max_relative_ligand_step_A=0.,
                max_independent_distance_error_A=0.,max_independent_angle_error_deg=0.)
    audit_frames=set(np.linspace(0,n-1,41,dtype=int)); audit=[]; previous=None; count=0
    event_fields=['frame','time_ns','feature','triplet','donor','hydrogen','acceptor','DA_A','HA_A','DHA_deg']
    salt_columns=sorted(set(sfi.tolist()))
    salt_distances=np.full((n,len(salt_columns)),np.inf,dtype=np.float32)
    with gzip.open(out/'hbond_events.csv.gz','wt',newline='') as event_handle, XTCFile(case['trajectory'],'r') as trajectory:
        writer=csv.DictWriter(event_handle,fieldnames=event_fields);writer.writeheader()
        for i,frame in enumerate(trajectory):
            assert i<n
            x=(frame.x*10).astype(float)
            assert x.shape==(len(u.atoms),3) and np.isfinite(x).all()
            times[i]=frame.time/1000.
            da,ha,angle,valid=short_linear_geometry(x,triplets,screen['hbond']['primary'])
            sd=np.linalg.norm(x[pairs[:,0]]-x[pairs[:,1]],axis=1)
            minimum=np.full(nf,np.inf)
            np.minimum.at(minimum,sfi,sd)
            salt_distances[i]=minimum[salt_columns]
            for mode in matrices:
                limits=screen['hbond'][mode]
                hv=((da<=limits['donor_acceptor_max_A']) & (ha<=limits['H_acceptor_max_A']) & (angle>=limits['DHA_min_deg']))
                matrices[mode][i]=aggregate_flags(hv,hfi,nf) | (minimum<=screen['salt_cutoffs_A'][mode])
            for t in np.flatnonzero(valid):
                meta=defs['triplet_mapping'][int(t)]
                writer.writerow(dict(frame=i,time_ns=float(times[i]),feature=meta['feature'],triplet=int(t),
                    donor=meta['donor'],hydrogen=meta['hydrogen'],acceptor=meta['acceptor'],
                    DA_A=float(da[t]),HA_A=float(ha[t]),DHA_deg=float(angle[t])))
            checks['max_heavy_bond_A']=max(checks['max_heavy_bond_A'],float(np.linalg.norm(x[bonds[:,0]]-x[bonds[:,1]],axis=1).max()))
            checks['max_adjacent_CA_A']=max(checks['max_adjacent_CA_A'],float(np.linalg.norm(np.diff(x[ca],axis=0),axis=1).max()))
            relative=x[lh].mean(axis=0)-x[ca].mean(axis=0)
            if previous is not None:
                checks['max_relative_ligand_step_A']=max(checks['max_relative_ligand_step_A'],float(np.linalg.norm(relative-previous)))
            previous=relative
            if i in audit_frames:
                d,h,a=triplets.T
                oda=calc_bonds(x[d],x[a]);oha=calc_bonds(x[h],x[a]);oa=np.degrees(calc_angles(x[d],x[h],x[a]))
                osd=calc_bonds(x[pairs[:,0]],x[pairs[:,1]])
                de=max(float(np.max(np.abs(da-oda))),float(np.max(np.abs(ha-oha))),float(np.max(np.abs(sd-osd))))
                ae=float(np.max(np.abs(angle-oa)))
                checks['max_independent_distance_error_A']=max(checks['max_independent_distance_error_A'],de)
                checks['max_independent_angle_error_deg']=max(checks['max_independent_angle_error_deg'],ae)
                assert de<1e-4 and ae<.01
                for mode in matrices:
                    lim=screen['hbond'][mode]
                    ohv=(oda<=lim['donor_acceptor_max_A'])&(oha<=lim['H_acceptor_max_A'])&(oa>=lim['DHA_min_deg'])
                    other=np.zeros(nf,dtype=bool)
                    for fi in set(hfi[ohv]):other[fi]=True
                    for fi in set(sfi[osd<=screen['salt_cutoffs_A'][mode]]):other[fi]=True
                    assert np.array_equal(other,matrices[mode][i]),(case['id'],i,mode,'flag mismatch')
                audit.append(dict(frame=i,time_ns=float(times[i]),triplets=len(triplets),salt_pairs=len(pairs),
                    distance_error_A=de,angle_error_deg=ae,all_feature_flags_match=True))
            count+=1
            if i%1000==0:print(case['id'],'frame',i,'time_ns',times[i],flush=True)
    assert count==n and np.allclose(times,np.arange(n)*cfg['expected_step_ps']/1000.,atol=1e-6)
    assert checks['max_heavy_bond_A']<3 and checks['max_adjacent_CA_A']<4.5 and checks['max_relative_ligand_step_A']<10
    assert np.all(matrices['tight']<=matrices['primary']) and np.all(matrices['primary']<=matrices['relaxed'])
    np.savez_compressed(out/'polar_features.npz',times_ns=times,features=np.array([f['feature'] for f in features]),
        salt_feature_indices=np.array(salt_columns),salt_min_NO_A=salt_distances,**matrices)
    table(out/'independent_geometry_samples.csv',audit)
    dump(out/'input_identity.json',identity)
    dump(out/'all_polar_validation.json',dict(status='PASS',frames=count,atoms=len(u.atoms),
        first_ns=float(times[0]),last_ns=float(times[-1]),protein_residues=len(protein.residues),
        ligand_atoms=len(ligand),features=nf,triplets=len(triplets),salt_pairs=len(pairs),
        n_hbond_features=sum(f['kind']=='HB' for f in features),
        n_salt_features=sum(f['kind']=='SB' for f in features),independent_frames=len(audit),**checks))
    dump(out/'all_polar_completion.json',dict(status='TERMINAL_NORMAL',at=datetime.now(timezone.utc).isoformat(),
        config_sha256=sha(cfg['_path']),screen_sha256=sha(screen_path),source_sha256=sha(__file__),
        files={p.name:sha(p) for p in out.iterdir() if p.is_file()}))


def compare_all_polar(cfg):
    import pandas as pd
    from .analyze_ligand_contact_patterns import sha,dump,table,frequent_itemsets,runs
    from .polar_contacts import summarize_mask
    screen=json.loads(Path(cfg['screen_path']).read_text())
    assert sha(cfg['screen_path'])==cfg['screen_sha256']
    out=ensure_dir(Path(cfg['shared_output']))
    data={}; metadata={}; local_defs={}; times=None
    for case in cfg['cases']:
        base=Path(case['output']);completion=json.loads((base/'all_polar_completion.json').read_text())
        assert completion['status']=='TERMINAL_NORMAL'
        for name,digest in completion['files'].items():assert sha(base/name)==digest,(case['id'],name)
        data[case['id']]=np.load(base/'polar_features.npz')
        ds=json.loads((base/'atom_definitions.json').read_text())['features']
        local_defs[case['id']]={d['feature']:d for d in ds}
        metadata.update({d['feature']:d for d in ds})
        if times is None:times=data[case['id']]['times_ns']
        else:assert np.array_equal(times,data[case['id']]['times_ns'])
    names=sorted(metadata); definitions=[metadata[k] for k in names]
    lookup={k:i for i,k in enumerate(names)}
    aligned={};available={}
    for cid,d in data.items():
        ids=np.array([lookup[str(f)] for f in d['features']])
        available[cid]=set(str(f) for f in d['features'])
        aligned[cid]={}
        for mode in ('primary','tight','relaxed'):
            matrix=np.zeros((len(times),len(names)),dtype=bool);matrix[:,ids]=d[mode]
            aligned[cid][mode]=matrix
    summary_screen=dict(primary_window_ns=screen['primary_window_ns'],
        secondary_windows_ns=[[0,100]]+[[n,n+10] for n in range(0,100,10)])
    singles=[]
    for cid,modes in aligned.items():
        for j,name in enumerate(names):
            if name not in available[cid]:
                singles.append(dict(case_id=cid,feature=name,mode='primary',window='primary',
                    status='CHEMICALLY_ABSENT',resid=definitions[j]['resid'],kind=definitions[j]['kind']))
                continue
            for mode,matrix in modes.items():
                rows=summarize_mask(matrix[:,j],times,summary_screen,screen['gap_tolerance_ns'])
                if mode!='primary':rows=[r for r in rows if r['window'] in ('primary','full')]
                singles.extend(dict(case_id=cid,feature=name,mode=mode,status='APPLICABLE',
                    resid=definitions[j]['resid'],kind=definitions[j]['kind'],**r) for r in rows)
    table(out/'single_feature_summary.csv',singles)
    primary_rows={ (r['case_id'],r['feature']):r for r in singles if r['mode']=='primary' and r['window']=='primary'}
    candidates=[]; coverage=[]
    windows=[('primary',50.,100.)]+[(f'exploratory_{i}_{i+20}',float(i),float(i+20)) for i in range(0,81,10)]
    controls=cfg['control_cases']
    for target in cfg['target_cases']:
        for mode in ('primary','tight','relaxed'):
            for window_name,lo,hi in windows:
                # Exploratory intervals have the same half-open convention as blocks.
                selected=(times>=lo)&(times<=hi if hi==100 else times<hi)
                minimum=math.ceil(screen['target_min']*int(selected.sum())-1e-9)
                counts=Counter();hits=0;frequent=0
                for itemset,joint in frequent_itemsets(aligned[target][mode][selected],minimum):
                    if not distinct_residues(itemset,definitions):continue
                    counts[len(itemset)]+=1;frequent+=1
                    keys=[names[j] for j in itemset]
                    stats={};block_ok=False
                    for cid,modes in aligned.items():
                        mask=modes[mode][:,itemset].all(axis=1)
                        values=mask[selected]
                        blocks=[]
                        for a in np.arange(lo,hi,10.):
                            b=min(a+10,hi);sel=(times>=a)&(times<=b if b==100 else times<b)
                            blocks.append(float(mask[sel].mean()))
                        spans=runs(values,times[selected])
                        stats[cid]=dict(occupancy=float(values.mean()),
                            longest_uninterrupted_ns=max((r[2] for r in spans),default=0.),
                            blocks=blocks,chemically_absent=[k for k in keys if k not in available[cid]])
                    block_ok=any(a>=screen['target_min'] and b>=screen['target_min']
                        for a,b in zip(stats[target]['blocks'],stats[target]['blocks'][1:]))
                    low_controls=all(stats[c]['occupancy']<=screen['control_max'] for c in controls)
                    beyond_345=any(definitions[j]['resid']!=345 for j in itemset)
                    passes=block_ok and low_controls and beyond_345
                    hits+=int(passes)
                    candidates.append(dict(target_case=target,mode=mode,window=window_name,start_ns=lo,end_ns=hi,
                        order=len(itemset),features=' & '.join(keys),resids=';'.join(str(definitions[j]['resid']) for j in itemset),
                        target_occupancy=stats[target]['occupancy'],
                        max_control_occupancy=max(stats[c]['occupancy'] for c in controls),
                        target_longest_uninterrupted_ns=stats[target]['longest_uninterrupted_ns'],
                        two_adjacent_high_blocks=block_ok,beyond_arg345=beyond_345,qualifies=passes,
                        stats=json.dumps(stats,separators=(',',':'))))
                coverage.append(dict(target_case=target,mode=mode,window=window_name,frames=int(selected.sum()),
                    target_frequent_patterns=frequent,orders=dict(counts),qualifying_patterns=hits))
    frame=pd.DataFrame(candidates)
    frame.to_csv(out/'all_frequent_patterns.csv',index=False)
    frame[frame.qualifies].to_csv(out/'qualifying_patterns.csv',index=False)
    # Save every primary single, including intermittent contacts below the screen.
    display=[]
    for target in cfg['target_cases']:
        for name in names:
            r=primary_rows[(target,name)]
            if r['status']!='APPLICABLE' or r.get('occupancy',0)==0:continue
            row=dict(target_case=target,feature=name,resid=r['resid'],kind=r['kind'],
                target_occupancy=r['occupancy'],longest_uninterrupted_ns=r['longest_uninterrupted_ns'])
            for c in controls:
                cr=primary_rows[(c,name)];row[c]=cr.get('occupancy')
                row[c+'_status']=cr['status']
            display.append(row)
    pd.DataFrame(display).sort_values(['target_case','target_occupancy'],ascending=[True,False]).to_csv(out/'observed_primary_contacts.csv',index=False)
    dump(out/'screen_coverage.json',dict(feature_union=len(names),cases=list(data),targets=cfg['target_cases'],
        controls=controls,windows=coverage,combination_definition='AND across distinct receptor residues; no order cap; exact frequency pruning',
        chemical_absence='No eligible atom group in topology, not a measured 0% interaction. Such a pattern cannot exist in that control.',
        primary_hits=int(((frame['mode']=='primary')&(frame.window=='primary')&frame.qualifies).sum()),
        exploratory_primary_geometry_hits=int(((frame['mode']=='primary')&(frame.window!='primary')&frame.qualifies).sum())))
    dump(out/'comparison_completion.json',dict(status='TERMINAL_NORMAL',at=datetime.now(timezone.utc).isoformat(),
        config_sha256=sha(cfg['_path']),screen_sha256=sha(cfg['screen_path']),source_sha256=sha(__file__),
        files={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
    print('COMPARISON',json.dumps(json.loads((out/'screen_coverage.json').read_text()),separators=(',',':')),flush=True)
