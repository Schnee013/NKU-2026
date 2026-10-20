"""Reconstruct supplied MD statistics and check independent snapshot geometry.

This reads scientific time series, fixed selections and six-frame structures.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .common import REPO_ROOT, ensure_dir
from .independent_geometry import independent_fit
from .analyze_local_mechanism import stretches, first_run

STUDY = REPO_ROOT / 'studies/receptor-engineering/molecular-dynamics'
# Fixed in advance: CSV calculations preserve source float64 precision. PDB
# coordinates are rounded to 0.001 A, so coordinate-derived checks allow 0.003 A.
CSV_ATOL = 1e-10
PDB_ATOL = 0.003


def pdb_frames(path):
    frames, xyz, names, resids = [], [], [], []
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                xyz.append([float(line[s:s+8]) for s in (30,38,46)])
                if not frames:
                    names.append(line[12:16].strip());resids.append(int(line[22:26]))
            elif line.startswith('ENDMDL'):
                frames.append(np.array(xyz));xyz=[]
    if len(frames) != 6:
        raise ValueError('Expected actual 0,10,25,50,75,100 ns PDB models')
    return np.array(frames), np.array(names), np.array(resids)


def load_case_summary(path, cases, keys, expected_keys=None):
    """Read one row per case and scientific key from a complete summary table."""
    table = pd.read_csv(path)
    identity = ['case', *keys]
    if not set(identity) <= set(table.columns):
        raise ValueError(f'{path}: missing summary identity fields')
    if table[identity].isna().any().any() or table.duplicated(identity).any():
        raise ValueError(f'{path}: missing or duplicate summary key')
    if set(table['case']) != set(cases):
        raise ValueError(f'{path}: summary must contain exactly the configured cases')
    if expected_keys is not None:
        actual = set(table[identity].itertuples(index=False, name=None))
        expected = {(case, *key) for case in cases for key in expected_keys}
        if actual != expected:
            raise ValueError(f'{path}: summary keys differ from the configured definitions')
    return table


def verify(study: Path):
    inputs = json.loads((study/'inputs.json').read_text())
    selections = json.loads((study/'configs/selections.json').read_text())
    mechanism = json.loads((study/'configs/mechanism.json').read_text())
    if len(inputs) != 8 or len({c['dataset_id'] for c in inputs}) != 8:
        raise ValueError('All eight predefined systems are required')
    cases = [case['dataset_id'] for case in inputs]
    systems = json.loads((study/'configs/systems.json').read_text())['systems']
    if set(cases) != set(systems) or set(cases) != {c['id'] for c in mechanism['cases']}:
        raise ValueError('Input, basic-analysis and local-analysis case identities differ')
    tables = study/'results/local-summary'
    anchors_all = load_case_summary(
        tables/'anchor_sensitivity.csv', cases, ['resid', 'cutoff_A'],
        [(resid, cutoff) for resid in mechanism['sites'] + mechanism['anchors']
         for cutoff in mechanism['contact_cutoffs_A']],
    )
    pairs_all = load_case_summary(
        tables/'local_pair_summary.csv', cases, ['site', 'neighbor'],
        mechanism['receptor_pairs'],
    )
    rotamers_all = load_case_summary(
        tables/'sidechain_rotamers.csv', cases, ['site', 'chi1_bin_center_deg'],
        [(site, center) for site in mechanism['sites'] for center in (-60, 60, 180)],
    )
    bridges_all = load_case_summary(
        tables/'proximity_bridges.csv', cases, ['site', 'neighbor', 'anchor'],
    )
    helix_sets=[set(r['helix_resids']) & set(r['membrane_resids']) for r in selections['records']]
    core=stretches(set.intersection(*helix_sets),7)
    pocket=sorted(set.union(*(set(r['initial_pocket_resids']) for r in selections['records']))-set(selections['lid_resids']))
    if core!=selections['core_resids'] or pocket!=selections['pocket_resids'] or len(core)!=123 or len(pocket)!=37:
        raise ValueError('Fixed selection differs from recorded frame0 definition')
    check_count = 0
    rebuilt = []
    def check(case, name, error=0., tolerance=CSV_ATOL):
        nonlocal check_count
        error=float(error)
        if not np.isfinite(error) or error > tolerance:
            raise ValueError(f'{case}: {name}: error {error} exceeds {tolerance}')
        check_count += 1
    base_summary=pd.read_csv(study/'results/base-summary/case_summary.csv').set_index('case_id')
    local_summary=pd.read_csv(study/'results/local-summary/case_summary.csv').set_index('case')
    if (base_summary.index.has_duplicates or local_summary.index.has_duplicates or
            set(base_summary.index) != set(cases) or set(local_summary.index) != set(cases)):
        raise ValueError('Basic and local summaries require one row for each configured system')
    paths = pd.read_csv(study/'results/base-summary/ligand-paths/centroid_path.csv')
    path_summary = pd.read_csv(study/'results/base-summary/ligand_path_primary_summary.csv')
    if (set(paths.system) != set(cases) or set(path_summary.case_id) != set(cases) or
            path_summary.case_id.duplicated().any()):
        raise ValueError('Common-frame paths require all configured systems')
    for row in path_summary.itertuples():
        path = paths[(paths.arm == row.arm) & (paths.system == row.case_id)]
        if len(path) != 201 or path.time_ns.duplicated().any():
            raise ValueError('Invalid common-frame sample identities')
        late_path = path[(path.time_ns >= row.window_start_ns) & (path.time_ns <= row.window_end_ns)]
        check(row.case_id, 'path_denominator', abs(len(late_path)-row.n_samples))
        check(row.case_id, 'path_displacement', abs(late_path.centroid_displacement_from_arm_initial_in_common_frame_A.mean()-row.centroid_displacement_mean_A))
        check(row.case_id, 'path_contact_count', abs(late_path.contacting_receptor_residues.mean()-row.dynamic_contact_residues_mean))
    for case in inputs:
        cid=case['dataset_id'];base=study/'results/base'/cid;local=study/'results/local'/cid
        f=pd.read_csv(local/'frame_metrics.csv');times=f.time_ns.to_numpy();late=times>=50
        check(cid,'100ps_grid',np.max(abs(times-np.arange(1001)*.1)),1e-7)
        if len(f)!=1001 or int(late.sum())!=501 or not np.isfinite(f.drop(columns='case').to_numpy()).all():
            raise ValueError(f'{cid}: invalid full local sample set')
        metrics=pd.read_csv(local/'metric_summary.csv').set_index('metric')
        for metric in metrics.index:
            arr=f[metric].to_numpy()
            expected={'initial':arr[0],'early_mean':arr[~late].mean(),'late_mean':arr[late].mean(),'late_min':arr[late].min(),'late_max':arr[late].max()}
            for stat,value in expected.items():check(cid,metric+'_'+stat,abs(value-metrics.loc[metric,stat]))
            if metric in local_summary.columns:check(cid,'case_'+metric,abs(arr[late].mean()-local_summary.loc[cid,metric]))
            rebuilt.append(dict(case=cid,metric=metric,**expected))
        rotamers=rotamers_all.loc[rotamers_all.case==cid]
        for row in rotamers.itertuples():
            angles=f[f'site{row.site}_chi1_deg'].to_numpy()
            centers=np.array([-60,60,180]); delta=np.abs((angles[:,None]-centers+180)%360-180)
            state=centers[delta.argmin(1)]
            check(cid,f'rotamer_{row.site}_{row.chi1_bin_center_deg}',abs(np.mean(state[late]==row.chi1_bin_center_deg)-row.late_occupancy))
        hbond_events=pd.read_csv(local/'hbond_events.csv')
        union=pd.read_csv(study/'results/local-summary/anchor_hbond_union.csv')
        for row in union[union.case==cid].itertuples():
            event_times=hbond_events.loc[(hbond_events.resid==row.resid)&(hbond_events.time_ns>=50),'time_ns'].unique()
            check(cid,f'union_hbond_{row.resid}',abs(len(event_times)/501-row.late_any_hbond_occupancy))
        anchors=anchors_all.loc[anchors_all.case==cid]
        for row in anchors.itertuples():
            arr=f[f'res{row.resid}_ligand_distance_A'].to_numpy()
            check(cid,f'anchor_{row.resid}_{row.cutoff_A}',abs(np.mean(arr[late]<=row.cutoff_A)-row.late_occupancy))
        bridges=bridges_all.loc[bridges_all.case==cid]
        series=pd.read_csv(local/'proximity_bridge_timeseries.csv')
        bridge_keys=['site','neighbor','anchor']
        summary_keys=set(bridges[bridge_keys].itertuples(index=False,name=None))
        series_keys=set(series[bridge_keys].itertuples(index=False,name=None))
        if summary_keys!=series_keys:
            raise ValueError(f'{cid}: bridge summary and time-series keys differ')
        for row in bridges.itertuples():
            q=series[(series.site==row.site)&(series.neighbor==row.neighbor)&(series.anchor==row.anchor)]
            if len(q)!=1001 or not np.allclose(q.time_ns,times,atol=1e-7,rtol=0):raise ValueError('Missing bridge timepoints')
            joint=(q.site_neighbor_A<=4.5)&(q.neighbor_anchor_A<=4.5)&(q.anchor_ligand_A<=4)
            check(cid,f'bridge_{row.site}_{row.neighbor}_{row.anchor}',abs(joint.to_numpy()[late].mean()-row.late_bridge_and_anchor_ligand_occupancy))
        # Compare original 10ps base series to 100ps local all-receptor fits.
        old=pd.read_csv(base/'selection_rmsd/selection_rmsd.csv')
        check(cid,'base_vs_local_global_fit',np.max(abs(old.selection_rmsd_A.to_numpy()[::10]-f.ligand_rmsd_global_A.to_numpy())),1e-4)
        specs=[('convergence/rmsd_summary.csv','backbone_rmsd_A','receptor_backbone_rmsd_A',10001),('convergence/rg_summary.csv','rg_A','receptor_CA_Rg_A',10001),('selection_rmsd/selection_rmsd.csv','selection_rmsd_A','ligand_heavy_rmsd_A',10001),('bsa/bsa.csv','bsa_A2','receptor_ligand_BSA_A2',201),('bsa/bsa.csv','bsa_fraction_of_initial','BSA_fraction_of_own_initial',201),('contacts/contact_persistence.csv','retained_fraction','own_initial_contact_retained_fraction',201)]
        for path,col,label,n in specs:
            q=pd.read_csv(base/path)
            if len(q)!=n or not np.isfinite(q[col]).all():raise ValueError('Invalid base sample set: '+path)
            check(cid,label+'_time_grid',np.max(abs(q.time_ns-np.linspace(0,100,n))),1e-7)
            v=q.loc[q.time_ns>=50,col].to_numpy();t=q.loc[q.time_ns>=50,'time_ns'].to_numpy()
            for suffix,value in [('mean',v.mean()),('max',v.max()),('slope_per_ns',np.polyfit(t,v,1)[0])]:
                check(cid,label+'_'+suffix,abs(value-base_summary.loc[cid,label+'_'+suffix]),1e-8)
        cp=pd.read_csv(base/'contacts/contact_persistence.csv')
        check(cid,'initial_contact_fraction',np.max(abs(cp.retained_fraction-cp.retained_contacts/cp.initial_contacts)))
        bs=pd.read_csv(base/'bsa/bsa.csv')
        check(cid,'BSA_normalization',np.max(abs(bs.bsa_fraction_of_initial-bs.bsa_A2/bs.bsa_A2.iloc[0])),1e-6)
        com=pd.read_csv(base/'com_distance/com.csv')
        for cluster,q in com.groupby('cluster'):
            if len(q)!=10001:raise ValueError('Invalid COM frame count')
            check(cid,'COM_'+cluster+'_grid',np.max(abs(q.time_ns-np.arange(10001)*.01)),1e-7)
            v=q.loc[q.time_ns>=50,'distance_A'].to_numpy();ts=q.loc[q.time_ns>=50,'time_ns'].to_numpy()
            for suffix,value in [('mean',v.mean()),('max',v.max()),('slope_per_ns',np.polyfit(ts,v,1)[0])]:
                check(cid,'COM_'+cluster+'_'+suffix,abs(value-base_summary.loc[cid,'ligand_COM_to_'+cluster+'_A_'+suffix]),1e-8)
        rf=pd.read_csv(base/'rmsf/rmsf.csv')
        if len(rf)!=406:raise ValueError('Invalid RMSF residue count')
        for suffix,value in [('mean',rf.rmsf_A.mean()),('max',rf.rmsf_A.max())]:check(cid,'RMSF_'+suffix,abs(value-base_summary.loc[cid,'receptor_CA_RMSF_A_'+suffix]),1e-8)
        hb=pd.read_csv(base/'hbond/hbond_summary.csv')
        if len(hb):
            check(cid,'base_hbond_denominator',np.max(abs(hb.occupancy-hb.events/101)))
            check(cid,'base_top_hbond',abs(hb.occupancy.max()-base_summary.loc[cid,'top_geometric_hbond_occupancy']))
        pair_series=pd.read_csv(local/'pair_details/receptor_pair_atoms.csv')
        pair_stats=pairs_all.loc[pairs_all.case==cid]
        for row in pair_stats.itertuples():
            q=pair_series[(pair_series.site==row.site)&(pair_series.neighbor==row.neighbor)]
            if len(q)!=1001:raise ValueError('Invalid atom-pair timepoint count')
            dist=q.distance_A.to_numpy()
            check(cid,f'pair_{row.site}_{row.neighbor}_mean',abs(dist[late].mean()-row.late_mean))
            check(cid,f'pair_{row.site}_{row.neighbor}_occupancy',abs(np.mean(dist[late]<=4)-row.occupancy_4A_late))
        charged=pd.read_csv(local/'pair_details/Arg368_Asp70_geometry.csv')
        charged_stats=pd.read_csv(study/'results/local-summary/Arg368_Asp70_summary.csv')
        if len(charged):
            r=charged_stats[charged_stats.case==cid].iloc[0]
            check(cid,'charged_group_mean',abs(charged.loc[charged.time_ns>=50,'guanidinium_Asp70_distance_A'].mean()-r.late_guanidinium_Asp70_mean_A))
            check(cid,'charged_group_hbond_union',abs(charged.loc[charged.time_ns>=50,'Arg368_Asp70_hbond'].mean()-r.late_any_hbond_occupancy))
        # Independent SVD and broadcast distances use actual rounded PDB examples.
        X,names,resids=pdb_frames(study/'structures'/f'{cid}_six_frames.pdb')
        protein=resids<=406;heavy=np.array([not n.startswith('H') for n in names]);lig=np.flatnonzero(~protein & heavy)
        fits={'global':np.flatnonzero(protein & np.isin(names,['N','CA','C','O']))}
        for label in ['core','pocket']:fits[label]=np.flatnonzero(protein & np.isin(names,['N','CA','C','O']) & np.isin(resids,selections[label+'_resids']))
        for i,t in enumerate([0,10,25,50,75,100]):
            expected=f.iloc[int(t*10)]
            for label,ids in fits.items():
                aligned=independent_fit(X[i,ids],X[0,ids],X[i,lig])
                value=np.sqrt(np.mean(np.sum((aligned-X[0,lig])**2,axis=1)))
                check(cid,f'PDB_{label}_{t}ns',abs(value-expected[f'ligand_rmsd_{label}_A']),PDB_ATOL)
            for resid in [84,368,102,186,189,345]:
                ids=np.flatnonzero(protein & heavy & (resids==resid))
                dist=np.sqrt(np.sum((X[i,ids,None,:]-X[i,lig][None,:,:])**2,axis=-1)).min()
                check(cid,f'PDB_anchor_{resid}_{t}ns',abs(dist-expected[f'res{resid}_ligand_distance_A']),PDB_ATOL)
    from .verify_polar_tables import verify_polar_tables
    check_count += verify_polar_tables(study)
    return check_count, rebuilt


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path,default=STUDY);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=ensure_dir(args.output)
    check_count, rebuilt = verify(args.study)
    pd.DataFrame(rebuilt).to_csv(out/'local-statistics.csv',index=False)
    print(f"{check_count} numerical comparisons; statistics: {out / 'local-statistics.csv'}")


if __name__=='__main__':main()
