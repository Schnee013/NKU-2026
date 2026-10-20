"""Recalculate complete Gpr2 comparisons from saved per-pose scores and contacts.

Only the supplied study directory is read. A fresh output directory is required;
reference results and inputs are never overwritten. No docking jobs are submitted.
"""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal as D
import json
from pathlib import Path
import re
import statistics
if __package__ in {None, ''}:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from src.modeling.docking.summarize import read_pose_models


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text())


POSE_COLUMNS = ['run_id','mode','affinity_kcal_per_mol','rmsd_lb','rmsd_ub',
                'pdbqt_affinity_kcal_per_mol','pdbqt_rmsd_lb','pdbqt_rmsd_ub']


def construct_residues(document):
    definition = document['construct_definition']
    retained = [n for start,end in definition['retained_original_ranges'] for n in range(start,end+1)]
    deleted = {n for start,end in definition['deleted_original_ranges'] for n in range(start,end+1)}
    require(len(retained) == len(set(retained)) == definition['docking_residue_count'], 'Invalid construct residue count')
    require(not (set(retained) & deleted) and set(retained) | deleted == set(range(1,definition['prediction_residue_count']+1)), 'Invalid retained/deleted residue partition')
    require(definition['chain'] == 'A' and definition['renumbered'] is False, 'Construct must retain chain A and original residue numbers')
    return set(retained)


def load_runs(study):
    document = read_json(study/'inputs/runs.json')
    defaults = document.get('defaults', {})
    return [{**defaults, **run,
             'configuration': {**defaults.get('configuration', {}), **run.get('configuration', {})},
             'contact_parameters': {**document['contact_parameters'], **run.get('contact_parameters', {})},
             'selected_pose': {**document['pose_selection'], **run['selected_pose']}}
            for run in document['runs']]


def read_csv(path, columns):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        require(set(reader.fieldnames or []) == set(columns), f'Unexpected columns: {path.name}')
        return list(reader)


def entries(contact):
    return {(int(number), kind) for number, value in contact.items() for kind in value['types']}


def changes(before, after):
    a, b = entries(before), entries(after)
    def labels(values, contacts):
        return [contacts[str(number)]['residue'] + ':' + kind for number, kind in sorted(values)]
    return {
        'gained': labels(b-a, after), 'lost': labels(a-b, before),
        'gained_residues': [after[n]['residue'] for n in sorted(after.keys()-before.keys(), key=int)],
        'lost_residues': [before[n]['residue'] for n in sorted(before.keys()-after.keys(), key=int)],
        'identity_changes': [{'number': int(n), 'wt': before[n]['residue'], 'mutant': after[n]['residue']}
                             for n in sorted(before.keys() & after.keys(), key=int)
                             if before[n]['residue'] != after[n]['residue']],
    }


def read_contacts(folder, run):
    rows = read_csv(folder/'contact_residues.csv', ['residue','original_number','interaction_types'])
    contacts = {r['original_number']: {'residue': r['residue'], 'types': sorted(r['interaction_types'].split(';'))} for r in rows}
    require(len(contacts) == len(rows), f'Duplicate residue: {folder}')
    with (folder/'ifp_presence.csv').open(newline='') as h:
        table = list(csv.reader(h))
    require(len(table) == 5, f'Expected one analyzed pose: {folder}')
    require(len({len(r) for r in table}) == 1, f'Ragged fingerprint: {folder}')
    observed = {}
    for i, cell in enumerate(table[4][1:], 1):
        require(cell.lower() in ['true','false'], f'Missing fingerprint value: {folder}')
        if cell.lower() == 'true':
            observed.setdefault(table[1][i], set()).add(table[2][i])
    require(observed == {v['residue']: set(v['types']) for v in contacts.values()}, f'Contact join mismatch: {folder}')
    parameters = run['contact_parameters']
    pose_count = len(table)-4
    require(pose_count == run['selected_pose']['analyzed_pose_count'] == 1, f'Contact denominator mismatch: {folder}')
    require(int(table[4][0]) == run['selected_pose']['pose_index_in_ifp'], f'Contact pose index mismatch: {folder}')
    keys = [(table[0][i],table[1][i],table[2][i]) for i in range(1,len(table[0]))]
    require(len(set(keys)) == len(keys), f'Duplicate fingerprint entry: {folder}')
    require(all(k[2] in parameters['interactions'] for k in keys), f'Unknown interaction type: {folder}')
    require(all(re.fullmatch(r'[A-Z]+'+re.escape(number)+r'\.A', value['residue']) for number,value in contacts.items()), f'Original residue number mismatch: {folder}')
    return contacts


def metrics(record, baseline):
    a, b = D(record['S3']) - D(baseline['S3']), D(record['S18']) - D(baseline['S18'])
    return {'deltaS3': str(a), 'deltaS18': str(b), 'D': str(b-a)}


def equal_number(a, b, label):
    require(a is not None and b is not None, f'Missing numerical value: {label}')
    require(D(str(a)) == D(str(b)), f'Numeric mismatch: {label}: {a} != {b}')


def load_pose_scores(study, runs):
    by_id = {run['run_id']: run for run in runs}
    require(len(by_id) == len(runs), 'Duplicate run ID')
    rows = read_csv(study/'reference/pose-scores.csv', POSE_COLUMNS)
    by_pose = {}
    for row in rows:
        key = (row['run_id'], int(row['mode']))
        require(key not in by_pose and row['run_id'] in by_id, f'Duplicate or unbound pose ID: {key}')
        for field in POSE_COLUMNS[2:]:
            require(D(row[field]).is_finite(), f'Missing/nonfinite {field}: {key}')
        by_pose[key] = row
    for run in runs:
        modes = sorted(mode for run_id, mode in by_pose if run_id == run['run_id'])
        require(run['pose_count'] > 0 and modes == list(range(1, run['pose_count']+1)),
                f'Score-row mode set mismatch: {run["run_id"]}')
    return rows, by_pose


def select_minimum(run, modes):
    minimum = min(D(row['affinity_kcal_per_mol']) for row in modes)
    tied = sorted(int(row['mode']) for row in modes if D(row['affinity_kcal_per_mol']) == minimum)
    require(tied == run['selected_pose']['tied_minimum_modes'], f'Minimum-score tie set mismatch: {run["run_id"]}')
    require(tied[0] == run['selected_mode'], f'Minimum-score selection mismatch: {run["run_id"]}')
    return next(row for row in modes if int(row['mode']) == tied[0])


def verify_scientific_bodies(study, runs):
    """Compare retained ligand poses, scores and contact-table identities."""
    _, saved_scores = load_pose_scores(study, runs)
    coordinate_poses = 0
    for run in runs:
        tag, alias, lig = run['batch'], run['alias'], run['ligand']
        folder = study/'reference/structures'/('batch'+tag)/alias/lig
        ligand = study.parent/'receptor-comparison/inputs/prepared/ligands/deprotonated-minus1/pdbqt'/(lig+'_deprotonated.pdbqt')
        for path in [ligand,folder/'poses.pdbqt',folder/'top1.sdf']:
            if not path.is_file():
                raise FileNotFoundError(f'Missing scientific input: {path}')
        read_contacts(study/run['contacts_path'], run)
        models = read_pose_models(folder/'poses.pdbqt')
        require(sorted(models) == list(range(1,run['pose_count']+1)), 'Coordinate pose count/identity mismatch: '+run['run_id'])
        modes = [saved_scores[(run['run_id'], mode)] for mode in sorted(models)]
        selected = select_minimum(run, modes)
        for mode, model in models.items():
            row = saved_scores[(run['run_id'],mode)]
            for column, value in zip(['affinity_kcal_per_mol','rmsd_lb','rmsd_ub'], model['scores']):
                require(value == D(row['pdbqt_'+column]), 'PDBQT per-pose value mismatch: '+run['run_id']+'/'+str(mode)+'/'+column)
            require(model['scores'][0] == D(row['affinity_kcal_per_mol']), 'Cross-source pose affinity mismatch: '+run['run_id'])
        coordinate_minimum = min(model['scores'][0] for model in models.values())
        coordinate_ties = sorted(mode for mode, model in models.items() if model['scores'][0] == coordinate_minimum)
        require(coordinate_ties == run['selected_pose']['tied_minimum_modes'] and coordinate_ties[0] == run['selected_mode'],
                'Coordinate minimum-score pose mismatch: '+run['run_id'])
        equal_number(coordinate_minimum, selected['affinity_kcal_per_mol'], run['run_id'])
        coordinate_poses += len(models)
    return {'ligand_runs': len(runs), 'coordinate_poses': coordinate_poses}


def reconstruct(study):
    study = study.resolve()
    runs = load_runs(study)
    receptor_document = read_json(study/'inputs/receptors.json')
    construct_residues(receptor_document)
    receptors = receptor_document['receptors']
    identities = {r['record_id']:r for r in receptors}
    comparison_sets = receptor_document['comparison_sets']
    for key in ['single_record_ids','saturation_record_ids','combination_record_ids']:
        ids = comparison_sets[key]
        require(len(ids) == len(set(ids)) and set(ids) <= identities.keys(), f'Invalid comparison collection: {key}')
    for batch, rid in comparison_sets['single_wt_record_ids'].items():
        require(rid in identities and identities[rid]['batch'] == batch and identities[rid]['variant'] == 'baseline', 'Unknown single-mutant WT')
    require(all(r['construct_id'] == receptor_document['construct_definition']['id'] for r in receptors), 'Unknown construct definition')
    by_id = {r['run_id']: r for r in runs}
    require(len(by_id) == len(runs) and len({r['record_id'] for r in receptors}) == len(receptors), 'Duplicate run or receptor ID')
    poses, by_pose = load_pose_scores(study, runs)
    records = {}
    for identity in receptors:
        rid = identity['record_id']
        record = {k: identity[k] for k in ['record_id','batch','variant','alias','af_model','seed']}
        record.update({k:identity[k] for k in ['module','role','mutations'] if k in identity})
        record['rmsd_A'] = identity.get('alignment', {}).get('ca_rmsd_after_angstrom')
        for ligand, score_key in [('ascr3','S3'),('ascr18','S18')]:
            run = by_id.get(rid+'/'+ligand)
            require(run is not None, f'Missing ligand arm: {rid}/{ligand}')
            require(run['ligand'] == ligand, f'Ligand identity mismatch: {rid}/{ligand}')
            for field in ['record_id','batch','variant','alias','af_model']:
                require(run[field] == identity[field], f'Identity mismatch: {rid}/{ligand}/{field}')
            modes = [p for (runid, _), p in by_pose.items() if runid == run['run_id']]
            best = select_minimum(run, modes)
            require(int(run['configuration']['seed']) == identity['seed'], f'Seed mismatch: {rid}')
            record[score_key] = best['affinity_kcal_per_mol']
            contact = read_contacts(study/run['contacts_path'], run)
            record[ligand] = {'score': record[score_key], 'model': int(best['mode']), 'score_rows': len(modes),
                              'contacts': contact, 'count': len(contact), 'entry_count': len(entries(contact))}
        record['G'] = str(D(record['S18'])-D(record['S3']))
        records[rid] = record
    expected_run_ids = {r['record_id']+'/'+l for r in receptors for l in ['ascr3','ascr18']}
    require(set(by_id) == expected_run_ids, 'Unbound run inventory entries')
    wt = {r['batch']: r for r in records.values() if r['variant'] == 'baseline'}
    for r in records.values():
        baseline_id = identities[r['record_id']]['baseline_record_id']
        require(all(baseline_id == by_id[r['record_id']+'/'+lig]['baseline_record_id'] for lig in ['ascr3','ascr18']), 'Ligand arms use different WT baselines')
        require(baseline_id in records and records[baseline_id]['variant'] == 'baseline', 'Unknown WT baseline')
        baseline = records[baseline_id]
        r['baseline_record_id'] = baseline['record_id']
        r['reference_sensitivity'] = {b: metrics(r, w) for b, w in wt.items()}
        r.update(metrics(r, baseline))
        for lig in ['ascr3','ascr18']:
            c, w = r[lig], baseline[lig]
            a, b = entries(w['contacts']), entries(c['contacts'])
            c['changes'] = changes(w['contacts'], c['contacts'])
            c['jaccard_to_WT'] = len(a&b)/len(a|b) if a|b else 1.0
            c['gained_entries'], c['lost_entries'] = len(b-a), len(a-b)
    combinations = []
    for rid in comparison_sets['combination_record_ids']:
        identity = identities[rid]
        if identity['batch'] != '08':
            continue
        r = records[rid]
        r['parents'], r['epsilon'] = {}, {}
        parent_records = []
        for name, parent_id in identity['parent_record_ids'].items():
            require(parent_id in records, f'Unknown parent: {parent_id}')
            parent = records[parent_id]
            require(parent['variant'] == name, f'Wrong parent design: {name}')
            require(parent['batch'] == ('07' if name.startswith('W338') else '05'), f'Wrong parent batch: {name}')
            parent_records.append(parent)
            dx, dy = D(r['S3'])-D(parent['S3']), D(r['S18'])-D(parent['S18'])
            values = {'score_change3': str(dx), 'score_change18': str(dy), 'D_change': str(dy-dx)}
            values['record_id'], values['contacts'] = parent['record_id'], {}
            for lig in ['ascr3','ascr18']:
                a, b = entries(parent[lig]['contacts']), entries(r[lig]['contacts'])
                c = {**changes(parent[lig]['contacts'],r[lig]['contacts']), 'retained_entries': len(a&b),
                     'parent_entries': len(a), 'combination_entries': len(b), 'jaccard': len(a&b)/len(a|b)}
                values['contacts'][lig] = c
            r['parents'][name] = values
        require(len(parent_records) == 2, f'Expected two parents: {r["record_id"]}')
        for batch, wt_id in identity['epsilon_wt_record_ids'].items():
            require(batch in wt and wt[batch]['record_id'] == wt_id, 'Unknown epsilon WT')
            a, b = parent_records
            x = D(r['S3'])-D(a['S3'])-D(b['S3'])+D(wt[batch]['S3'])
            y = D(r['S18'])-D(a['S18'])-D(b['S18'])+D(wt[batch]['S18'])
            r['epsilon'][batch] = {'S3': str(x), 'S18': str(y), 'D': str(y-x)}
        combinations.append(r)
    singles = [records[rid] for rid in comparison_sets['single_record_ids']]
    require(all(re.fullmatch(r'[A-Z]\d+[A-Z]',r['variant']) and r['batch'] in ['04','05','06','07'] for r in singles),'Non-single in single-mutant comparison')
    actual_singles = {r['record_id'] for r in records.values() if r['batch'] in ['04','05','06','07'] and re.fullmatch(r'[A-Z]\d+[A-Z]',r['variant'])}
    require({r['record_id'] for r in singles} == actual_singles,'Single-mutant collection incomplete')
    panel = [records[rid] for rid in comparison_sets['saturation_record_ids']]
    require(len({r['record_id'] for r in panel}) == len(panel), 'Duplicate saturation design')
    stats = {}
    for site in ['N84','V110','V368','W338']:
        rr = [r for r in panel if r['variant'][:-1] == site]
        require(len(rr) == 19 and len({r['variant'][-1] for r in rr}) == 19, f'Incomplete saturation set: {site}')
        ds = [D(r['D']) for r in rr]
        stats[site] = {'n':len(rr),'D_positive':sum(d>0 for d in ds),'D_negative':sum(d<0 for d in ds),'D_zero':sum(d==0 for d in ds),
                       'D_median':str(statistics.median(ds)),'D_min':str(min(ds)),'D_max':str(max(ds)),
                       'both_improve':sum(D(r['deltaS3'])<0 and D(r['deltaS18'])<0 for r in rr),
                       'only3_improves':sum(D(r['deltaS3'])<0 and D(r['deltaS18'])>0 for r in rr),
                       'only18_improves':sum(D(r['deltaS3'])>0 and D(r['deltaS18'])<0 for r in rr),
                       'both_worse':sum(D(r['deltaS3'])>0 and D(r['deltaS18'])>0 for r in rr)}
        stats[site]['contacts']={}
        for lig in ['ascr3','ascr18']:
            cs=[r[lig] for r in rr]
            def has(c,n,kind=None):
                return str(n) in c['contacts'] and (kind is None or kind in c['contacts'][str(n)]['types'])
            cstats={'median_residues': statistics.median(c['count'] for c in cs),
                    'range_residues':[min(c['count'] for c in cs),max(c['count'] for c in cs)],
                    'median_jaccard':statistics.median(c['jaccard_to_WT'] for c in cs),
                    'direct_site_count':sum(has(c,int(site[1:])) for c in cs),
                    'R189_Anionic':sum(has(c,189,'Anionic') for c in cs),
                    '338_Hydrophobic':sum(has(c,338,'Hydrophobic') for c in cs),
                    'R189_Anionic_and_338_Hydrophobic':sum(has(c,189,'Anionic') and has(c,338,'Hydrophobic') for c in cs),
                    'residue_frequency':{str(n):sum(has(c,n) for c in cs) for n in sorted({int(n) for c in cs for n in c['contacts']})}}
            stats[site]['contacts'][lig]=cstats
    rankings={}
    for batch in sorted(wt):
        rankings[batch]=[r['record_id'] for r in sorted((r for r in records.values() if r['batch']==batch and r['variant']!='baseline'),key=lambda r:(-D(r['D']),r['variant'],r['af_model']))]
    rankings['08']=[r['record_id'] for r in sorted(combinations,key=lambda r:(-D(r['D']),r['variant']))]
    return {'records':list(records.values()),'within_batch_D_rankings':rankings,'saturation_statistics':stats,
            'counts':{'receptors':len(receptors),'ligand_runs':len(runs),'score_rows':len(poses),'single_structure_records':len(singles),
                      'unique_single_designs':len({r['variant'] for r in singles}),'source_WT_single_batches':len(comparison_sets['single_wt_record_ids']),
                      'saturation_designs':len(panel),'batch08_combinations':len(combinations),
                      'batch08_parents':len({p['record_id'] for r in combinations for p in r['parents'].values()})}}


def validate_output(study, output):
    study, output = study.resolve(), output.resolve()
    require(study != output and study not in output.parents and output not in study.parents, 'Output must be outside the study inputs and references')
    require(not output.exists(), 'Output directory already exists; select a fresh directory')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    try:
        validate_output(a.study, a.output)
        result = reconstruct(a.study)
        a.output.mkdir(parents=True)
        (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
        rows = [{k:r[k] for k in ['record_id','batch','variant','af_model','seed','baseline_record_id','S3','S18','deltaS3','deltaS18','D']} for r in result['records']]
        with (a.output/'comparisons.csv').open('w', newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        print(json.dumps(result['counts']))
    except (ValueError,KeyError,OSError,ArithmeticError) as exc:
        raise SystemExit(f'Cannot reconstruct mutation-screen summary: {exc}')


if __name__ == '__main__':
    main()
