#!/usr/bin/env python3
"""Reconstruct paired docking scores from the run and per-mode score tables."""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal
import json
import math
import re
from collections import defaultdict
from pathlib import Path
try:
    from .output_safety import protect_output
except ImportError:
    from output_safety import protect_output


SCORE_FIELDS = ['affinity_kcal_mol', 'rmsd_lb_A', 'rmsd_ub_A']
POSE_FIELDS = ['pose_affinity_kcal_mol', 'pose_rmsd_lb_A', 'pose_rmsd_ub_A']
RECEPTOR_FASTA = {
    'daf37_wt': 'DAF-37_Caenorhabditis_elegans_N2_Q18321.fasta',
    'gpr2_twf154_wt': 'Gpr2_TWF154_TGJ68930.1.fasta',
    'gpr3_twf154_wt': 'Gpr3_TWF154_TGJ67963.1.fasta',
    'gprc_wt': 'GprC_Arthrobotrys_flagrans_A0A437AGY2.fasta',
    'srbc64_wt': 'SRBC-64_Caenorhabditis_elegans_N2_G5ECV4.fasta',
}
RESIDUE_NAMES = dict(zip('ACDEFGHIKLMNPQRSTVWY',
    'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split()))


def read_csv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def finite_number(value, label):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f'Nonfinite numerical value: {label}')
    return number


def load_score_tables(dataset):
    dataset = Path(dataset)
    runs = read_csv(dataset/'reference/runs.csv')
    scores = read_csv(dataset/'reference/pose-scores.csv')
    if not runs or not scores:
        raise ValueError('Nonempty run and score tables are required')
    identities, pairs = set(), set()
    for run in runs:
        key = run['result_id']
        pair = (run['model_rel'], run['ligand'])
        if key in identities or pair in pairs:
            raise ValueError(f'Duplicate docking run identity: {key}')
        if key != '/'.join(pair):
            raise ValueError(f'Run identity mismatch: {key}')
        identities.add(key)
        pairs.add(pair)
    grouped = defaultdict(list)
    for row in scores:
        key = row['result_id']
        if key not in identities:
            raise ValueError(f'Unknown result identity in score table: {key}')
        for field in SCORE_FIELDS:
            finite_number(row[field], f'{key}/{row["pose"]}/{field}')
        if row['pose_present'] not in {'true', 'false'}:
            raise ValueError(f'Invalid pose presence value: {key}')
        for field in POSE_FIELDS:
            if row['pose_present'] == 'true':
                finite_number(row[field], f'{key}/{row["pose"]}/{field}')
            elif row[field] != '':
                raise ValueError(f'Unexpected value for an unstored pose: {key}/{field}')
        grouped[key].append(row)
    for run in runs:
        key = run['result_id']
        count = int(run['logged_scores'])
        if count <= 0 or sorted(int(row['pose']) for row in grouped[key]) != list(range(1, count+1)):
            raise ValueError(f'Incomplete/duplicate score keys: {key}')
    return runs, scores, grouped


def read_pose_models(path):
    """Read complete PDBQT MODEL blocks, scores and their atom identities."""
    models, current, body = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('MODEL'):
            match = re.fullmatch(r'MODEL\s+(\d+)\s*', line)
            if current is not None or match is None:
                raise ValueError(f'Invalid or nested MODEL: {path}')
            current = int(match[1])
            if current < 1 or current in models:
                raise ValueError(f'Duplicate/invalid coordinate mode: {path}/{current}')
            body = []
        elif line.startswith('ENDMDL'):
            if current is None or line.strip() != 'ENDMDL':
                raise ValueError(f'Unmatched ENDMDL: {path}')
            triples = re.findall(r'^REMARK VINA RESULT:\s+(\S+)\s+(\S+)\s+(\S+)\s*$', '\n'.join(body), re.M)
            if len(triples) != 1:
                raise ValueError(f'Expected one score per coordinate mode: {path}/{current}')
            for value in triples[0]:
                finite_number(value, f'{path}/{current}/score')
            atoms = {}
            for atom in body:
                if not atom.startswith(('ATOM  ', 'HETATM')):
                    continue
                serial = int(atom[6:11])
                if serial in atoms:
                    raise ValueError(f'Duplicate coordinate atom: {path}/{current}/{serial}')
                for start in [30, 38, 46]:
                    finite_number(atom[start:start+8], f'{path}/{current}/coordinate')
                atoms[serial] = (atom[12:16].strip(), atom[17:20].strip(), atom[21:22].strip(), int(atom[22:26]), atom.split()[-1])
            if not atoms:
                raise ValueError(f'Coordinate mode has no atoms: {path}/{current}')
            if models and atoms != next(iter(models.values()))['atoms']:
                raise ValueError(f'Atom identity differs between coordinate modes: {path}/{current}')
            models[current] = {'scores': tuple(Decimal(value) for value in triples[0]), 'atoms': atoms}
            current = None
        elif current is not None:
            body.append(line)
        elif line.startswith(('ATOM  ', 'HETATM', 'REMARK VINA RESULT:')):
            raise ValueError(f'Pose data outside a MODEL: {path}')
    if current is not None or not models:
        raise ValueError(f'Missing or incomplete coordinate MODEL: {path}')
    return models


def verify_scientific_bodies(dataset):
    """Check coordinate modes and saved score/contact tables against one another."""
    dataset = Path(dataset)
    runs, _, grouped = load_score_tables(dataset)
    construct_rows = read_csv(dataset/'inputs/constructs.csv')
    constructs = {row['model_rel']: row for row in construct_rows}
    if len(constructs) != len(construct_rows) or set(constructs) != {run['model_rel'] for run in runs}:
        raise ValueError('Construct/model identity mismatch')
    sequences = {}
    for receptor in {row['receptor'] for row in construct_rows}:
        lines = (dataset/'inputs/receptors/fasta'/RECEPTOR_FASTA[receptor]).read_text().splitlines()
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        if sum(line.startswith('>') for line in lines) != 1 or not sequence or not set(sequence) <= RESIDUE_NAMES.keys():
            raise ValueError(f'Invalid receptor sequence: {receptor}')
        sequences[receptor] = sequence
    stored_poses, contacts = 0, 0
    for run in runs:
        key = run['result_id']
        models = read_pose_models(dataset/'reference/vina'/run['model_rel']/run['ligand']/'poses.pdbqt')
        delivered = {int(row['pose']): row for row in grouped[key] if row['pose_present'] == 'true'}
        if list(models) != sorted(delivered) or len(models) != int(run['actual_poses']):
            raise ValueError(f'Coordinate pose count/identity mismatch: {key}')
        for mode, model in models.items():
            if model['scores'] != tuple(Decimal(delivered[mode][field]) for field in POSE_FIELDS):
                raise ValueError(f'Pose score mismatch: {key}/{mode}')
        stored_poses += len(models)
        construct = constructs[run['model_rel']]
        residues = {number for bounds in construct['present_original_residue_ranges'].split(';')
                    for start, end in [map(int, bounds.split('-'))] for number in range(start, end+1)}
        sequence = sequences[construct['receptor']]
        if len(residues) != int(construct['residues']) or not residues <= set(range(1,len(sequence)+1)):
            raise ValueError(f'Construct residue count/range mismatch: {key}')
        residue_labels = {f'{RESIDUE_NAMES[sequence[number-1]]}{number}.A' for number in residues}
        ligand_residues = {f'{name}{number}' + (f'.{chain}' if chain else '')
                           for _, name, chain, number, _ in next(iter(models.values()))['atoms'].values()}
        base = dataset/'reference/prolif'/run['model_rel']/run['ligand']/'ifp'
        with (base/'ifp_counts.csv').open(newline='') as handle:
            table = list(csv.reader(handle))
        with (base/'ifp_presence.csv').open(newline='') as handle:
            presence = list(csv.reader(handle))
        if len(table) != len(models)+4 or len(presence) != len(table) or table[:4] != presence[:4]:
            raise ValueError(f'IFP pose denominator/header mismatch: {key}')
        if len({len(row) for row in table+presence}) != 1:
            raise ValueError(f'Ragged IFP table: {key}')
        if [row[0] for row in table[:3]] != ['ligand', 'protein', 'interaction']:
            raise ValueError(f'Unknown IFP identity headers: {key}')
        expected_indices = list(range(len(models)))
        if any([int(row[0]) for row in values[4:]] != expected_indices for values in [table, presence]):
            raise ValueError(f'IFP pose identity mismatch: {key}')
        features = list(zip(table[0][1:], table[1][1:], table[2][1:]))
        if len(set(features)) != len(features):
            raise ValueError(f'Repeated ligand/protein/interaction feature: {key}')
        if any(ligand not in ligand_residues or protein not in residue_labels or not interaction
               for ligand, protein, interaction in features):
            raise ValueError(f'IFP atom/residue identity mismatch: {key}')
        expected = {}
        for j, (_, protein, interaction) in enumerate(features, 1):
            feature = (protein, interaction)
            if feature in expected:
                raise ValueError(f'Ambiguous ligand contact frequency key: {key}')
            values = [int(row[j]) for row in table[4:]]
            if any(value < 0 for value in values):
                raise ValueError(f'Negative contact geometry count: {key}')
            flags = [row[j].lower() for row in presence[4:]]
            if any(flag not in {'true', 'false'} for flag in flags):
                raise ValueError(f'Missing/invalid contact presence: {key}')
            if any((flag == 'true') != (value > 0) for flag, value in zip(flags, values)):
                raise ValueError(f'Presence/count disagreement: {key}')
            count = sum(value > 0 for value in values)
            expected[feature] = (count, count/len(models), sum(values))
        observed = {}
        for row in read_csv(base/'interaction_frequency.csv'):
            feature = (row['protein'], row['interaction'])
            if feature in observed:
                raise ValueError(f'Duplicate contact frequency key: {key}')
            observed[feature] = (int(row['poses_with_interaction']), finite_number(row['pose_frequency'], key), int(row['total_geometries']))
        if observed.keys() != expected.keys():
            raise ValueError(f'Contact key mismatch: {key}')
        for feature, values in expected.items():
            actual = observed[feature]
            if values[0] != actual[0] or not math.isclose(values[1], actual[1], rel_tol=0, abs_tol=1e-14) or values[2] != actual[2]:
                raise ValueError(f'Contact summary mismatch: {key} {feature}')
        contacts += len(expected)
    return {'docking_runs': len(runs), 'coordinate_poses': stored_poses, 'contact_features': contacts}


def reconstruct(dataset):
    runs, scores, grouped = load_score_tables(dataset)
    paired = defaultdict(dict)
    for run in runs:
        paired[run['model_rel']][run['ligand']] = min(float(row['affinity_kcal_mol']) for row in grouped[run['result_id']])
    rows = []
    for model, ligands in sorted(paired.items()):
        if set(ligands) != {'ascr3_deprotonated', 'ascr18_deprotonated'}:
            raise ValueError(f'Incomplete dual ligand pair: {model}')
        a, b = ligands['ascr3_deprotonated'], ligands['ascr18_deprotonated']
        difference = finite_number(b-a, f'{model}/score difference')
        rows.append({'model_rel': model, 'ascr3_best_kcal_mol': a, 'ascr18_best_kcal_mol': b, 'ascr18_minus_ascr3_kcal_mol': difference})
    return rows, {'docking_runs': len(runs), 'receptor_models': len(paired), 'score_rows': len(scores)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    output = protect_output(args.output_dir, empty_directory=True)
    rows, stats = reconstruct(args.dataset)
    output.mkdir(parents=True, exist_ok=True)
    with (output/'score-comparison.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['model_rel', 'ascr3_best_kcal_mol', 'ascr18_best_kcal_mol', 'ascr18_minus_ascr3_kcal_mol'])
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(stats))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
