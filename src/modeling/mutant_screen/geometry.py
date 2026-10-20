"""Recalculate all batch-08 parent/combination pose displacements in a shared frame."""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
from .summarize import require, validate_output


def geometry(a,b):
    from rdkit import Chem
    molecules=[]
    for path in [a,b]:
        supplied=list(Chem.SDMolSupplier(str(path),removeHs=False))
        require(len(supplied)==1 and supplied[0] is not None,'Expected exactly one valid selected SDF pose')
        molecules.append(Chem.RemoveHs(supplied[0]))
    require(Chem.MolToSmiles(molecules[0],isomericSmiles=True)==Chem.MolToSmiles(molecules[1],isomericSmiles=True),'Ligand stereochemical graph mismatch')
    matches=molecules[1].GetSubstructMatches(molecules[0],uniquify=False,useChirality=True,maxMatches=10000)
    require(bool(matches) and len(matches)<10000,'Atom-map search absent or truncated')
    aa=molecules[0].GetConformer().GetPositions();bb=molecules[1].GetConformer().GetPositions()
    rmsd=min(math.sqrt(sum(sum((aa[i][k]-bb[j][k])**2 for k in range(3)) for i,j in enumerate(match))/len(aa)) for match in matches)
    return {'heavy_atom_rmsd_no_fit_A':rmsd,'centroid_shift_A':float(sum((aa.mean(axis=0)-bb.mean(axis=0))**2)**0.5),'atom_count':len(aa)}


def reconstruct(study):
    data=json.loads((study/'inputs/receptors.json').read_text())
    records={r['record_id']:r for r in data['receptors']}
    require(len(records)==len(data['receptors']),'Duplicate receptor ID')
    def path(record,ligand):
        return study/'reference/structures'/('batch'+record['batch'])/record['alias']/ligand/'top1.sdf'
    rows=[]
    for rid in data['comparison_sets']['combination_record_ids']:
        require(rid in records,'Unknown combination ID')
        record=records[rid]
        if record['batch']!='08':continue
        require(len(record['parent_record_ids'])==2,'Expected two parents')
        for name,parent_id in record['parent_record_ids'].items():
            require(parent_id in records,'Unknown parent ID')
            parent=records[parent_id]
            require(parent['variant']==name,'Parent identity mismatch')
            require(parent['batch']==('07' if name.startswith('W338') else '05'),'Parent batch mismatch')
            for ligand in ['ascr3','ascr18']:
                values=geometry(path(parent,ligand),path(record,ligand))
                rows.append({'record_id':record['record_id'],'parent_record_id':parent['record_id'],'ligand':ligand,**values})
    require(len({(r['record_id'],r['parent_record_id'],r['ligand']) for r in rows})==len(rows),'Duplicate geometry comparison')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:
        validate_output(a.study,a.output)
        rows=reconstruct(a.study)
        a.output.mkdir(parents=True)
        with (a.output/'parent-pose-geometry.csv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        (a.output/'geometry.json').write_text(json.dumps({'comparisons':len(rows),'definition':'stereochemical graph; minimum over equivalent atom maps; shared receptor frame; no ligand fitting'},indent=2)+'\n')
        print(a.output/'parent-pose-geometry.csv')
    except (ValueError,OSError,KeyError) as exc:raise SystemExit(f'Cannot compare ligand geometries: {exc}')


if __name__=='__main__':main()
