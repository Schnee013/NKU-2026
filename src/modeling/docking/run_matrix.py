#!/usr/bin/env python3
"""Run a 35-model by two-ligand matrix from a user-prepared input directory."""
from pathlib import Path
import argparse,csv,json
try:
    from .run_vina import run
    from .output_safety import protect_output
except ImportError:
    from run_vina import run
    from output_safety import protect_output

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset',type=Path,required=True,help='User-prepared directory containing inputs/constructs.csv, prepared receptors/ligands and configs/boxes.tsv');p.add_argument('--output-dir',required=True,type=Path);p.add_argument('--cpu',type=int,default=16);p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    out=protect_output(a.output_dir,empty_directory=True)
    with (a.dataset/'inputs/constructs.csv').open() as f:models=list(csv.DictReader(f))
    boxes={r[0]:r[1:] for r in (line.split() for line in (a.dataset/'configs/boxes.tsv').read_text().splitlines()) if r and not r[0].startswith('#')}
    pairs=[]
    for model in models:
        receptor=a.dataset/'inputs/prepared/receptors/af3/pdbqt'/(model['model_rel']+'.pdbqt')
        for ligand in ['ascr3_deprotonated','ascr18_deprotonated']:
            input_ligand=a.dataset/'inputs/prepared/ligands/deprotonated-minus1/pdbqt'/(ligand+'.pdbqt')
            if not receptor.is_file() or not input_ligand.is_file():raise FileNotFoundError('Missing formal matrix input')
            box=boxes[model['receptor']+'_model_0'];pairs.append((model['model_rel'],ligand,receptor,input_ligand,box))
    if len(models)!=35 or len(pairs)!=70 or len({(r[0],r[1]) for r in pairs})!=70:raise ValueError('Formal matrix identity incomplete or duplicated')
    if a.dry_run:print(json.dumps({'receptor_models':len(models),'ligands':2,'jobs':len(pairs),'mode':'input-discovery-only'}));return 0
    out.mkdir(parents=True,exist_ok=True)
    for model,ligand,receptor,input_ligand,box in pairs:
        target=out/model/ligand;target.parent.mkdir(parents=True,exist_ok=True);config=target.parent/(ligand+'.txt')
        values={'receptor':str(receptor.resolve()),'ligand':str(input_ligand.resolve()),**dict(zip(['center_x','center_y','center_z','size_x','size_y','size_z'],box)),'exhaustiveness':8,'num_modes':9,'energy_range':3,'seed':1919,'cpu':a.cpu,'spacing':0.375,'min_rmsd':1,'max_evals':0}
        config.write_text(''.join(f'{k} = {v}\n' for k,v in values.items()));run(config,target)
    return 0
if __name__=='__main__':raise SystemExit(main())
