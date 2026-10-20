#!/usr/bin/env python3
"""Run one saved Vina configuration into a new directory using official bindings.

Config input paths are relative to the configuration file, not the caller cwd.
The Python backend uses the installed Vina version and records its search settings.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,time
from pathlib import Path
try:
    from .output_safety import protect_output
except ImportError:
    from output_safety import protect_output

def sha256(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def parse_config(path):
    values={}
    for raw in Path(path).read_text().splitlines():
        line=raw.split('#',1)[0].strip()
        if not line:continue
        if '=' not in line: raise ValueError(f'Invalid config line: {raw}')
        key,value=map(str.strip,line.split('=',1))
        if key in values: raise ValueError(f'Duplicate config key: {key}')
        values[key]=value
    return values

def run(config,output_dir,*,cpu=None):
    from vina import Vina,__version__
    config=Path(config).resolve();v=parse_config(config)
    output_dir=protect_output(output_dir,empty_directory=True)
    required={'receptor','ligand','center_x','center_y','center_z','size_x','size_y','size_z','seed','exhaustiveness','num_modes','energy_range'}
    missing=required-v.keys()
    if missing:raise ValueError(f'Missing configuration fields: {sorted(missing)}')
    inputs={k:(config.parent/Path(v[k])).resolve() for k in ['receptor','ligand']}
    for key,path in inputs.items():
        if not path.is_file(): raise FileNotFoundError(f'Missing {key}: {path}')
        if not any(x.startswith(('ATOM  ','HETATM')) for x in path.read_text().splitlines()):raise ValueError(f'No PDBQT atoms in {key}')
    center=[float(v[f'center_{a}']) for a in 'xyz'];size=[float(v[f'size_{a}']) for a in 'xyz']
    numeric=center+size
    if not all(math.isfinite(x) for x in numeric) or any(x<=0 for x in size):raise ValueError('Invalid search box')
    settings={'seed':int(v['seed']),'cpu':int(cpu if cpu is not None else v.get('cpu',4)),
       'exhaustiveness':int(v['exhaustiveness']),'num_modes':int(v['num_modes']),'energy_range':float(v['energy_range']),
       'spacing':float(v.get('spacing',0.375)),'min_rmsd':float(v.get('min_rmsd',1.0)),'max_evals':int(v.get('max_evals',0))}
    if any(settings[k]<1 for k in ['cpu','exhaustiveness','num_modes']) or settings['max_evals']<0 or any(settings[k]<=0 for k in ['spacing','min_rmsd','energy_range']):raise ValueError('Invalid Vina search parameter')
    output_dir.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    engine=Vina(sf_name='vina',cpu=settings['cpu'],seed=settings['seed'])
    engine.set_receptor(str(inputs['receptor']));engine.set_ligand_from_file(str(inputs['ligand']))
    engine.compute_vina_maps(center=center,box_size=size,spacing=settings['spacing'])
    engine.dock(exhaustiveness=settings['exhaustiveness'],n_poses=settings['num_modes'],min_rmsd=settings['min_rmsd'],max_evals=settings['max_evals'])
    poses=engine.poses(n_poses=settings['num_modes'],energy_range=settings['energy_range'])
    scores=re.findall(r'^REMARK VINA RESULT:\s+(\S+)\s+(\S+)\s+(\S+)',poses,re.M)
    if not scores:raise RuntimeError('Vina produced no poses')
    (output_dir/'poses.pdbqt').write_text(poses)
    with (output_dir/'scores.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['pose','affinity_kcal_mol','rmsd_lb_A','rmsd_ub_A']);w.writerows((i+1,*row) for i,row in enumerate(scores))
    report={'vina_version':__version__,'backend':'official-python-bindings','inputs':{k:{'name':p.name,'sha256':sha256(p)} for k,p in inputs.items()},'settings':settings,'center_A':center,'size_A':size,'actual_poses':len(scores),'config_sha256':sha256(config),'elapsed_seconds':time.monotonic()-started}
    (output_dir/'metadata.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--cpu',type=int)
    a=p.parse_args()
    try: report=run(a.config,a.output_dir,cpu=a.cpu)
    except (ValueError,FileNotFoundError,FileExistsError) as e:p.error(str(e))
    print(json.dumps(report,indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
