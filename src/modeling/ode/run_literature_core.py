"""Run configurable constant-input scenarios with the six-state model."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import sys
import scipy
if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from src.modeling.ode.literature_core import (
    LiteratureCoreModel, ROOT, PROJECT_PATH, SOURCE_PATH, STATES, RECONSTRUCTED,
)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parameters',type=Path,default=PROJECT_PATH)
    parser.add_argument('--source',type=Path,default=SOURCE_PATH)
    parser.add_argument('--ligand',type=float,nargs='+',default=[0,.1,1,10])
    parser.add_argument('--hours',type=float,default=8)
    parser.add_argument('--output',type=Path,default=ROOT/'temp/ode-literature-core')
    args=parser.parse_args()
    if not np.isfinite(args.hours) or args.hours<=0:
        parser.error('--hours must be finite and positive')
    model=LiteratureCoreModel(args.parameters,args.source)
    for ligand in args.ligand:
        model.drive(ligand)
    # Resolve aliases before checking protected inputs and packaged references.
    destination = args.output.resolve()
    protected = [ROOT/'studies/signal-model/pathway',
                 model.parameters_path.parent, model.source_path.parent]
    if any(destination == p.resolve() or destination.is_relative_to(p.resolve())
           for p in protected):
        parser.error('--output must be outside input and reference directories')
    # Atomically reserve a fresh directory before writing results.
    try:
        args.output.mkdir(parents=True,exist_ok=False)
    except FileExistsError:
        parser.error('--output already exists; choose a new directory for each run')
    end=args.hours*3600
    times=np.unique(np.r_[0,np.geomspace(end*1e-8,end,200),np.linspace(0,end,201)])
    rows=[]
    for ligand in args.ligand:
        sol=model.solve(ligand,times)
        expanded=model.reconstruct(sol.y)
        if not np.isfinite(expanded).all() or expanded.min() < -1e-8:
            raise RuntimeError('Nonphysical trajectory; output not saved')
        rows.extend([ligand,time,*expanded[:,i],expanded[7,i]/model.pools[3]] for i,time in enumerate(times))
    with (args.output/'timecourses.csv').open('w') as f:
        writer=csv.writer(f)
        writer.writerow(['ligand_normalized','time_s',*RECONSTRUCTED,'MAPK_fraction'])
        writer.writerows(rows)
    record={'environment':{'python':sys.version.split()[0],'numpy':np.__version__,'scipy':scipy.__version__},
            'model_states':STATES,'hours':args.hours,'ligand_normalized':args.ligand,
            'resolved_endpoint_parameters':model.parameters,'fixed_source_rates':model.rates,
            'fixed_source_pools_uM':model.pools.tolist()}
    (args.output/'run.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'output':str(args.output),'rows':len(rows)}))


if __name__=='__main__':
    main()
