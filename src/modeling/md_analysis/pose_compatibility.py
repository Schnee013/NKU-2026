"""Static, fixed-scaffold ligand compatibility in sampled receptor frames.

Scores combine internal MMFF energy and a steric penalty. Exported alternate
ligands are constructed models in the sampled receptor frame.
"""
from pathlib import Path
import json
import warnings
import numpy as np
from .common import ensure_dir
from .polar_io import load_config
from .polar_transition import dump, table, sha, fit


def mapped_names(mol, contract):
    names = contract['heavy_atoms'].copy()
    assert all(a.GetAtomicNum() > 1 for a in list(mol.GetAtoms())[:len(names)])
    for atom in list(mol.GetAtoms())[len(names):]:
        assert atom.GetAtomicNum() == 1 and atom.GetDegree() == 1
        parent = names[atom.GetNeighbors()[0].GetIdx()]
        options = [b if a == parent else a for a,b,_ in contract['bonds']
                   if (a==parent and b.startswith('H')) or (b==parent and a.startswith('H'))]
        names.append(next(n for n in options if n not in names))
    assert sorted(names) == sorted(contract['atom_names'])
    assert {tuple(sorted((names[b.GetBeginAtomIdx()],names[b.GetEndAtomIdx()]))) for b in mol.GetBonds()} == {tuple(sorted((a,b))) for a,b,_ in contract['bonds']}
    return names


def stereo_valid(x, names, contract):
    for center in contract['stereocenters']:
        a,b,c,d = x[[names.index(n) for n in center['neighbors']]]
        volume = np.linalg.det(np.array([a-d,b-d,c-d]))
        if volume * center['source_volume_A3'] <= 0: return False
    return True


def energy_gradient(x, ff, free, moving_heavy, environment, limits, weight):
    flat = tuple(x.ravel())
    energy = ff.CalcEnergy(flat)
    gradient = np.array(ff.CalcGrad(flat)).reshape((-1,3))
    if environment is not None:
        diff = x[moving_heavy,None,:] - environment[None,:,:]
        distances = np.linalg.norm(diff,axis=2)
        overlap = np.maximum(limits-distances,0.)
        energy += weight * (overlap**2).sum()
        gradient[moving_heavy] += (-2*weight*overlap[:,:,None]*diff / np.maximum(distances[:,:,None],1e-8)).sum(1)
    return float(energy), gradient[free].ravel()


def run_compatibility(cfg):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    from scipy.optimize import minimize
    from scipy.spatial import cKDTree
    import MDAnalysis as mda

    root=Path(cfg['output_root']);out=ensure_dir(root/'compatibility')
    ensure_dir(root/'structures')
    if (out/'completion.json').exists():raise RuntimeError('Completed compatibility output exists')
    settings=cfg['compatibility']
    contract_path=Path(settings['chemistry_contract'])
    assert sha(contract_path)==settings['chemistry_contract_sha256']
    contracts=json.loads(contract_path.read_text())
    source_structures=Path(settings['representatives_root'])
    reps=json.loads((source_structures/'representatives.json').read_text())['records']
    frames=[r for r in reps if r['state'].startswith('compat_')]
    source_case=next(c for c in cfg['cases'] if c['id']==cfg['transition']['target_case'])
    source_contract=contracts[settings['source_contract_key']]
    radius={'C':1.70,'N':1.55,'O':1.52,'S':1.80,'H':1.20}
    all_results=[];selections=[];atom_maps=[];pools=[]
    for ligand_spec in settings['ligands']:
        name=ligand_spec['name'];contract=contracts[ligand_spec['contract_key']]
        assert sha(ligand_spec['sdf'])==ligand_spec['sdf_sha256']
        mol=Chem.SDMolSupplier(ligand_spec['sdf'],removeHs=False)[0]
        assert mol is not None and Chem.GetFormalCharge(mol)==-1
        names=mapped_names(mol,contract)
        assert stereo_valid(mol.GetConformer().GetPositions(),names,contract)
        props=AllChem.MMFFGetMoleculeProperties(mol,mmffVariant='MMFF94');assert props is not None
        heavy=np.array([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1])
        radii=np.array([radius[a.GetSymbol()] for a in mol.GetAtoms()])
        mapping=dict(ligand_spec['scaffold_map'])
        # Map scaffold H atoms through their actual bonded parent, never their labels alone.
        for to_name, from_name in list(mapping.items()):
            target_h=[names[a.GetIdx()] for a in mol.GetAtomWithIdx(names.index(to_name)).GetNeighbors() if a.GetAtomicNum()==1]
            source_h=[b if a==from_name else a for a,b,_ in source_contract['bonds'] if (a==from_name and b.startswith('H')) or (b==from_name and a.startswith('H'))]
            assert len(target_h)==len(source_h)
            mapping.update(zip(target_h,source_h))
        frozen=np.array([names.index(n) for n in mapping]);free=np.array([i for i in range(len(names)) if i not in frozen])
        moving_heavy=np.array([i for i in heavy if i in free])
        assert {'O3','O4','O2'}.issubset(mapping)
        for to_name,from_name in mapping.items():atom_maps.append(dict(ligand=name,model_atom=to_name,actual_ascr3_atom=from_name,role='fixed_scaffold'))
        params=AllChem.ETKDGv3();params.randomSeed=settings['seed'];params.numThreads=1;params.pruneRmsThresh=.3
        conf_ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=settings['conformers'],params=params))
        opt=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=1000)
        pool=[(i,float(status[1]),mol.GetConformer(i).GetPositions().copy()) for i,status in zip(conf_ids,opt) if status[0]==0]
        assert pool
        pools.append(dict(ligand=name,requested=settings['conformers'],embedded=len(conf_ids),converged=len(pool)))
        model=Chem.Mol(mol);model.RemoveAllConformers();model.AddConformer(Chem.Conformer(len(names)))
        template_case=next(c for c in cfg['cases'] if c['id']==ligand_spec['template_case'])
        if sha(template_case['topology']) != template_case['expected_topology_sha256']:
            raise ValueError('Template topology differs from configured input')
        template=mda.Universe(template_case['topology'])
        template_ligand=template.select_atoms('not protein')
        assert set(template_ligand.names)==set(names)
        assert {tuple(sorted((template.atoms[int(a)].name,template.atoms[int(b)].name))) for a,b in template_ligand.bonds.indices} == {tuple(sorted((a,b))) for a,b,_ in contract['bonds']}
        for frame in frames:
            path=source_structures/frame['filename'];assert sha(path)==frame['sha256']
            actual=mda.Universe(path);actual_ligand=actual.select_atoms('not protein')
            actual_xyz={a.name:a.position.astype(float) for a in actual_ligand}
            frozen_xyz=np.array([actual_xyz[mapping[names[i]]] for i in frozen])
            environment_group=actual.select_atoms('protein and not name H*')
            env_xyz=environment_group.positions.astype(float)
            env_r=np.array([radius.get(a.name[0],1.7) for a in environment_group])
            # Complete spatial envelope: any ligand atom can move within its molecular reach.
            center=frozen_xyz.mean(0);near=np.linalg.norm(env_xyz-center,axis=1)<settings['environment_radius_A']
            environment,er=env_xyz[near],env_r[near]
            limits=radii[moving_heavy,None]+er[None,:]-settings['soft_overlap_A']
            ff=AllChem.MMFFGetMoleculeForceField(model,props,confId=0)
            def relax(start, use_environment):
                fixed=start.copy();fixed[frozen]=frozen_xyz
                def fun(values):
                    x=fixed.copy();x[free]=values.reshape((-1,3))
                    return energy_gradient(x,ff,free,moving_heavy,environment if use_environment else None,
                                           limits,settings['repulsion_weight'])
                result=minimize(fun,fixed[free].ravel(),jac=True,method='L-BFGS-B',
                                options={'maxiter':settings['max_iterations'],'ftol':1e-9,'gtol':1e-3,'maxls':40})
                x=fixed.copy();x[free]=result.x.reshape((-1,3))
                return x,result
            aligned=[]
            for cid,energy,x in pool:
                start=fit(x,frozen,frozen_xyz)
                start[frozen]=frozen_xyz
                free_pose,result=relax(start,False)
                if not stereo_valid(free_pose,names,contract):continue
                energy=ff.CalcEnergy(tuple(free_pose.ravel()))
                clash=energy_gradient(free_pose,ff,free,moving_heavy,environment,limits,settings['repulsion_weight'])[0]-energy
                aligned.append((clash,cid,energy,free_pose,result.success))
            assert aligned
            reference_energy=min(r[2] for r in aligned)
            starts=sorted(aligned,key=lambda r:r[0])[:settings['refine_candidates']]
            if name==settings['source_ligand']:
                starts.append((0,-1,0,np.array([actual_xyz[n] for n in names]),True))
            candidates=[]
            for _,cid,_,start,_ in starts:
                x,result=relax(start,True)
                energy=ff.CalcEnergy(tuple(x.ravel()))
                d=np.linalg.norm(x[heavy,None,:]-env_xyz[None,:,:],axis=2)
                overlaps=radii[heavy,None]+env_r[None,:]-d
                bond_errors=[]
                for bond in mol.GetBonds():
                    i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
                    expected=props.GetMMFFBondStretchParams(mol,i,j)[2]
                    bond_errors.append(abs(np.linalg.norm(x[i]-x[j])-expected))
                stereo=stereo_valid(x,names,contract)
                # E/Z double-bond identity is checked from regenerated 3D stereo.
                check=Chem.Mol(model)
                for k,v in enumerate(x):check.GetConformer().SetAtomPosition(k,v)
                Chem.AssignStereochemistryFrom3D(check,confId=0,replaceExistingTags=True)
                cip_ok=all(check.GetAtomWithIdx(names.index(s['center'])).GetProp('_CIPCode')==s['source_cip'] for s in contract['stereocenters'])
                source_double=[b.GetStereo() for b in mol.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE and all(mol.GetAtomWithIdx(i).GetAtomicNum()==6 for i in [b.GetBeginAtomIdx(),b.GetEndAtomIdx()])]
                model_double=[b.GetStereo() for b in check.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE and all(check.GetAtomWithIdx(i).GetAtomicNum()==6 for i in [b.GetBeginAtomIdx(),b.GetEndAtomIdx()])]
                double_ok=source_double==model_double
                rec=dict(ligand=name,state=frame['state'],time_ns=frame['time_ns'],conformer=int(cid),
                         converged=bool(result.success),iterations=int(result.nit),
                         minimum_heavy_distance_A=float(d.min()),max_vdw_overlap_A=float(overlaps.max()),
                         pairs_below_2p2_A=int((d<2.2).sum()),scaffold_max_displacement_A=float(np.max(np.abs(x[frozen]-frozen_xyz))),
                         max_bond_error_A=float(max(bond_errors)),stereochemistry_valid=bool(stereo and cip_ok and double_ok),
                         MMFF_internal_kcal=float(energy),fixed_scaffold_reference_kcal=float(reference_energy),
                         internal_strain_kcal=float(energy-reference_energy),
                         objective=float(result.fun))
                rec['geometry_pass']=bool(rec['pairs_below_2p2_A']==0 and rec['max_vdw_overlap_A']<=settings['max_vdw_overlap_A'] and
                    rec['max_bond_error_A']<=settings['max_bond_error_A'] and rec['stereochemistry_valid'] and
                    rec['internal_strain_kcal']<=settings['max_internal_strain_kcal'] and rec['converged'])
                all_results.append(rec);candidates.append((rec,x))
            accepted=[item for item in candidates if item[0]['geometry_pass']]
            selected=min(accepted or candidates,key=lambda item:item[0]['objective'])
            rec,x=selected;selection=dict(rec,accepted_candidates=len(accepted),refined_candidates=len(candidates))
            filename=f"{name}_{frame['state']}_constructed.pdb"
            protein=actual.select_atoms('protein')
            merged=mda.Merge(protein,template_ligand)
            merged.atoms.positions=np.concatenate([protein.positions,np.array([x[names.index(n)] for n in template_ligand.names])])
            with warnings.catch_warnings():
                warnings.simplefilter('ignore');merged.atoms.write(str(root/'structures'/filename))
            reread=mda.Universe(root/'structures'/filename)
            assert len(reread.atoms)==len(protein)+len(names)
            assert np.max(np.abs(reread.select_atoms('protein').positions-protein.positions))<.0006
            selection.update(filename=filename,sha256=sha(root/'structures'/filename),
                             model_type='Constructed ligand compatibility model in an unchanged actual MD receptor frame')
            selections.append(selection)
            # Save chemically explicit SDF as well as the full receptor+ligand PDB.
            for k,v in enumerate(x):model.GetConformer().SetAtomPosition(k,v)
            writer=Chem.SDWriter(str(root/'structures'/filename.replace('.pdb','.sdf')));writer.write(model);writer.close()
            table(out/'candidates.csv',all_results);table(out/'selected_models.csv',selections)
            print('COMPATIBILITY_FRAME',name,frame['state'],'accepted',len(accepted),'of',len(candidates),flush=True)
    table(out/'atom_mapping.csv',atom_maps)
    dump(out/'completion.json',dict(status='complete',config_sha256=cfg['_config_sha256'],module_sha256=sha(__file__),
        rdkit_version=rdBase.rdkitVersion,conformer_pools=pools,selected_models=selections,
        interpretation='Fixed-scaffold compatibility scores combine internal MMFF energy and a steric-overlap penalty.'))


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True)
    args=parser.parse_args();cfg=load_config(Path(args.config))
    run_compatibility(cfg)


if __name__=='__main__':main()
