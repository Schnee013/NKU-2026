#!/usr/bin/env python3
"""Restore receptor chemistry and pose atom mapping, then compute ProLIF.

Exactly regenerating the saved receptor PDBQT proves that Meeko restoration uses
its docking coordinate frame and parameterization. Pose recovery uses serialized
Meeko SMILES/index mapping, never bond inference from PDBQT coordinates.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re,subprocess,sys
from pathlib import Path
try:
    from .output_safety import protect_output
    from .prolif_ifp import load_receptor,load_poses,summarize_interactions,hydrogen_count
except ImportError:
    from output_safety import protect_output
    from prolif_ifp import load_receptor,load_poses,summarize_interactions,hydrogen_count

def main():
    import prolif as plf
    import pandas as pd
    from meeko import PDBQTMolecule,RDKitMolCreate
    from rdkit import Chem,rdBase
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receptor-pdb',type=Path)
    p.add_argument('--receptor-json',type=Path)
    p.add_argument('--receptor-pdbqt',type=Path,required=True)
    p.add_argument('--poses-pdbqt',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--profile',choices=['cross-receptor','mutant-top1'],default='cross-receptor')
    a=p.parse_args()
    if bool(a.receptor_pdb)==bool(a.receptor_json):p.error('Provide exactly one of --receptor-pdb or --receptor-json')
    out=protect_output(a.output_dir,empty_directory=True)
    for source in [a.receptor_pdb or a.receptor_json,a.receptor_pdbqt,a.poses_pdbqt]:
        if not source.is_file():p.error(f'Missing required input {source.name}')
    out.mkdir(parents=True,exist_ok=True)
    if a.receptor_pdb:
        receptor_json=out/'receptor.json'; regenerated=out/'receptor.pdbqt'
        command=[sys.executable,'-m','meeko.cli.mk_prepare_receptor','--read_pdb',str(a.receptor_pdb),'-p',str(regenerated),'-j',str(receptor_json)]
        subprocess.run(command,check=True)
        if regenerated.read_bytes()!=a.receptor_pdbqt.read_bytes():raise ValueError('Regenerated receptor differs from exact docking input; stop before contact calculation')
        route='PDB to byte-identical Meeko PDBQT and JSON'
    else:
        from meeko import Polymer,PDBQTWriterLegacy
        receptor_json=a.receptor_json
        polymer=Polymer.from_json(receptor_json.read_text())
        rigid,flex=PDBQTWriterLegacy.write_from_polymer(polymer)
        if flex or rigid.encode()!=a.receptor_pdbqt.read_bytes():raise ValueError('Saved Meeko object differs from exact docking PDBQT')
        route='saved Meeko JSON object with byte-identical PDBQT export'
    receptor=load_receptor(receptor_json)
    pdbqt=PDBQTMolecule.from_file(str(a.poses_pdbqt),skip_typing=True)
    sdf,failures=RDKitMolCreate.write_sd_string(pdbqt)
    if failures or not sdf.strip():raise ValueError(f'Failed chemical pose recovery: {failures}')
    posefile=out/'poses.sdf';posefile.write_text(sdf)
    poses=load_poses(posefile)
    models=re.split(r'^MODEL\s+\d+\s*$',a.poses_pdbqt.read_text(),flags=re.M)[1:]
    if not models:models=[a.poses_pdbqt.read_text()]
    if len(poses)!=len(models):raise ValueError('PDBQT and recovered SDF pose count mismatch')
    # Confirm every serialized nonhydrogen coordinate is retained with element identity.
    maxdelta=0.0
    for i,(block,mol) in enumerate(zip(models,poses,strict=True)):
        serialized=[]
        for line in block.splitlines():
            if line.startswith(('ATOM  ','HETATM')) and line.split()[-1] not in ['H','HD','HS']:
                typ=line.split()[-1];element={'A':'C','NA':'N','OA':'O','SA':'S'}.get(typ,typ)
                serialized.append((element,tuple(float(line[x:x+8]) for x in [30,38,46])))
        coords=[(atom.GetSymbol(),tuple(mol.GetConformer().GetAtomPosition(atom.GetIdx()))) for atom in mol.GetAtoms() if atom.GetAtomicNum()!=1]
        if len(serialized)!=len(coords):raise ValueError(f'Heavy atom count mismatch at pose {i+1}')
        for element,coord in serialized:
            choices=[(sum((x-y)**2 for x,y in zip(coord,xyz))**0.5,j) for j,(el,xyz) in enumerate(coords) if el==element]
            if not choices:raise ValueError('Recovered pose element identity mismatch')
            delta,j=min(choices);maxdelta=max(maxdelta,delta);coords.pop(j)
    if maxdelta>0.001:raise ValueError('Chemical recovery changed docking frame')
    expected_charge=sum(Chem.MolFromSmiles(re.search(r'^REMARK SMILES (.+)$',a.poses_pdbqt.read_text(),re.M)[1]).GetAtomWithIdx(i).GetFormalCharge() for i in range(Chem.MolFromSmiles(re.search(r'^REMARK SMILES (.+)$',a.poses_pdbqt.read_text(),re.M)[1]).GetNumAtoms()))
    if any(sum(x.GetFormalCharge() for x in pose.GetAtoms())!=expected_charge for pose in poses):raise ValueError('Recovered pose formal charge mismatch')
    interactions=['Anionic','CationPi','Cationic','HBAcceptor','HBDonor','Hydrophobic','PiCation','PiStacking','VdWContact']
    if a.profile=='mutant-top1':
        if len(receptor.residues)!=370:raise ValueError('Mutant contact profile requires the 370-residue Gpr2 construct')
        energies=[float(x) for x in re.findall(r'^REMARK VINA RESULT:\s+(\S+)',a.poses_pdbqt.read_text(),re.M)]
        selected=min(range(len(energies)),key=lambda i:energies[i])
        poses=[poses[selected]]
    fp=plf.Fingerprint(interactions=interactions if a.profile=='mutant-top1' else None,count=a.profile=='cross-receptor',vicinity_cutoff=6.0)
    fp.run_from_iterable(poses,receptor,n_jobs=1,progress=True)
    if a.profile=='cross-receptor':
        counts=fp.to_dataframe(index_col='Pose')
        counts.to_csv(out/'ifp_counts.csv');counts.astype(bool).to_csv(out/'ifp_presence.csv');summarize_interactions(counts).to_csv(out/'interaction_frequency.csv',index=False)
    else:
        presence=fp.to_dataframe()
        presence.to_csv(out/'ifp_presence.csv')
        residue_types={}
        for (_,residue,interaction),present in presence.any(axis=0).items():
            if present:residue_types.setdefault(str(residue),set()).add(str(interaction))
        with (out/'contact_residues.csv').open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['residue','original_number','interaction_types'])
            for residue in sorted(residue_types,key=lambda r:(int(re.search(r'-?\d+',r)[0]),r)):
                writer.writerow([residue,int(re.search(r'-?\d+',residue)[0]),';'.join(sorted(residue_types[residue]))])
    meta={'preparation_route':route,'contact_profile':a.profile,'interactions':interactions if a.profile=='mutant-top1' else sorted(fp.interactions),'selected_pose_1based':selected+1 if a.profile=='mutant-top1' else None,'prolif_version':plf.__version__,'rdkit_version':rdBase.rdkitVersion,'pandas_version':pd.__version__,'n_poses':len(poses),'protein_residue_count':len(receptor.residues),'ligand_formal_charge':expected_charge,'max_heavy_atom_coordinate_delta_A':maxdelta,'vicinity_cutoff_A':6.0,'inputs':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in {'receptor_pdbqt':a.receptor_pdbqt,'poses_pdbqt':a.poses_pdbqt,'receptor_chemistry':a.receptor_pdb or a.receptor_json}.items()}}
    (out/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
