"""Prepare a selected Gpr2 docking construct.

The AF3 route retains original residue numbers 22–163 and 179–406, aligns the
153 specified transmembrane C-alpha atoms, and prepares Gasteiger PDBQT/Meeko
chemistry. The fixed-backbone route removes the designated alanine side chains
without moving retained atoms. Predictions are supplied as local input files.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import shutil
import sys
from .summarize import require, validate_output, construct_residues
from ..af3 import single_protein_sequence

EXPECTED = set(range(22,164)) | set(range(179,407))
TM_RANGES = '30-50,62-83,101-126,144-163,182-201,326-346,358-380'
FIXED_GROUPS = {'baseline':{}, 'N84A_L90A_R345A':{84:'ASN',90:'LEU',345:'ARG'},
                'V110A_R181A_W338A_V368A':{110:'VAL',181:'ARG',338:'TRP',368:'VAL'}}
RESIDUE_NAMES = dict(zip('ACDEFGHIKLMNPQRSTVWY',
    'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split()))


def prediction_sequence(study,identity):
    """Read the retained single-chain request and any matching FASTA records."""
    prediction=identity['prediction']
    sequences=[]
    if prediction.get('job_request_path'):
        request=json.loads((study/prediction['job_request_path']).read_text())
        sequences.append(single_protein_sequence(request))
    for relative in prediction.get('fasta_paths',[]):
        lines=(study/relative).read_text().splitlines()
        require(sum(line.startswith('>') for line in lines)==1,'Expected one FASTA sequence')
        sequences.append(''.join(line.strip() for line in lines if not line.startswith('>')))
    require(sequences and len(set(sequences))==1,'Missing or inconsistent prediction sequences')
    require(set(sequences[0])<=RESIDUE_NAMES.keys(),'Unsupported prediction residue identity')
    return sequences[0]


def construct_sequence(study,identity,document):
    """Bind original residue numbers to the selected record's retained sequence."""
    if identity['batch']=='01':
        baselines=[prediction_sequence(study,r) for r in document['receptors']
                   if r['variant']=='baseline' and r.get('prediction')]
        require(baselines and len(set(baselines))==1,'Missing or inconsistent WT prediction sequences')
        sequence=baselines[0]
    else:
        sequence=prediction_sequence(study,identity)
    require(len(sequence)==document['construct_definition']['prediction_residue_count'],
            'Prediction sequence length mismatch')
    residues={n:RESIDUE_NAMES[sequence[n-1]] for n in construct_residues(document)}
    if identity['batch']=='01':
        require(identity['variant'] in FIXED_GROUPS,'Unsupported fixed-backbone design')
        for number,original in FIXED_GROUPS[identity['variant']].items():
            require(residues.get(number)==original,'WT sequence differs at fixed-backbone mutation site')
            residues[number]='ALA'
    return residues


def atoms(path):
    result={}
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM  '): continue
        key=(line[21],int(line[22:26]),line[26],line[12:16].strip())
        if line[76:78].strip() in ['H','D'] or key[3].lstrip('0123456789').startswith(('H','D')):continue
        require(line[16]==' ' and key not in result,'Alternate or duplicate heavy atom')
        result[key]=line
    require(result,'No protein heavy atoms')
    return result


def trim_full_length(source,destination,expected=EXPECTED,prediction_residue_count=406):
    source_atoms=atoms(source)
    require({k[:3] for k in source_atoms}=={('A',n,' ') for n in range(1,prediction_residue_count+1)},'Unexpected full-length chain A residue inventory')
    retained=[line for key,line in source_atoms.items() if key[1] in expected]
    destination.write_text('\n'.join(retained)+'\nEND\n')
    require({k[:3] for k in atoms(destination)}=={('A',n,' ') for n in expected},'Trimmed residue inventory changed')


def build_fixed_backbone(source,destination,variant,expected=EXPECTED):
    require(variant in FIXED_GROUPS,'Unsupported fixed-backbone design')
    original=atoms(source)
    require({k[:3] for k in original}=={('A',n,' ') for n in expected},'Reference residue inventory mismatch')
    mutations=FIXED_GROUPS[variant]
    for n,resname in mutations.items():
        require(original[('A',n,' ','CA')][17:20]==resname,'Reference sequence mismatch')
    lines=[]
    for key,line in original.items():
        if key[1] in mutations:
            if key[3] not in {'N','CA','C','O','CB','OXT'}:continue
            line=line[:17]+'ALA'+line[20:]
        lines.append(line)
    destination.write_text('\n'.join(lines)+'\nEND\n')
    current=atoms(destination)
    expected={k for k in original if k[1] not in mutations or k[3] in {'N','CA','C','O','CB','OXT'}}
    require(set(current)==expected,'Fixed-backbone atom inventory changed')
    require(all(line[30:54]==original[k][30:54] for k,line in current.items()),'Fixed-backbone retained coordinates moved')


def validate_construct(reference,generated,tolerance_A=0.0018):
    a,b=atoms(reference),atoms(generated)
    require(set(a)==set(b),'Rebuilt heavy-atom inventory differs from reference')
    maximum=0.0
    for key in a:
        require(a[key][17:20]==b[key][17:20],'Rebuilt sequence differs from reference')
        delta=sum((float(a[key][i:i+8])-float(b[key][i:i+8]))**2 for i in [30,38,46])**0.5
        maximum=max(maximum,delta)
    require(maximum<=tolerance_A,f'Rebuilt coordinate frame differs: max {maximum} A')
    return {'heavy_atoms':len(a),'residues':len({k[:3] for k in a}),'max_heavy_atom_delta_A':maximum,'tolerance_A':tolerance_A}


def prepare(study,record_id,reference_pdb,output,cif=None,route="full-prediction",*,construct_pdb,construct_pdbqt):
    study,output=study.resolve(),output.resolve()
    validate_output(study,output)
    receptor_document=json.loads((study/'inputs/receptors.json').read_text())
    expected=construct_residues(receptor_document)
    identities=receptor_document['receptors']
    chosen=[r for r in identities if r['record_id']==record_id]
    require(len(chosen)==1,'Unknown or ambiguous record ID')
    identity=chosen[0]
    require(identity['construct_id']==receptor_document['construct_definition']['id'],'Unknown construct definition')
    require(route in {'full-prediction','saved-construct'},'Unknown preparation route')
    require(construct_pdb.is_file(),'Missing user-supplied construct PDB')
    require(construct_pdbqt.is_file(),'Missing user-supplied construct PDBQT')
    supplied_atoms=atoms(construct_pdb)
    require({k[:3] for k in supplied_atoms}=={('A',n,' ') for n in expected},'Construct residue inventory mismatch')
    sequence=construct_sequence(study,identity,receptor_document)
    require(all(line[17:20]==sequence[key[1]] for key,line in supplied_atoms.items()),
            'Construct sequence does not match selected record: '+record_id)
    if route=='full-prediction':
        require(reference_pdb is not None and reference_pdb.is_file(),'Missing common alignment reference')
    if identity['batch']!='01' and route=='full-prediction':require(cif is not None and cif.is_file(),'AF3 route requires a supplied CIF')
    output.mkdir(parents=True)
    if route=='saved-construct':
        shutil.copyfile(construct_pdb,output/'input.pdb')
    elif identity['batch']=='01':
        build_fixed_backbone(reference_pdb,output/'input.pdb',identity['variant'],expected)
    else:
        from ..docking.cif_to_pdb import convert_one
        convert_one(cif,output/'full-length.pdb')
        trim_full_length(output/'full-length.pdb',output/'trimmed.pdb',expected,receptor_document['construct_definition']['prediction_residue_count'])
        command=[sys.executable,'-m','src.modeling.docking.align_pdb_to_reference',
                 '--moving-pdb',str(output/'trimmed.pdb'),'--reference-pdb',str(reference_pdb),
                 '--output-pdb',str(output/'input.pdb'),'--moving-ranges',TM_RANGES,
                 '--reference-ranges',TM_RANGES,'--report-json',str(output/'alignment.json')]
        subprocess.run(command,check=True)
    coordinates=validate_construct(construct_pdb,output/'input.pdb')
    command=[sys.executable,'-m','src.modeling.docking.prepare_receptor_pdbqt',
             '--receptor-pdb',str(output/'input.pdb'),'--output-pdbqt',str(output/'receptor.pdbqt'),
             '--meeko-executable',sys.executable,'--charge-model','gasteiger',
             '--meeko-arg=-j','--meeko-arg='+str(output/'receptor.json')]
    subprocess.run(command,check=True)
    from meeko import Polymer
    polymer=Polymer.from_json((output/'receptor.json').read_text())
    require(set(polymer.get_valid_monomers())=={f'A:{n}' for n in expected},'Meeko chemical residue inventory mismatch')
    prepared_coordinates=validate_construct(construct_pdbqt,output/'receptor.pdbqt')
    old_atoms,new_atoms=atoms(construct_pdbqt),atoms(output/'receptor.pdbqt')
    require(all(float(old_atoms[k][70:76])==float(new_atoms[k][70:76]) and
                old_atoms[k][77:].strip()==new_atoms[k][77:].strip() for k in old_atoms),
            'Rebuilt heavy-atom charge or AutoDock type differs')
    result={'record_id':record_id,'route':route if route=='saved-construct' else identity['preparation_route'],
            'coordinate_validation':coordinates,'prepared_coordinate_validation':prepared_coordinates,
            'charge_model':'gasteiger','coordinate_precision_A':0.001}
    (output/'preparation.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True);p.add_argument('--record-id',required=True)
    p.add_argument('--construct-pdb',type=Path,required=True,help='User-supplied construct PDB; also the coordinate comparison target')
    p.add_argument('--construct-pdbqt',type=Path,required=True,help='Matching user-supplied PDBQT for coordinate and chemistry comparison')
    p.add_argument('--reference-pdb',type=Path,help='Common alignment/fixed-backbone reference; required for full-prediction')
    p.add_argument('--cif',type=Path,help='User-supplied full-length prediction; required for batches 02–08 with full-prediction')
    p.add_argument('--route',choices=['full-prediction','saved-construct'],default='full-prediction')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:print(json.dumps(prepare(a.study,a.record_id,a.reference_pdb,a.output,a.cif,a.route,
                                construct_pdb=a.construct_pdb,construct_pdbqt=a.construct_pdbqt)))
    except (ValueError,OSError,KeyError,subprocess.CalledProcessError) as exc:raise SystemExit(f'Cannot prepare Gpr2 construct: {exc}')


if __name__=='__main__':main()
