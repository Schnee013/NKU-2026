"""Resolve the public receptor comparison dataset independently of caller cwd."""
from pathlib import Path
REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
def resolve_paths(root=None):
    root=Path(root or REPOSITORY_ROOT).resolve()
    data=root/'studies/receptor-engineering/receptor-comparison'
    return {'ligand_inputs':data/'inputs/ligands',
      'ligand_pdbqt':data/'inputs/prepared/ligands/deprotonated-minus1/pdbqt','vina_results':data/'reference/vina',
      'prolif_results':data/'reference/prolif','box_config':data/'configs/boxes.tsv'}
