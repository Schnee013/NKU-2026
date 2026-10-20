"""Generate invented docking records for testing parsers and score arithmetic.

Ligand names and batch labels exercise the software schema using synthetic
ethanol coordinates.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def pose_text(score, shift=0.0):
    atom = (f"HETATM{1:5d}  C1  UNL B{1:4d}    "
            f"{shift:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00     0.000 C\n")
    return f"MODEL 1\nREMARK VINA RESULT: {score} 0 0\n{atom}ENDMDL\n"


def sdf_text(shift=0.0):
    atoms = [(shift, 0., 0., "C"), (shift + 1.5, 0., 0., "C"),
             (shift + 2.5, 1., 0., "O")]
    body = "".join(f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3s}"
                   " 0  0  0  0  0  0  0  0  0  0  0  0\n"
                   for x, y, z, symbol in atoms)
    return ("Synthetic ethanol example\n  NKU example\n\n"
            "  3  2  0  0  0  0  0  0  0  0999 V2000\n" + body +
            "  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n$$$$\n")


def write_docking_dataset(root):
    """Two invented one-pose searches, scores -6 and -7 kcal/mol."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    runs, scores = [], []
    for index, ligand in enumerate(("ascr3_deprotonated", "ascr18_deprotonated")):
        score = str(-6 - index)
        identity = "toy/" + ligand
        runs.append(dict(result_id=identity, model_rel="toy", ligand=ligand,
                         logged_scores=1, actual_poses=1))
        scores.append(dict(result_id=identity, pose=1, affinity_kcal_mol=score,
                           rmsd_lb_A=0, rmsd_ub_A=0, pose_affinity_kcal_mol=score,
                           pose_rmsd_lb_A=0, pose_rmsd_ub_A=0, pose_present="true"))
        folder = root / "reference/vina/toy" / ligand
        folder.mkdir(parents=True)
        (folder / "poses.pdbqt").write_text(pose_text(score))
        folder = root / "reference/prolif/toy" / ligand / "ifp"
        folder.mkdir(parents=True)
        header = "ligand,UNL1.B\nprotein,ALA1.A\ninteraction,Hydrophobic\nFrame,\n"
        (folder / "ifp_counts.csv").write_text(header + "0,2\n")
        (folder / "ifp_presence.csv").write_text(header + "0,True\n")
        write_csv(folder / "interaction_frequency.csv", [
            dict(protein="ALA1.A", interaction="Hydrophobic",
                 poses_with_interaction=1, pose_frequency=1., total_geometries=2)])
    write_csv(root / "reference/runs.csv", runs)
    write_csv(root / "reference/pose-scores.csv", scores)
    write_csv(root / "inputs/constructs.csv", [
        dict(model_rel="toy", receptor="daf37_wt", residues=1,
             present_original_residue_ranges="1-1")])
    fasta = root / "inputs/receptors/fasta/DAF-37_Caenorhabditis_elegans_N2_Q18321.fasta"
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">synthetic one-residue alanine sequence\nA\n")
    return root


def write_mutant_dataset(root):
    """Invent the complete four-site schema required by the study summarizer.

All single mutants use S3=-7 and S18=-6; WT uses -6,-6; the double uses
-8,-5. Thus singles have D=1, the double D=3, and its epsilon is (0,1,1).
Coordinates for the double are translated exactly 2 A from both parents.
"""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    receptors, singles, saturation = [], [], []
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    def add(batch, variant):
        rid = f"{batch}/synthetic/{variant}/af0/seed1919"
        receptors.append(dict(record_id=rid, batch=batch, variant=variant,
                              alias=variant, af_model=0, seed=1919,
                              construct_id="toy", baseline_record_id=
                              f"{batch if batch != '08' else '05'}/synthetic/baseline/af0/seed1919"))
        return rid

    baselines = {batch: add(batch, "baseline") for batch in ("05", "06", "07")}
    for site in ("N84", "V110", "V368", "W338"):
        for amino_acid in amino_acids:
            if amino_acid != site[0]:
                rid = add("07" if site == "W338" else "05", site + amino_acid)
                singles.append(rid)
                saturation.append(rid)
    singles.append(add("06", "W338L"))
    combination = add("08", "N84A_V368A")
    receptors[-1].update(
        parent_record_ids={name: f"05/synthetic/{name}/af0/seed1919"
                           for name in ("N84A", "V368A")},
        epsilon_wt_record_ids={"05": baselines["05"], "07": baselines["07"]})
    document = dict(
        construct_definition=dict(id="toy", chain="A", prediction_residue_count=1,
                                  docking_residue_count=1, retained_original_ranges=[[1, 1]],
                                  deleted_original_ranges=[], renumbered=False),
        receptors=receptors,
        comparison_sets=dict(single_record_ids=singles, saturation_record_ids=saturation,
                             combination_record_ids=[combination],
                             single_wt_record_ids=baselines))
    runs, scores = [], []
    fields = ["run_id", "mode", "affinity_kcal_per_mol", "rmsd_lb", "rmsd_ub",
              "pdbqt_affinity_kcal_per_mol", "pdbqt_rmsd_lb", "pdbqt_rmsd_ub"]
    for record in receptors:
        values = ("-6", "-6") if record["variant"] == "baseline" else (
            ("-8", "-5") if record["batch"] == "08" else ("-7", "-6"))
        for ligand, score in zip(("ascr3", "ascr18"), values):
            run_id = record["record_id"] + "/" + ligand
            contact_path = f"results/contacts/batch{record['batch']}/{record['alias']}/{ligand}"
            runs.append(dict(
                **{key: record[key] for key in ("record_id", "batch", "variant", "alias",
                                               "af_model", "baseline_record_id")},
                run_id=run_id, ligand=ligand, pose_count=1, selected_mode=1,
                selected_pose={"tied_minimum_modes": [1]}, contacts_path=contact_path))
            scores.append(dict(zip(fields, [run_id, 1, score, 0, 0, score, 0, 0])))
            folder = root / contact_path
            write_csv(folder / "contact_residues.csv", [
                dict(residue="ALA1.A", original_number=1, interaction_types="VdWContact")])
            (folder / "ifp_presence.csv").write_text(
                "ligand,UNL1.B\nprotein,ALA1.A\ninteraction,VdWContact\nFrame,\n0,True\n")
            folder = root / "reference/structures" / ("batch" + record["batch"]) / record["alias"] / ligand
            folder.mkdir(parents=True)
            shift = 2. if record["batch"] == "08" else 0.
            (folder / "poses.pdbqt").write_text(pose_text(score, shift))
            (folder / "top1.sdf").write_text(sdf_text(shift))
    write_json(root / "inputs/receptors.json", document)
    write_json(root / "inputs/runs.json", dict(
        defaults={"configuration": {"seed": 1919}},
        contact_parameters={"interactions": ["VdWContact"]},
        pose_selection={"analyzed_pose_count": 1, "pose_index_in_ifp": 0},
        runs=runs))
    write_csv(root / "reference/pose-scores.csv", scores)
    ligands = root.parent / "receptor-comparison/inputs/prepared/ligands/deprotonated-minus1/pdbqt"
    ligands.mkdir(parents=True, exist_ok=True)
    for ligand in ("ascr3", "ascr18"):
        with (ligands / (ligand + "_deprotonated.pdbqt")).open("x") as handle:
            handle.write(pose_text("-6"))
    return root


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    from src.modeling.docking.output_safety import protect_output
    output = protect_output(args.output)
    output.mkdir(parents=True, exist_ok=False)
    write_docking_dataset(output / "docking")
    write_mutant_dataset(output / "mutants")
    print(f"Synthetic inputs only: {output}")


if __name__ == "__main__":
    main()
