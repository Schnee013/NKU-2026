"""Generate a toy trajectory and run the seven MD analysis programs.

Four stationary alanines and a three-atom artificial partner are placed by
coordinates. The partner moves 0, 1 and 2 A. The three frames provide known
distances for checking trajectory analysis.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DRIVERS = ("convergence", "selection_rmsd", "hbond", "com_distance", "rmsf", "bsa", "contacts")


def create_inputs(destination):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    processed = destination / "processed/rep1"
    processed.mkdir(parents=True)
    lines, coordinates = [], []
    for resid, (y, z) in enumerate(((0, 0), (5, 0), (0, 5), (5, 5)), 1):
        for name, x, dy, element in (("N", 0, 0, "N"), ("H", -1, 0, "H"),
                                     ("CA", 1.5, 0, "C"), ("C", 2.5, 1, "C"),
                                     ("O", 3.5, 1, "O"), ("CB", 1.5, -1, "C")):
            xyz = (x, y + dy, z)
            coordinates.append(xyz)
            lines.append(f"ATOM  {len(coordinates):5d} {name:^4s} ALA A{resid:4d}    "
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}\n")
    for name, x, element in (("O1", -2.8, "O"), ("C1", -3.8, "C"), ("O2", -4.8, "O")):
        coordinates.append((x, 0., 0.))
        lines.append(f"HETATM{len(coordinates):5d} {name:^4s} LIG B{5:4d}    "
                     f"{x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {element:>2s}\n")
    topology = processed / "protein_ref.pdb"
    topology.write_text("".join(lines) + "END\n")
    universe = mda.Universe(str(topology))
    original = np.array(coordinates)
    with XTCWriter(str(processed / "prod_protein.xtc"), n_atoms=len(original)) as writer:
        for index in range(3):
            xyz = original.copy()
            xyz[-3:, 0] -= index
            universe.atoms.positions = xyz
            universe.trajectory.ts.time = index * 50000.
            writer.write(universe.atoms)
    config = dict(default="synthetic", systems={"synthetic": dict(
        label="Synthetic geometry", result_dir=str(destination), interaction_partner_role="ligand",
        backbone_selection="protein and backbone", ca_selection="protein and name CA",
        summary_window_ns=[50, 100], occupancy_window_ns=[0, 100], rmsf_window_ns=[50, 100],
        selection_rmsd_selection="resname LIG", plot_segids=["A"],
        receptor=dict(name="Toy receptor", selection="protein", interface_selection="protein",
                      interface_clusters={"whole_receptor": "protein"}),
        partner=dict(name="Toy partner", selection="resname LIG", com_selection="resname LIG"))})
    path = destination / "systems.json"
    path.write_text(json.dumps(config, indent=2) + "\n")
    return path


def run_drivers(config, output, cache):
    env = dict(os.environ, PYTHONPATH=str(ROOT), MPLCONFIGDIR=str(cache),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MD_CONFIG=str(config), MD_SYSTEM="synthetic",
               MD_OUTPUT_ROOT=str(output), MD_REPS="1")
    # Do not let a caller's real-data location override the generated inputs.
    env.pop("MD_PAYLOAD_ROOT", None)
    for name in DRIVERS:
        extra = ["--step", "1"] if name == "hbond" else (
            ["--stride", "1"] if name in ("bsa", "contacts") else [])
        subprocess.run([sys.executable, "-m", "src.modeling.md_analysis.analyze_" + name, *extra],
                       env=env, check=True)


def check_outputs(output):
    folder = Path(output) / "synthetic/analysis"
    rmsd = pd.read_csv(folder / "selection_rmsd/selection_rmsd.csv")
    np.testing.assert_allclose(rmsd.selection_rmsd_A, [0, 1, 2], atol=1e-4)
    backbone = pd.read_csv(folder / "convergence/rmsd_summary.csv")
    np.testing.assert_allclose(backbone.backbone_rmsd_A, 0, atol=1e-6)
    rmsf = pd.read_csv(folder / "rmsf/rmsf.csv")
    np.testing.assert_allclose(rmsf.rmsf_A, 0, atol=1e-6)
    contacts = pd.read_csv(folder / "contacts/contact_persistence.csv")
    np.testing.assert_allclose(contacts.retained_fraction, [1, 1, 0])
    hbond = pd.read_csv(folder / "hbond/hbond_summary.csv")
    if not ((hbond.events == 1) & np.isclose(hbond.occupancy, 1 / 3)).any():
        raise AssertionError("Expected one H-bond event in three sampled frames")
    com = pd.read_csv(folder / "com_distance/com.csv")
    if not np.all(np.diff(com.distance_A) > 0):
        raise AssertionError("Toy partner COM distance must increase")
    bsa = pd.read_csv(folder / "bsa/bsa.csv")
    if not ((bsa.bsa_A2 >= -1e-3).all() and bsa.bsa_A2.iloc[0] > bsa.bsa_A2.iloc[-1]):
        raise AssertionError("Expected decreasing, nonnegative buried area")
    return dict(scope="synthetic software checks only", drivers=list(DRIVERS),
                expected_partner_rmsd_A=[0, 1, 2], expected_retained_fraction=[1, 1, 0],
                expected_hbond_occupancy=1 / 3, status="passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    from src.modeling.docking.output_safety import protect_output
    output = protect_output(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    config = create_inputs(output / "input")
    run_drivers(config, output / "analysis", output / "matplotlib-cache")
    checks = check_outputs(output / "analysis")
    (output / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks))


if __name__ == "__main__":
    main()
