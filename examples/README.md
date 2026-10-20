# Synthetic examples

Run from the repository root after installing the relevant [dependencies](../environments/README.md). Every output directory below must be new. Examples run locally, without downloading project datasets.

## MD geometry

```sh
python -m examples.md_geometry --output runs/examples/md-geometry
```

The generator creates four stationary alanines and an artificial three-atom partner at three coordinate frames. Frame labels are 0, 50 and 100 ns; no physical time integration takes place. The partner translates by 0, 1 and 2 Å. The example runs convergence, selection RMSD, H-bond, COM-distance, RMSF, BSA and contact analysis.

`checks.json` is written only after numerical checks pass: partner RMSD `[0, 1, 2]` Å, stationary-receptor RMSD/RMSF near zero, initial-contact retention `[1, 1, 0]`, and one H-bond event over three sampled frames. COM distance must increase and BSA must decrease without becoming materially negative. Generated PDB/XTC files and analysis CSVs remain in the chosen output directory.

## Docking and mutation arithmetic

The generator uses only the Python standard library:

```sh
python -m examples.synthetic_data --output runs/examples/scores
python -m src.modeling.docking.summarize --dataset runs/examples/scores/docking --output-dir runs/examples/docking-summary
python -m src.modeling.mutant_screen.summarize --study runs/examples/scores/mutants --output runs/examples/mutant-summary
```

The paired docking example uses invented scores −6 and −7 kcal/mol, so the ASCR18-minus-ASCR3 contrast is −1. The ligand names follow the study schema; the coordinates describe synthetic ethanol.

The mutation summarizer requires all 19 alternatives at each of four sites. The generator constructs that schema in temporary output files, plus WT controls, a repeated design in another batch and one double mutant. Invented singles have `(S3, S18) = (−7, −6)`, WT `(−6, −6)`, and the double `(−8, −5)`. Thus single-mutant `D = 1`, double-mutant `D = 3`, and its `(epsilon3, epsilon18, epsilonD) = (0, 1, 1)`.

With RDKit installed, also run:

```sh
python -m src.modeling.mutant_screen.geometry --study runs/examples/scores/mutants --output runs/examples/mutant-geometry
```

The artificial ethanol coordinates of the double are translated 2 Å from each parent, so both no-fit RMSD and centroid shift are 2 Å in all four comparisons.

## ODE-generated paired readings

```sh
python -m src.modeling.dual_channel.generate_inputs --input-dir studies/signal-model/dual-channel/inputs --source studies/signal-model/pathway/configs/literature_core_source.json --output runs/examples/dual-inputs
python -m src.modeling.dual_channel.replay --input-dir runs/examples/dual-inputs --output runs/examples/dual-replay --eco-config studies/economic-model/configs/eco_model_example_model.json --protocol-id SIM_V4_MI1TO10_CE10TO1_Q1_8H_A480
```

The generator solves the supplied ODE, producing 18 channel endpoints, nine calibration pairs and 36 evaluation pairs. All are marked `synthetic_ODE`; calibration and evaluation use distinct sample identities and target/background combinations. Calibration and evaluation share the same ODE generator.

For a standalone pathway run, follow the [pathway guide](../studies/signal-model/pathway/README.md). The [economic guide](../studies/economic-model/README.md) provides a standard-library-only scenario.

## Automated checks

[Tests](../tests/README.md) create these inputs during execution and verify independent arithmetic/geometry expectations, conservation, regression metrics, malformed-input rejection and output protection.
