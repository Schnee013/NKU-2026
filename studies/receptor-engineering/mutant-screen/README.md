# Mutation-screen tools

Programs for checking construct inputs, comparing paired ligand scores and contacts, and measuring parent/combination ligand-pose displacement. Sequence requests, construct definitions and run parameters specify the preparation workflow. Use the synthetic example below to try the analysis, or supply your own scores, poses and contact tables. Scientific findings are on the [docking Wiki](https://2026.igem.wiki/nku-china/model-docking/).

## Methods and interpretation

For ligand scores S3 and S18 in kcal/mol and a declared WT baseline:

```text
deltaS3 = S3(variant) - S3(WT)
deltaS18 = S18(variant) - S18(WT)
D = deltaS18 - deltaS3
```

A positive D is a relative score shift toward ascr#3 under the declared protocol. Inspect both ligand scores: a contrast can reflect worsening of the other ligand. These values are not measured binding affinities.

Each batch must identify its own baseline. Batch 08 parent comparisons use batch-05 N84/V110/V368 parents and batch-07 W338 parents, with explicit WT references for sensitivity comparisons. For each ligand, `epsilon = S(AB)-S(A)-S(B)+S(WT)`. Parent and combination structures are separate predictions; these residuals do not establish experimental epistasis.

The [construct definition](inputs/receptors.json) deletes residues 1–21 and 164–178 while keeping original chain-A numbering. Fixed-backbone alanine substitution and full-prediction routes are distinct. Full-prediction preparation fits C-alpha residues 30–50, 62–83, 101–126, 144–163, 182–201, 326–346 and 358–380 to a supplied common reference. Native prediction requests specify sequences and integer seeds. Numbered model labels identify exports; the [AF3 workflow guide](../../../REPRODUCIBILITY.md#alphafold-3-predictions) describes their association with source predictions.

[Run parameters](inputs/runs.json) define search settings and the contact-selection rule. To summarize completed runs, provide a result inventory with observed pose counts, selected modes and contact-file paths, as described under [Your own data](#your-own-data).

Vina settings include center (20, −20, 5) Å, size (26, 26, 26) Å, exhaustiveness 8, maximum nine modes, energy range 3 kcal/mol and seed 1919. Each run specifies its CPU count. Contact selection uses the minimum score, breaking ties by original mode number. Contacts match original residue number plus interaction type; Jaccard uses entry intersection/union. Other poses from one search are not independent repeats.

## Synthetic comparison

```sh
python -m examples.synthetic_data --output runs/examples/scores
python -m src.modeling.mutant_screen.summarize --study runs/examples/scores/mutants --output runs/examples/mutant-summary
python -m src.modeling.mutant_screen.geometry --study runs/examples/scores/mutants --output runs/examples/mutant-geometry
```

The first two commands use the standard library; geometry requires RDKit. [Examples](../../../examples/README.md#docking-and-mutation-arithmetic) define the invented scores and exact 2 Å coordinate translation. The generated study deliberately has the complete saturation schema required by this summarizer. [Tests](../../../tests/README.md) check score arithmetic, parent links, cross-batch identity, contact denominators and coordinate calculations.

## Your own data

`--study` must point to a complete external result collection for result summaries or pose geometry. `inputs/receptors.json` identifies constructs, baselines and comparison sets. `inputs/runs.json` must additionally contain actual `pose_count`, `selected_mode`, `selected_pose.tied_minimum_modes` and `contacts_path` for every ligand run. `reference/pose-scores.csv` has one row per run/mode; matching contact folders contain `contact_residues.csv` and `ifp_presence.csv`. The synthetic generator supplies a complete schema example.

Summary outputs are `summary.json` and `comparisons.csv`. Geometry reads chemically consistent `top1.sdf` files and reports heavy-atom RMSD over equivalent stereochemical atom mappings and centroid displacement, without fitting away ligand motion. All outputs require a new directory outside the input study.

For new docking, install the [docking environment](../../../environments/README.md#docking), copy an `inputs/representative-*.conf` template and supply receptor/ligand paths in a common frame:

```sh
python -m src.modeling.docking.run_vina --config /absolute/path/to/your/vina.conf --output-dir runs/mutants/docking
python -m src.modeling.docking.analyze_contacts --receptor-pdb /absolute/path/to/your/input.pdb --receptor-pdbqt /absolute/path/to/your/receptor.pdbqt --poses-pdbqt runs/mutants/docking/poses.pdbqt --profile mutant-top1 --output-dir runs/mutants/contacts
```

The preparation helper checks a supplied construct PDB and matching PDBQT against its sequence request:

```sh
python -m src.modeling.mutant_screen.prepare --study studies/receptor-engineering/mutant-screen --record-id 05/run-01/N84V/af0/seed1919 --route saved-construct --construct-pdb /absolute/path/to/your/input.pdb --construct-pdbqt /absolute/path/to/your/receptor.pdbqt --output runs/mutants/preparation
```

`--route full-prediction` also needs `--reference-pdb` and, for prediction-derived constructs, `--cif`. The helper checks residue inventory, sequence, coordinates, charge and AutoDock types. [Reproducibility](../../../REPRODUCIBILITY.md) describes input preparation, software setup and execution steps. [Attributions](../../../ATTRIBUTIONS.md) and [AF3 terms](../../../third-party/AlphaFold3OutputTerms.md) describe external software and materials.
