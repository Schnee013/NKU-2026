# Cross-receptor docking tools

These tools prepare user-supplied receptors, run Vina, analyze pose contacts and summarize paired ligand scores. Source sequences, construct definitions, ligand inputs and search-box configurations support new calculations; the synthetic example demonstrates summary arithmetic. See the [docking Wiki](https://2026.igem.wiki/nku-china/model-docking/) for scientific results.

## Synthetic example

```sh
python -m examples.synthetic_data --output runs/examples/scores
python -m src.modeling.docking.summarize --dataset runs/examples/scores/docking --output-dir runs/examples/docking-summary
```

The invented score pair gives an ASCR18-minus-ASCR3 contrast of −1 kcal/mol. [Tests](../../../tests/README.md) check score calculations, coordinate/contact matching, invalid inputs and output protection.

## New calculations

Install the [docking environment](../../../environments/README.md#docking). Copy [the Vina configuration](configs/wt-gpr2-ascr18.txt) outside the checkout's study inputs, supply absolute receptor/ligand paths and set the search box in the receptor's coordinate frame. Prepared PDBQT ligands are available under [ligand inputs](inputs/prepared/ligands/deprotonated-minus1/pdbqt/).

```sh
python -m src.modeling.docking.run_vina --config /absolute/path/to/your/vina.conf --output-dir runs/docking/new-search
python -m src.modeling.docking.analyze_contacts --receptor-pdb /absolute/path/to/your/receptor.pdb --receptor-pdbqt /absolute/path/to/your/receptor.pdbqt --poses-pdbqt runs/docking/new-search/poses.pdbqt --output-dir runs/docking/new-contacts
```

The contact runner prepares a Meeko receptor object from the supplied PDB, or reads one with `--receptor-json`, and checks its PDBQT export against the docking receptor. Ligand restoration checks the serialized chemical map, charge, atom identities and coordinates. `--profile mutant-top1` selects the mutation-screen definition.

The matrix runner requires a user-prepared dataset with [construct definitions](inputs/constructs.csv), [search boxes](configs/boxes.tsv), receptors at `inputs/prepared/receptors/af3/pdbqt/<model_rel>.pdbqt` and ligands at `inputs/prepared/ligands/deprotonated-minus1/pdbqt/<ligand>.pdbqt`. The dataset must contain 35 prepared receptor models and two ligands:

```sh
python -m src.modeling.docking.run_matrix --dataset /absolute/path/to/your/dataset --output-dir runs/docking/matrix --dry-run
```

Remove `--dry-run` only when intending to execute all searches. Defaults are exhaustiveness 8, seed 1919, at most nine modes, energy range 3 kcal/mol, spacing 0.375 Å and 16 CPUs (`--cpu`). Outputs must be new or empty and outside `studies/`.

## Summaries and interfaces

For your own result collection, `summarize --dataset ...` expects `reference/runs.csv` and `reference/pose-scores.csv`. The [synthetic generator](../../../examples/synthetic_data.py) demonstrates the schema. Paired scores use the lowest score per ligand and preserve the signed ASCR18-minus-ASCR3 difference. Missing values are not measured zeroes.

The optional `verify_scientific_bodies()` function checks supplied pose MODEL blocks, scores, residue identities, fingerprint counts/presence and their sampled-pose denominators. It requires the corresponding external coordinate and contact files. `merge_prolif_csv` accepts an explicit run table and batch root:

```sh
python -m src.modeling.docking.merge_prolif_csv --runs /absolute/path/to/runs.csv --batch-root /absolute/path/to/prolif --output-dir runs/docking/merged --layout both --kinds frequency presence counts
```

[Methods](METHODS.md) define preparation, scores, contact denominators and interpretation. [Reproducibility](../../../REPRODUCIBILITY.md) describes workflow inputs, execution steps and outputs. See [attributions](../../../ATTRIBUTIONS.md) for software and scientific sources.
