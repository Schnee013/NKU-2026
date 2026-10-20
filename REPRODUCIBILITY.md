# Reproduction guide

Run commands from the repository root. Choose a workflow below, install its dependencies, prepare the listed inputs and use a new output directory for each run.

## Set up an environment

The [environment guide](environments/README.md) lists dependency files and installation commands for each workflow. For MD analysis and its synthetic examples:

```sh
python3.12 -m venv /tmp/nku-analysis-env
. /tmp/nku-analysis-env/bin/activate
python -m pip install -r environments/md-compatibility-requirements.txt
```

Use separate environments for dependency sets with different pinned versions, including the paired-channel and docking workflows.

## Run the examples

The [example guide](examples/README.md) provides commands and output descriptions for synthetic docking, mutation geometry, MD analysis and paired-channel readings. For the MD example:

```sh
python -m examples.md_geometry --output runs/examples/md-geometry
```

This creates artificial PDB/XTC inputs, runs seven analysis drivers and writes CSV outputs and `checks.json` under the chosen directory. For automated checks, use the commands in the [test guide](tests/README.md).

## Prepare workflow inputs

| Workflow | Inputs and configuration | Run instructions and outputs |
|---|---|---|
| Docking summaries | Completed-run inventory, per-pose scores and matching model/ligand identifiers | [Receptor comparison](studies/receptor-engineering/receptor-comparison/README.md): paired score tables |
| Receptor and ligand preparation | Receptor structures, ligand SDFs, atom maps and search-box settings | [Preparation methods](studies/receptor-engineering/receptor-comparison/METHODS.md): prepared structures and contact-table schemas |
| Mutation comparisons | Construct definitions, completed-run inventory, selected pose coordinates and contact tables | [Mutation tools](studies/receptor-engineering/mutant-screen/README.md#your-own-data): score contrasts and pose-geometry tables |
| MD analysis | Processed trajectories, matching topology, atom order, selections and analysis settings; bonded topology for hydrogen-bond analysis | [MD guide](studies/receptor-engineering/molecular-dynamics/README.md): time series, summaries and plots |
| Gpa1 sequence design and confidence parsing | Supplied source sequences and graft definitions; AF3 output directories for confidence parsing | [Gpa1 guide](studies/gprotein-coupling/README.md): FASTA designs and confidence tables |
| Pathway and paired-channel models | Parameter JSONs, protocol settings and paired readings | [Pathway](studies/signal-model/pathway/README.md) and [paired-channel](studies/signal-model/dual-channel/README.md) guides: time courses, predictions and figures |
| Economic scenarios | Density inputs, cost assumptions and model configuration | [Economic guide](studies/economic-model/README.md): scenario reports and comparison tables |

Update configuration paths and expected SHA-256 values to match your files. Keep atom order, residue numbering, construct IDs and run IDs consistent across inputs.

## Generate MD summaries and figures

Organize analysis outputs using the filenames and fields in the [MD data dictionary](studies/receptor-engineering/molecular-dynamics/DATA.md). Set each program's `--study`, `--data-root`, `--summary-dir` or `--table` argument to the corresponding input directory or file.

For a complete study containing all eight configured systems, compare the per-system tables with combined summaries:

```sh
python -m src.modeling.md_analysis.verify_tables --study /absolute/path/to/complete-md-study --output runs/md/verified-tables
```

To draw the structural summary from `results/local-summary/`:

```sh
python -m src.modeling.md_analysis.plot_structural_summary --study /absolute/path/to/complete-md-study --language en --output runs/md/structural-summary
```

The plotting command writes PNG and SVG files. For Chinese labels, use `--language zh` and provide a Chinese-capable font with `--font /absolute/path/to/font.otf`.

## AlphaFold 3 predictions

1. Set up the [AF3 environment](environments/alphafold3.md), model parameters and databases. Select the request JSON for the desired construct from the [mutation](studies/receptor-engineering/mutant-screen/README.md#your-own-data) or [Gpa1](studies/gprotein-coupling/README.md#sequence-design-and-job-correspondence) workflow.
2. Set the protein sequences, chain IDs and integer seeds in the request, then run prediction with the configured AF3 installation.
3. Keep each job's [native output directory](environments/alphafold3.md#native-prediction-files) together: mmCIF structures, confidence files, processed input JSON, ranking CSV and accompanying terms.
4. Select structures by their `(job, seed, sample)` identities. If using numbered exports, map each model number to its selected source file and use that identity in downstream construct and run tables.

Convert a selected structure for PDB-based preparation with the docking environment:

```sh
python -m src.modeling.docking.cif_to_pdb --input-file /absolute/path/to/model.cif --output-file runs/structures/model.pdb
```

Keep the source mmCIF and accompanying terms alongside the converted structure.
