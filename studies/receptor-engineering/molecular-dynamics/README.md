# Molecular-dynamics analysis tools

These tools analyze receptor/partner motion, geometry, contacts, local fits and polar interactions. Use the source code, method definitions and eight-system configurations with your processed trajectories and matching topologies, or start with the synthetic example below. Scientific findings are presented on the [MD Wiki page](https://2026.igem.wiki/nku-china/model-md/).

| Material | Entry |
|---|---|
| Analysis dependencies and original computational environment | [Environment](ENVIRONMENT.md) |
| Preparation schedule, fitting, geometry, sampling and interpretation | [Methods](METHODS.md) |
| System identities and receptor sequences | [Inputs](inputs.json), [FASTA](sequences.fasta) |
| Integration and equilibration parameters | [Simulation settings](configs/simulation/systems.json) |
| Selections and local definitions | [Basic selections](configs/systems.json), [local parameters](configs/mechanism.json), [fixed selections](configs/selections.json) |
| Polar definitions | [Targeted anchors](configs/polar-anchor-definition.json), [all-residue screen](configs/polar-screen-definition.json) |
| Generated output fields and units | [Data dictionary](DATA.md) |
| Implementations | [Source code](../../../src/modeling/md_analysis/) |

## Synthetic geometry example

```sh
python -m examples.md_geometry --output runs/examples/md-geometry
```

This generates an artificial three-frame trajectory, runs all seven basic drivers and checks known numerical answers. Install the [MD dependencies](../../../environments/md-requirements.txt) before running it. See the [example description](../../../examples/README.md#md-geometry) and [tests](../../../tests/README.md). The toy geometry is not a physical MD simulation.

## Analyze your own trajectories

Run commands from the repository root in the [MD analysis environment](ENVIRONMENT.md). For each dataset, provide `processed/rep1/protein_ref.pdb` and `processed/rep1/prod_protein.xtc` under its input directory. The example project selections expect receptor and ligand coordinates at 10 ps intervals. Match selections, atom order, time spacing and analysis windows to your own data.

```sh
export MD_PAYLOAD_ROOT=/absolute/path/to/native_wt_ascr18
export MD_OUTPUT_ROOT="$PWD/runs/md/base"
export MD_SYSTEM=native_wt_ascr18
python -m src.modeling.md_analysis.analyze_convergence
python -m src.modeling.md_analysis.analyze_selection_rmsd
python -m src.modeling.md_analysis.analyze_com_distance
python -m src.modeling.md_analysis.analyze_rmsf
python -m src.modeling.md_analysis.analyze_bsa --stride 50 --n-sphere-points 240
python -m src.modeling.md_analysis.analyze_contacts --stride 50 --cutoff 4
python -m src.modeling.md_analysis.analyze_hbond --step 50
```

These seven commands use the processed PDB and XTC inputs above. Their atom selections and CSV outputs are:

| Driver | Atom selection / measurement | Output under `analysis/` |
|---|---|---|
| `analyze_convergence` | Receptor backbone RMSD to frame 0; receptor CA radius of gyration | `convergence/rmsd_summary.csv`, `convergence/rg_summary.csv` |
| `analyze_selection_rmsd` | Ligand heavy atoms after equal-weight receptor backbone fitting | `selection_rmsd/selection_rmsd.csv` |
| `analyze_com_distance` | Mass-weighted ligand heavy-atom COM to whole receptor and its own initial-contact residue cluster | `com_distance/com.csv` |
| `analyze_rmsf` | Receptor CA fluctuations over 50–100 ns on the processed coordinates | `rmsf/rmsf.csv` |
| `analyze_bsa` | Receptor and ligand SASA, individually and together | `bsa/bsa.csv` |
| `analyze_contacts` | Receptor/ligand residue heavy-atom distances; initial-contact retention and late occupancy | `contacts/contact_persistence.csv`, `contacts/initial_contact_occupancy.csv` |
| `analyze_hbond` | Candidate N/O/S donor–H and acceptor geometry between receptor and ligand | `hbond/hbond_events.csv`, `hbond/hbond_summary.csv` |

[Methods](METHODS.md#basic-metrics) give the geometric definitions, units and sampling denominators. Each driver writes its CSVs and figures under `MD_OUTPUT_ROOT/<system>/analysis/<metric>` and prints numerical summaries. Use a new output location for each calculation. `--system` overrides `MD_SYSTEM`; without either, the configuration selects WT–ascr18. `MD_CONFIG` selects another systems file (default `studies/receptor-engineering/molecular-dynamics/configs/systems.json`), and `MD_REPS` selects replicate numbers (default `1`). Relative input paths in configurations are resolved from the repository root. BSA/contact stride and H-bond step can be changed for other frame intervals; the published 500 ps sampling uses 50 frames from these 10 ps inputs.

For a common-frame comparison, pair each arm label with its configuration, dataset ID and input directory. This example compares WT–ascr18 and N84V–ascr18:

```sh
MD_OUTPUT_ROOT=runs/md/path-context python -m src.modeling.md_analysis.compare_ligand_paths \
  --arm-config WT=studies/receptor-engineering/molecular-dynamics/configs/systems.json \
  --arm-system WT=native_wt_ascr18 \
  --arm-payload WT=/absolute/path/to/native_wt_ascr18 \
  --arm-config N84V=studies/receptor-engineering/molecular-dynamics/configs/systems.json \
  --arm-system N84V=raw_n84v_ascr18 \
  --arm-payload N84V=/absolute/path/to/raw_n84v_ascr18 \
  --reference-arm WT --replicate 1 --fit-selection 'protein and backbone' \
  --contact-cutoff-A 4 --window-start-ns 50 --window-end-ns 100 --stride 50 \
  --output-dir runs/md/ligand-paths
```

The comparison writes `centroid_path.csv`, `residue_contact_occupancy.csv`, `representative_frames.csv`, figures and input/selection metadata to `--output-dir`. The actual sampling denominators are recorded per residue; [Methods](METHODS.md#basic-metrics) define the reference frame and early/late comparison. Fit, contact cutoff, time window and stride are explicit command-line options.

For local and atom-pair analyses, [mechanism.json](configs/mechanism.json) names each dataset's matching `solution.tpr`, `protein_ref.pdb` and `prod_protein.xtc`. Put these inputs at the configured paths or pass a configuration with their actual paths. The full-system `initial_full.tpr` and first-frame `initial_full.xtc` are needed only when deriving membrane-based selections with `--prepare`.

```sh
python -m src.modeling.md_analysis.analyze_local_mechanism --config studies/receptor-engineering/molecular-dynamics/configs/mechanism.json --output-root runs/md/local
python -m src.modeling.md_analysis.analyze_local_mechanism --config studies/receptor-engineering/molecular-dynamics/configs/mechanism.json --output-root runs/md/local --pair-detail
```

Local analysis uses the [fixed backbone fits and ligand heavy-atom selections](METHODS.md#local-fits-and-selection-definitions), with topology bonds defining ligand groups and candidate hydrogen bonds. It writes per-frame metrics, contact events, `anchor_summary.csv`, `proximity_bridges.csv`, `sidechain_rotamers.csv` and metric summaries under `--output-root/<dataset>`. `--pair-detail` measures the configured receptor heavy-atom pairs and Arg368 charged-end geometry, writing atom-pair time series and `pair_summary.csv` under that dataset's `pair_details/`. Output fields and units are listed in [DATA.md](DATA.md). `--case DATASET_ID` selects one dataset. To derive selections from the current initial structures, run with `--prepare` and an independent output root, then pass its `selections/selections.json` through `--selections` for local analysis. The default selections are [configs/selections.json](configs/selections.json). The published eight-system comparison uses the original six-system core and pocket selections unchanged; deriving selections from eight starting structures defines a different comparison.

## Polar analyses and external result plotting

The following interfaces analyze supplied processed trajectories. Copy the relevant configuration, set its input paths and expected SHA256 identities to the actual inputs, and choose new output directories outside the published study. Relative configuration paths resolve from the repository root. Each calculation checks its own inputs, atom correspondence and time grid.

```sh
python -m src.modeling.md_analysis.analyze_ligand_contact_patterns \
  --config studies/receptor-engineering/molecular-dynamics/configs/polar-anchors.json \
  --case v368r_ascr3_arg368_start --polar
python -m src.modeling.md_analysis.analyze_ligand_contact_patterns \
  --config studies/receptor-engineering/molecular-dynamics/configs/polar-screen.json \
  --case v368r_ascr3_arg368_start --all-polar
```

Repeat the chosen mode for each case in its configuration. Targeted anchors use five cases: both directed starts, the original V368R–ascr3, WT–ascr3 and V368R–ascr18. The all-residue screen uses four: both directed starts, WT–ascr3 and V368R–ascr18. Once all four scans exist under the configured case outputs, `--compare-all-polar` with the same screen configuration writes all single features, all frequent combinations, threshold results and coverage into `shared_output`.

`python -m src.modeling.md_analysis.polar_transition --config studies/receptor-engineering/molecular-dynamics/configs/polar-transition.json` reads each trajectory and its prior polar feature matrix, then writes timelines, phase statistics, selected-atom geometry and actual representative frames. `--case` restricts it to one configured case. Set the configuration's matrix paths to the matching polar-scan outputs.

`python -m src.modeling.md_analysis.pose_compatibility --config studies/receptor-engineering/molecular-dynamics/configs/pose-compatibility.json` requires external pocket frames and matching ligand topologies, together with the supplied chemically explicit SDF templates. It writes all candidate scores and selected constructed PDB/SDF models to a separate output root. This static calculation requires the [MD compatibility dependencies](../../../environments/md-compatibility-requirements.txt). [Methods](METHODS.md#polar-contacts-and-transition-geometry) define the thresholds, windows and interpretation.

Plotting existing results requires external result tables. Supply a complete study layout explicitly, for example:

```sh
python -m src.modeling.md_analysis.plot_md_comparison --study /absolute/path/to/complete-md-study --output runs/md/comparison-figures
python -m src.modeling.md_analysis.plot_structural_summary --study /absolute/path/to/complete-md-study --language en --output runs/md/summary-figures
```

Other plotting interfaces accept `--data-root`, `--summary-dir` or `--table` inputs. Set these arguments to your analysis outputs, using the file layout and field definitions in [DATA.md](DATA.md). See [reproducibility](../../../REPRODUCIBILITY.md) for the inputs needed to reproduce a scientific result.
