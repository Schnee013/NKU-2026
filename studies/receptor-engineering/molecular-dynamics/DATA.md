# MD inputs and generated output fields

This dictionary defines configuration fields, output filenames, identity keys and units used by the analysis and plotting programs. Use these schemas to organize your inputs and generated results. See [reproducibility](../../../REPRODUCIBILITY.md) for workflow requirements.

In [inputs.json](inputs.json), `receptor_prediction.method` names the prediction method, `execution` specifies the execution mode, and `environment_path` links the deployment configuration relative to this study directory. See [AF3 workflow guide](../../../REPRODUCIBILITY.md#alphafold-3-predictions) for output files and downstream structure mapping.

`inputs.json` records dataset IDs, genotype, ligand, replicate number, `production_ns` duration, receptor / ligand atom counts and preparation differences. `protein_atoms` counts all receptor atoms, including hydrogens; `protein_charge_e` is the receptor net charge in elementary-charge units. `processed_atoms` counts receptor plus ligand atoms, while `full_atoms` also includes membrane, water and ions. `sequences.fasta` gives residues 1–406 from each system's first representative structure. The two ligand arms retain their separate identities when receptor sequences coincide. [Simulation configurations](configs/simulation/systems.json) map system composition to parameter files; [Methods](METHODS.md) define the preparation schedule, selections and numerical thresholds.

| Data | Primary identity | Fields / units |
|---|---|---|
| Base per-system CSVs | Dataset directory, `rep`, `time_ns`; RMSF adds residue / segment | `backbone_rmsd_A`, `rg_A`, `selection_rmsd_A`, `distance_A`, `rmsf_A` are Å; `bsa_A2` is Å²; fractions / occupancy dimensionless |
| Base H-bond events | `rep`, `time_ns`, direction, donor–H–acceptor pair | Distance Å and angle degrees; `events` counts observations, `n_frames` is the sampled-frame denominator, and `occupancy = events / n_frames` |
| Base summaries | `case_id`, metric / cluster, time block where present | Means, maxima, slopes per ns, initial contact counts, complete eight-system comparison; `unique_frames` counts deduplicated production timepoints; `max_heavy_bond_A` is the maximum heavy-atom covalent bond length in Å; `max_relative_ligand_COM_step_A` is the maximum change in ligand COM relative to the receptor between consecutive 10 ps frames, in Å |
| Local `frame_metrics.csv` | `case`, `time_ns` | Global / core / pocket fits, ligand RMSD and displacement, lid RMSD, distances Å; sidechain chi1 degrees |
| Local `anchor_distances.csv` | Time, residue | Minimum distance and nearest atom names; Å |
| Local `atom_contact_events.csv` | Time, residue, protein atom, ligand atom | Atom-pair distances and chemical ligand group |
| Local `hbond_events.csv` | Time, residue, direction, donor, H, acceptor | Geometric events; distance Å, angle degrees; union occupancy counts each timepoint once |
| Local `site_neighborhoods.csv` | Site, neighbor, cutoff | Sequence gap, distance and occupancy; cutoffs Å |
| Local rotamer summary | Case, site, circular bin center | Initial state and early / late occupancy; centers degrees |
| Local spatial adjacency | Site, neighbor, anchor, time where present | Both edge distances and anchor–ligand distance Å; simultaneous occupancy dimensionless |
| Local `pair_details/` | Site, neighbor, time | Actual closest atom names, distances, Arg368 charged-end geometry and chi2/3/4 degrees |
| Local summaries | Case plus metric, residue / threshold or pair | Full comparison collection; late means / minima / maxima and descriptive first-event time ns |
| Representative PDBs | Dataset, MODEL 1–6 | Full 406-residue receptor and ligand coordinates at 0, 10, 25, 50, 75, 100 ns; Å units and 0.001 Å output precision; equilibrated production-start and later configurations for inspecting ligand position, lid geometry and contacts |
| Fixed selections | `core_resids`, `pocket_resids`, `lid_resids` | One-based receptor residues and the selection definition; per-system derivation observations are external |
| Simulation settings | Common MDP profiles and `configs/simulation/systems.json` | Integration time ps; cutoff nm; temperature K; pressure bar; molecule counts and per-system atom counts |

Aggregated local summary tables use `case` plus the keys below. New trajectory analysis writes the corresponding per-case files under its chosen output directory.

| Aggregated table | Additional key fields | New per-case output |
|---|---|---|
| `results/local-summary/anchor_sensitivity.csv` | `resid`, `cutoff_A` | `anchor_summary.csv` |
| `results/local-summary/local_pair_summary.csv` | `site`, `neighbor` | `pair_details/pair_summary.csv` |
| `results/local-summary/proximity_bridges.csv` | `site`, `neighbor`, `anchor` | `proximity_bridges.csv` |
| `results/local-summary/sidechain_rotamers.csv` | `site`, `chi1_bin_center_deg` | `sidechain_rotamers.csv` |

A complete external collection places initial-contact occupancies and CA RMSF at `results/base/<dataset>/contacts/initial_contact_occupancy.csv` and `results/base/<dataset>/rmsf/rmsf.csv`. The directory's dataset ID links each table to [inputs.json](inputs.json).

`early_mean` in local output uses t<50 ns; `late_mean` includes both 50 and 100 ns. Occupancies use the sampling grid and eligible atoms in [Methods](METHODS.md#basic-metrics). Empty event tables indicate that the defined event was absent. Empty Arg368–Asp70 tables for WT/N84V indicate that Arg368-specific geometry is inapplicable to those genotypes.

Base and local geometric H-bond definitions are given in [Methods](METHODS.md). Tables retain their numeric precision. Plotting programs produce PNG and, where supported, SVG.

`plot_structural_summary` reads seven files from `results/local-summary/`: `local_pair_summary.csv`, `anchor_hbond_union.csv`, `anchor_group_occupancy.csv`, `Arg368_Asp70_summary.csv`, `proximity_bridges.csv`, `case_summary.csv` and `event_timeline.csv`. Supply separate WT–ascr3 and V368R–ascr3 rows for the Arg345 carboxylate contact comparison. An empty sustained-event onset is displayed as a dash; zero occupancy is displayed as zero. The command writes `mechanism_evidence_map.png` and `.svg`.

## Polar contacts and transition geometry

| Output layout | Identity, values and use |
|---|---|
| `results/polar/<case>/polar/` | Five cases; 10 ps `polar_timeseries.csv.gz`, per-frame nearest and best qualifying pair geometry, atom definitions/mapping, uninterrupted episodes, per-case windows and anchor competition. `structures.json` identifies actual initial, modal late-state and final PDB frames. |
| `results/polar-summary/` | `case_id`, feature, metric and window identify occupancy, present-frame count, denominator and contiguous spans. `primary_polar_summary.csv` and `all_windows_polar_summary.csv` retain controls and chemically absent entries. |
| `results/polar-screen/<case>/` | Four cases; `features.csv` identifies residue, group, direction and applicability. `polar_features.npz` stores 10,001 `times_ns`, feature names, boolean primary/tight/relaxed matrices and float32 minimum salt-bridge distances. Hydrogen-bond event CSVs retain atom-pair geometry at source precision. |
| `results/polar-screen/comparison/` | Complete `single_feature_summary.csv` and `all_frequent_patterns.csv`, including nonqualifying patterns. Keys are case/feature/mode/window and target/mode/window/pattern, respectively. |
| `results/polar-screen-summary/` | Qualifying patterns, observed primary features, full screen coverage, Asp89 window comparisons and 1 ns timeline values. Missing chemical groups have missing occupancy; a measured absence has zero occupancy. |
| `results/polar-transition/<case>/` | Four-case timelines at 10 ps, phase and neighbor statistics; `geometry.npz` retains selected-atom positions under core/local fits, atom identities and 100 ps all-residue neighbor distances. Selected coordinates are float32; time and neighbor distances retain float64. These derived arrays contain local analysis atoms, not the full raw trajectory. |
| `results/polar-transition-summary/` | Full phase summaries, motion decomposition, proximity checks and mid-state interruptions; candidate scores, selected models, chemical atom mapping and conformer-pool counts. Energies/strain are kcal/mol, distances Å, angles degrees and time ns. |
| `results/polar-transition/structures/` | Actual phase/compatibility frames indexed in `representatives.json`; `*_constructed.pdb` and SDF are static ligand models in fixed actual pockets. `geometry_pass` and accepted/refined candidate counts must accompany their interpretation. |
| `structures/polar-screen/` | Actual all-residue-screen representative frames and a later-state frame; filenames identify the case and time. |

For the coordinate checks in `verify_tables`, provide `structures/<dataset>_six_frames.pdb` at 0/10/25/50/75/100 ns with full receptor and ligand atom order. PDB coordinate output is rounded to 0.001 Å; numerical arrays can retain higher precision. All eight configured systems use the same fixed fit residue sets.
