# Methods and data interfaces

This guide defines receptor preparation, score and contact calculations, and the input/output schemas for analyzing docking runs. Scientific results are described on the [docking Wiki](https://2026.igem.wiki/nku-china/model-docking/).

## Input identity and constructs

[Constructs](inputs/constructs.csv) and [source sequences](inputs/receptors/fasta/) define receptor identities. `model_index` is a numbered prediction-export identity; `model_rel` identifies a prepared docking construct. Use the [AF3 workflow guide](../../../REPRODUCIBILITY.md#alphafold-3-predictions) to prepare source predictions and map selected structures to these construct identities.

Molecular coordinates are in Å. Gpr2 full-length spans residues 1–406; its open-pocket construct retains 22–163 and 179–406 (370 residues), deleting N-terminal residues 1–21 and ECL2 lid 164–178. Gpr3 full-length spans 1–401; its docking construct retains 8–149 and 164–401 (380 residues), deleting 1–7 and 150–163. Original residue numbering is retained. Separately predicted delta-N inputs and trimmed full-length predictions have distinct identities.

The lid deletion is a preparation choice for rigid-receptor search. Full-length MD preparation is described in the [MD methods](../molecular-dynamics/METHODS.md). New calculations require receptor coordinates matching the intended construct, numbering and search-box frame.

Five family boxes in [boxes.tsv](configs/boxes.tsv) are shared by the five aligned models and, where present, full-length/truncated forms of each family. `align_pdb_to_reference.py` fits explicitly matched C-alpha residue ranges by a rigid transformation; `cif_to_pdb.py` converts user-supplied mmCIF coordinates through Biopython.

Ligand inputs include the PubChem structures (ascr#3 CID 16066476; ascr#18 CID 71514784) and selected deprotonated microstates, each with formal charge −1. Meeko 0.7.1 prepares Gasteiger partial charges. Explicit atom/bond maps and `check_sdf_pdbqt_mapping.py` check chemical identity, atom correspondence and coordinate conservation. The serialized Meeko map supplies ligand bond orders during pose restoration.

## Scores and pose counts

The runner uses official Vina 1.2.7 Python bindings. Defaults are exhaustiveness 8, seed 1919, at most nine modes, energy range 3 kcal/mol and grid spacing 0.375 Å. Record the actual version, CPU count and box with each new run. Versions and platforms can affect poses even with a fixed seed.

Affinity is the Vina score in kcal/mol; lower scores are more favorable under that scoring model. Each receptor model is paired with itself across ligands. `ascr18_minus_ascr3_kcal_mol` is the lowest logged ASCR18 score minus the lowest logged ASCR3 score. The full `model_rel` identity is the join key.

The number of logged search scores can differ from the number of stored coordinate models. Mark absent coordinate records with `pose_present=false`. Separate `pose_*` columns preserve PDBQT values at their own precision; stdout may use coarser precision. Missing coordinate values remain distinct from zero.

## Contacts

The cross-receptor profile uses ProLIF 2.1.0, `count=True`, `vicinity_cutoff=6 Å` and that version's default interaction set. Meeko per-residue chemical templates supply explicit hydrogens; residues are combined without inferring bonds from spatial proximity. Each feature is a protein residue–interaction pair. Geometry counts can exceed one per pose; presence is `count > 0`. Frequency is the number of stored poses with the feature divided by the actual number of stored poses for that run.

The mutant-screen profile uses only the selected minimum-scoring pose, nine explicit interaction classes and presence semantics (`count=False`). `analyze_contacts.py --profile mutant-top1` implements this definition. Single-pose contacts, pose frequencies, recurrence across models and MD time occupancy have different denominators.

## Residue comparisons

The original comparison design uses truncated Gpr2/Gpr3 and the other three WT constructs, with five models and two ligands per construct; full-length Gpr2/Gpr3 are additional comparisons. A recurrence analysis must preserve those construct/model identities and the actual available-model denominator.

Recognition recurrence counts models containing the exact residue–interaction feature. A median over available models includes zero for an absent interaction. Unavailable models remain missing. Contact-only VdW, hydrophobic and directional/charge-related interactions must retain their separate meanings.

The descriptive contact-frequency rule requires preferred-ligand recurrence, a greater preferred-ligand frequency in at least three of five paired models and a median paired difference of at least 0.20.

Sequence alignments, topology predictions, profile coverage and structural mappings are external evidence for cross-receptor residue correspondence. Similar contacts or a common low-identity fold-template position alone do not establish homology. Distinguish alignment gaps, construct-deleted positions and unavailable evidence.

## Data fields

The following paths describe a complete external result collection. The [synthetic generator](../../../examples/synthetic_data.py) supplies a small invented instance for software checks.

| File | Identity and fields |
|---|---|
| `inputs/constructs.csv` | `model_rel`; receptor, construct, model index, residue count/ranges and original numbering |
| `reference/runs.csv` | Model–ligand `result_id`; executable version, seed, CPU, grid spacing/center/size, exhaustiveness, logged/stored pose counts, maximum modes and energy range |
| `reference/pose-scores.csv` | `(result_id, pose)` with 1-based pose numbers; stdout affinity/RMSDs, separate PDBQT values and `pose_present` |
| `reference/vina/<model_rel>/<ligand>/` | Stored pose MODEL blocks and their chemical maps |
| `reference/prolif/<model_rel>/<ligand>/` | Matching SDF poses and per-pose contact fingerprints |
| Generated paired scores | Lowest ASCR3/ASCR18 scores and signed ASCR18-minus-ASCR3 difference for each model |

Contact `ifp_counts.csv` and `ifp_presence.csv` use `ligand`, `protein`, `interaction` header levels and zero-based `Pose` rows. `interaction_frequency.csv` uses a `(protein, interaction)` key, poses with contact, pose frequency and total geometries. `metadata.json` records the actual pose denominator, software versions, residue count, charge and explicit hydrogen counts.

To join a residue comparison back to a pose, resolve `(receptor, construct, model_index)` to `model_rel`, then resolve its ligand to `result_id = model_rel/ligand`. Match the exact residue and interaction in the corresponding frequency table. Zero means an observed absence in available poses; missing input is not zero.

## Scientific sources

The following studies describe receptor function and provide the biological context for receptor selection.

- DAF-37/GprC and fungal receptor function: [10.1038/s41564-024-01731-9](https://doi.org/10.1038/s41564-024-01731-9).
- Gpr2/Gpr3 receptor responses: [10.1038/s41564-024-01679-w](https://doi.org/10.1038/s41564-024-01679-w).
- SRBC-64/66 function: [10.1126/science.1176331](https://doi.org/10.1126/science.1176331).
