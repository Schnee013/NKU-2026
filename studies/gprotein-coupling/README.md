# Preliminary design exploration of chimeric Gpa1

This module designs terminal-segment Gpa1 chimeras and parses prediction confidence fields. Ten chimeras replace the terminal segment of the 472-residue *Saccharomyces cerevisiae* Gpa1 backbone with residues from *Arthrobotrys oligospora* Gpa2 or *Arthrobotrys flagrans* GasA, alongside an unmodified Gpa1 WT control.

Use the source sequences for graft generation and supply AF3 outputs to the confidence parser. See the [modeling Wiki](https://2026.igem.wiki/nku-china/model/) for project context and [reproducibility](../../REPRODUCIBILITY.md) for input requirements.

## Materials

| Material | Complete record |
|---|---|
| Backbone, donors and receptor sequences | [Source FASTA records](sequences/) |
| Ten graft designs and WT control | [11 individual FASTA files including WT](chimeras/), [combined FASTA of the ten chimeras](chimeras/all_chimeras.fasta) |
| 42 prediction input specifications | [Job/score correspondence](prediction-jobs.csv), [native-format input specifications](prediction-inputs/) |

## Sequence design and job correspondence

For graft length `n`, the design is `Gpa1[0:472-n] + donor[-n:]`, with the one-based retained backbone ending at `472-n`. Each sequence has 472 residues. The graft lengths are 5, 10, 15, 25 and 37 for each of the AoGpa2 and AfGasA donors, plus unmodified Gpa1 WT. Sequence generation uses the source FASTA records directly and does not require a computed alignment.

[Batch 1](prediction-inputs/batch1/) contains 26 single-job JSON files and [batch 2](prediction-inputs/batch2/) contains 16. These input specifications express the recorded sequences and seeds in the native `alphafold3` dialect, schema version 4, with integer seeds and explicit protein chain IDs. Omitted MSA/template fields enable database searches by the [local AF3 deployment](../../environments/alphafold3.md); the resulting processed `<job>_data.json` is a separate file containing the actual MSA and template data. The `input_path` column in the job table identifies each specification. The request plan comprises 11 monomers, 30 fungal receptor–chimera complexes and one Ste2–WT control. The input name `Gpa1-WT_Ste2_complex` maps to the parser's `Ste2_Gpa1-WT_complex` score ID. Both names appear in [prediction-jobs.csv](prediction-jobs.csv) and identify the same sequence pair.

The `cognate` group comprises Gpr2/AoGpa2, Gpr3/AoGpa2 and GprC/AfGasA designs; the other fungal receptor/donor combinations are `cross`. Ste2–WT is `control`. These labels record the intended pairings.

## Confidence fields

For new calculations, the [confidence parser](../../src/modeling/gprotein/parse_af3_results.py) supports two input layouts. With native outputs it selects the largest full-precision score from `ranking_scores.csv`, retaining the first CSV row on an exact tie, and writes `best_seed` and `best_sample` with `best_model` empty. With numbered exports it writes `best_model` and leaves the native identifiers empty. The confidence fields are:

- `plddt_global`: mean Cα B-factor from the selected monomer CIF, rounded to 2 decimals on the pLDDT scale. `plddt_cterminal` is the mean of the last `n` Cα values for a graft, and `N/A` for WT.
- `ptm`, `iptm` and `ranking_score`: selected model-summary confidence values, rounded to 4 decimals and dimensionless. `fraction_disordered` is the supplied fraction rounded to 3 decimals; `has_clash` is the supplied flag.
- `ipae_mean`: mean of all **directed token pairs from different chains** in the selected full PAE matrix, rounded to 2 decimals in Å. The mean includes non-contacting cross-chain pairs.
- `chain_pair_pae_min`: a separate AF3 summary field read from `[0][1]`, rounded to 2 decimals in Å.
- `score_range` and `iptm_range`: minimum/maximum supplied model-summary values over the predictions. They describe variation in model confidence. Native model selection uses the ranking CSV; these reported confidence fields retain the summary JSON values and their precision.

## Use

Use Python 3.12 and the standard library from the repository root, with a new or empty output directory outside the study:

```sh
python -m src.modeling.gprotein.design --output /tmp/gpa1-designs
```

The [sequence generator](../../src/modeling/gprotein/design.py) writes 11 FASTA files and `designs.csv`, with `design_id`, `donor`, `graft_length`, `backbone_end_1based` and `length`. See [test commands](../../tests/README.md) for sequence and prediction-input comparisons.

To parse predictions, supply the local AF3 output directories for all 42 jobs as `--downloads`. Native jobs are matched by the `name` in `<job>_data.json`, including the original `Gpa1-WT_Ste2_complex` name. Preserve `<job>_ranking_scores.csv` and each `seed-<seed>_sample-<sample>/` directory with its model CIF, confidence JSON and summary JSON. The [native output layout](../../environments/alphafold3.md#native-prediction-files) gives the filenames. Multiple seeds and different numbers of diffusion samples are supported.

Numbered exports remain readable: directory names use the job ID or score ID from `prediction-jobs.csv`, lowercased with hyphens replaced by underscores. Supply `*summary_confidences_<index>.json` for indices 0–4, the selected monomer's `*model_<index>.cif`, and the selected complex's `*full_data_<index>.json`. PAE dimensions follow `token_chain_ids`; monomer pLDDT is read from Cα records in the CIF. Incomplete or ambiguous inputs fail before score tables are written.

```sh
python -m src.modeling.gprotein.parse_af3_results --downloads /path/to/af3-outputs --output /tmp/gpa1-scores
```

The parser writes `monomer_scores.csv` and `complex_scores.csv`. Its output directory must be new or empty and outside both the study and the AF3 input directory.

## Design references

The design set draws on prior chimeric Gα studies:

- Hara et al. (2012), *Chimeric Yeast G-Protein alpha Subunit Harboring a 37-Residue C-Terminal Gustducin-Specific Sequence Is Functional in Saccharomyces cerevisiae*, DOI [10.1271/bbb.110820](https://doi.org/10.1271/bbb.110820).
- Ueda et al. (2003), *Functional Interaction between T2R Taste Receptors and G-Protein α Subunits Expressed in Taste Receptor Cells*, DOI [10.1523/JNEUROSCI.23-19-07376.2003](https://doi.org/10.1523/JNEUROSCI.23-19-07376.2003).

Sequence accession identifiers appear in the FASTA filenames and headers. See [attributions](../../ATTRIBUTIONS.md) for software sources and redistribution terms, and [AF3 workflow guide](../../REPRODUCIBILITY.md#alphafold-3-predictions) for prediction inputs and output handling.
