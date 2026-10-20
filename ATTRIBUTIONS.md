# Scientific references, software and licenses

## Scientific references

Study-specific references, sequence and chemical accessions, and sources for numerical assumptions are listed in the [study documentation](README.md).

AlphaFold 3: Abramson et al., *Nature* (2024), DOI [10.1038/s41586-024-07487-w](https://doi.org/10.1038/s41586-024-07487-w).

## Software sources

AlphaFold 3 source code is provided by [Google DeepMind](https://github.com/google-deepmind/alphafold3/tree/85c4d20505fd5cef05eac22b534d4e793971ae69). The analysis tools use separately installed dependencies including NumPy, SciPy, Matplotlib, pandas, MDAnalysis, MDTraj, RDKit, Biopython, Meeko, ProLIF and AutoDock Vina. Dependency sets and version specifications are listed in the [software environments](environments/README.md).

## Licensing

The [Apache License 2.0](LICENSE) covers the original source code, analysis scripts and tests that the NKU-China team has the right to license. Third-party materials, scientific data and separately installed dependencies retain their applicable licenses and attribution.

AlphaFold 3 model parameters are subject to the official [model parameter terms](https://github.com/google-deepmind/alphafold3/blob/85c4d20505fd5cef05eac22b534d4e793971ae69/WEIGHTS_TERMS_OF_USE.md). The included [AlphaFold 3 Output Terms](third-party/AlphaFold3OutputTerms.md) are reproduced from a native prediction output's `TERMS_OF_USE.md`. Each prediction output and its derivatives retain the terms supplied with that output.
