# AlphaFold 3 prediction environment

This file records the hardware, installed software and process settings of the local AlphaFold 3 deployment using official model parameters and sequence/template databases. Prediction input specifications and model-selection conventions are described in the [mutation screen](../studies/receptor-engineering/mutant-screen/README.md#your-own-data) and [Gpa1 study](../studies/gprotein-coupling/README.md#sequence-design-and-job-correspondence). The [AF3 workflow guide](../REPRODUCIBILITY.md#alphafold-3-predictions) describes prediction inputs and output handling. Docking preparation and molecular-dynamics environments are documented separately in the [environment index](README.md).

## Hardware and operating system

| Component | Configuration |
|---|---|
| GPU | 1 × NVIDIA A100-SXM4-80GB |
| GPU memory | 81,920 MiB |
| GPU partitioning | MIG disabled |
| NVIDIA driver | 590.48.01 |
| Reported CPU model | AMD EPYC 7763 64-Core Processor |
| Prediction process CPU affinity | 32 logical CPUs, IDs 0–31 |
| Host-visible system memory | 1,861,362,676 KiB (approximately 1.73 TiB) |
| Architecture | x86_64 |
| Operating system | AlmaLinux 10.0 (Purple Lion) |
| Linux kernel | 6.12.0-124.27.1.el10_1.x86_64 |

## AlphaFold source and software

| Component | Version or setting |
|---|---|
| AlphaFold 3 release tag | [v3.0.4](https://github.com/google-deepmind/alphafold3/tree/v3.0.4) |
| AlphaFold 3 source commit | [85c4d20505fd5cef05eac22b534d4e793971ae69](https://github.com/google-deepmind/alphafold3/tree/85c4d20505fd5cef05eac22b534d4e793971ae69) |
| Package-internal version | 3.0.2 |
| Python | 3.12.12 |
| uv | 0.9.24 |
| HMMER | 3.4, with the AF3 `jackhmmer_seq_limit.patch` |
| JAX / jaxlib / JAX CUDA plugin | 0.10.2 / 0.10.2 / 0.10.2 |
| Haiku | 0.0.16 |
| Tokamax | 0.0.12 |
| RDKit | 2025.9.4 |
| NumPy | 2.4.1 |
| CUDA runtime package (`nvidia-cuda-runtime-cu12`) | 12.9.79 |
| cuDNN package (`nvidia-cudnn-cu12`) | 9.17.1.4 |
| System CUDA compiler | nvcc 13.1.115 |
| GCC | 14.3.1; Red Hat 14.3.1-2 |
| glibc | 2.39 |

The v3.0.4 tag points to the source commit above. That commit's [`version.py`](https://github.com/google-deepmind/alphafold3/blob/85c4d20505fd5cef05eac22b534d4e793971ae69/src/alphafold3/version.py) defines the internal version as `3.0.2`, which appears in the native CIF metadata. The release tag and internal version are recorded separately.

The complete installed third-party Python package versions are listed in [alphafold3-python-packages.txt](alphafold3-python-packages.txt). AlphaFold 3 itself is identified by the source commit above. The package list records the installed environment; source installation uses the AF3 lockfile and build steps below.

## Installation and runtime settings

| Setting | Value |
|---|---|
| Dependency installation | `uv sync --frozen --all-groups --no-editable --offline`, using previously obtained dependencies |
| Chemical component dictionary | Built with `uv run --no-sync build_data` |
| HMMER and CMake build parallelism | 16 workers |
| CMake compatibility option | `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` |
| Optional Python bytecode precompilation | Disabled |
| Process open-file soft limit | 8,192 |
| JAX backend | `gpu` |
| GPU device index / `CUDA_VISIBLE_DEVICES` | `0` / `0` |
| `XLA_FLAGS` | `--xla_gpu_enable_triton_gemm=false` |
| `XLA_PYTHON_CLIENT_PREALLOCATE` | `false` |
| `XLA_CLIENT_MEM_FRACTION` | `0.5` |
| `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` | `8` / `8` / `8` |

Source-file identifiers:

| File | SHA-256 |
|---|---|
| AF3 `uv.lock` | `ff4275f82dcd554c3244ab92703ab7097f5c3844f403900ae08a3a375513bdd9` |
| `jackhmmer_seq_limit.patch` | `df9e3ae35ad1659921d96ebfca67a9616a7a467ddde2be18a56f9bd3edb38c41` |

## Native prediction files

The native output layout for a job with one seed and five diffusion samples contains 21 files:

| Location | Files |
|---|---|
| Each of five `seed-<seed>_sample-<sample>/` directories | `<job>_seed-<seed>_sample-<sample>_model.cif`, `_confidences.json` and `_summary_confidences.json` |
| Job directory | `<job>_model.cif`, `<job>_confidences.json` and `<job>_summary_confidences.json` for the top-ranked sample |
| Job directory | `<job>_ranking_scores.csv`, `<job>_data.json` and `TERMS_OF_USE.md` |

The processed input in `<job>_data.json` uses the `alphafold3` dialect, schema version 4, integer `modelSeeds`, and `protein` entities with explicit chain IDs. It contains the MSA and template data used for inference. Prediction identities are the `(seed, sample)` pairs in the ranking CSV; the CSV retains the ranking-score precision used for model selection. Summary JSON scores are rounded. The top-level structure and confidence files represent one of the sampled predictions.

Native structures are mmCIF files, accompanied by `TERMS_OF_USE.md`. PDB files are conversion products. The [AlphaFold 3 output terms](https://github.com/google-deepmind/alphafold3/blob/85c4d20505fd5cef05eac22b534d4e793971ae69/OUTPUT_TERMS_OF_USE.md) and [model parameter terms](https://github.com/google-deepmind/alphafold3/blob/85c4d20505fd5cef05eac22b534d4e793971ae69/WEIGHTS_TERMS_OF_USE.md) describe the corresponding materials.

## Subsequent CIF-to-PDB conversion

The recorded separate conversion environment used Biopython 1.85 and NumPy 2.3.5. For the [conversion tool supplied in this repository](../src/modeling/docking/cif_to_pdb.py), use the [docking environment](README.md#docking). The tool reads the supplied mmCIF with `MMCIFParser` and writes PDB coordinates and B-factor values with `PDBIO`. See [AF3 workflow guide](../REPRODUCIBILITY.md#alphafold-3-predictions) for source-file mapping and accompanying terms.
