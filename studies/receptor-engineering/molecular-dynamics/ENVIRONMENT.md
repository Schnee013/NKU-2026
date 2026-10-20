# Computational environments

## Molecular-dynamics production

The original six 100 ns production trajectories used the following environment.

| Component | Specification |
|---|---|
| CPU | Intel Xeon Platinum 8470Q; 104 physical cores / 208 logical processors visible on the host; instance CPU quota equivalent to 25 processors |
| RAM | 90 GiB instance memory limit (96,636,764,160 bytes) |
| GPU | NVIDIA GeForce RTX 5090; 32 GB class, 32,607 MiB visible device memory; compute capability 12.0 |
| Operating-system platform | Linux x86_64 |
| NVIDIA driver | 580.105.08 |
| GROMACS | 2025.4, conda-forge build; mixed precision; 64-bit memory model; AVX2_256 |
| Parallel libraries | Open MPI 5.0.8 and OpenMP; one MPI process per trajectory |
| CUDA | Toolkit/compiler 12.9.86; runtime 12.9 |
| FFT libraries | FFTW 3.3.11 on CPU and cuFFT on GPU |

Production used six OpenMP threads per trajectory, with eight threads in the initial WT–ascr18 segments. Two independent trajectories shared the GPU during paired execution. Nonbonded and PME calculations, supported bonded terms and coordinate updates were assigned to the GPU; CMAP terms were evaluated on the CPU. The simulation MDP files specify `nstlist = 20`; optimized production execution used `-nstlist 200`. The [simulation methods and parameters](METHODS.md#preparation-and-production) define the physical model and integration settings.

The two directed V368R–ascr3 trajectories use the same CPU model, 25-processor quota, 90 GiB memory limit, RTX 5090 and GROMACS/CUDA builds, with NVIDIA driver **595.71.05**. Each uses six OpenMP threads and runtime `-nstlist 100`; PME tuning is disabled. The runtime Verlet list adapts from 1.216 to 1.387 nm. These runtime settings are distinct from the original six-system `-nstlist 200` setting and from the requested MDP `nstlist = 20`.

## Trajectory processing and analysis

Initial preparation also used GROMACS 2025.4. The four predicted/native-lid systems used a double-precision, thread-MPI build for their selected energy minimization, with one MPI process and 16 OpenMP threads on CPU. The two original rebuilt-lid systems and the directed-start systems used the mixed-precision production build. [Per-system preparation parameters](configs/simulation/systems.json) specify the minimization profiles.

| Component | Specification |
|---|---|
| CPU | AMD EPYC 7763 64-Core Processor |
| RAM | Approximately 1.73 TiB visible system memory |
| GPU installed on the analysis host | NVIDIA A100-SXM4-80GB, 81,920 MiB; driver 590.48.01 |
| Operating system | AlmaLinux 10.0; Linux kernel 6.12.0-124.27.1.el10_1.x86_64; glibc 2.39 |
| GROMACS trajectory tools | 2025.4, conda-forge MPI/CUDA build |
| Python | CPython 3.11.17 |
| MDAnalysis | 2.10.0 |
| MDTraj | 1.11.1 |
| NumPy | 2.4.6 |
| pandas | 2.3.3 |
| SciPy | 1.17.1 |
| Matplotlib | 3.11.2 |

The trajectory metrics were calculated on the CPU. The original basic and local analyses used two CPU threads per job. The two directed-start local analyses and common-frame path extraction use the same numerical stack with one CPU/BLAS thread per job. Polar screening and static compatibility also use one CPU thread; the latter uses RDKit 2025.09.5. The [original analysis dependency specification](../../../environments/md-analysis-original.yml) records this numerical stack.

## Portable analysis dependencies

The repository also supplies a CPython 3.12 environment through [fixed requirements](../../../environments/md-requirements.txt): MDAnalysis 2.9.0, MDTraj 1.11.0, NumPy 2.2.6, pandas 2.2.3, SciPy 1.15.3 and Matplotlib 3.10.3.

Static compatibility and the complete MD test collection additionally require RDKit 2025.09.5; [md-compatibility-requirements.txt](../../../environments/md-compatibility-requirements.txt) extends the portable MD dependencies with that package.
