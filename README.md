# NKU-2026

Dry-lab software for the NKU-China 2026 iGEM project: receptor docking and mutation comparisons, molecular-dynamics analysis, sequence design, signaling models and conditional economic analysis.

This repository provides source code, dependencies, methods, sequence and chemical inputs, parameter configurations and synthetic examples. Use the examples to explore each workflow, or supply your own data for analysis. Scientific results and their interpretation are presented on the [team Wiki](https://2026.igem.wiki/nku-china/model/).

| Tool | Documentation |
|---|---|
| Cross-receptor docking and contact analysis | [Receptor comparison](studies/receptor-engineering/receptor-comparison/README.md) |
| Mutation-score, contact and parent comparisons | [Mutant screening](studies/receptor-engineering/mutant-screen/README.md) |
| Trajectory geometry, contacts and local analyses | [Molecular dynamics](studies/receptor-engineering/molecular-dynamics/README.md) |
| Terminal-graft sequence design and prediction parsing | [Gpa1 tools](studies/gprotein-coupling/README.md) |
| Six-state signaling and effective color output | [ODE pathway](studies/signal-model/pathway/README.md) |
| Synthetic paired-sensor generation and regression | [Dual-channel model](studies/signal-model/dual-channel/README.md) |
| Conditional economic scenarios | [Economic model](studies/economic-model/README.md) |

## Get started

Use the checkout supplied through the team's Software deliverable, or the [GitHub source repository](https://github.com/Schnee013/NKU-2026). Commands below run from its root:

```sh
git clone https://github.com/Schnee013/NKU-2026.git
cd NKU-2026
```

Choose the appropriate [environment](environments/README.md). For a standard-library-only example:

```sh
python -m src.modeling.economic.eco_model 100 --report-dir runs/examples/economics
```

With the MD dependencies installed, generate a small artificial trajectory, run all seven analysis drivers and check their numerical outputs:

```sh
python -m examples.md_geometry --output runs/examples/md-geometry
```

The [examples guide](examples/README.md) also covers synthetic docking records, mutation comparisons and ODE-generated paired readings. Generated outputs go to new directories under `runs/`; they are excluded from version control. A synthetic example demonstrates software behavior and does not validate a biological claim.

See [tests](tests/README.md) for executable checks and [reproducibility](REPRODUCIBILITY.md) for each workflow's input requirements. Reproducing a scientific MD result requires the original trajectories, matching topology and analysis settings.

Team-owned original source code, analysis scripts and tests use the [Apache License 2.0](LICENSE). [Attributions](ATTRIBUTIONS.md) describe scientific sources and the separate terms for third-party materials and dependencies.
