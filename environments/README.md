# Software environments

Programs run from the repository root. The pinned portable analysis environments target CPython 3.12.14 on Linux x86_64. The documented AlphaFold 3 deployment uses CPython 3.12.12, and the original MD analysis environment uses CPython 3.11.17. Create a separate environment for each dependency set.

| Workflow | Environment and dependencies |
|---|---|
| AlphaFold 3 prediction | [Hardware, software and runtime settings](alphafold3.md); [installed Python package versions](alphafold3-python-packages.txt) |
| Six-state pathway | [pathway-requirements.txt](pathway-requirements.txt) |
| Paired-channel generation and regression | [dual-channel-requirements.txt](dual-channel-requirements.txt) |
| Docking and chemical contacts | [docking-requirements.txt](docking-requirements.txt) |
| MD analysis | [Python 3.12 requirements](md-requirements.txt), [static compatibility and MD software tests](md-compatibility-requirements.txt), [original Python 3.11 environment](md-analysis-original.yml); [hardware and software specifications](../studies/receptor-engineering/molecular-dynamics/ENVIRONMENT.md) |
| Gpa1 sequence/confidence tools | Python 3.12 standard library |
| Economic accounting | Python 3.12 standard library |

Docking requires `libXrender.so.1`, `libXext.so.6`, `libX11.so.6` and `libexpat.so.1` for ProLIF's RDKit drawing support. On Debian these are supplied by `libxrender1`, `libxext6`, `libx11-6` and `libexpat1`; the Dockerfile includes them for the docking environment.

For example, from the repository root:

```sh
python3.12 -m venv /tmp/nku-pathway-env
/tmp/nku-pathway-env/bin/python -m pip install -r environments/pathway-requirements.txt
/tmp/nku-pathway-env/bin/python -m src.modeling.ode.run_literature_core --output /tmp/nku-pathway-default
```

This command writes `timecourses.csv` and `run.json` to the new directory `/tmp/nku-pathway-default`. The [pathway study](../studies/signal-model/pathway/README.md) defines its parameters, units, results and plotting command. Use the dependency file for the selected workflow and choose a new output directory for each analysis run. [Tests](../tests/README.md) lists the test and numerical-check commands.

## Docking

```sh
python3.12 -m venv /tmp/nku-docking-env
/tmp/nku-docking-env/bin/python -m pip install -r environments/docking-requirements.txt
```

Summary arithmetic for both receptor modules uses the Python standard library; the synthetic generator supplies runnable inputs. SDF geometry requires RDKit. New docking, preparation and contact calculations require user-supplied receptor structures or predictions. Use this environment's Python for the [receptor comparison](../studies/receptor-engineering/receptor-comparison/README.md) and [mutant screen](../studies/receptor-engineering/mutant-screen/README.md). The pinned packages include Meeko 0.7.1, RDKit 2025.09.5, ProLIF 2.1.0, pandas 3.0.2, Biopython 1.87 and Vina 1.2.7. Vina runs through its Python bindings. Commands accept absolute script and configuration paths; use fresh output directories outside the study inputs.

## Containers

The [Dockerfile](Dockerfile) pins its Python base-image digest and takes the requirement file as a build argument:

```sh
docker build --build-arg REQUIREMENTS=pathway-requirements.txt \
  -t nku-pathway -f environments/Dockerfile environments

mkdir -p /tmp/nku-pathway-results

docker run --rm --network none --read-only \
  --tmpfs /tmp:rw,exec,size=1g \
  --mount "type=bind,src=$PWD,dst=/app,readonly" \
  --mount "type=bind,src=/tmp/nku-pathway-results,dst=/results" \
  nku-pathway python -m src.modeling.ode.run_literature_core --output /results/default
```

The container writes the output to `/tmp/nku-pathway-results/default` on the host; use a new output path for another run.
