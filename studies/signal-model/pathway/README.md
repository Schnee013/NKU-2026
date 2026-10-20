# Six-state pathway and effective color model

This model connects a constant normalized ligand input to a conserved yeast
G-protein/scaffold/Fus3/Sst2 reaction core and an effective color output.
Simulations start from the supplied reaction rates, initial concentrations and
input/output parameters. Reaction rates come from the literature model described
below; input coupling and color parameters are assumed for this computational
scenario.

The independent states are `G_trimer`, `Ga_GTP`, `scaffold_bound`,
`Fus3_active`, `Sst2_active`, and `color_au`. Five further molecular states are
calculated algebraically from conserved pools. Concentrations are in micromolar
and solver time is in seconds. `ligand_normalized` is dimensionless and
`color_au` is an arbitrary color unit. Conversion to ascaroside concentration,
absorbance, nematode density or economic benefit requires a separate calibration.
The [dual-channel study](../dual-channel/README.md) specifies its own simulated
optical conditions and channel parameters.

## Inputs and generated outputs

| Material | Purpose |
|---|---|
| [Source configuration](configs/literature_core_source.json) | Literature rates, initial concentrations, Delft–Supady identifiers, reaction names, source locations and units |
| [Project configuration](configs/literature_core_project.json) | Assumed input and effective output parameters |
| [Alternative configuration](configs/alternate_project.json) | Assumed endpoint scenario with zero color loss and a fractional Hill exponent |

The commands below generate time courses and figures from the supplied parameters. Scientific results are on the [signaling Wiki](https://2026.igem.wiki/nku-china/model-ode/); [tests](../../../tests/README.md) check the solver against independent conservation and analytic expectations.

CSV columns are `ligand_normalized`, `time_s`, the 11 states listed below, and
`MAPK_fraction = Fus3_active / total_Fus3`. Molecular states are in micromolar;
`color_au` has arbitrary color units and `MAPK_fraction` is dimensionless.
The time grid combines dense early logarithmic sampling with uniform sampling
across the full duration.

## Run

Use CPython 3.12 on Linux x86_64 and install the pinned dependencies from
the repository root:

```sh
python -m venv /tmp/pathway-env
. /tmp/pathway-env/bin/activate
python -m pip install -r environments/pathway-requirements.txt
python -m src.modeling.ode.run_literature_core --output /tmp/pathway-default
python -m src.modeling.ode.plot_literature_core /tmp/pathway-default/timecourses.csv --output /tmp/pathway-default/timecourses.png
```

The runner requires a new output directory, including when a previous directory
is empty. The plotting command requires a new output file. Both commands protect
input and packaged reference locations from writes. A failed calculation may
leave a reserved directory; restart with a fresh output path.

To run the alternative parameter scenario from outside the checkout:

```sh
pathway_repo=$(pwd)
cd /tmp
python "$pathway_repo/src/modeling/ode/run_literature_core.py" --parameters "$pathway_repo/studies/signal-model/pathway/configs/alternate_project.json" --ligand 0 .2 1 --output /tmp/pathway-alternate
python "$pathway_repo/src/modeling/ode/plot_literature_core.py" /tmp/pathway-alternate/timecourses.csv --output /tmp/pathway-alternate/timecourses.png
```

The runner accepts `--parameters`, `--source`, `--ligand`, `--hours` and
`--output`. Defaults are the project/source configurations, ligand inputs
`0 .1 1 10`, and eight hours. It saves `timecourses.csv` and `run.json`; the latter
contains the resolved numerical inputs, fixed rates and pools, model-state order,
duration, ligand values and software versions.

## Equations and conservation

Let `H = G_trimer`, `A = Ga_GTP`, `Q = scaffold_bound`,
`F* = Fus3_active`, `Z* = Sst2_active`, and `C = color_au`.
The five conserved pools, calculated from the source initial state, are:

| Pool | Definition | Default total (micromolar) |
|---|---|---:|
| Ta | H + A + D | 1.667 |
| Tb | H + B + Q | 1.667 |
| Ts | S + Q | 0.036 |
| Tf | F + F* | 0.686 |
| Tz | Z + Z* | 0.505 |

Reconstruction uses `D = Ta-H-A`, `B = Tb-H-Q`, `S = Ts-Q`,
`F = Tf-F*`, and `Z = Tz-Z*`. Output ordering is
`H, A, D, B, S, Q, F, F*, Z, Z*, C`, with CSV names
`G_trimer, Ga_GTP, Ga_GDP, Gbg_free, scaffold_free, scaffold_bound,
Fus3, Fus3_active, Sst2, Sst2_active, color_au`.

For constant ligand `L`, define
`u = u_max L/(K_L+L)` and `m = F*/Tf`. Derivatives with respect to seconds are:

```text
dH/dt  = -u H + k5 D B
dA/dt  =  u H - (k4 + k41 Z*) A
dQ/dt  =  k6 S B - (k7 + k71 Z*) Q
dF*/dt =  k8 Q F - (k9 + k91 Z) F*
dZ*/dt =  k10 F* Z - k11 Z*
dC/dt  = [b + a m^n/(K^n + m^n) - delta C] / 3600
```

The Fus3 derivative uses activation minus inactivation (`v11-v12` in Delft
Table 4.1), conserving total Fus3. The printed active-state derivative in Delft
Table 4.3 uses their sum. Delft K12 corresponds to `k91`: its flux multiplies Sst2
and active Fus3 concentrations, so the rate has bimolecular units
`uM^-1 s^-1`, with the supplied numerical value 0.017. The printed Delft unit is
`s^-1`.

The eleven internal rate entries are `k4, k41, k5, k6, k7, k71, k8,
k9, k91, k10, k11`, including the source value `k9 = 0`. First-order rates have
units `s^-1`; bimolecular rates have units `uM^-1 s^-1`.
The source receptor rates `k1, k2, k12, k3` remain identified in the
configuration. This model replaces receptor dynamics with the sustained drive
`u`; `k3 × receptor_initial_uM` sets the default demonstration scale for
`u_max`.

The six-state representation eliminates five dependent states by conservation.
It introduces no fast-equilibrium approximation. Constant protein pools over
hours and the sustained receptor drive are model assumptions. SciPy `solve_ivp`
uses Radau with relative tolerance `1e-9` and absolute tolerance `1e-12`.
A negative numerical intermediate `m` is replaced by zero only when evaluating
the fractional Hill power, keeping intermediate evaluations real. The runner
rejects nonfinite trajectories and concentrations below `-1e-8` micromolar.

## Effective output and units

The seven endpoint parameters are `u_max`, `K_L`, `b`, `a`, `K`, `n` and
`delta`; all default values are assumed. The default output Hill exponent is
`n = 2`.

Write the transcription activity as
`TF = bT + E m^n/(KT^n+m^n)`, output production as
`Prod = l + N g TF/(D0+TF)`, and color dynamics in hours as
`dC/dh = p Prod - delta C`. Exact algebraic merging gives:

```text
b = p [l + N g bT/(D0+bT)]
a = p N g D0 E / [(D0+bT)(D0+bT+E)]
K = KT [(D0+bT)/(D0+bT+E)]^(1/n)
```

The exponent `n` and loss constant `delta` retain their values.
`b` and `a` are color units per hour; `delta` is per hour. Division by
3600 converts their contribution to solver seconds. At zero ligand the core
stays at its initial state, and color has the independent solution
`C(t) = b/delta [1-exp(-delta t/3600)]`; for `delta = 0`,
`C(t) = b t/3600`. The mapping of transcriptional activity to the Fus3 fraction
is assumed. This effective output combines transcription, enzyme accumulation
and pigment maturation into a single color state.

## Literature definition and applicability

The reaction core follows Supady's 2009 Model B
([thesis](https://rumo.biologie.hu-berlin.de/tbp/images/Theses/bsc_supady.pdf),
PDF pages 48–49, sections 5.1 and 6.1) and the Delft 2012
[single-cell model](https://2012.igem.org/Team:TU-Delft/Modeling/SingleCellModel)
([details](https://static.igem.org/mediawiki/2012/0/05/ReducedModelDetails.pdf),
Tables 4.1–4.3). The source configuration maps Delft K1–K14 to the
14 shared nonzero Supady rates and their reactions; Supady additionally
specifies zero basal Fus3 inactivation `k9`.
Initial concentrations are supplied in the source configuration.

The Sst2/Fus3 terms represent the source model's effective feedback mechanism.
Applying the model assumes functional Sst2 feedback and transfer of the source
rates and conserved pools. A strain lacking SST2, or altered GTP hydrolysis or
G-protein reassembly, requires corresponding changes to the equations or rates.
The [Gpa1 study](../../gprotein-coupling/README.md) covers preliminary sequence
design and prediction comparison. The color parameters represent the effective
project output; native Ste12 and Delft GFP expression parameters are separate
from this output model.
