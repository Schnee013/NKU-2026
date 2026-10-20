# Economic model

This module estimates single-season conditional cash flows for a Spanish protected-tomato example with *Meloidogyne incognita*. A separate framework demonstrates how public loss tables and explicit costs enter the accounting.

## Generated scenarios

The [density input](configs/eco_model_example_input.json) and [Spanish model configuration](configs/eco_model_example_model.json) provide example parameters. The commands below generate reports and comparison tables from these inputs. See the [economic Wiki](https://2026.igem.wiki/nku-china/model-eco/) for project results and [reproducibility](../../REPRODUCIBILITY.md) for input requirements.

## Model, units and inputs

The loss function is Seinhorst: relative yield `y=1` for `d<=T`; otherwise `y=m+(1-m)*z**(d/T-1)`. The packaged parameters are `T=6`, `m=0.34`, `z=0.95`. Loss fraction is `1-y`. Spanish density inputs must be finite and within 0–869.6 **J2/100 cm^3 soil**. Costs and net income are **EUR/ha**. Yield is kg/ha and price EUR/kg. Input/config units are checked exactly.

Let `R=yield*price`, `L=R*(1-y)`, `B=production cost`, `C=test cost`, and `(K,a)` be treatment cost and assumed fraction of loss avoided. Guided net income is `R-B-L*(1-a)-K-C`. No-control income is `R-B-L`; uniform-control income is `R-B-L*(1-a_uniform)-K_uniform`. Each reported net benefit is the difference between the appropriate incomes.

[Model configuration](configs/eco_model_example_model.json) supplies all costs, treatments, English descriptions, laboratory fees, sample count and one-way sensitivity ranges. [Density input](configs/eco_model_example_input.json) supplies `density` and `density_unit`. The preset tiers are low below 10, medium from 10 through below 80, and high from 80. Stable API category tokens are `低`, `中`, `高`, displayed as Low, Medium, High in English HTML. These tiers are model-scenario strategies, distinct from the option giving the highest modelled income under the current parameters.

Treatment population reductions 0.637 and 0.7213 are converted through the loss curve at reference density 100 to fixed avoided-loss fractions 0.586805 and 0.691234. Reusing these fractions at other densities is an assumption. Efficacy-source plots started above 50; lower-density use is flagged, including comparator and unselected alternatives. Inputs above 759 receive an additional economic-management support warning. The 5 EUR/ha test cost is a team pricing assumption. The regional yield proxy is not an uninfested control measurement.

[Parameter sources and classifications](sources/parameter-sources.csv) identify the numerical input families. The Spanish scenario combines the spring ungrafted-tomato density–yield fit of Expósito et al. (2020), DOI [10.3389/fpls.2020.560024](https://doi.org/10.3389/fpls.2020.560024), with tomato economic inputs, treatment costs, population reductions and economic-management density estimates from Talavera-Rubia et al. (2022), DOI [10.3390/plants11202774](https://doi.org/10.3390/plants11202774). These studies use different fields. `literature_direct` records source values; `literature_converted` and `literature_derived` identify unit conversions and calculations; `demo_derived`, `team_pricing_assumption` and `scenario_assumption` identify model conversions or assumptions.

## Run

Use Python 3.12 on Linux x86_64, with the standard library only. Run directly from source. From the repository root:

```bash
python -m src.modeling.economic.eco_model 100 --report-dir /tmp/nku-economic-example
python -m src.modeling.economic.run_mvp --out /tmp/nku-economic-framework/results.csv
```

The Spanish command writes four files to `--report-dir` (default: `economic-results/` under the current working directory): `diagnostic_report.txt` and `diagnostic_report.html` for Chinese reading, `diagnostic_report_en.html` for English reading, and `diagnostic_report.json` for programmatic reuse. The browser stays closed unless `--open` is supplied. `--no-report` prints results without writing files. `--break-even-only` prints strict roots and threshold jumps without generating a report. The framework command writes all 21 comparison rows to the specified CSV path. Use a new or empty output directory outside packaged `studies/`; a nonempty directory is rejected.

From another working directory, set `PYTHONPATH` to the absolute repository root before the same `python -m ...` command. Default example configuration paths are resolved relative to the module location, so they do not depend on the current working directory. All scientific computation is local and needs no network.

## Generic framework

[Framework configuration](configs/public_mvp_framework.json) contains three separate examples: a 22-point spring processing-tomato table (spring J2/kg soil, USD/acre), a normalized potato-PCN linear approximation (eggs/g soil), and a threshold-only soybean-SCN example (eggs/100 cc soil). Tomato yield 30 source tons/acre and price 54 USD/source ton are illustrative public values. Potato quantities are normalized, not local monetary income. Soybean monetary fields stay null. The density–yield and threshold sources are the UC IPM spring processing-tomato table, the EPPO potato-PCN approximation and the UMN SCN management bands. Spring J2/kg soil cannot be converted to the Spanish or dual-channel J2/100 cm³ soil interface without additional soil information. Effect and early-program-cost sweeps are scenario assumptions. Piecewise endpoints and intervention effects are clamped as specified by the framework; this behavior does not relax the Spanish model's density rejection.

## Sensitivity and paired-channel inputs

Sensitivity is one parameter changed at a time, **not a confidence interval**. Input density is treated as known and accurate; sampling uncertainty, sensor sensitivity/specificity and misclassification are absent. The uniform comparison is conditional on density, not an ex-ante expected value of information. [Dual-channel calibration](../signal-model/dual-channel/README.md) may pass an accepted, unit-compatible density through its explicit quality and range checks before this model is evaluated.
