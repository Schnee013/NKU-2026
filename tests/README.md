# Tests

Run commands from the repository root using the appropriate [environment](../environments/README.md). Tests generate synthetic inputs in temporary directories and check their outputs against known arithmetic, geometric and analytical expectations.

| Workflow | Command | Main checks |
|---|---|---|
| AF3 input and output parsing | `python -m unittest discover -s tests -p 'test_af3.py' -v` | Fabricated native and numbered outputs, ranking precision, sequence identity, PAE and pLDDT |
| Six-state pathway | `python -m unittest discover -s tests -p 'test_ode*.py' -v` | Independent full-state integration, conservation, zero-input analytic solution, CLI output |
| Paired-channel generation and regression | `python -m unittest discover -s tests -p 'test_dual*.py' -v` | Actual ODE generation, independent regression equations, units, interval/QC gates and economic adapter |
| Economic accounting | `python -m unittest discover -s tests -p 'test_eco*.py' -v` | Independent cash-flow equations, thresholds, report serialization and output protection |
| Gpa1 sequence design | `python -m unittest discover -s tests -p 'test_gprotein.py' -v` | Graft sequences, prediction input identities and parser boundaries |
| MD analysis | `python -m unittest discover -s tests -p 'test_md*.py' -v` | Seven real drivers on generated geometry, independent fits, contact denominators, polar primitives and bilingual summary figures from synthetic tables |
| Cross-receptor tools | `python -m unittest discover -s tests -p 'test_docking*.py' -v` | Synthetic score/coordinate/contact joins, malformed inputs and protected outputs |
| Mutation tools | `python -m unittest discover -s tests -p 'test_mutant*.py' -v` | Invented complete saturation schema, known parent contrasts, 2 Å pose displacement and preparation contracts |

The full MD collection needs [MD compatibility requirements](../environments/md-compatibility-requirements.txt), including RDKit. Mutation geometry also requires RDKit. Preparation contract tests substitute the external preparation process; they do not run a full docking search. AF3 parser tests fabricate output files and do not execute AlphaFold.

In an environment containing all required modules:

```sh
python -m unittest discover -s tests -v
```

To check the consistency of your own analysis tables, supply the complete dataset to the result-table checker described in [reproducibility](../REPRODUCIBILITY.md). Synthetic tests establish software behavior; biological and field validation require independent measurements.
