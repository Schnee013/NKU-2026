"""Fit paired-signal regression and evaluate conditional economic scenarios."""
import argparse
import csv
import json
from collections import Counter
from copy import deepcopy
from pathlib import Path

import numpy as np

from .calibration import fit, predict, eco_scenarios
from ._output import new_output

PROTOCOL_ID = "SIM_V4_MI1TO10_CE10TO1_Q1_8H_A480"
from src.modeling.economic.eco_model_pkg.config import load_model


def read_pairs(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        required = {'sample_id', 'target', 'background', 'x3', 'x18', 'data_origin'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('missing paired columns: ' + ', '.join(sorted(required-set(reader.fieldnames or []))))
        rows = list(reader)
    if not rows or len({r['sample_id'] for r in rows}) != len(rows):
        raise ValueError('nonempty unique paired samples required')
    for row in rows:
        if row['data_origin'] != 'synthetic_ODE':
            raise ValueError('replay accepts synthetic ODE fixtures only')
        for key in ('target', 'background', 'x3', 'x18'):
            row[key] = float(row[key])
            if not np.isfinite(row[key]):
                raise ValueError('nonfinite paired value')
        if min(row['target'], row['background']) < 0:
            raise ValueError('reference densities must be nonnegative')
    return rows


def error_metrics(rows, estimates):
    residual = np.asarray(estimates)-np.array([r['target'] for r in rows])
    return dict(n=len(rows), rmse=float(np.sqrt(np.mean(residual**2))),
                mae=float(np.abs(residual).mean()), max_absolute_error=float(np.abs(residual).max()))


def serialize_predictions(results):
    """Group shared metadata when saving complete prediction/economic results."""
    if not results:
        raise ValueError('nonempty prediction batch required')
    samples = deepcopy(results)
    fields = ('protocol_id', 'evidence_status', 'reference_basis',
              'calibration_density_range', 'density_unit')
    first = samples[0]['prediction']
    if any(key not in first for key in fields):
        raise ValueError('complete prediction metadata required')
    metadata = {key: deepcopy(first[key]) for key in fields}
    for sample in samples:
        prediction = sample['prediction']
        for key, value in metadata.items():
            if key not in prediction or prediction[key] != value:
                raise ValueError(f'inconsistent prediction metadata: {key}')
            del prediction[key]
        eco = sample['eco']
        for key in ('evidence_status', 'reference_basis'):
            if key in eco:
                if eco[key] != metadata[key]:
                    raise ValueError(f'inconsistent economic metadata: {key}')
                del eco[key]
        for scenario in eco['scenarios']:
            if scenario.get('density_unit') != metadata['density_unit']:
                raise ValueError('inconsistent economic density unit')
            del scenario['density_unit']
    return dict(metadata=metadata, samples=samples)


def run(input_dir, output, eco_config, protocol_id):
    input_dir, output, eco_config = Path(input_dir), Path(output), Path(eco_config)
    if protocol_id != PROTOCOL_ID:
        raise ValueError('protocol does not match the synthetic signal inputs')
    training = read_pairs(input_dir/'calibration.csv')
    evaluation = read_pairs(input_dir/'evaluation.csv')
    if {r['sample_id'] for r in training} & {r['sample_id'] for r in evaluation}:
        raise ValueError('calibration/evaluation sample overlap')
    if {(r['target'], r['background']) for r in training} & {(r['target'], r['background']) for r in evaluation}:
        raise ValueError('evaluation combinations must differ from calibration')
    x = np.array([[r['x3'], r['x18']] for r in training])
    y = np.array([r['target'] for r in training])
    model = fit(x, y, sample_ids=[r['sample_id'] for r in training], protocol_id=protocol_id,
                evidence_status='simulated', reference_basis='spike_count')
    eco = load_model(eco_config)
    variants = [('channel3', [0]), ('channel18', [1]), ('dual', [0, 1])]
    betas = {name: np.linalg.lstsq(np.column_stack((np.ones(len(x)), x[:, cols])), y, rcond=None)[0]
             for name, cols in variants}
    results = []
    for row in evaluation:
        signal = np.array([row['x3'], row['x18']])
        prediction = predict(model, signal, protocol_id=protocol_id, saturated=False, qc_pass=True)
        # Saturation and QC are explicit assumptions of the simulated generator.
        estimates = {name: float(np.r_[1., signal[cols]] @ betas[name]) for name, cols in variants}
        results.append(dict(sample_id=row['sample_id'], ratio=row['mi_ascr3_to_ascr18'],
                            reference_density=row['target'], background=row['background'],
                            estimates=estimates, prediction=prediction, eco=eco_scenarios(prediction, eco)))
    metrics = {}
    for ratio in sorted({r['mi_ascr3_to_ascr18'] for r in evaluation}):
        indices = [i for i, r in enumerate(evaluation) if r['mi_ascr3_to_ascr18'] == ratio]
        metrics[ratio] = {name: error_metrics([evaluation[i] for i in indices],
                                             [results[i]['estimates'][name] for i in indices]) for name, _ in variants}
        selected = [results[i] for i in indices]
        def band(density):
            if density < 0 or density > eco.max_density:
                return 'invalid_physical_density'
            return '低' if density < eco.middle_from_density else ('中' if density < eco.high_from_density else '高')
        metrics[ratio].update(
            valid_point_n=sum(r['prediction']['within_calibration_range'] for r in selected),
            eco_conditional_n=sum(r['eco']['status'] == 'conditional_demo_only' for r in selected),
            raw_point_band_match_n=sum(band(r['reference_density']) == band(r['estimates']['dual']) for r in selected),
            note='Band matches use nonnegative points within economic-model support; all 12 points enter the error metrics.')
    economic_parameters = {name: getattr(eco, name) for name in (
        'density_unit', 'area_unit', 'currency', 'tolerance_density', 'minimum_relative_yield',
        'seinhorst_z', 'max_density', 'baseline_yield_kg_per_ha', 'price_per_kg',
        'baseline_production_cost_per_ha', 'test_cost_per_ha',
        'middle_from_density', 'high_from_density')}
    economic_parameters['treatments'] = {
        name: dict(cost_per_ha=getattr(eco, name).cost_per_ha,
                   avoided_loss_fraction=getattr(eco, name).avoided_loss_fraction)
        for name in ('without_test', 'low', 'middle', 'high')}
    summary = dict(evidence_status='simulated', protocol_id=protocol_id,
                   density_unit=eco.density_unit,
                   calibration_metrics=error_metrics(training, np.c_[np.ones(len(x)), x] @ model.coefficients),
                   scenario_metrics=metrics,
                   input_files=[dict(path=str(path), role=role) for path, role in (
                       (input_dir/'calibration.csv', 'paired regression calibration'),
                       (input_dir/'evaluation.csv', 'ratio comparisons'),
                       (eco_config, 'conditional economic scenarios'))],
                   economic_parameters=economic_parameters,
                   counts=dict(evaluation=len(results), valid_point=sum(r['prediction']['within_calibration_range'] for r in results),
                               eco_conditional=sum(r['eco']['status'] == 'conditional_demo_only' for r in results)),
                   rejection_reasons=dict(Counter(reason for r in results for reason in r['prediction']['reasons'])),
                   thresholds=dict(middle=eco.middle_from_density, high=eco.high_from_density, unit=eco.density_unit),
                   frozen_during_ratio_test=['both ODE configurations', 'dual and single regression coefficients',
                                             'Ce ratio', 'common Q'])
    model_data = dict(equation='d=beta0+beta3*x3+beta18*x18',
                      coefficients=model.coefficients.tolist(), covariance=model.covariance.tolist(),
                      residual_variance=model.residual_variance, df=model.df,
                      hull_equations=model.hull_equations.tolist(),
                      calibration_density_range=list(model.density_range),
                      condition_number=model.condition_number,
                      input_unit='blank-corrected A480, 1 cm effective path; both channels paired',
                      density_unit=eco.density_unit, protocol_id=model.protocol_id,
                      evidence_status=model.evidence_status, reference_basis=model.reference_basis,
                      interval_interpretation='OLS formal residual-based interval only; deterministic simulated data, no empirical 95% coverage')
    batch = serialize_predictions(results)
    with new_output(output, input_dir, eco_config) as staging:
        for filename, data in [('summary.json', summary), ('model.json', model_data), ('predictions.json', batch)]:
            (staging/filename).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        plot(results, staging/'dual-channel-check.png')
    return summary


def plot(results, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    baseline = [r for r in results if r['ratio'] == '1:10']
    for name in ('channel3', 'channel18', 'dual'):
        axes[0].scatter([r['reference_density'] for r in baseline],
                        [r['estimates'][name] for r in baseline], label=name)
    axes[0].plot([-10, 100], [-10, 100], '--', color='grey')
    axes[0].set(xlabel='Assumed Mi density (J2/100 cm3)', ylabel='Estimated Mi density (J2/100 cm3)', title='New combinations: Mi 1:10')
    for ratio in sorted({r['ratio'] for r in results}):
        rows = [r for r in results if r['ratio'] == ratio]
        axes[1].scatter([r['reference_density'] for r in rows],
                        [r['estimates']['dual']-r['reference_density'] for r in rows], label=ratio)
    axes[1].axhline(0, linestyle='--', color='grey')
    axes[1].set(xlabel='Assumed Mi density (J2/100 cm3)', ylabel='Estimate minus reference (J2/100 cm3)', title='Fixed coefficients; changed Mi ratio')
    for ax in axes:
        ax.legend()
        ax.grid(alpha=.2)
    fig.suptitle('Deterministic paired-sensor ODE simulations')
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--eco-config', type=Path, required=True)
    parser.add_argument('--protocol-id', required=True)
    args = parser.parse_args()
    run(args.input_dir, args.output, args.eco_config, args.protocol_id)
