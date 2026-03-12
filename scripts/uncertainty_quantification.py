#!/usr/bin/env python3
"""
Uncertainty Quantification for ELM Out-of-Sample Predictions

Latin Hypercube Sampling over the top-15 OAT-sensitive frozen priors,
each within its documented plausible range. For each sample:
  1. Perturb all 15 parameters simultaneously
  2. Re-calibrate: root-find the 6 ITP compound scale factors so male
     predictions still match observed ITP targets exactly
  3. Record all 6 female predictions and the rapa+acarb combo prediction

Output: percentile bands (5th, 25th, 50th, 75th, 95th) on each prediction.

Usage:
    python scripts/uncertainty_quantification.py              # 500 samples
    python scripts/uncertainty_quantification.py --samples 1000
    python scripts/uncertainty_quantification.py --verbose
"""

import sys
import time
import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
from scipy.optimize import brentq
from scipy.stats.qmc import LatinHypercube

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from elm.model import simulate, run_control, calculate_lifespan_extension
from elm.compounds import (
    get_compound, get_itp_start_time, ITP_VALIDATION, COMPOUNDS,
)
from elm.pathways import ALL_PATHWAY_PARAMS
from elm.sex_mechanisms import apply_sex_modifier

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
T_MAX = 2.5
DT = 0.01  # coarse dt for speed (matches OAT)

ITP_COMPOUNDS = ['rapamycin', 'acarbose', 'canagliflozin',
                 '17_alpha_estradiol', 'aspirin', 'glycine']

# Pathway ratios for each ITP compound
COMPOUND_RATIOS: Dict[str, Dict[str, float]] = {}
for cname in ITP_COMPOUNDS:
    inj = COMPOUNDS[cname]
    total = sum(inj.values())
    COMPOUND_RATIOS[cname] = {k: v / total for k, v in inj.items()}

# Male targets
MALE_TARGETS = {cname: ITP_VALIDATION[cname].target_male
                for cname in ITP_COMPOUNDS}

# ---------------------------------------------------------------------------
# Top-15 OAT-sensitive parameters with documented plausible ranges
# Format: (subsystem, key, lo, hi)
# Ranges from PARAM_AUDIT2.md
# ---------------------------------------------------------------------------
TOP15_PARAMS = [
    ('heteroplasmy', 'transcription_advantage', 0.15, 0.30),
    ('heteroplasmy', 'time_scale',              8.0,  18.0),
    ('heteroplasmy', 'k_selection_foxo3',       0.05, 0.20),
    ('nad',          'k_nampt_base',            0.40, 0.70),
    ('mtorc1',       'k_mtorc1_damage_reduction', 0.15, 0.40),
    ('sirtuin',      'k_sirt1_foxo3',           0.60, 1.00),
    ('tsc2',         'k_tsc2_mtorc1_max',       0.50, 0.80),
    ('nad',          'k_consumption_basal',      0.25, 0.55),
    ('mtorc1',       'k_mtorc1_senescence',     0.20, 0.50),
    ('dna_repair',   'k_damage_ros',            0.40, 0.80),
    ('sirtuin',      'Km_nad_sirt1',            0.30, 0.60),
    ('tsc2',         'Km_tsc2_mtorc1',          0.30, 0.70),
    ('nad',          'k_cd38_base',             0.25, 0.50),
    ('nad',          'k_parp_consumption',       0.50, 0.90),
    ('methylation',  'k_ezh2_base',             1.50, 2.20),
]


# ---------------------------------------------------------------------------
# Core functions (reused from OAT)
# ---------------------------------------------------------------------------

def build_injections(compound: str, alpha: float) -> Dict[str, float]:
    """Build injection dict from pathway ratios and scale factor."""
    ratios = COMPOUND_RATIOS[compound]
    inj = {k: v * alpha for k, v in ratios.items()}
    inj['start_time'] = get_itp_start_time(compound)
    return inj


def male_extension(compound: str, alpha: float, control_m) -> float:
    """Compute male lifespan extension for a given scale factor."""
    inj = build_injections(compound, alpha)
    inj_sex = apply_sex_modifier(compound, inj, 'M')
    result = simulate(interventions=inj_sex, sex='M', t_max=T_MAX, dt=DT)
    return calculate_lifespan_extension(result, control_m)


def calibrate_compound(compound: str, control_m, target: float,
                       alpha_hint: float = None) -> float:
    """Root-find scale factor so male prediction matches target."""
    def objective(alpha):
        return male_extension(compound, alpha, control_m) - target

    if alpha_hint is not None and alpha_hint > 0:
        alpha_lo = alpha_hint * 0.3
        alpha_hi = alpha_hint * 3.0
    else:
        alpha_lo = 0.01
        alpha_hi = 3.0

    f_lo = objective(alpha_lo)
    f_hi = objective(alpha_hi)
    if f_lo * f_hi > 0:
        alpha_lo = 0.01
        alpha_hi = 5.0
        f_lo = objective(alpha_lo)
        f_hi = objective(alpha_hi)
        if f_lo * f_hi > 0:
            total = sum(COMPOUNDS[compound].values())
            return total

    return brentq(objective, alpha_lo, alpha_hi, xtol=5e-3)


_BASELINE_ALPHAS: Dict[str, float] = {}


def calibrate_all(control_m) -> Dict[str, float]:
    """Root-find all 6 scale factors against male targets."""
    alphas = {}
    for cname in ITP_COMPOUNDS:
        hint = _BASELINE_ALPHAS.get(cname)
        alphas[cname] = calibrate_compound(
            cname, control_m, MALE_TARGETS[cname], alpha_hint=hint)
    return alphas


def female_predictions(alphas: Dict[str, float],
                       control_f) -> Dict[str, float]:
    """Compute all 6 female predictions with given scale factors."""
    preds = {}
    for cname in ITP_COMPOUNDS:
        inj = build_injections(cname, alphas[cname])
        inj_sex = apply_sex_modifier(cname, inj, 'F')
        result = simulate(interventions=inj_sex, sex='F', t_max=T_MAX, dt=DT)
        preds[cname] = calculate_lifespan_extension(result, control_f)
    return preds


def combo_prediction(alphas: Dict[str, float], control_m) -> float:
    """Compute rapa+acarb combination prediction (male, 9-month start)."""
    rapa_inj = build_injections('rapamycin', alphas['rapamycin'])
    rapa_inj_sex = apply_sex_modifier('rapamycin', rapa_inj, 'M')
    acarb_inj = build_injections('acarbose', alphas['acarbose'])
    acarb_inj_sex = apply_sex_modifier('acarbose', acarb_inj, 'M')

    combo = {}
    for k, v in rapa_inj_sex.items():
        if k != 'start_time':
            combo[k] = v
    for k, v in acarb_inj_sex.items():
        if k != 'start_time':
            combo[k] = combo.get(k, 0) + v
    combo['start_time'] = 9.0 / 30.0

    result = simulate(interventions=combo, sex='M', t_max=T_MAX, dt=DT)
    return calculate_lifespan_extension(result, control_m)


# ---------------------------------------------------------------------------
# LHS Sampling
# ---------------------------------------------------------------------------

def set_params(values: List[float]):
    """Set all 15 parameters to the given values."""
    for (subsystem, key, _, _), val in zip(TOP15_PARAMS, values):
        ALL_PATHWAY_PARAMS[subsystem][key] = val


def get_original_values() -> List[float]:
    """Get current (original) values of all 15 parameters."""
    return [ALL_PATHWAY_PARAMS[sub][key] for sub, key, _, _ in TOP15_PARAMS]


def restore_params(originals: List[float]):
    """Restore all 15 parameters to their original values."""
    set_params(originals)


def generate_lhs_samples(n_samples: int, seed: int = 42) -> np.ndarray:
    """Generate LHS samples in [0,1]^15, then map to parameter ranges.

    Returns array of shape (n_samples, 15) with actual parameter values.
    """
    sampler = LatinHypercube(d=len(TOP15_PARAMS), seed=seed)
    unit_samples = sampler.random(n=n_samples)  # shape (n, 15) in [0,1]

    # Map to parameter ranges
    samples = np.empty_like(unit_samples)
    for j, (_, _, lo, hi) in enumerate(TOP15_PARAMS):
        samples[:, j] = lo + unit_samples[:, j] * (hi - lo)

    return samples


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_uq(n_samples: int = 500, verbose: bool = False, seed: int = 42):
    """Run LHS uncertainty quantification."""

    print("Uncertainty Quantification for ELM Out-of-Sample Predictions")
    print(f"Samples: {n_samples}, Parameters: {len(TOP15_PARAMS)}")
    print()

    t_start = time.time()

    # Save original values
    originals = get_original_values()

    # Baseline
    print("Computing baseline (unperturbed)...")
    control_m = run_control(sex='M', t_max=T_MAX, dt=DT)
    control_f = run_control(sex='F', t_max=T_MAX, dt=DT)
    base_alphas = calibrate_all(control_m)
    _BASELINE_ALPHAS.update(base_alphas)
    base_f_preds = female_predictions(base_alphas, control_f)
    base_combo = combo_prediction(base_alphas, control_m)

    print("  Baseline female predictions:")
    for cname in ITP_COMPOUNDS:
        target = ITP_VALIDATION[cname].target_female
        pred = base_f_preds[cname]
        print(f"    {cname:25s}  pred={pred:+.1f}%  target={target:+.1f}%")
    print(f"  Rapa+acarb combo: {base_combo:+.1f}% (observed: +34.0%)")
    print()

    # Generate LHS samples
    print("Generating Latin Hypercube samples...")
    samples = generate_lhs_samples(n_samples, seed=seed)

    # Results storage
    pred_keys = [f"F_{c}" for c in ITP_COMPOUNDS] + ['combo_rapa_acarb']
    all_results = {k: [] for k in pred_keys}
    n_failures = 0

    # Sweep
    print(f"Running {n_samples} samples...")
    t_sweep = time.time()

    for i in range(n_samples):
        if (i + 1) % 25 == 0 or verbose:
            elapsed = time.time() - t_sweep
            rate = elapsed / (i + 1)
            eta = rate * (n_samples - i - 1)
            print(f"  [{i+1}/{n_samples}] "
                  f"elapsed={elapsed:.0f}s, ETA={eta:.0f}s")

        # Set perturbed parameter values
        set_params(samples[i])

        try:
            cm = run_control(sex='M', t_max=T_MAX, dt=DT)
            cf = run_control(sex='F', t_max=T_MAX, dt=DT)
            alphas = calibrate_all(cm)
            f_preds = female_predictions(alphas, cf)
            combo = combo_prediction(alphas, cm)

            for cname in ITP_COMPOUNDS:
                all_results[f"F_{cname}"].append(f_preds[cname])
            all_results['combo_rapa_acarb'].append(combo)

        except Exception as e:
            n_failures += 1
            if verbose:
                print(f"    WARNING: sample {i+1} failed: {e}")
        finally:
            restore_params(originals)

    elapsed_total = time.time() - t_start
    n_valid = n_samples - n_failures
    print(f"\nCompleted in {elapsed_total:.0f}s "
          f"({elapsed_total/n_samples:.2f}s per sample)")
    if n_failures > 0:
        print(f"  {n_failures} failures ({n_failures/n_samples*100:.1f}%)")
    print()

    # Compute percentiles
    percentiles = [5, 25, 50, 75, 95]
    print("=" * 95)
    print(f"{'Prediction':<28s}  {'Baseline':>8s}  "
          + "  ".join(f"{'p'+str(p):>7s}" for p in percentiles)
          + f"  {'90% CI':>14s}")
    print("=" * 95)

    summary = {}
    for key in pred_keys:
        vals = np.array(all_results[key])
        if len(vals) == 0:
            continue
        pcts = np.percentile(vals, percentiles)

        if key == 'combo_rapa_acarb':
            display = 'Rapa+Acarb combo (M)'
            baseline = base_combo
        else:
            cname = key[2:]  # strip "F_"
            display = f"F_{cname}"
            baseline = base_f_preds[cname]

        ci90 = f"[{pcts[0]:+.1f}, {pcts[4]:+.1f}]"
        print(f"{display:<28s}  {baseline:>+7.1f}%  "
              + "  ".join(f"{p:>+6.1f}%" for p in pcts)
              + f"  {ci90:>14s}")

        summary[key] = {
            'baseline': round(baseline, 2),
            'p5': round(float(pcts[0]), 2),
            'p25': round(float(pcts[1]), 2),
            'p50': round(float(pcts[2]), 2),
            'p75': round(float(pcts[3]), 2),
            'p95': round(float(pcts[4]), 2),
            'n_valid': n_valid,
        }

    print()

    # Compact ± format for paper
    print("Paper format (median [90% CI]):")
    for key in pred_keys:
        if key not in summary:
            continue
        s = summary[key]
        if key == 'combo_rapa_acarb':
            display = 'Rapa+Acarb combo'
        else:
            display = key[2:]
        lo = s['p5'] - s['p50']
        hi = s['p95'] - s['p50']
        print(f"  {display:<25s}  {s['p50']:+.1f}% "
              f"[{lo:+.1f}, {hi:+.1f}]")

    # Write JSON results
    out_dir = Path(__file__).resolve().parent.parent / 'docs'
    json_path = out_dir / 'uq_results.json'
    with open(json_path, 'w') as f:
        json.dump({
            'n_samples': n_samples,
            'n_valid': n_valid,
            'seed': seed,
            'parameters': [
                {'subsystem': sub, 'key': key, 'lo': lo, 'hi': hi,
                 'current': originals[i]}
                for i, (sub, key, lo, hi) in enumerate(TOP15_PARAMS)
            ],
            'predictions': summary,
            'raw': {k: [round(v, 3) for v in vals]
                    for k, vals in all_results.items()},
        }, f, indent=2)
    print(f"\nResults written to {json_path}")

    return summary


# ---------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Uncertainty quantification for ELM predictions')
    parser.add_argument('--samples', type=int, default=500,
                        help='Number of LHS samples (default: 500)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print progress per sample')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for LHS (default: 42)')
    args = parser.parse_args()

    run_uq(n_samples=args.samples, verbose=args.verbose, seed=args.seed)
