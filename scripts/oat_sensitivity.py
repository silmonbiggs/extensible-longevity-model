#!/usr/bin/env python3
"""
One-at-a-Time (OAT) Sensitivity Analysis for ELM Frozen Priors

For each frozen parameter:
  1. Perturb by ±10%
  2. Re-calibrate: root-find the 6 ITP compound scale factors so male
     predictions still match observed ITP targets exactly
  3. Record shift in all 6 female predictions and rapa+acarb combo prediction
  4. Report max |shift| across the 7 out-of-sample predictions

This answers: "Which frozen priors most affect the paper's central claims
(female predictions and combination prediction)?"

Usage:
    python scripts/oat_sensitivity.py              # full sweep, TSV output
    python scripts/oat_sensitivity.py --top 20     # show only top 20
    python scripts/oat_sensitivity.py --verbose     # print progress per param
"""

import sys
import os
import time
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from scipy.optimize import brentq

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
DT = 0.01            # coarse dt for speed (~2x faster than DT_FAST)
PERTURB = 0.10       # ±10%

ITP_COMPOUNDS = ['rapamycin', 'acarbose', 'canagliflozin',
                 '17_alpha_estradiol', 'aspirin', 'glycine']

# Pathway ratios for each ITP compound (extracted from COMPOUNDS dict)
# ratio_i = injection_i / sum(injections)
COMPOUND_RATIOS: Dict[str, Dict[str, float]] = {}
for cname in ITP_COMPOUNDS:
    inj = COMPOUNDS[cname]
    total = sum(inj.values())
    COMPOUND_RATIOS[cname] = {k: v / total for k, v in inj.items()}

# Male targets
MALE_TARGETS = {cname: ITP_VALIDATION[cname].target_male
                for cname in ITP_COMPOUNDS}

# Parameters to skip (non-numeric, structural, or normalization)
SKIP_KEYS = {
    # Normalization constants (derived from control, not independent)
    'meth_norm', 'damage_norm', 'H_norm', 'sen_norm',
    'death_threshold',
    # BioAge weights (swept separately in weight sensitivity analysis)
    'w_meth', 'w_damage', 'w_hetero', 'w_sen',
    # Bounds and floors that are 0.0 (can't perturb multiplicatively)
    # — handled by the numeric filter below
}


# ---------------------------------------------------------------------------
# Core functions
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
    """Root-find scale factor so male prediction matches target.

    alpha_hint seeds the bracket near the expected value for speed.
    """
    def objective(alpha):
        return male_extension(compound, alpha, control_m) - target

    # Tight bracket around hint (for ±10% parameter perturbations,
    # optimal alpha changes by much less than ±50%)
    if alpha_hint is not None and alpha_hint > 0:
        alpha_lo = alpha_hint * 0.3
        alpha_hi = alpha_hint * 3.0
    else:
        alpha_lo = 0.01
        alpha_hi = 3.0

    # Check bracket
    f_lo = objective(alpha_lo)
    f_hi = objective(alpha_hi)
    if f_lo * f_hi > 0:
        # Widen
        alpha_lo = 0.01
        alpha_hi = 5.0
        f_lo = objective(alpha_lo)
        f_hi = objective(alpha_hi)
        if f_lo * f_hi > 0:
            # Fallback
            total = sum(COMPOUNDS[compound].values())
            return total

    return brentq(objective, alpha_lo, alpha_hi, xtol=5e-3)


# Baseline alphas — filled by run_baseline(), used to seed re-calibration
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

    # Merge (sum shared pathways)
    combo = {}
    for k, v in rapa_inj_sex.items():
        if k != 'start_time':
            combo[k] = v
    for k, v in acarb_inj_sex.items():
        if k != 'start_time':
            combo[k] = combo.get(k, 0) + v
    combo['start_time'] = 9.0 / 30.0  # 9-month start

    result = simulate(interventions=combo, sex='M', t_max=T_MAX, dt=DT)
    return calculate_lifespan_extension(result, control_m)


def run_baseline():
    """Compute baseline: calibrate and predict with unperturbed parameters."""
    control_m = run_control(sex='M', t_max=T_MAX, dt=DT)
    control_f = run_control(sex='F', t_max=T_MAX, dt=DT)
    alphas = calibrate_all(control_m)
    # Cache baseline alphas for seeding perturbed re-calibrations
    _BASELINE_ALPHAS.update(alphas)
    f_preds = female_predictions(alphas, control_f)
    combo = combo_prediction(alphas, control_m)
    return control_m, control_f, alphas, f_preds, combo


def collect_parameters() -> List[Tuple[dict, str, str, float]]:
    """Collect all numeric, perturbable parameters from ALL_PATHWAY_PARAMS.

    Returns list of (param_dict, key, display_name, current_value).
    """
    params = []
    for system_name, pdict in ALL_PATHWAY_PARAMS.items():
        for key, value in pdict.items():
            # Skip non-numeric
            if not isinstance(value, (int, float)):
                continue
            # Skip zeros (can't perturb multiplicatively)
            if value == 0.0:
                continue
            # Skip explicitly excluded keys
            if key in SKIP_KEYS:
                continue
            display = f"{system_name}.{key}"
            params.append((pdict, key, display, float(value)))
    return params


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def _raw_screen(all_params, base_lifespan_m, verbose=False):
    """Pass 1: fast control-only sensitivity screen.

    For each param, perturb ±10%, run male control only, measure lifespan
    shift.  This is a proxy: if the control doesn't move, the parameter is
    insensitive (compound-specific pathways are always also active in
    control since they represent endogenous biology).  Returns list sorted
    by max |lifespan shift|.
    """
    results = []
    for i, (pdict, key, display, orig_value) in enumerate(all_params):
        if verbose and (i + 1) % 20 == 0:
            print(f"    Screening [{i+1}/{len(all_params)}]...")

        max_shift = 0.0
        for factor in [1 - PERTURB, 1 + PERTURB]:
            pdict[key] = orig_value * factor
            try:
                cm = run_control(sex='M', t_max=T_MAX, dt=DT)
                shift = abs(cm.t_death - base_lifespan_m)
                if shift > max_shift:
                    max_shift = shift
            except Exception:
                pass
            finally:
                pdict[key] = orig_value

        results.append((pdict, key, display, orig_value, max_shift))

    results.sort(key=lambda r: r[4], reverse=True)
    return results


def run_oat(verbose: bool = False, top_n: int = 0,
            recal_cutoff: float = 0.001, max_recal: int = 60):
    """Run two-pass OAT sweep.

    Pass 1: Fast raw screen (no re-calibration) for all parameters.
    Pass 2: Re-calibrated sensitivity for top parameters (above recal_cutoff
            or up to max_recal, whichever is larger).
    """

    print("OAT Sensitivity Analysis for ELM Frozen Priors")
    print(f"Perturbation: ±{PERTURB*100:.0f}%")
    print()

    t_start = time.time()

    # Collect parameters
    all_params = collect_parameters()
    print(f"Found {len(all_params)} perturbable parameters across "
          f"{len(ALL_PATHWAY_PARAMS)} subsystems")
    print()

    # Baseline
    print("Computing baseline (unperturbed)...")
    t0 = time.time()
    ctrl_m, ctrl_f, base_alphas, base_f_preds, base_combo = run_baseline()
    t_base = time.time() - t0

    base_lifespan_m = ctrl_m.t_death
    print(f"  Baseline computed in {t_base:.1f}s")
    print("  Baseline female predictions:")
    for cname in ITP_COMPOUNDS:
        target = ITP_VALIDATION[cname].target_female
        pred = base_f_preds[cname]
        print(f"    {cname:25s}  pred={pred:+.1f}%  target={target:+.1f}%")
    print(f"  Rapa+acarb combo: {base_combo:+.1f}% (observed: +34.0%)")
    print()

    # --- Pass 1: Fast control-only screen ---
    print("Pass 1: Control-lifespan screen (2 sims per parameter)...")
    t1 = time.time()
    raw_results = _raw_screen(all_params, base_lifespan_m, verbose=verbose)
    print(f"  Screened {len(raw_results)} parameters in {time.time()-t1:.0f}s")

    # Determine which params need re-calibration
    n_above = sum(1 for r in raw_results if r[4] >= recal_cutoff)
    n_recal = max(n_above, min(max_recal, len(raw_results)))
    recal_params = raw_results[:n_recal]
    print(f"  {n_above} params above {recal_cutoff:.4f} lifespan-shift threshold")
    print(f"  Re-calibrating top {n_recal} parameters...")
    print()

    # --- Pass 2: Re-calibrated sensitivity ---
    print("Pass 2: Re-calibrated sensitivity (root-finding scale factors)...")
    t2 = time.time()
    results = []

    for i, (pdict, key, display, orig_value, raw_shift) in enumerate(recal_params):
        if verbose or (i + 1) % 10 == 0:
            eta = (time.time() - t2) / max(i, 1) * (n_recal - i) if i > 0 else 0
            print(f"  [{i+1}/{n_recal}] {display} = {orig_value:.4g} "
                  f"(raw: {raw_shift:.6f}) ETA: {eta:.0f}s")

        shifts_lo = {}
        shifts_hi = {}

        for factor, shifts in [(1 - PERTURB, shifts_lo),
                               (1 + PERTURB, shifts_hi)]:
            pdict[key] = orig_value * factor

            try:
                cm = run_control(sex='M', t_max=T_MAX, dt=DT)
                cf = run_control(sex='F', t_max=T_MAX, dt=DT)
                alphas = calibrate_all(cm)
                f_preds = female_predictions(alphas, cf)
                combo = combo_prediction(alphas, cm)

                for cname in ITP_COMPOUNDS:
                    shifts[f"F_{cname}"] = f_preds[cname] - base_f_preds[cname]
                shifts["combo_rapa_acarb"] = combo - base_combo

            except Exception as e:
                if verbose:
                    print(f"    WARNING: {display} at {factor:.2f}x: {e}")
                for cname in ITP_COMPOUNDS:
                    shifts[f"F_{cname}"] = 0.0
                shifts["combo_rapa_acarb"] = 0.0

            finally:
                pdict[key] = orig_value

        all_shifts = list(shifts_lo.values()) + list(shifts_hi.values())
        max_shift = max(abs(s) for s in all_shifts) if all_shifts else 0.0

        per_pred_max = {}
        for pred_key in shifts_lo:
            per_pred_max[pred_key] = max(abs(shifts_lo[pred_key]),
                                         abs(shifts_hi[pred_key]))

        results.append((display, max_shift, per_pred_max, orig_value))

    # Add insensitive params (those below cutoff, not re-calibrated)
    # Convert raw lifespan shift to approximate pp for consistent display
    for pdict, key, display, orig_value, raw_shift in raw_results[n_recal:]:
        approx_pp = raw_shift / base_lifespan_m * 100.0 if base_lifespan_m > 0 else 0.0
        results.append((display, approx_pp, {}, orig_value))

    elapsed = time.time() - t_start
    print(f"\nSweep completed in {elapsed:.1f}s "
          f"({elapsed/len(all_params):.2f}s per parameter)")

    # Sort by max shift (descending)
    results.sort(key=lambda r: r[1], reverse=True)

    # Print results
    print()
    print("=" * 90)
    print(f"{'Parameter':<45s}  {'Value':>8s}  {'Max |shift|':>11s}  "
          f"{'Label':>20s}")
    print("=" * 90)

    n_show = top_n if top_n > 0 else len(results)
    for display, max_shift, per_pred, value in results[:n_show]:
        if max_shift >= 1.0:
            label = f"±{max_shift:.1f} pp"
        elif max_shift >= 0.1:
            label = f"±{max_shift:.2f} pp"
        else:
            label = "insensitive"
        print(f"{display:<45s}  {value:>8.4g}  {max_shift:>10.3f} pp  "
              f"{label:>20s}")

    # Summary statistics
    n_sensitive = sum(1 for _, ms, _, _ in results if ms >= 1.0)
    n_moderate = sum(1 for _, ms, _, _ in results if 0.1 <= ms < 1.0)
    n_insensitive = sum(1 for _, ms, _, _ in results if ms < 0.1)
    print()
    print(f"Summary: {n_sensitive} sensitive (>=1 pp), "
          f"{n_moderate} moderate (0.1-1 pp), "
          f"{n_insensitive} insensitive (<0.1 pp)")
    print(f"Total: {len(results)} parameters swept")

    # Write TSV
    tsv_path = Path(__file__).resolve().parent.parent / 'docs' / 'oat_sensitivity.tsv'
    with open(tsv_path, 'w') as f:
        headers = ['parameter', 'value', 'max_shift_pp']
        headers += [f"F_{c}" for c in ITP_COMPOUNDS]
        headers += ['combo_rapa_acarb']
        f.write('\t'.join(headers) + '\n')

        for display, max_shift, per_pred, value in results:
            row = [display, f"{value:.6g}", f"{max_shift:.4f}"]
            for cname in ITP_COMPOUNDS:
                key = f"F_{cname}"
                row.append(f"{per_pred.get(key, 0):.4f}")
            row.append(f"{per_pred.get('combo_rapa_acarb', 0):.4f}")
            f.write('\t'.join(row) + '\n')

    print(f"\nResults written to {tsv_path}")
    return results


# ---------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='OAT sensitivity analysis for ELM frozen priors')
    parser.add_argument('--top', type=int, default=0,
                        help='Show only top N most sensitive parameters')
    parser.add_argument('--verbose', action='store_true',
                        help='Print progress per parameter')
    args = parser.parse_args()

    run_oat(verbose=args.verbose, top_n=args.top)
