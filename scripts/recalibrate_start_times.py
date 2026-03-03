#!/usr/bin/env python3
"""
Recalibrate compound pathway injections for correct ITP start times.

Each ITP compound was previously calibrated assuming start_time=0.56.
Now that we use the actual ITP protocol start times, the injection
values need to be re-tuned to match the observed lifespan extensions.

Strategy: For each compound, uniformly scale all pathway injection values
by a single multiplier until the male extension matches the ITP target.
Then verify the female prediction (sex mechanisms should still work).
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from scipy.optimize import brentq

from elm.model import simulate, run_control, calculate_lifespan_extension
from elm.compounds import ITP_VALIDATION, COMPOUNDS, get_compound, get_itp_start_time
from elm.sex_mechanisms import apply_sex_modifier


DT = 0.002
T_MAX = 2.5


def run_scaled(compound, scale, sex='M'):
    """Run simulation with uniformly scaled pathway injections."""
    base = get_compound(compound)
    scaled = {k: v * scale for k, v in base.items()}
    scaled = apply_sex_modifier(compound, scaled, sex)
    scaled['start_time'] = get_itp_start_time(compound)
    return simulate(interventions=scaled, sex=sex, t_max=T_MAX, dt=DT)


def extension_error(scale, compound, target, sex, control):
    """Returns (model_extension - target) for a given scale factor."""
    result = run_scaled(compound, scale, sex)
    ext = calculate_lifespan_extension(result, control)
    return ext - target


def find_scale(compound, target, sex='M'):
    """Find scale factor that matches target extension."""
    control = run_control(sex=sex, t_max=T_MAX, dt=DT)

    # Check if target is 0 (male-only compounds for females)
    if target == 0:
        return 1.0  # Scale doesn't matter for zero-effect compounds

    # Binary search
    try:
        scale = brentq(extension_error, 0.01, 5.0,
                       args=(compound, target, sex, control),
                       xtol=0.001)
    except ValueError:
        # If brentq fails, try a wider range
        print(f"  WARNING: brentq failed for {compound} {sex}, trying wider range")
        try:
            scale = brentq(extension_error, 0.001, 20.0,
                           args=(compound, target, sex, control),
                           xtol=0.001)
        except ValueError:
            print(f"  ERROR: Could not find scale for {compound} {sex}")
            return 1.0

    return scale


def main():
    print("=" * 70)
    print("  RECALIBRATING COMPOUNDS FOR ITP START TIMES")
    print("=" * 70)

    compounds = ['rapamycin', 'acarbose', 'canagliflozin',
                 '17_alpha_estradiol', 'aspirin', 'glycine']

    results = {}

    for cname in compounds:
        val = ITP_VALIDATION[cname]
        start_time = get_itp_start_time(cname)
        print(f"\n--- {val.name} (start={start_time:.2f}, was 0.56) ---")
        print(f"  Targets: M={val.target_male}%, F={val.target_female}%")
        print(f"  Current injections: {COMPOUNDS[cname]}")

        # Calibrate to male target
        scale = find_scale(cname, val.target_male, 'M')
        print(f"  Scale factor (to match male): {scale:.4f}")

        # Compute new injection values
        base = get_compound(cname)
        new_injections = {k: round(v * scale, 4) for k, v in base.items()}
        print(f"  New injections: {new_injections}")

        # Verify male
        control_m = run_control(sex='M', t_max=T_MAX, dt=DT)
        scaled_m = {k: v * scale for k, v in base.items()}
        scaled_m = apply_sex_modifier(cname, scaled_m, 'M')
        scaled_m['start_time'] = start_time
        result_m = simulate(interventions=scaled_m, sex='M', t_max=T_MAX, dt=DT)
        ext_m = calculate_lifespan_extension(result_m, control_m)

        # Verify female
        control_f = run_control(sex='F', t_max=T_MAX, dt=DT)
        scaled_f = {k: v * scale for k, v in base.items()}
        scaled_f = apply_sex_modifier(cname, scaled_f, 'F')
        scaled_f['start_time'] = start_time
        result_f = simulate(interventions=scaled_f, sex='F', t_max=T_MAX, dt=DT)
        ext_f = calculate_lifespan_extension(result_f, control_f)

        print(f"  Verification: M={ext_m:.1f}% (target {val.target_male}%), "
              f"F={ext_f:.1f}% (target {val.target_female}%)")

        results[cname] = {
            'scale': scale,
            'new_injections': new_injections,
            'ext_m': ext_m,
            'ext_f': ext_f,
        }

    # Summary
    print("\n" + "=" * 70)
    print("  SUMMARY: New COMPOUNDS dict entries")
    print("=" * 70)
    for cname in compounds:
        r = results[cname]
        val = ITP_VALIDATION[cname]
        err_m = abs(r['ext_m'] - val.target_male)
        err_f = abs(r['ext_f'] - val.target_female)
        print(f"\n    '{cname}': {{")
        for k, v in r['new_injections'].items():
            print(f"        '{k}': {v},")
        print(f"    }},")
        print(f"    # M: {r['ext_m']:.1f}% (target {val.target_male}%, err={err_m:.2f}%) "
              f"F: {r['ext_f']:.1f}% (target {val.target_female}%, err={err_f:.2f}%)")


if __name__ == '__main__':
    main()
