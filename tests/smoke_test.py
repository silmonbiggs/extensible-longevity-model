#!/usr/bin/env python3
"""
Smoke test — verifies the model runs and reproduces reported predictions.

Run: python smoke_test.py
Pass: all 12 ITP outputs within tolerance of targets.

Male predictions are calibrated by root-finding (tight tolerance).
Female predictions are genuine out-of-sample: the model uses four
sex-specific mechanisms from prior literature with no female fitting.
The paper reports 1.3 pp mean absolute error for females with the
pure-AMPK model (0.65 pp with dual-action canagliflozin extension).
The wider female tolerance verifies that the code reproduces the
paper's reported values, not that the model is perfect.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from elm.model import simulate, run_control, calculate_lifespan_extension
from elm.compounds import ITP_VALIDATION

TOLERANCE_MALE = 1.0    # pp — calibrated by root-finding, should be near-exact
TOLERANCE_FEMALE = 5.0  # pp — out-of-sample predictions, not fitted

def main():
    results = {}
    failures = []
    
    for sex in ['M', 'F']:
        control = run_control(sex=sex, t_max=1.5, dt=0.002)
        
        for cname, val in ITP_VALIDATION.items():
            result = simulate(compound=cname, sex=sex, t_max=1.5, dt=0.002)
            ext = calculate_lifespan_extension(result, control)
            
            target = val.target_male if sex == 'M' else val.target_female
            tol = TOLERANCE_MALE if sex == 'M' else TOLERANCE_FEMALE
            error = abs(ext - target)
            status = "PASS" if error <= tol else "FAIL"
            
            key = f"{cname}_{sex}"
            results[key] = {
                'extension': ext,
                'target': target,
                'error': error,
                'status': status,
            }
            
            if status == "FAIL":
                failures.append(key)
    
    # Report
    print("=" * 65)
    print(f"  ELM rev_00 Smoke Test")
    print("=" * 65)
    print(f"  {'Compound':30s} {'Got':>7s} {'Target':>7s} {'Error':>7s}  Status")
    print("-" * 65)
    
    for key, r in sorted(results.items()):
        marker = "+" if r['status'] == "PASS" else "X"
        print(f"  {key:30s} {r['extension']:6.1f}% {r['target']:6.1f}% "
              f"{r['error']:6.2f}%  {marker} {r['status']}")
    
    print("-" * 65)
    
    if failures:
        print(f"  FAILED: {', '.join(failures)}")
        sys.exit(1)
    else:
        print(f"  All {len(results)} targets within tolerance "
              f"(M: ±{TOLERANCE_MALE}pp, F: ±{TOLERANCE_FEMALE}pp). PASS.")
        sys.exit(0)


if __name__ == '__main__':
    main()
