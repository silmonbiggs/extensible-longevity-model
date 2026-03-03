#!/usr/bin/env python3
"""
Smoke test for rev_00 — verifies the model runs and calibration holds.

Run: python smoke_test.py
Pass: all 12 ITP outputs within tolerance of targets.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from elm.model import simulate, run_control, calculate_lifespan_extension
from elm.compounds import ITP_VALIDATION

TOLERANCE = 1.0  # percentage points

def main():
    results = {}
    failures = []
    
    for sex in ['M', 'F']:
        control = run_control(sex=sex, t_max=1.5, dt=0.002)
        
        for cname, val in ITP_VALIDATION.items():
            result = simulate(compound=cname, sex=sex, t_max=1.5, dt=0.002)
            ext = calculate_lifespan_extension(result, control)
            
            target = val.target_male if sex == 'M' else val.target_female
            error = abs(ext - target)
            status = "PASS" if error <= TOLERANCE else "FAIL"
            
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
        print(f"  All {len(results)} targets within ±{TOLERANCE}%. PASS.")
        sys.exit(0)


if __name__ == '__main__':
    main()
