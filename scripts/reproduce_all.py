#!/usr/bin/env python3
"""
Reproduce all paper results in one command.

Runs the smoke test, regenerates all figures, and runs the OAT
sensitivity analysis. Exit code is non-zero if any step fails.

Usage:
    python scripts/reproduce_all.py
"""

import subprocess
import sys
import time

STEPS = [
    ("Smoke test (12 ITP targets)", [sys.executable, "tests/smoke_test.py"]),
    ("Figure generation", [sys.executable, "scripts/generate_figures.py"]),
    ("OAT sensitivity analysis", [sys.executable, "scripts/oat_sensitivity.py"]),
]


def main():
    print("=" * 65)
    print("  ELM — Reproduce All Paper Results")
    print("=" * 65)

    t0 = time.time()
    failed = []

    for name, cmd in STEPS:
        print(f"\n{'─' * 65}")
        print(f"  Running: {name}")
        print(f"{'─' * 65}\n")

        result = subprocess.run(cmd)
        if result.returncode != 0:
            failed.append(name)
            print(f"\n  ✗ FAILED: {name}")
        else:
            print(f"\n  ✓ Done: {name}")

    elapsed = time.time() - t0
    print(f"\n{'=' * 65}")
    if failed:
        print(f"  FAILED: {', '.join(failed)}")
        print(f"  Elapsed: {elapsed:.0f}s")
        sys.exit(1)
    else:
        print(f"  All steps passed. Elapsed: {elapsed:.0f}s")
        print(f"  Outputs: docs/figures/*.png, oat_sensitivity.tsv")
    print(f"{'=' * 65}")


if __name__ == '__main__':
    main()
