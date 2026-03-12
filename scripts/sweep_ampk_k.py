#!/usr/bin/env python3
"""
Sensitivity sweep of the AMPK saturation curvature parameter k.

This parameter controls how much higher female baseline AMPK attenuates
the effect of AMPK-activating compounds (acarbose, canagliflozin, glycine).
The attenuation factor is f = 1 / (1 + k * (baseline - 1)), with baseline=1.3.

The sweep varies k across its plausible range [3, 20] and plots:
  Panel A: Female prediction (%) vs observed for each AMPK compound
  Panel B: Mean absolute female error (all 6 compounds) vs k

Male predictions are unaffected (k only enters the female sex modifier).

Usage:
    python scripts/sweep_ampk_k.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from contextlib import contextmanager

from elm.model import simulate, run_control, calculate_lifespan_extension
from elm.compounds import ITP_VALIDATION
from elm.sex_mechanisms import SEX_MECHANISMS

FIGURES_DIR = Path(os.path.dirname(__file__)) / '..' / 'docs' / 'figures'
DT = 0.002
T_MAX = 2.5
DPI = 150

AMPK_COMPOUNDS = ['acarbose', 'canagliflozin', 'glycine']
ALL_COMPOUNDS = list(ITP_VALIDATION.keys())

# Plausible range from the audit: predecessor compound-specific values were
# k in {5, 11, 15}, so the range spans roughly 3 to 20.
K_VALUES = np.array([3, 4, 5, 6, 7, 7.5, 8, 9, 10, 12, 15, 20])
K_CURRENT = 8.0


@contextmanager
def set_ampk_k(k_value):
    """Temporarily set AMPK saturation curvature for all AMPK compounds."""
    orig = {c: SEX_MECHANISMS[c].ampk_k_diminishing for c in AMPK_COMPOUNDS}
    try:
        for c in AMPK_COMPOUNDS:
            SEX_MECHANISMS[c].ampk_k_diminishing = k_value
        yield
    finally:
        for c, v in orig.items():
            SEX_MECHANISMS[c].ampk_k_diminishing = v


def run_sweep():
    """Run the full k-sweep and return results."""
    # Cache controls (independent of k)
    ctrl_m = run_control(sex='M', t_max=T_MAX, dt=DT)
    ctrl_f = run_control(sex='F', t_max=T_MAX, dt=DT)

    # Cache male predictions (independent of k)
    male_preds = {}
    for cn in ALL_COMPOUNDS:
        r = simulate(compound=cn, sex='M', t_max=T_MAX, dt=DT)
        male_preds[cn] = calculate_lifespan_extension(r, ctrl_m)

    # Cache female predictions for non-AMPK compounds (independent of k)
    fixed_female_preds = {}
    for cn in ALL_COMPOUNDS:
        if cn not in AMPK_COMPOUNDS:
            r = simulate(compound=cn, sex='F', t_max=T_MAX, dt=DT)
            fixed_female_preds[cn] = calculate_lifespan_extension(r, ctrl_f)

    # Sweep k
    results = {cn: [] for cn in AMPK_COMPOUNDS}
    mean_female_errors = []
    max_female_errors = []

    for k_val in K_VALUES:
        with set_ampk_k(k_val):
            # Female predictions for AMPK compounds at this k
            ampk_female = {}
            for cn in AMPK_COMPOUNDS:
                r = simulate(compound=cn, sex='F', t_max=T_MAX, dt=DT)
                pred = calculate_lifespan_extension(r, ctrl_f)
                ampk_female[cn] = pred
                results[cn].append(pred)

            # Compute female errors across all 6 compounds
            female_errors = []
            for cn in ALL_COMPOUNDS:
                target_f = ITP_VALIDATION[cn].target_female
                if cn in AMPK_COMPOUNDS:
                    pred_f = ampk_female[cn]
                else:
                    pred_f = fixed_female_preds[cn]
                female_errors.append(abs(pred_f - target_f))

            mean_female_errors.append(np.mean(female_errors))
            max_female_errors.append(np.max(female_errors))

    # Convert to arrays
    for cn in AMPK_COMPOUNDS:
        results[cn] = np.array(results[cn])

    return results, np.array(mean_female_errors), np.array(max_female_errors)


def plot_sweep(results, mean_errors, max_errors):
    """Create the 2-panel sweep figure."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Observed female targets
    targets = {cn: ITP_VALIDATION[cn].target_female for cn in AMPK_COMPOUNDS}
    colors = {'acarbose': '#2196F3', 'canagliflozin': '#FF9800', 'glycine': '#4CAF50'}
    labels = {'acarbose': 'Acarbose', 'canagliflozin': 'Canagliflozin', 'glycine': 'Glycine'}

    # --- Panel A: Female predictions vs k ---
    for cn in AMPK_COMPOUNDS:
        ax1.plot(K_VALUES, results[cn], '-o', color=colors[cn],
                 label=f'{labels[cn]} (predicted)', linewidth=2, markersize=5)
        ax1.axhline(y=targets[cn], color=colors[cn], linestyle='--', alpha=0.5,
                    label=f'{labels[cn]} (observed: {targets[cn]:.0f}%)')

    ax1.axvline(x=K_CURRENT, color='gray', linestyle=':', alpha=0.7, linewidth=1.5,
                label=f'Current k = {K_CURRENT}')
    ax1.set_xlabel('AMPK saturation curvature k', fontsize=12)
    ax1.set_ylabel('Female lifespan extension (%)', fontsize=12)
    ax1.set_title('A. Female predictions across plausible k', fontsize=13)
    ax1.legend(fontsize=8, loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(2, 21)

    # --- Panel B: Mean female error vs k ---
    ax2.plot(K_VALUES, mean_errors, 'k-s', linewidth=2, markersize=7,
             label='Mean |error| (all 6 compounds)')
    ax2.plot(K_VALUES, max_errors, 'r-^', linewidth=1.5, markersize=6, alpha=0.7,
             label='Max |error| (worst compound)')

    # Mark current k
    idx_current = np.argmin(np.abs(K_VALUES - K_CURRENT))
    ax2.plot(K_CURRENT, mean_errors[idx_current], 'ks', markersize=12, markerfacecolor='gold',
             markeredgewidth=2, zorder=5, label=f'k = {K_CURRENT}: {mean_errors[idx_current]:.1f} pp')

    ax2.axvline(x=K_CURRENT, color='gray', linestyle=':', alpha=0.7, linewidth=1.5)
    ax2.set_xlabel('AMPK saturation curvature k', fontsize=12)
    ax2.set_ylabel('Female prediction error (pp)', fontsize=12)
    ax2.set_title('B. Model error across plausible k', fontsize=13)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(2, 21)

    fig.tight_layout()
    return fig


def main():
    print('Running AMPK curvature k sweep...')
    print(f'  k values: {K_VALUES.tolist()}')
    print(f'  Compounds affected: {AMPK_COMPOUNDS}')

    results, mean_errors, max_errors = run_sweep()

    print('\nResults:')
    print(f'  {"k":>5}  {"Acarb F%":>9}  {"Cana F%":>9}  {"Gly F%":>9}  {"Mean err":>9}  {"Max err":>9}')
    for i, k_val in enumerate(K_VALUES):
        marker = ' <--' if k_val == K_CURRENT else ''
        print(f'  {k_val:5.1f}  {results["acarbose"][i]:9.1f}  '
              f'{results["canagliflozin"][i]:9.1f}  {results["glycine"][i]:9.1f}  '
              f'{mean_errors[i]:9.1f}  {max_errors[i]:9.1f}{marker}')

    fig = plot_sweep(results, mean_errors, max_errors)
    outpath = FIGURES_DIR / 'k_sweep.png'
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f'\nFigure saved: {outpath}')


if __name__ == '__main__':
    main()
