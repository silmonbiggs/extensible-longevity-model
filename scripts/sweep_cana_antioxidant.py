#!/usr/bin/env python3
"""
Sweep of canagliflozin antioxidant fraction.

Canagliflozin is currently modeled as a pure AMPK activator. This sweep
tests a dual-action hypothesis: part AMPK activation (subject to female
AMPK saturation) plus part antioxidant/FGF21 (sex-independent).

At each fraction f, the base ratios are {'ampk': 1-f, 'antioxidant': f},
with a scale factor s root-found to match the male ITP target (+14%).
Female predictions then follow from the AMPK saturation sex mechanism,
which diminishes only the AMPK portion.

Two-phase adaptive sweep:
  Phase 1: Coarse grid (f in {0.0, 0.2, 0.4, 0.6, 0.8, 1.0})
  Phase 2: Fine grid (5 points around the Phase 1 optimum)

Output: 2-panel figure docs/figures/cana_antioxidant_sweep.png

Usage:
    python scripts/sweep_cana_antioxidant.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq

from elm.model import simulate, run_control, calculate_lifespan_extension
from elm.compounds import ITP_VALIDATION, get_itp_start_time
from elm.sex_mechanisms import apply_sex_modifier

FIGURES_DIR = Path(os.path.dirname(__file__)) / '..' / 'docs' / 'figures'
DT = 0.002
T_MAX = 2.5
DPI = 300

CANA_TARGET_M = ITP_VALIDATION['canagliflozin'].target_male   # 14.0%
CANA_TARGET_F = ITP_VALIDATION['canagliflozin'].target_female  # 9.0%
CANA_START = get_itp_start_time('canagliflozin')

ALL_COMPOUNDS = list(ITP_VALIDATION.keys())


def run_cana(f, s, sex, ctrl):
    """Run canagliflozin simulation with antioxidant fraction f and scale s."""
    if f <= 0:
        base = {'ampk': 1.0}
    elif f >= 1:
        base = {'antioxidant': 1.0}
    else:
        base = {'ampk': 1.0 - f, 'antioxidant': f}

    scaled = {k: v * s for k, v in base.items()}
    scaled = apply_sex_modifier('canagliflozin', scaled, sex)
    scaled['start_time'] = CANA_START
    result = simulate(interventions=scaled, sex=sex, t_max=T_MAX, dt=DT)
    return calculate_lifespan_extension(result, ctrl)


def find_male_scale(f, ctrl_m):
    """Find scale factor s so male canagliflozin extension = target."""
    def error(s):
        return run_cana(f, s, 'M', ctrl_m) - CANA_TARGET_M

    try:
        return brentq(error, 0.01, 5.0, xtol=0.001)
    except ValueError:
        try:
            return brentq(error, 0.001, 20.0, xtol=0.001)
        except ValueError:
            print(f"  WARNING: brentq failed for f={f:.2f}")
            return None


def run_sweep():
    """Run the two-phase adaptive sweep."""
    ctrl_m = run_control(sex='M', t_max=T_MAX, dt=DT)
    ctrl_f = run_control(sex='F', t_max=T_MAX, dt=DT)

    # Pre-compute fixed female errors for the other 5 compounds
    fixed_errors = []
    for cn in ALL_COMPOUNDS:
        if cn == 'canagliflozin':
            continue
        result = simulate(compound=cn, sex='F', t_max=T_MAX, dt=DT)
        pred_f = calculate_lifespan_extension(result, ctrl_f)
        target_f = ITP_VALIDATION[cn].target_female
        fixed_errors.append(abs(pred_f - target_f))

    sum_fixed = sum(fixed_errors)
    max_fixed = max(fixed_errors)

    def sweep_points(f_values):
        results = []
        for f in f_values:
            s = find_male_scale(f, ctrl_m)
            if s is None:
                continue
            ext_m = run_cana(f, s, 'M', ctrl_m)
            ext_f = run_cana(f, s, 'F', ctrl_f)
            cana_err = abs(ext_f - CANA_TARGET_F)
            mean_err = (sum_fixed + cana_err) / 6
            max_err = max(max_fixed, cana_err)
            results.append({
                'f': f, 's': s,
                'ext_m': ext_m, 'ext_f': ext_f,
                'cana_err': cana_err,
                'mean_err': mean_err, 'max_err': max_err,
            })
        return results

    # Phase 1: Coarse sweep
    print('Phase 1: Coarse sweep (6 points)...')
    coarse_f = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    coarse = sweep_points(coarse_f)

    best_coarse = min(coarse, key=lambda r: r['mean_err'])
    gamma_star = best_coarse['f']
    print(f'  Phase 1 optimum: f={gamma_star:.2f} '
          f'(mean err={best_coarse["mean_err"]:.2f} pp)')

    # Phase 2: Fine sweep (+-10 in steps of 5)
    print('Phase 2: Fine sweep (around f={:.2f})...'.format(gamma_star))
    fine_f = [gamma_star + d for d in [-0.10, -0.05, 0.0, 0.05, 0.10]]
    fine_f = [round(f, 2) for f in fine_f if 0.0 <= f <= 1.0]
    already = set(r['f'] for r in coarse)
    fine_f = [f for f in fine_f if f not in already]
    fine = sweep_points(fine_f)

    phase2_all = coarse + fine
    best_fine = min(phase2_all, key=lambda r: r['mean_err'])
    gamma2 = best_fine['f']
    print(f'  Phase 2 optimum: f={gamma2:.2f} '
          f'(mean err={best_fine["mean_err"]:.2f} pp)')

    # Phase 3: Superfine sweep (+-8 in steps of 2)
    print('Phase 3: Superfine sweep (around f={:.2f})...'.format(gamma2))
    superfine_f = [gamma2 + d / 100.0 for d in [-8, -6, -4, -2, 0, 2, 4, 6, 8]]
    superfine_f = [round(f, 2) for f in superfine_f if 0.0 <= f <= 1.0]
    already2 = set(r['f'] for r in phase2_all)
    superfine_f = [f for f in superfine_f if f not in already2]
    superfine = sweep_points(superfine_f)

    all_results = sorted(phase2_all + superfine, key=lambda r: r['f'])
    best = min(all_results, key=lambda r: r['mean_err'])

    return all_results, best, coarse_f


def plot_sweep(results, best, coarse_f):
    """Create the 2-panel sweep figure (Aging Cell style)."""
    # Journal style: Arial, 9 pt base, 300 DPI, minimal chrome
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'lines.linewidth': 1.2,
        'lines.markersize': 4,
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

    f_vals = [r['f'] for r in results]
    ext_f_vals = [r['ext_f'] for r in results]
    mean_errs = [r['mean_err'] for r in results]
    max_errs = [r['max_err'] for r in results]

    r0 = next(r for r in results if r['f'] == 0.0)
    f_star = best['f']

    # --- Panel A: Canagliflozin female prediction vs f ---
    ax1.plot(f_vals, ext_f_vals, '-o', color='#333333', markerfacecolor='#333333',
             markeredgecolor='#333333', zorder=2)
    ax1.axhline(y=CANA_TARGET_F, color='#CC0000', linestyle='--', linewidth=0.8,
                label=f'ITP observed ({CANA_TARGET_F:.0f}%)')
    ax1.axvline(x=f_star, color='#666666', linestyle=':', linewidth=0.8,
                label=f'$f^*$ = {f_star:.2f}')

    ax1.annotate(f'f = 0\n{r0["ext_f"]:.1f}%',
                 xy=(0.0, r0['ext_f']), xytext=(0.10, r0['ext_f'] + 0.8),
                 fontsize=7, ha='left',
                 arrowprops=dict(arrowstyle='-', color='#999999', lw=0.6))
    ax1.annotate(f'$f^*$ = {f_star:.2f}\n{best["ext_f"]:.1f}%',
                 xy=(f_star, best['ext_f']), xytext=(f_star - 0.22, best['ext_f'] + 1.8),
                 fontsize=7, ha='left',
                 arrowprops=dict(arrowstyle='-', color='#999999', lw=0.6))

    ax1.set_xlabel('Antioxidant fraction ($f$)')
    ax1.set_ylabel('Female lifespan extension (%)')
    ax1.legend(loc='upper left', frameon=False)
    ax1.set_xlim(-0.03, 1.03)
    ax1.set_ylim(0, None)
    ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes,
             fontsize=12, fontweight='bold', va='bottom')

    # --- Panel B: Mean and max female error vs f ---
    ax2.plot(f_vals, mean_errs, '-s', color='#333333', markerfacecolor='#333333',
             markeredgecolor='#333333', label='Mean |error| (6 compounds)')
    ax2.plot(f_vals, max_errs, '-^', color='#CC0000', markerfacecolor='white',
             markeredgecolor='#CC0000', label='Max |error| (worst compound)')
    ax2.axvline(x=f_star, color='#666666', linestyle=':', linewidth=0.8)

    # Noise-floor references
    itp_se = 2.5
    expected_mae = itp_se * np.sqrt(2 / np.pi)  # ~2.0 pp
    ax2.axhline(y=itp_se, color='#CC0000', linestyle=':', linewidth=0.6,
                alpha=0.7, label=f'ITP measurement SE ({itp_se} pp)')
    ax2.axhline(y=expected_mae, color='#333333', linestyle=':', linewidth=0.6,
                alpha=0.7, label=f'Expected MAE under noise ({expected_mae:.1f} pp)')

    ax2.set_xlabel('Antioxidant fraction ($f$)')
    ax2.set_ylabel('Female prediction error (pp)')
    ax2.legend(loc='upper right', frameon=False)
    ax2.set_xlim(-0.03, 1.03)
    ax2.set_ylim(0, 8)
    ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes,
             fontsize=12, fontweight='bold', va='bottom')

    for ax in (ax1, ax2):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(direction='out')

    fig.tight_layout(w_pad=2.5)
    return fig


def main():
    print('=' * 60)
    print('  CANAGLIFLOZIN ANTIOXIDANT FRACTION SWEEP')
    print('=' * 60)
    print(f'  Male target: {CANA_TARGET_M:.0f}%')
    print(f'  Female target: {CANA_TARGET_F:.0f}%')
    print()

    results, best, coarse_f = run_sweep()

    # Console table
    print(f'\n{"=" * 80}')
    print(f'  RESULTS (sorted by f)')
    print(f'{"=" * 80}')
    print(f'  {"f":>5}  {"scale":>8}  {"M ext%":>7}  {"F ext%":>7}  '
          f'{"F err":>6}  {"Mean err":>9}  {"Max err":>9}  {"Phase":>6}')
    print(f'  {"-" * 72}')
    for r in results:
        phase = 'coarse' if r['f'] in coarse_f else 'fine'
        marker = ' <-- optimal' if r['f'] == best['f'] else ''
        print(f'  {r["f"]:5.2f}  {r["s"]:8.4f}  {r["ext_m"]:7.1f}  {r["ext_f"]:7.1f}  '
              f'{r["cana_err"]:6.2f}  {r["mean_err"]:9.2f}  {r["max_err"]:9.2f}  '
              f'{phase:>6}{marker}')

    print(f'\n  Optimal: f* = {best["f"]:.2f}, scale = {best["s"]:.4f}')
    print(f'  At f*: M = {best["ext_m"]:.1f}%, F = {best["ext_f"]:.1f}%')
    print(f'  Mean female error at f*: {best["mean_err"]:.2f} pp')

    f_star = best['f']
    s_star = best['s']
    ampk_val = round(s_star * (1.0 - f_star), 4)
    antioxidant_val = round(s_star * f_star, 4)
    print(f'\n  New COMPOUNDS entry:')
    print(f"    'canagliflozin': {{")
    if ampk_val > 0:
        print(f"        'ampk': {ampk_val},")
    if antioxidant_val > 0:
        print(f"        'antioxidant': {antioxidant_val},")
    print(f"    }},")

    fig = plot_sweep(results, best, coarse_f)
    outpath = FIGURES_DIR / 'cana_antioxidant_sweep.png'
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f'\n  Figure saved: {outpath}')


if __name__ == '__main__':
    main()
