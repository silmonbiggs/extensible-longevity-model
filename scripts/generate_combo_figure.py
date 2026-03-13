#!/usr/bin/env python3
"""
Generate rapa+acarbose combination figure (journal style, with error bars).

Usage:
    python scripts/generate_combo_figure.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from elm.model import simulate, run_control, calculate_lifespan_extension, run_combination_extension
from elm.compounds import ITP_VALIDATION

FIGURES_DIR = Path(os.path.dirname(__file__)) / '..' / 'docs' / 'figures'
DT = 0.002
T_MAX = 2.5
DPI = 300

# ITP combo target (Strong 2022): +34% male, 9-month start
COMBO_TARGET = 34.0
COMBO_ITP_SE = 2.5   # estimated, same basis as individual compounds
COMBO_MODEL_SE = 1.9  # from OAT RSS (see oat_sensitivity.tsv)


def main():
    ctrl_m = run_control(sex='M', t_max=T_MAX, dt=DT)

    # Individual predictions
    rapa_r = simulate(compound='rapamycin', sex='M', t_max=T_MAX, dt=DT)
    rapa_ext = calculate_lifespan_extension(rapa_r, ctrl_m)

    acarb_r = simulate(compound='acarbose', sex='M', t_max=T_MAX, dt=DT)
    acarb_ext = calculate_lifespan_extension(acarb_r, ctrl_m)

    naive_sum = rapa_ext + acarb_ext

    # Combination prediction (9-month start)
    combo_ext = run_combination_extension(
        {'rapamycin': 1.0, 'acarbose': 1.0},
        sex='M', t_max=T_MAX, dt=DT, start_time=9.0/30.0)

    print(f'Rapa alone:       {rapa_ext:.1f}%')
    print(f'Acarbose alone:   {acarb_ext:.1f}%')
    print(f'Naive sum:        {naive_sum:.1f}%')
    print(f'ELM combination:  {combo_ext:.1f}%')
    print(f'ITP observed:     {COMBO_TARGET:.1f}%')

    # --- Figure ---
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
    })

    fig, ax = plt.subplots(figsize=(4.5, 4.0))

    labels = ['Rapa\nalone', 'Acarb\nalone', 'Naive\nsum', 'ELM', 'ITP\nobserved']
    values = [rapa_ext, acarb_ext, naive_sum, combo_ext, COMBO_TARGET]
    errors = [0, 0, 0, COMBO_MODEL_SE, COMBO_ITP_SE]

    # Colors: ITP data = saturated blue (male measurement), naive = gray,
    # ELM = light blue (modeled male). Consistent with sex-diff figure hue code.
    colors = ['#1565C0', '#1565C0', '#BDBDBD', '#90CAF9', '#1565C0']
    edge_colors = ['#0D47A1', '#0D47A1', '#757575', '#64B5F6', '#0D47A1']

    bw = 0.55
    x = np.arange(len(labels))

    bars = ax.bar(x, values, bw, color=colors, edgecolor=edge_colors,
                  linewidth=0.6, yerr=errors, error_kw=dict(
                      elinewidth=0.8, capsize=3, capthick=0.6, ecolor='#333333'),
                  zorder=2)

    # Observed line
    ax.axhline(y=COMBO_TARGET, color='#1565C0', linestyle='--', linewidth=0.7,
               alpha=0.4, zorder=1)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Male lifespan extension (%)')
    ax.set_ylim(0, 52)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out')

    fig.tight_layout()

    outpath = FIGURES_DIR / 'combo_rapa_acarb.png'
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f'\nFigure saved: {outpath}')


if __name__ == '__main__':
    main()
