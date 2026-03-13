#!/usr/bin/env python3
"""
Generate rapamycin survival curves (journal style) with median-gap arrows.

Male panel: arrow showing control→treated median gap.
Female panel: arrow showing female gap, plus a ghost arrow below showing
the male gap for direct visual comparison.

Usage:
    python scripts/generate_rapa_survival.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from elm.model import simulate, run_control, calculate_lifespan_extension

FIGURES_DIR = Path(os.path.dirname(__file__)) / '..' / 'docs' / 'figures'
DT = 0.002
T_MAX = 2.5
DPI = 300


def draw_median_drop(ax, t_median, survival, t_array, color, alpha=0.5):
    """Draw a vertical drop line from the survival curve down to the 50% line."""
    ax.plot([t_median, t_median], [0, 50], color=color, linestyle=':',
            linewidth=0.7, alpha=alpha, zorder=1)


def main():
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2))

    medians = {}

    for ax, sex, panel_label in [(ax1, 'M', 'Male'), (ax2, 'F', 'Female')]:
        c = run_control(sex=sex, t_max=T_MAX, dt=DT)
        r = simulate(compound='rapamycin', sex=sex, t_max=T_MAX, dt=DT)
        ext = calculate_lifespan_extension(r, c)

        ax.plot(c.t, c.Survival * 100, 'k-', linewidth=1.4, label='Control')
        ax.plot(r.t, r.Survival * 100, '-', color='#CC0000', linewidth=1.4,
                label='Rapamycin')
        ax.axhline(50, color='#999999', linestyle='--', linewidth=0.5, alpha=0.5)

        # Vertical drop lines at each median
        draw_median_drop(ax, c.t_death, c.Survival, c.t, '#333333')
        draw_median_drop(ax, r.t_death, r.Survival, r.t, '#CC0000')

        ax.set_xlabel('Normalized age')
        if sex == 'M':
            ax.set_ylabel('Survival (%)')
        ax.set_xlim(0, 2.0)
        ax.set_ylim(0, 105)
        ax.legend(fontsize=7, loc='lower left', frameon=False)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(direction='out')

        # Panel letter + title
        ax.text(-0.12, 1.05, 'A' if sex == 'M' else 'B',
                transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')
        ax.set_title(panel_label, fontsize=10, fontweight='bold', pad=8)

        medians[sex] = {'ctrl': c.t_death, 'rapa': r.t_death, 'ext': ext}

    # --- Arrow on male panel ---
    mc = medians['M']['ctrl']
    mr = medians['M']['rapa']
    y_arrow = 50
    ax1.annotate('', xy=(mr, y_arrow), xytext=(mc, y_arrow),
                 arrowprops=dict(arrowstyle='->', color='#1565C0', lw=1.4,
                                 shrinkA=0, shrinkB=0))
    ax1.text((mc + mr) / 2, y_arrow + 3,
             f'+{medians["M"]["ext"]:.0f}%', fontsize=7.5, ha='center',
             color='#1565C0', fontweight='bold')

    # --- Arrows on female panel ---
    fc = medians['F']['ctrl']
    fr = medians['F']['rapa']

    # Female gap arrow (primary, red)
    ax2.annotate('', xy=(fr, y_arrow), xytext=(fc, y_arrow),
                 arrowprops=dict(arrowstyle='->', color='#CC0000', lw=1.4,
                                 shrinkA=0, shrinkB=0))
    ax2.text((fc + fr) / 2, y_arrow + 3,
             f'+{medians["F"]["ext"]:.0f}%', fontsize=7.5, ha='center',
             color='#CC0000', fontweight='bold')

    # Male gap reference arrow (ghost, below — same start, male arrow length)
    y_ghost = 45
    male_gap = mr - mc
    ax2.annotate('', xy=(fc + male_gap, y_ghost), xytext=(fc, y_ghost),
                 arrowprops=dict(arrowstyle='->', color='#1565C0', lw=1.0,
                                 linestyle='dashed', shrinkA=0, shrinkB=0))
    ax2.text(fc - 0.02, y_ghost,
             f'male +{medians["M"]["ext"]:.0f}%', fontsize=6.5, ha='right',
             va='center', color='#1565C0', style='italic')

    fig.tight_layout(w_pad=2.0)

    outpath = FIGURES_DIR / 'rapa_survival_mf.png'
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    plt.close(fig)

    print(f'Male:   ctrl median={mc:.3f}, rapa median={mr:.3f}, ext=+{medians["M"]["ext"]:.1f}%')
    print(f'Female: ctrl median={fc:.3f}, rapa median={fr:.3f}, ext=+{medians["F"]["ext"]:.1f}%')
    print(f'Figure saved: {outpath}')


if __name__ == '__main__':
    main()
