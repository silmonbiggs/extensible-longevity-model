#!/usr/bin/env python3
"""
Generate all figures for the ELM interactive diagram.

Usage:
    python generate_figures.py              # regenerate all figures
    python generate_figures.py stage_01     # regenerate just stage_01 figures
    python generate_figures.py w_meth       # regenerate just the w_meth sweep

All figures are saved to figures/ and referenced by the HTML slide deck.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches
from itertools import combinations, product
from pathlib import Path

from elm.model import (
    simulate, run_control, calculate_lifespan_extension,
    run_combination, run_combination_extension,
    derive_normalization_constants, michaelis_menten,
)
from elm.compounds import (
    ITP_VALIDATION, COMPOUNDS, get_compound,
    TAGUCHI_8_STACK, get_itp_start_time,
)
from elm.pathways import (
    BIOAGE_PARAMS, HETEROPLASMY_PARAMS, SENESCENCE_PARAMS,
    SASP_PARAMS, AMPK_PARAMS, NAD_PARAMS, SIRTUIN_PARAMS,
    AUTOPHAGY_PARAMS, TSC2_PARAMS, MTORC1_PARAMS,
    DNA_REPAIR_PARAMS, METHYLATION_PARAMS,
)
from elm.sex_mechanisms import apply_sex_modifier
from elm.dose_response import effective_dose, DOSE_RESPONSE_PARAMS

FIGURES_DIR = Path(os.path.dirname(__file__)) / '..' / 'docs' / 'figures'
DT = 0.002
DT_FAST = 0.005  # Faster dt for mass enumeration (stage 14/15)
DT_FINE = 0.001  # Fine dt for whatif/audit scenarios (long t_max)
T_MAX = 2.5
DPI = 150


# =============================================================================
# HELPERS
# =============================================================================

def ctrl(sex='M'):
    return run_control(sex=sex, t_max=T_MAX, dt=DT)

def treated(compound, sex='M'):
    return simulate(compound=compound, sex=sex, t_max=T_MAX, dt=DT)

def ext(compound, sex='M'):
    c = ctrl(sex)
    t = treated(compound, sex)
    return calculate_lifespan_extension(t, c)

def save(fig, name):
    fig.savefig(FIGURES_DIR / name, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f'  {name}')


# =============================================================================
# CONSTANTS
# =============================================================================

COMPOUND_COLORS = {
    'rapamycin': '#e74c3c', 'acarbose': '#6495ed', 'glycine': '#1abc9c',
    'canagliflozin': '#2ecc71', '17_alpha_estradiol': '#9b59b6', 'aspirin': '#f39c12',
}
COMPOUND_NAMES = {
    'rapamycin': 'Rapamycin', 'acarbose': 'Acarbose', 'glycine': 'Glycine',
    'canagliflozin': 'Canagliflozin', '17_alpha_estradiol': '17\u03b1-Estradiol',
    'aspirin': 'Aspirin',
}
COMPOUND_FILE_KEYS = {
    'rapamycin': 'rapa', 'acarbose': 'acarb', 'glycine': 'glycine',
    'canagliflozin': 'cana', '17_alpha_estradiol': '17ae2', 'aspirin': 'aspirin',
}
COMPOUND_STATE_COLORS = {
    'rapamycin': '#fa8072', 'acarbose': '#6495ed', 'glycine': '#1abc9c',
    'canagliflozin': '#2ecc71', '17_alpha_estradiol': '#9b59b6', 'aspirin': '#f39c12',
}
# Short labels for calibration bar chart
SHORT_NAMES = {
    'rapamycin': 'Rapa', 'acarbose': 'Acarb', 'canagliflozin': 'Cana',
    '17_alpha_estradiol': '17\u03b1E2', 'aspirin': 'Aspirin', 'glycine': 'Glycine',
}
ITP_COMPOUND_ORDER = ['rapamycin', 'acarbose', 'canagliflozin',
                       '17_alpha_estradiol', 'aspirin', 'glycine']

BLUE = '#6495ed'
SALMON = '#fa8072'


# =============================================================================
# STAGE 00: The Thesis
# =============================================================================

def gen_stage_00():
    c = ctrl()
    idx_end = int(1.05 / DT) + 1
    t = c.t[:idx_end]
    H = c.Heteroplasmy[:idx_end]
    BA = c.BioAge[:idx_end]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 5), sharex=True,
                                     gridspec_kw={'hspace': 0.15})
    fig.suptitle('Core Driver of Aging', fontsize=13, fontweight='bold')

    ax1.plot(t, H, '-', color='#e74c3c', linewidth=2.5)
    ax1.set_ylabel('Heteroplasmy (H)')
    ax1.set_ylim(0, 0.20)
    ax1.grid(True, alpha=0.3)

    ax2.plot(t[1:], BA[1:], '-', color='#2c3e50', linewidth=2.5)
    ax2.axhline(1.0, color='gray', linestyle='--', alpha=0.7, label='Death threshold')
    ax2.set_ylabel('Phenotypic Age')
    ax2.set_xlabel('Normalized Age')
    ax2.legend(loc='lower right', fontsize=9)
    ax2.grid(True, alpha=0.3)

    save(fig, 'stage_00_thesis.png')


# =============================================================================
# STAGE 01: Heteroplasmy + NAD decline
# =============================================================================

def gen_stage_01():
    c = ctrl()
    n = int(1.0 / DT) + 1
    t = c.t[:n]
    H = c.Heteroplasmy[:n]
    foxo3 = c.FOXO3[:n]

    p = HETEROPLASMY_PARAMS
    tau = p['time_scale']
    Neff = p['Neff']
    mu = p['transcription_advantage']
    s0 = p['k_selection_base']
    sf = p['k_selection_foxo3']

    drift = H * (1 - H) / Neff * tau
    mutation_pressure = mu * H * (1 - H) * tau
    selection = (s0 + sf * foxo3) * H * tau
    net = drift + mutation_pressure - selection

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 5), sharex=True,
                                     gridspec_kw={'hspace': 0.15})
    fig.suptitle('Kowald-Kirkwood mtDNA Dynamics', fontsize=13, fontweight='bold')

    ax1.plot(t, H, '-', color='#e74c3c', linewidth=2.5)
    ax1.set_ylabel('Heteroplasmy (H)')
    ax1.set_ylim(0, 0.20)
    ax1.grid(True, alpha=0.3)

    ax2.plot(t, drift, '-', color='#ff9900', linewidth=1.5, label='Drift')
    ax2.plot(t, mutation_pressure, '-', color='#e74c3c', linewidth=2, label='Mutation pressure')
    ax2.plot(t, -selection, '-', color='#2ecc71', linewidth=2, label='\u2212Selection (FOXO3)')
    ax2.plot(t, net, '--', color='#333333', linewidth=2.5, label='Net dH/dt')
    ax2.set_ylabel('Rate (per unit time)')
    ax2.set_xlabel('Normalized Age')
    ax2.set_ylim(-0.2, 0.6)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left', fontsize=9, framealpha=0.9)

    save(fig, 'stage_01_heteroplasmy.png')

    # NAD decline
    fig2, ax = plt.subplots(figsize=(5, 3.2))
    ax.plot(t, c.NAD[:n], '-', color='#7b4ea3', linewidth=2.5)
    ax.set_ylabel('NAD+ (normalized)')
    ax.set_xlabel('Normalized Age')
    ax.set_title('NAD+ Decline', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    save(fig2, 'stage_01_nad_decline.png')


# =============================================================================
# STAGE 02: DNA damage
# =============================================================================

def gen_stage_02():
    c = ctrl()
    n = int(1.0 / DT) + 1
    t = c.t[:n]
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.plot(t, c.DNA_damage[:n], '-', color='#f0a500', linewidth=2.5)
    ax.set_ylabel('DNA Damage')
    ax.set_xlabel('Normalized Age')
    ax.set_title('DNA Damage Accumulation', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    save(fig, 'stage_02_dna_damage.png')


# =============================================================================
# STAGE 03: Autophagy (control vs rapamycin)
# =============================================================================

def gen_stage_03():
    c = ctrl()
    r = treated('rapamycin')
    n = len(c.t)  # full length to 1.5
    t = c.t[:n]
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.plot(t, c.Autophagy[:n], '--', color='gray', linewidth=2, label='Control')
    ax.plot(t, r.Autophagy[:n], '-', color='#00897b', linewidth=2, label='+ Rapamycin')
    ax.set_ylabel('Autophagy Activity')
    ax.set_xlabel('Normalized Age')
    ax.set_title('Autophagy Activity', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 1.0)
    ax.set_xlim(0, 1.5)
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(True, alpha=0.3)
    save(fig, 'stage_03_autophagy.png')


# =============================================================================
# STAGE 04: Senescence + SASP
# =============================================================================

def gen_stage_04():
    c = ctrl()
    n = int(1.0 / DT) + 1
    t = c.t[:n]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 5), sharex=True,
                                     gridspec_kw={'hspace': 0.15})
    fig.suptitle('Senescence & SASP', fontsize=13, fontweight='bold')

    ax1.plot(t, c.SenCells[:n], '-', color='#f0a500', linewidth=2.5)
    ax1.set_ylabel('Senescent Cells')
    ax1.grid(True, alpha=0.3)

    ax2.plot(t, c.SASP[:n], '-', color='#c0392b', linewidth=2.5)
    ax2.set_ylabel('SASP')
    ax2.set_xlabel('Normalized Age')
    ax2.grid(True, alpha=0.3)

    save(fig, 'stage_04_senescence.png')


# =============================================================================
# STAGE 05: All intermediates (6-panel)
# =============================================================================

def gen_stage_05():
    c = ctrl()
    n = int(1.0 / DT) + 1
    t = c.t[:n]

    fig, axes = plt.subplots(3, 2, figsize=(5, 6.5))
    vars_info = [
        ('NAD', 'NAD+', '#7b4ea3'),
        ('SenCells', 'SenCells', '#f0a500'),
        ('SASP', 'SASP', '#c0392b'),
        ('Methylation', 'Methylation', '#27ae60'),
        ('DNA_damage', 'DNA Damage', '#f0a500'),
        ('BioAge', 'Phenotypic Age', '#2c3e50'),
    ]
    for idx, (attr, label, color) in enumerate(vars_info):
        ax = axes.flatten()[idx]
        # Skip first data point in bottom-right panel to avoid jump at t=0
        if idx == 5:
            ax.plot(t[1:], getattr(c, attr)[1:n], '-', color=color, linewidth=2)
        else:
            ax.plot(t, getattr(c, attr)[:n], '-', color=color, linewidth=2)
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        if attr == 'BioAge':
            ax.axhline(1.0, color='gray', linestyle='--', alpha=0.7)
        # x-label only on bottom row
        if idx >= 4:
            ax.set_xlabel('Normalized Age')

    plt.tight_layout()
    save(fig, 'stage_05_intermediates.png')


# =============================================================================
# STAGE 06: Compound bar + states (12 figures)
# =============================================================================

def gen_stage_06():
    c = ctrl()
    for cname in ITP_COMPOUND_ORDER:
        val = ITP_VALIDATION[cname]
        t_result = treated(cname)
        e = calculate_lifespan_extension(t_result, c)
        target = val.target_male
        fkey = COMPOUND_FILE_KEYS[cname]
        color = COMPOUND_STATE_COLORS[cname]
        display = COMPOUND_NAMES[cname]

        # Bar chart
        fig, ax = plt.subplots(figsize=(4.5, 4))
        bars = ax.bar(['ITP\nObserved', 'PALM\nModel'], [target, e],
                       color=[BLUE, SALMON], edgecolor='black', linewidth=0.5)
        ax.set_ylabel('Male Lifespan Extension (%)')
        ax.set_title(display, fontsize=13, fontweight='bold')
        ax.set_ylim(0, max(target, e) * 1.3 + 1)
        for i, v in enumerate([target, e]):
            ax.text(i, v + 0.3, f'{v:.1f}%', ha='center', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        save(fig, f'stage_06_{fkey}_bar.png')

        # State variables (3x2)
        n_pts = len(c.t)  # full to 1.5
        tt = c.t[:n_pts]
        fig, axes = plt.subplots(3, 2, figsize=(6, 6))
        fig.suptitle(f'{display}: State Variables', fontsize=13, fontweight='bold')
        state_vars = [
            ('Heteroplasmy', 'Heteroplasmy (H)'),
            ('NAD', 'NAD+'),
            ('DNA_damage', 'DNA Damage'),
            ('SenCells', 'Senescent Cells'),
            ('Methylation', 'Methylation'),
            ('BioAge', 'BioAge'),
        ]
        for idx, (attr, label) in enumerate(state_vars):
            ax = axes.flatten()[idx]
            ax.plot(tt, getattr(c, attr)[:n_pts], '--', color='gray', linewidth=1.5, label='Control')
            ax.plot(tt, getattr(t_result, attr)[:n_pts], '-', color=color, linewidth=2, label=display)
            ax.set_title(label, fontsize=10)
            ax.grid(True, alpha=0.3)
            if attr == 'BioAge':
                ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
            if idx == 0:
                ax.legend(fontsize=8)
            ax.set_xlim(0, 1.5)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        save(fig, f'stage_06_{fkey}_states.png')


# =============================================================================
# STAGE 07: Male calibration bar chart + survival curves
# =============================================================================

def gen_stage_07():
    c = ctrl()
    compounds = ITP_COMPOUND_ORDER
    targets = [ITP_VALIDATION[cn].target_male for cn in compounds]
    preds = [calculate_lifespan_extension(treated(cn), c) for cn in compounds]
    names = [SHORT_NAMES[cn] for cn in compounds]

    # Bar chart
    fig, ax = plt.subplots(figsize=(5, 4))
    x = np.arange(len(names))
    w = 0.35
    ax.bar(x - w/2, targets, w, label='ITP Observed', color=BLUE)
    ax.bar(x + w/2, preds, w, label='ELM Model', color=SALMON)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('Male Lifespan Extension (%)')
    ax.set_title('ITP Calibration (Male)', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 30)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    save(fig, 'stage_07_male_calibration.png')

    # Survival curves
    fig2, ax2 = plt.subplots(figsize=(5, 4))
    ax2.plot(c.t, c.Survival, 'k-', linewidth=2.5, label='Control')
    for cn in compounds:
        r = treated(cn)
        ax2.plot(r.t, r.Survival, '-', color=COMPOUND_COLORS[cn],
                 linewidth=1.5, label=COMPOUND_NAMES[cn])
    ax2.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
    ax2.text(0.05, 0.52, '50%', fontsize=9, color='gray')
    ax2.set_xlabel('Normalized Age')
    ax2.set_ylabel('Survival')
    ax2.set_title('Male Survival Curves', fontsize=12, fontweight='bold')
    ax2.set_xlim(0, 1.5)
    ax2.set_ylim(0, 1.0)
    ax2.legend(fontsize=8, loc='lower left', ncol=2)
    ax2.grid(True, alpha=0.3)
    save(fig2, 'stage_07_male_survival.png')


# =============================================================================
# STAGE 08: Rapa+Acarbose validation
# =============================================================================

def gen_stage_08():
    c = ctrl()
    rapa_ext = ext('rapamycin')
    acarb_ext = ext('acarbose')
    # ITP combo study: both drugs started at 9 months (Strong 2022, PMID 36179270)
    RAPA_ACARB_START = 9.0 / 30.0
    combo_ext = run_combination_extension({'rapamycin': 1.0, 'acarbose': 1.0}, sex='M',
                                          t_max=T_MAX, dt=DT, start_time=RAPA_ACARB_START)
    naive = rapa_ext + acarb_ext
    observed = 34.0

    # Bar chart
    fig, ax = plt.subplots(figsize=(5, 4.5))
    labels = ['Rapa\nalone', 'Acarb\nalone', 'Naive\nsum', 'Model\nprediction', 'ITP\nobserved']
    vals = [rapa_ext, acarb_ext, naive, combo_ext, observed]
    colors = [BLUE, BLUE, '#aaaaaa', SALMON, '#2ecc71']
    bars = ax.bar(labels, vals, color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(34.0, color='#2ecc71', linestyle='--', linewidth=1.5, alpha=0.7)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, v + 0.5, f'{v:.1f}%',
                ha='center', fontsize=10, fontweight='bold')
    ax.set_ylabel('Lifespan Extension (%)')
    ax.set_title('Rapa + Acarbose Combo', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 60)
    ax.grid(True, alpha=0.3, axis='y')
    save(fig, 'stage_08_mixing_bars.png')

    # Pathway saturation curve
    fig2, ax2 = plt.subplots(figsize=(5, 4))
    x = np.linspace(0, 2.0, 200)
    y = np.array([michaelis_menten(xi, AMPK_PARAMS['Km_ampk'], AMPK_PARAMS['Vmax_ampk']) for xi in x])
    ax2.plot(x, y, '-', color='#1a3a5c', linewidth=2.5)

    rapa_ampk = 0.2991
    acarb_ampk = 0.1927
    combined = rapa_ampk + acarb_ampk
    markers = [(rapa_ampk, 'Rapa', '#e74c3c', (10, -20)),
               (acarb_ampk, 'Acarb', '#6495ed', (-40, 10)),
               (combined, 'R+A', '#9b59b6', (10, 10))]
    for inp, label, clr, offset in markers:
        yv = michaelis_menten(inp, AMPK_PARAMS['Km_ampk'], AMPK_PARAMS['Vmax_ampk'])
        ax2.plot(inp, yv, 'o', color=clr, markersize=10, zorder=5)
        ax2.annotate(label, (inp, yv), textcoords='offset points',
                     xytext=offset, fontsize=10, color=clr, fontweight='bold')
    ax2.axhline(1.0, color='gray', linestyle=':', alpha=0.5)
    ax2.text(0.05, 1.02, 'Vmax', fontsize=9, color='gray')
    ax2.set_xlabel('Total AMPK Input')
    ax2.set_ylabel('Effective AMPK Activity')
    ax2.set_title('AMPK Saturation (Hill)', fontsize=13, fontweight='bold')
    ax2.set_xlim(0, 2.0)
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)
    save(fig2, 'stage_08_pathway_saturation.png')

    # Survival curves
    fig3, ax3 = plt.subplots(figsize=(6, 4.5))
    ax3.plot(c.t, c.Survival * 100, 'k-', linewidth=2.5, label='Control')
    for cn in ITP_COMPOUND_ORDER:
        r = treated(cn)
        ax3.plot(r.t, r.Survival * 100, '-', color=COMPOUND_COLORS[cn],
                 linewidth=1.5, label=COMPOUND_NAMES[cn])
    ax3.axhline(50, color='gray', linestyle='--', alpha=0.5)
    ax3.text(0.05, 52, '50%', fontsize=9, color='gray')
    ax3.set_xlabel('Normalized Age')
    ax3.set_ylabel('Surviving (%)')
    ax3.set_title('Male Survival Curves', fontsize=13, fontweight='bold')
    ax3.set_xlim(0, 1.5)
    ax3.set_ylim(0, 105)
    ax3.legend(fontsize=8, loc='lower left', ncol=2)
    ax3.grid(True, alpha=0.3)
    save(fig3, 'stage_08_survival_curves.png')


# =============================================================================
# STAGE 09: Sex differences -- AMPK compounds
# =============================================================================

def gen_stage_09():
    ampk_compounds = ['acarbose', 'canagliflozin', 'glycine']
    names = [COMPOUND_NAMES[c_] for c_ in ampk_compounds]

    fig, ax = plt.subplots(figsize=(5, 4))
    x = np.arange(len(names))
    w = 0.35
    males = [ext(c_, 'M') for c_ in ampk_compounds]
    females = [ext(c_, 'F') for c_ in ampk_compounds]
    targets_m = [ITP_VALIDATION[c_].target_male for c_ in ampk_compounds]
    targets_f = [ITP_VALIDATION[c_].target_female for c_ in ampk_compounds]

    ax.bar(x - w/2, males, w, label='Male', color=BLUE)
    ax.bar(x + w/2, females, w, label='Female', color=SALMON)
    ax.scatter(x - w/2, targets_m, marker='D', color='#1a1a1a', s=50, zorder=5)
    ax.scatter(x + w/2, targets_f, marker='D', color='#1a1a1a', s=50, zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('Lifespan Extension (%)')
    ax.set_title('AMPK Compounds: M vs F', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 25)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    save(fig, 'stage_09_ampk_mf_bars.png')

    # Saturation curve with M/F
    fig2, ax2 = plt.subplots(figsize=(5, 4))
    x_range = np.linspace(0, 1.5, 200)
    y_range = np.array([michaelis_menten(xi, AMPK_PARAMS['Km_ampk']) for xi in x_range])
    ax2.plot(x_range, y_range, '-', color='#1a3a5c', linewidth=2)

    m_input = 0.55
    f_input = 0.15
    m_y = michaelis_menten(m_input, AMPK_PARAMS['Km_ampk'])
    f_y = michaelis_menten(f_input, AMPK_PARAMS['Km_ampk'])
    ax2.plot(m_input, m_y, 's', color=BLUE, markersize=12, label='Male', zorder=5)
    ax2.plot(f_input, f_y, 's', color='#e74c3c', markersize=12, label='Female', zorder=5)
    ax2.annotate('Dim. factor\n= 0.14', (f_input, f_y), textcoords='offset points',
                 xytext=(20, -15), fontsize=9, fontstyle='italic', color='#e74c3c')
    ax2.axhline(1.0, color='gray', linestyle=':', alpha=0.5)
    ax2.text(0.05, 1.02, 'Vmax', fontsize=9, color='gray')
    ax2.set_xlabel('Effective AMPK Input')
    ax2.set_ylabel('AMPK Activity')
    ax2.set_title('AMPK Saturation: M vs F', fontsize=13, fontweight='bold')
    ax2.set_xlim(0, 1.5)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    save(fig2, 'stage_09_ampk_saturation.png')


# =============================================================================
# STAGE 10: Rapamycin PK + full M/F comparison
# =============================================================================

def gen_stage_10():
    # Rapamycin dose-response
    fig, ax = plt.subplots(figsize=(5, 4))
    doses = np.linspace(0.01, 6, 200)
    # sigmoid dose-response
    base_inhib = 0.84
    ec50 = (1.0 - base_inhib) / base_inhib  # ~0.19
    male_inh = np.array([d / (ec50 + d) * 100 for d in doses])
    female_inh = np.array([3.48 * d / (ec50 + 3.48 * d) * 100 for d in doses])
    female_inh = np.clip(female_inh, 0, 100)

    ax.plot(doses, male_inh, '-', color='#1a3a5c', linewidth=2.5)
    # Mark M and F at dose=1
    m_val = 1.0 / (ec50 + 1.0) * 100
    f_val = min(100, 3.48 / (ec50 + 3.48) * 100)
    ax.plot(1.0, m_val, 's', color=BLUE, markersize=12, zorder=5)
    ax.plot(3.48, f_val, 's', color='#e74c3c', markersize=12, zorder=5)
    ax.annotate(f'Male ({m_val:.0f}%)', (1.0, m_val), textcoords='offset points',
                xytext=(-15, -20), fontsize=9, color=BLUE, fontweight='bold')
    ax.annotate(f'Female ({f_val:.0f}%)', (3.48, f_val), textcoords='offset points',
                xytext=(-15, 10), fontsize=9, color='#e74c3c', fontweight='bold')
    # Double-headed arrow between x=1 and x=3.48
    ax.annotate('', xy=(3.48, 70), xytext=(1.0, 70),
                arrowprops=dict(arrowstyle='<->', color='#555555', lw=1.5))
    ax.text(2.24, 56, '3.5\u00d7 blood\nlevel', ha='center', fontsize=9,
            fontstyle='italic', color='#555555')
    ax.axhline(100, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Relative Blood Concentration')
    ax.set_ylabel('mTORC1 Inhibition (%)')
    ax.set_title('Rapamycin Dose-Response', fontsize=13, fontweight='bold')
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)
    save(fig, 'stage_10_rapa_dose_response.png')

    # Full M/F comparison with PK + AMPK annotations
    compounds = ['rapamycin', 'acarbose', 'canagliflozin', 'glycine']
    names = [COMPOUND_NAMES[c_] for c_ in compounds]
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    x = np.arange(len(names))
    w = 0.35
    males = [ext(c_, 'M') for c_ in compounds]
    females = [ext(c_, 'F') for c_ in compounds]
    targets_m = [ITP_VALIDATION[c_].target_male for c_ in compounds]
    targets_f = [ITP_VALIDATION[c_].target_female for c_ in compounds]
    ax2.bar(x - w/2, males, w, label='Male', color=BLUE)
    ax2.bar(x + w/2, females, w, label='Female', color=SALMON)
    ax2.scatter(x - w/2, targets_m, marker='D', color='#1a1a1a', s=50, zorder=5)
    ax2.scatter(x + w/2, targets_f, marker='D', color='#1a1a1a', s=50, zorder=5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(names, rotation=15)
    ax2.set_ylabel('Lifespan Extension (%)')
    ax2.set_title('M vs F: + Rapamycin', fontsize=13, fontweight='bold')
    ax2.set_ylim(-10, 35)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')

    # Mechanism annotations below x-axis
    # PK bracket under Rapa
    y_ann = -8
    ax2.annotate('', xy=(0.3, y_ann), xytext=(-0.3, y_ann),
                 arrowprops=dict(arrowstyle='-', color='blue', lw=1.5),
                 annotation_clip=False)
    ax2.text(0, y_ann - 1.5, 'PK', ha='center', fontsize=9,
             fontweight='bold', color='blue')
    # AMPK bracket spanning Acarb-Glycine
    ax2.annotate('', xy=(3.3, y_ann), xytext=(0.7, y_ann),
                 arrowprops=dict(arrowstyle='-', color='black', lw=1.5),
                 annotation_clip=False)
    ax2.text(2.0, y_ann - 1.5, 'AMPK', ha='center', fontsize=9,
             fontweight='bold', color='black')

    plt.tight_layout()
    save(fig2, 'stage_10_mf_updated.png')


# =============================================================================
# STAGE 11: Testosterone + Estrogen mechanisms
# =============================================================================

def gen_stage_11():
    # 11a: Testosterone
    e2_ext_m = ext('17_alpha_estradiol', 'M')
    fig, ax = plt.subplots(figsize=(4.5, 4))
    bars = ax.bar(['Male', 'Female', 'Castrated\nMale'], [e2_ext_m, 0, 0],
                   color=[BLUE, '#f0f0f0', '#f0f0f0'], edgecolor='black', linewidth=0.5)
    ax.text(0, e2_ext_m + 0.5, f'+{e2_ext_m:.0f}%', ha='center', fontsize=11, fontweight='bold')
    ax.text(1, 0.5, '+0%', ha='center', fontsize=11, fontweight='bold')
    ax.text(2, 0.5, '+0%', ha='center', fontsize=11, fontweight='bold')
    ax.set_ylabel('Lifespan Extension (%)')
    ax.set_title('17\u03b1E2: Requires Testosterone', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 25)
    ax.text(0.5, 22, '5\u03b1-reductase requires\ntestosterone substrate',
            ha='center', fontsize=8, fontstyle='italic', transform=ax.transData)
    ax.grid(True, alpha=0.3, axis='y')
    save(fig, 'stage_11a_testosterone.png')

    # 11b: Aspirin estrogen interference
    asp_ext_m = ext('aspirin', 'M')
    fig2, ax2 = plt.subplots(figsize=(4.5, 4))
    ax2.bar(['Male', 'Female'], [asp_ext_m, 0],
            color=[BLUE, '#f0f0f0'], edgecolor='black', linewidth=0.5)
    ax2.text(0, asp_ext_m + 0.3, f'+{asp_ext_m:.0f}%', ha='center', fontsize=11, fontweight='bold')
    ax2.text(1, 0.3, '+0%', ha='center', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Lifespan Extension (%)')
    ax2.set_title('Aspirin: Estrogen Interference', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 12)
    ax2.text(1.0, max(asp_ext_m * 0.8, 6), 'COX inhibition blocks\nestrogen cardioprotection',
             ha='center', fontsize=8, fontstyle='italic')
    ax2.grid(True, alpha=0.3, axis='y')
    save(fig2, 'stage_11b_estrogen_intf.png')

    # 11a cumulative: M vs F + 17aE2
    _gen_mf_cumulative(
        ['17_alpha_estradiol', 'acarbose', 'canagliflozin', 'glycine', 'rapamycin'],
        'stage_11a_mf_cumulative.png', 'M vs F: + 17\u03b1E2',
        brackets=[('TEST', '#e74c3c', 0, 0), ('AMPK', 'black', 1, 3), ('PK', 'blue', 4, 4)]
    )

    # 11b cumulative: M vs F + Aspirin
    _gen_mf_cumulative(
        ['aspirin', '17_alpha_estradiol', 'acarbose', 'canagliflozin', 'glycine', 'rapamycin'],
        'stage_11b_mf_cumulative.png', 'M vs F: + Aspirin',
        brackets=[('E\u2082 INTF', '#cc00cc', 0, 0), ('TEST', '#e74c3c', 1, 1),
                  ('AMPK', 'black', 2, 4), ('PK', 'blue', 5, 5)]
    )

    # 11 complete: Four Sex Mechanisms (order matches rev_04)
    _gen_mf_cumulative(
        ['aspirin', '17_alpha_estradiol', 'canagliflozin', 'glycine', 'acarbose', 'rapamycin'],
        'stage_11_mf_complete.png', 'Four Sex Mechanisms',
        brackets=[('E\u2082 INTF', '#cc00cc', 0, 0), ('TEST', '#e74c3c', 1, 1),
                  ('AMPK', 'black', 2, 4), ('PK', 'blue', 5, 5)]
    )


def _gen_mf_cumulative(compounds, filename, title, brackets=None):
    names = [COMPOUND_NAMES[c_] for c_ in compounds]
    fig, ax = plt.subplots(figsize=(6, 5))
    x = np.arange(len(names))
    w = 0.35
    males = [ext(c_, 'M') for c_ in compounds]
    females = [ext(c_, 'F') for c_ in compounds]
    targets_m = [ITP_VALIDATION[c_].target_male for c_ in compounds]
    targets_f = [ITP_VALIDATION[c_].target_female for c_ in compounds]
    ax.bar(x - w/2, males, w, label='Male', color=BLUE)
    ax.bar(x + w/2, females, w, label='Female', color=SALMON)
    ax.scatter(x - w/2, targets_m, marker='D', color='#1a1a1a', s=50, zorder=5, label='ITP target')
    ax.scatter(x + w/2, targets_f, marker='D', color='#1a1a1a', s=50, zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('Lifespan Extension (%)')
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_ylim(-10, 35)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    # Add mechanism bracket annotations if provided
    # Two tiers: AMPK on lower tier spanning all columns; E2INF/TEST/PK on upper tier
    if brackets:
        n_items = len(x)
        tick_h = 1.0
        for label, color, i_start, i_end in brackets:
            if label == 'AMPK':
                # Lower tier: spans ALL columns
                y_ann = -7
                x0, x1 = -0.3, (n_items - 1) + 0.3
            else:
                # Upper tier
                y_ann = -3
                x0, x1 = i_start - 0.3, i_end + 0.3
            mid = (x0 + x1) / 2.0
            # Horizontal line
            ax.plot([x0, x1], [y_ann, y_ann], color=color, lw=2.0,
                    clip_on=False, solid_capstyle='butt')
            # Vertical end caps
            ax.plot([x0, x0], [y_ann + tick_h, y_ann], color=color, lw=2.0,
                    clip_on=False, solid_capstyle='butt')
            ax.plot([x1, x1], [y_ann + tick_h, y_ann], color=color, lw=2.0,
                    clip_on=False, solid_capstyle='butt')
            ax.text(mid, y_ann - 1.8, label, ha='center', fontsize=9,
                    fontweight='bold', color=color)

    plt.tight_layout()
    save(fig, filename)


# =============================================================================
# STAGE 14: Greedy subtraction + dose matrix
# =============================================================================

def gen_stage_14():
    compounds = ITP_COMPOUND_ORDER

    # --- Greedy subtraction at 1x ---
    for sex, suffix, bar_color in [('M', 'male', BLUE), ('F', 'female', SALMON)]:
        remaining = list(compounds)
        c = ctrl(sex)
        labels = []
        values = []
        remaining_names = []

        # Full cocktail
        full = {cn: 1.0 for cn in remaining}
        full_ext = run_combination_extension(full, sex=sex, t_max=T_MAX, dt=DT)
        labels.append(f'All {len(remaining)}')
        values.append(full_ext)
        remaining_names.append(', '.join([SHORT_NAMES[cn] for cn in remaining]))

        while len(remaining) > 1:
            min_loss = float('inf')
            min_cn = None
            for cn in remaining:
                reduced = {k: 1.0 for k in remaining if k != cn}
                r_ext = run_combination_extension(reduced, sex=sex, t_max=T_MAX, dt=DT)
                loss = full_ext - r_ext
                if loss < min_loss:
                    min_loss = loss
                    min_cn = cn
            remaining.remove(min_cn)
            combo = {cn: 1.0 for cn in remaining}
            step_ext = run_combination_extension(combo, sex=sex, t_max=T_MAX, dt=DT)
            label_count = len(remaining)
            labels.append(f'{label_count}: \u2212{COMPOUND_NAMES[min_cn]}')
            values.append(step_ext)
            remaining_names.append(', '.join([SHORT_NAMES[cn] for cn in remaining]))
            full_ext = step_ext

        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        y_pos = np.arange(len(values))
        bars = ax.barh(y_pos, values, color=bar_color, edgecolor='black', linewidth=0.5, height=0.7)
        for i, (bar, v) in enumerate(zip(bars, values)):
            ax.text(v * 0.5, i, f'{v:.1f}%', ha='center', va='center',
                    fontsize=10, fontweight='bold', color='white')
            ax.text(v + 1, i, remaining_names[i], ha='left', va='center',
                    fontsize=8, color='gray', fontstyle='italic')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlabel('Lifespan Extension (%)')
        sex_label = 'Male' if sex == 'M' else 'Female'
        ax.set_title(f'Predicted Benefit @ 1x Dose, {sex_label}', fontsize=13, fontweight='bold')
        ax.invert_yaxis()
        ref_line = 50 if sex == 'M' else 40
        ax.axvline(ref_line, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlim(0, 100)
        ax.grid(True, alpha=0.3, axis='x')
        save(fig, f'stage_14_greedy_{suffix}.png')

    # --- Greedy subtraction at 2x ---
    for sex, suffix, bar_color in [('M', 'male', BLUE), ('F', 'female', SALMON)]:
        remaining = list(compounds)
        labels = []
        values = []
        full = {cn: 2.0 for cn in remaining}
        full_ext = run_combination_extension(full, sex=sex, t_max=T_MAX, dt=DT)
        labels.append(f'All {len(remaining)} @2\u00d7')
        values.append(full_ext)

        while len(remaining) > 1:
            min_loss = float('inf')
            min_cn = None
            for cn in remaining:
                reduced = {k: 2.0 for k in remaining if k != cn}
                r_ext = run_combination_extension(reduced, sex=sex, t_max=T_MAX, dt=DT)
                loss = full_ext - r_ext
                if loss < min_loss:
                    min_loss = loss
                    min_cn = cn
            remaining.remove(min_cn)
            combo = {cn: 2.0 for cn in remaining}
            step_ext = run_combination_extension(combo, sex=sex, t_max=T_MAX, dt=DT)
            labels.append(f'\u2212{COMPOUND_NAMES[min_cn]}')
            values.append(step_ext)
            full_ext = step_ext

        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        y_pos = np.arange(len(values))
        bars = ax.barh(y_pos, values, color=bar_color, edgecolor='black', linewidth=0.5, height=0.7)
        for bar, v in zip(bars, values):
            ax.text(v * 0.5, bar.get_y() + bar.get_height()/2, f'{v:.1f}%',
                    ha='center', va='center', fontsize=10, fontweight='bold', color='white')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlabel('Lifespan Extension (%)')
        sex_label = 'Male' if sex == 'M' else 'Female'
        ax.set_title(f'Predicted Benefit @ 2x Dose, {sex_label}', fontsize=13, fontweight='bold')
        ax.invert_yaxis()
        ax.set_xlim(0, 70)
        ax.grid(True, alpha=0.3, axis='x')
        save(fig, f'stage_14_greedy_{suffix}_2x.png')

    # --- Greedy subtraction paired M/F ---
    # Run 1x greedy for both sexes for the paired chart
    paired_data = {}
    for sex in ['M', 'F']:
        remaining = list(compounds)
        labels_list = []
        vals_list = []
        full = {cn: 1.0 for cn in remaining}
        full_ext = run_combination_extension(full, sex=sex, t_max=T_MAX, dt=DT)
        labels_list.append(f'All {len(remaining)}')
        vals_list.append(full_ext)
        while len(remaining) > 1:
            min_loss = float('inf')
            min_cn = None
            for cn in remaining:
                reduced = {k: 1.0 for k in remaining if k != cn}
                r_ext = run_combination_extension(reduced, sex=sex, t_max=T_MAX, dt=DT)
                loss = full_ext - r_ext
                if loss < min_loss:
                    min_loss = loss
                    min_cn = cn
            remaining.remove(min_cn)
            combo = {cn: 1.0 for cn in remaining}
            step_ext = run_combination_extension(combo, sex=sex, t_max=T_MAX, dt=DT)
            labels_list.append(f'\u2212{COMPOUND_NAMES[min_cn]}')
            vals_list.append(step_ext)
            full_ext = step_ext
        paired_data[sex] = (labels_list, vals_list)

    # Use male labels (order may differ)
    n_rows = max(len(paired_data['M'][1]), len(paired_data['F'][1]))
    fig, ax = plt.subplots(figsize=(6.5, 5))
    y_pos = np.arange(n_rows)
    h = 0.35
    m_vals = paired_data['M'][1]
    f_vals = paired_data['F'][1]
    m_labels = paired_data['M'][0]

    ax.barh(y_pos - h/2, m_vals, h, color=BLUE, edgecolor='black', linewidth=0.5, label='Male')
    ax.barh(y_pos + h/2, f_vals, h, color=SALMON, edgecolor='black', linewidth=0.5, label='Female')
    for i in range(n_rows):
        ax.text(m_vals[i] * 0.5, i - h/2, f'{m_vals[i]:.1f}%', ha='center', va='center',
                fontsize=9, fontweight='bold', color='white')
        ax.text(f_vals[i] * 0.5, i + h/2, f'{f_vals[i]:.1f}%', ha='center', va='center',
                fontsize=9, fontweight='bold', color='white')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(m_labels)
    ax.set_xlabel('Lifespan Extension (%)')
    ax.set_title('Greedy Subtraction: Minimal Cocktail', fontsize=13, fontweight='bold')
    ax.invert_yaxis()
    ax.axvline(50, color='gray', linestyle='--', alpha=0.4)
    ax.axvline(40, color='gray', linestyle='--', alpha=0.4)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='x')
    save(fig, 'stage_14_greedy_subtraction.png')

    # --- Dose Matrix Heatmap ---
    _gen_dose_matrix()


def _gen_dose_matrix():
    compounds = ITP_COMPOUND_ORDER
    dose_levels = [0, 0.5, 1.0, 2.0]
    n_compounds = len(compounds)
    n_doses = len(dose_levels)
    dim = n_compounds * n_doses  # 24

    # Build labels
    axis_labels = []
    for cn in compounds:
        for d in dose_levels:
            axis_labels.append(f'{SHORT_NAMES[cn]} {d}\u00d7')

    # Compute matrix: lower-left = male, upper-right = female
    matrix = np.full((dim, dim), np.nan)
    for i in range(dim):
        for j in range(dim):
            ci, di = divmod(i, n_doses)
            cj, dj = divmod(j, n_doses)
            if ci == cj:
                continue  # skip same-compound
            if i == j:
                continue
            # Lower-left triangle: male (i > j)
            # Upper-right triangle: female (i < j)
            sex = 'M' if i > j else 'F'
            cn_i = compounds[ci]
            cn_j = compounds[cj]
            d_i = dose_levels[di]
            d_j = dose_levels[dj]
            doses = {}
            if d_i > 0:
                doses[cn_i] = d_i
            if d_j > 0:
                doses[cn_j] = d_j
            if len(doses) == 0:
                matrix[i, j] = 0.0
            else:
                matrix[i, j] = run_combination_extension(doses, sex=sex, t_max=T_MAX, dt=DT_FAST)

    fig, ax = plt.subplots(figsize=(18, 17))
    # Use green colormap
    cmap = plt.cm.Greens.copy()
    cmap.set_bad('lightgray')
    vmax = np.nanmax(matrix)
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=0, vmax=max(vmax, 50),
                   interpolation='nearest')

    # Cell text
    for i in range(dim):
        for j in range(dim):
            if not np.isnan(matrix[i, j]):
                val = matrix[i, j]
                color = 'white' if val > vmax * 0.6 else 'black'
                ax.text(j, i, f'{val:.0f}', ha='center', va='center', fontsize=12, fontweight='bold', color=color)

    ax.set_xticks(range(dim))
    ax.set_yticks(range(dim))
    ax.set_xticklabels(axis_labels, rotation=90, fontsize=12)
    ax.set_yticklabels(axis_labels, fontsize=12)

    # Color group labels on axes
    group_colors = [COMPOUND_COLORS[cn] for cn in compounds for _ in dose_levels]
    for i, (tick, clr) in enumerate(zip(ax.get_yticklabels(), group_colors)):
        tick.set_color(clr)
        tick.set_fontweight('bold')
    for i, (tick, clr) in enumerate(zip(ax.get_xticklabels(), group_colors)):
        tick.set_color(clr)
        tick.set_fontweight('bold')

    # Diagonal gray blocks
    for ci in range(n_compounds):
        start = ci * n_doses
        rect = plt.Rectangle((start - 0.5, start - 0.5), n_doses, n_doses,
                              facecolor='lightgray', edgecolor='gray', linewidth=0.5)
        ax.add_patch(rect)

    # Watermarks
    ax.text(dim * 0.25, dim * 0.75, 'MALE', fontsize=36, color='red', alpha=0.15,
            ha='center', va='center', fontweight='bold', rotation=45)
    ax.text(dim * 0.75, dim * 0.25, 'FEMALE', fontsize=36, color='red', alpha=0.15,
            ha='center', va='center', fontweight='bold', rotation=45)

    ax.set_title(
        'Pairwise Dose-Response Matrix  (Male lower-left, Female upper-right)\n'
        '6 compounds \u00d7 4 doses (0, 0.5x, 1x, 2x)',
        fontsize=16, fontweight='bold'
    )
    plt.colorbar(im, ax=ax, label='Lifespan Extension (%)', shrink=0.8)
    plt.tight_layout()
    save(fig, 'stage_14_dose_matrix.png')


# =============================================================================
# STAGE 15: Risk-benefit
# =============================================================================

def gen_stage_15():
    # Risk scores per compound per dose
    risk_data = {
        'rapamycin':          {'risks': [2.5, 4.0, 6.5], 'desc': 'Immunosuppression'},
        'acarbose':           {'risks': [1.5, 2.5, 4.0], 'desc': 'GI distress'},
        'canagliflozin':      {'risks': [1.0, 2.0, 3.5], 'desc': 'UTI / ketoacidosis'},
        '17_alpha_estradiol': {'risks': [1.5, 3.0, 5.0], 'desc': 'Hormonal effects'},
        'aspirin':            {'risks': [0.5, 1.0, 2.5], 'desc': 'GI bleeding'},
        'glycine':            {'risks': [0.3, 0.5, 1.0], 'desc': 'Minimal'},
    }

    compounds = ITP_COMPOUND_ORDER
    doses = [0.5, 1.0, 2.0]
    matrix = np.array([[risk_data[c_]['risks'][d] for d in range(3)] for c_ in compounds])
    descs = [risk_data[c_]['desc'] for c_ in compounds]

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=8)
    ax.set_xticks(range(len(doses)))
    ax.set_xticklabels([f'{d}\u00d7' for d in doses])
    ax.set_yticks(range(len(compounds)))
    ax.set_yticklabels([COMPOUND_NAMES[c_] for c_ in compounds], fontsize=11)
    for i in range(len(compounds)):
        for j in range(len(doses)):
            ax.text(j, i, f'{matrix[i,j]:.1f}', ha='center', va='center',
                    fontsize=10, fontweight='bold')
        # Risk description to the right
        ax.text(len(doses) - 0.3, i, descs[i], ha='left', va='center',
                fontsize=11, color='#555555')
    ax.set_title('Side-Effect Risk (0\u201310) by Intervention and Dose',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    save(fig, 'stage_15_risk_profile.png')

    # --- Precompute all cocktail extensions (one pass, both sexes) ---
    cocktails = _enumerate_cocktails(risk_data)
    ext_cache = _precompute_all_cocktail_extensions(cocktails)

    # --- Optimal cocktails ---
    _gen_optimal_cocktails(risk_data, cocktails, ext_cache)

    # --- Pareto ---
    _gen_pareto(risk_data, cocktails, ext_cache)


def _enumerate_cocktails(risk_data):
    """Enumerate all cocktails: each compound at dose 0, 0.5, 1.0, 2.0."""
    compounds = ITP_COMPOUND_ORDER
    dose_levels = [0, 0.5, 1.0, 2.0]
    risk_map = {}
    for c_ in compounds:
        for di, d in enumerate(dose_levels):
            if d == 0:
                risk_map[(c_, d)] = 0.0
            else:
                risk_map[(c_, d)] = risk_data[c_]['risks'][di - 1]

    results = []
    for combo in product(range(4), repeat=6):
        doses_dict = {}
        total_risk = 0.0
        n_active = 0
        for ci, di in enumerate(combo):
            c_ = compounds[ci]
            d = dose_levels[di]
            total_risk += risk_map[(c_, d)]
            if d > 0:
                doses_dict[c_] = d
                n_active += 1
        results.append({
            'doses': doses_dict,
            'risk': total_risk,
            'n_compounds': n_active,
            'combo': combo,
        })
    return results


def _precompute_all_cocktail_extensions(cocktails):
    """Precompute extensions for all cocktails, both sexes. Returns dict keyed by (combo_tuple, sex)."""
    cache = {}
    total = len(cocktails)
    for sex in ['M', 'F']:
        print(f'    Precomputing {total} cocktails for {sex}...')
        for i, c_ in enumerate(cocktails):
            key = (c_['combo'], sex)
            if len(c_['doses']) == 0:
                cache[key] = 0.0
            else:
                cache[key] = run_combination_extension(c_['doses'], sex=sex, t_max=T_MAX, dt=DT_FAST)
            if (i+1) % 500 == 0:
                print(f'      {i+1}/{total}')
    return cache


def _gen_optimal_cocktails(risk_data, cocktails, ext_cache):
    """Generate optimal cocktails at each risk budget using precomputed extensions."""
    risk_budgets = [5, 8, 12, 18, 25]

    fig, (ax_m, ax_f) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    for sex, ax, bar_color, title in [('M', ax_m, BLUE, 'Best Cocktail at Each Risk Budget: Male'),
                                       ('F', ax_f, SALMON, 'Best Cocktail at Each Risk Budget: Female')]:
        best_per_budget = []
        for budget in risk_budgets:
            eligible = [c for c in cocktails if c['risk'] <= budget]
            best = None
            best_ext = -999
            for c_ in eligible:
                ext_val = ext_cache[(c_['combo'], sex)]
                if ext_val > best_ext:
                    best_ext = ext_val
                    best = c_
            desc_parts = []
            if best and best['doses']:
                for cn in ITP_COMPOUND_ORDER:
                    if cn in best['doses']:
                        d = best['doses'][cn]
                        prefix = f'{d:g}\u00d7' if d != 1.0 else ''
                        desc_parts.append(f'{prefix}{SHORT_NAMES[cn]}')
            desc = ' + '.join(desc_parts)
            risk_used = best['risk'] if best else 0
            best_per_budget.append((budget, best_ext, f'{desc}  (risk={risk_used:.0f})'))

        y_pos = np.arange(len(risk_budgets))
        labels = [f'Risk \u2264 {b}' for b in risk_budgets]
        vals = [b[1] for b in best_per_budget]
        bars = ax.barh(y_pos, vals, color=bar_color, edgecolor='black', linewidth=0.5, height=0.6)
        for i, (budget, v, desc) in enumerate(best_per_budget):
            ax.text(v - 2, i, f'{v:.0f}%', ha='right', va='center',
                    fontsize=10, fontweight='bold', color='white')
            ax.text(v + 1, i, desc, ha='left', va='center', fontsize=8, color='gray')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlim(0, 100)
        ax.grid(True, alpha=0.3, axis='x')

    ax_f.set_xlabel('Lifespan Extension (%)')
    plt.tight_layout()
    save(fig, 'stage_15_optimal_cocktails.png')


def _gen_pareto(risk_data, cocktails, ext_cache):
    """Generate Pareto frontier scatter plots using precomputed extensions."""

    fig, (ax_m, ax_f) = plt.subplots(1, 2, figsize=(12.5, 5.5))
    fig.suptitle(
        "Risk-Benefit Pareto Frontier: What's the best cocktail at each risk budget?",
        fontsize=13, fontweight='bold'
    )

    for sex, ax, line_color, title in [('M', ax_m, 'blue', 'Male'),
                                         ('F', ax_f, 'red', 'Female')]:
        risks_all = []
        exts_all = []
        n_compounds_all = []

        for c_ in cocktails:
            ext_val = ext_cache[(c_['combo'], sex)]
            risks_all.append(c_['risk'])
            exts_all.append(ext_val)
            n_compounds_all.append(c_['n_compounds'])

        risks_arr = np.array(risks_all)
        exts_arr = np.array(exts_all)
        n_arr = np.array(n_compounds_all)

        # Scatter colored by # compounds
        sc = ax.scatter(risks_arr, exts_arr, c=n_arr, cmap='viridis', s=8, alpha=0.4,
                        vmin=0, vmax=6)

        # Compute Pareto frontier
        sorted_idx = np.argsort(risks_arr)
        frontier_risk = []
        frontier_ext = []
        best_so_far = -999
        for idx in sorted_idx:
            if exts_arr[idx] > best_so_far:
                best_so_far = exts_arr[idx]
                frontier_risk.append(risks_arr[idx])
                frontier_ext.append(exts_arr[idx])

        # Step function for frontier
        step_risk = []
        step_ext = []
        for i in range(len(frontier_risk)):
            if i > 0:
                step_risk.append(frontier_risk[i])
                step_ext.append(frontier_ext[i-1])
            step_risk.append(frontier_risk[i])
            step_ext.append(frontier_ext[i])

        ax.plot(step_risk, step_ext, '-', color=line_color, linewidth=2.5, alpha=0.8)

        ax.set_xlabel('Total Risk Score')
        ax.set_ylabel('Lifespan Extension (%)')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)

    plt.colorbar(sc, ax=ax_f, label='# Compounds', shrink=0.8)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, 'stage_15_pareto.png')


# =============================================================================
# STAGE 16: Novel targets + BioAge composition
# =============================================================================

def _build_whatif_interventions(base_compounds, novel_inj, sex):
    """Build merged interventions dict with per-compound sex modifiers applied.

    base_compounds: list of (compound_name, injections_dict) — each compound's
                    raw pathway injections (before sex adjustment).
    novel_inj:      optional dict of additional (non-sex-adjusted) pathway injections.
    sex:            'M' or 'F'
    """
    merged = {}
    earliest_start = 1.0
    for name, inj in base_compounds:
        start = get_itp_start_time(name)
        if start < earliest_start:
            earliest_start = start
        adjusted = apply_sex_modifier(name, dict(inj), sex)
        for k, v in adjusted.items():
            merged[k] = merged.get(k, 0) + v
    if novel_inj:
        for k, v in novel_inj.items():
            merged[k] = merged.get(k, 0) + v
    merged['start_time'] = earliest_start
    return merged


def gen_stage_16():
    # Base compounds: ITP unisex 4 at 1x calibrated dose
    base_compounds = [
        ('rapamycin', get_compound('rapamycin')),
        ('acarbose', get_compound('acarbose')),
        ('canagliflozin', get_compound('canagliflozin')),
        ('glycine', get_compound('glycine')),
    ]

    # Novel targets at 1x dose (Taguchi 8 audited values)
    novel_targets = {
        '+ NMN': get_compound('NMN'),
        '+ CD38 inhibitor': {'cd38_inhibitor': 0.50},
        '+ Mitophagy\n(Urolithin A)': get_compound('urolithin_A'),
        '+ Senolytic\n(D+Q)': get_compound('dasatinib_quercetin'),
    }

    # Compute base extension (use large t_max and fine dt for numerical accuracy)
    T_MAX_WHATIF = 5.0
    ctrl_m = run_control(sex='M', t_max=T_MAX_WHATIF, dt=DT_FINE)
    ctrl_f = run_control(sex='F', t_max=T_MAX_WHATIF, dt=DT_FINE)

    base_m_inj = _build_whatif_interventions(base_compounds, None, 'M')
    base_f_inj = _build_whatif_interventions(base_compounds, None, 'F')
    base_m_result = simulate(interventions=base_m_inj, sex='M', t_max=T_MAX_WHATIF, dt=DT_FINE)
    base_m_ext = calculate_lifespan_extension(base_m_result, ctrl_m)
    base_f_result = simulate(interventions=base_f_inj, sex='F', t_max=T_MAX_WHATIF, dt=DT_FINE)
    base_f_ext = calculate_lifespan_extension(base_f_result, ctrl_f)

    rows = [('ITP base\n(Rapa+Acarb+Cana+Gly)', base_m_ext, base_f_ext)]

    # Each novel target added individually
    for label, novel_inj in novel_targets.items():
        m_inj = _build_whatif_interventions(base_compounds, novel_inj, 'M')
        f_inj = _build_whatif_interventions(base_compounds, novel_inj, 'F')
        m_result = simulate(interventions=m_inj, sex='M', t_max=T_MAX_WHATIF, dt=DT_FINE)
        f_result = simulate(interventions=f_inj, sex='F', t_max=T_MAX_WHATIF, dt=DT_FINE)
        m_ext = calculate_lifespan_extension(m_result, ctrl_m)
        f_ext = calculate_lifespan_extension(f_result, ctrl_f)
        rows.append((label, m_ext, f_ext))

    # All 4 novel
    all_novel_inj = {}
    for novel_inj in novel_targets.values():
        for k, v in novel_inj.items():
            all_novel_inj[k] = all_novel_inj.get(k, 0) + v
    m_all_inj = _build_whatif_interventions(base_compounds, all_novel_inj, 'M')
    f_all_inj = _build_whatif_interventions(base_compounds, all_novel_inj, 'F')
    m_all = simulate(interventions=m_all_inj, sex='M', t_max=T_MAX_WHATIF, dt=DT_FINE)
    f_all = simulate(interventions=f_all_inj, sex='F', t_max=T_MAX_WHATIF, dt=DT_FINE)
    m_all_ext = calculate_lifespan_extension(m_all, ctrl_m)
    f_all_ext = calculate_lifespan_extension(f_all, ctrl_f)
    rows.append(('+ All 4 novel', m_all_ext, f_all_ext))

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6.5))
    n_rows = len(rows)
    y_pos = np.arange(n_rows)
    h = 0.35
    m_vals = [r[1] for r in rows]
    f_vals = [r[2] for r in rows]
    labels_list = [r[0] for r in rows]

    ax.barh(y_pos - h/2, m_vals, h, color=BLUE, edgecolor='black', linewidth=0.5, label='Male')
    ax.barh(y_pos + h/2, f_vals, h, color=SALMON, edgecolor='black', linewidth=0.5, label='Female')

    for i in range(n_rows):
        # Bold colored text inside bars
        mv, fv = m_vals[i], f_vals[i]
        if mv > 5:
            ax.text(mv * 0.5, i - h/2, f'+{mv:.0f}%', ha='center', va='center',
                    fontsize=9, fontweight='bold', color='white')
        if fv > 5:
            ax.text(fv * 0.5, i + h/2, f'+{fv:.0f}%', ha='center', va='center',
                    fontsize=9, fontweight='bold', color='white')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels_list, fontsize=9)
    ax.set_xlabel('Predicted Lifespan Extension (%)')
    ax.set_title('Adding Novel Targets to ITP Base Cocktail', fontsize=13, fontweight='bold')
    ax.invert_yaxis()
    ax.axvline(50, color='gray', linestyle='--', alpha=0.4)
    ax.axvline(40, color='gray', linestyle='--', alpha=0.4)
    max_val = max(max(m_vals), max(f_vals))
    ax.set_xlim(0, max_val * 1.15)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    save(fig, 'stage_16_whatif.png')

    # --- BioAge Decomposition ---
    _gen_bioage_composition(base_m_inj, m_all_inj, ctrl_m)


def _gen_bioage_composition(base_interventions, all_novel_interventions, ctrl_m):
    """Stacked bar showing BioAge decomposition at death for 3 scenarios."""
    c = ctrl_m

    # Scenario 1: Control
    idx_ctrl = int(round(c.t_death / DT_FINE))
    idx_ctrl = min(idx_ctrl, len(c.t) - 1)

    # Scenario 2: ITP base
    base_result = simulate(interventions=dict(base_interventions), sex='M', t_max=3.0, dt=DT_FINE)
    idx_base = int(round(base_result.t_death / DT_FINE))
    idx_base = min(idx_base, len(base_result.t) - 1)

    # Scenario 3: + Novel targets
    novel_result = simulate(interventions=dict(all_novel_interventions), sex='M', t_max=5.0, dt=DT_FINE)
    idx_novel = int(round(novel_result.t_death / DT_FINE))
    idx_novel = min(idx_novel, len(novel_result.t) - 1)

    def get_fracs(result, idx):
        meth = BIOAGE_PARAMS['w_meth'] * result.Methylation[idx] / BIOAGE_PARAMS['meth_norm']
        dmg = BIOAGE_PARAMS['w_damage'] * result.DNA_damage[idx] / BIOAGE_PARAMS['damage_norm']
        h = BIOAGE_PARAMS['w_hetero'] * result.Heteroplasmy[idx] / BIOAGE_PARAMS['H_norm']
        sen = BIOAGE_PARAMS['w_sen'] * result.SenCells[idx] / BIOAGE_PARAMS['sen_norm']
        return [meth, dmg, h, sen]

    fracs_ctrl = get_fracs(c, idx_ctrl)
    fracs_base = get_fracs(base_result, idx_base)
    fracs_novel = get_fracs(novel_result, idx_novel)

    bar_labels = ['Control', 'ITP base\n(Rapa+Acarb\n+Cana+Gly)', '+ Novel targets\n(NMN+CD38i\n+D+Q+UroA)']
    times = [c.t_death, base_result.t_death, novel_result.t_death]
    all_fracs = [fracs_ctrl, fracs_base, fracs_novel]
    component_names = ['Methylation', 'DNA Damage', 'Heteroplasmy', 'Senescent Cells']
    comp_colors = ['#27ae60', '#f39c12', '#e74c3c', '#9b59b6']

    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(3)
    for ci in range(4):
        bottoms = [sum(all_fracs[bi][:ci]) for bi in range(3)]
        heights = [all_fracs[bi][ci] for bi in range(3)]
        ax.bar(x, heights, bottom=bottoms, color=comp_colors[ci], label=component_names[ci],
               edgecolor='white', linewidth=0.5, width=0.6)
        # % labels inside each segment
        for bi in range(3):
            total = sum(all_fracs[bi])
            pct = all_fracs[bi][ci] / total * 100 if total > 0 else 0
            if pct > 5:
                y_pos = bottoms[bi] + heights[bi] / 2
                ax.text(bi, y_pos, f'{pct:.0f}%', ha='center', va='center',
                        fontsize=8, fontweight='bold', color='white')

    # Time annotations above each bar
    for bi in range(3):
        total = sum(all_fracs[bi])
        t_val = times[bi]
        ext_pct = (t_val - c.t_death) / c.t_death * 100
        if bi == 0:
            label = f't={t_val:.2f}'
        else:
            label = f't={t_val:.2f}/+{ext_pct:.0f}%'
        ax.text(bi, total + 0.03, label, ha='center', fontsize=8, fontweight='bold')

    ax.axhline(1.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.text(2.5, 1.02, 'BioAge = 1.0', fontsize=8, color='gray')

    # Arrow pointing to heteroplasmy segment in control
    h_bottom_ctrl = sum(fracs_ctrl[:2])
    h_mid_ctrl = h_bottom_ctrl + fracs_ctrl[2] / 2
    ax.annotate('Target H', xy=(0, h_mid_ctrl), xytext=(0.8, h_mid_ctrl + 0.15),
                fontsize=9, fontweight='bold', color='#e74c3c',
                arrowprops=dict(arrowstyle='->', color='#e74c3c', lw=1.5))

    ax.set_xticks(x)
    ax.set_xticklabels(bar_labels, fontsize=8)
    ax.set_ylabel('BioAge (fraction)')
    ax.set_title('BioAge Decomposition at Death', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 1.5)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    save(fig, 'stage_16_bioage_composition.png')


# =============================================================================
# STAGE 17: Feedback loops, node clamping, pathway sensitivity
# =============================================================================

def gen_stage_17():
    # Loop breaking (hardcoded)
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    loops = ['PARP-NAD\nVicious Cycle', 'SASP-NAD-Damage\nLoop', 'SASP Paracrine\nLoop']
    gains = [3.0, 0.2, 0.1]
    colors = ['#e74c3c', '#e67e22', '#9b59b6']
    bars = ax.barh(loops, gains, color=colors, edgecolor='black', linewidth=0.5, height=0.6)
    for bar, v in zip(bars, gains):
        ax.text(v + 0.05, bar.get_y() + bar.get_height()/2, f'+{v:.1f}%',
                ha='left', va='center', fontsize=10, fontweight='bold')
    ax.set_xlabel('Lifespan Extension from Breaking Loop (%)')
    ax.set_title('Breaking Feedback Loops\n(untreated control)', fontsize=12, fontweight='bold')
    ax.set_xlim(0, 3.5)
    ax.grid(True, alpha=0.3, axis='x')
    save(fig, 'stage_17_loop_breaking.png')

    # Node clamping (hardcoded)
    fig2, ax2 = plt.subplots(figsize=(5.5, 4.5))
    clamp_data = [
        ('DNA Damage', 68, '#e67e22'),
        ('Methylation', 36, '#27ae60'),
        ('Heteroplasmy', 15, '#e74c3c'),
        ('Senescent Cells', 14, '#9b59b6'),
        ('NAD+', 6, '#3498db'),
        ('SASP', 0, '#95a5a6'),
    ]
    labels_c = [d[0] for d in clamp_data]
    vals_c = [d[1] for d in clamp_data]
    colors_c = [d[2] for d in clamp_data]
    bars = ax2.barh(labels_c, vals_c, color=colors_c, edgecolor='black', linewidth=0.5, height=0.6)
    for bar, v in zip(bars, vals_c):
        ax2.text(v + 1, bar.get_y() + bar.get_height()/2, f'+{v}%',
                 ha='left', va='center', fontsize=10, fontweight='bold')
    ax2.set_xlabel('Lifespan Extension if Clamped (%)')
    ax2.set_title('"What If This Never Deteriorated?"', fontsize=12, fontweight='bold')
    ax2.set_xlim(0, 100)
    ax2.grid(True, alpha=0.3, axis='x')
    save(fig2, 'stage_17_node_clamping.png')

    # Pathway sensitivity (hardcoded)
    fig3, ax3 = plt.subplots(figsize=(6, 6))
    pathways = [
        ('AMPK', 62, True),
        ('Senolytic', 60, False),
        ('Gut microbiome', 22, True),
        ('Mitophagy', 13, False),
        ('mTORC1', 7, True),
        ('Antioxidant', 6, True),
        ('NAD+ precursor', 4, False),
        ('SIRT1 activator', 4, False),
        ('CD38 inhibitor', 3, False),
        ('Anti-inflammatory', 1, True),
        ('AKG', 1, False),
    ]
    p_labels = [p[0] for p in pathways]
    p_vals = [p[1] for p in pathways]
    p_colors = [BLUE if p[2] else '#17a2b8' for p in pathways]

    bars = ax3.barh(p_labels, p_vals, color=p_colors, edgecolor='black', linewidth=0.5, height=0.6)
    for bar, v in zip(bars, p_vals):
        ax3.text(v + 0.5, bar.get_y() + bar.get_height()/2, str(v),
                 ha='left', va='center', fontsize=9, fontweight='bold')
    ax3.invert_yaxis()
    ax3.set_xlabel('Sensitivity Score')
    ax3.set_title('Intervention Point Sensitivity', fontsize=13, fontweight='bold')
    ax3.set_xlim(0, 100)
    ax3.grid(True, alpha=0.3, axis='x')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=BLUE, label='ITP compound targets'),
                       Patch(facecolor='#17a2b8', label='Novel / untested')]
    ax3.legend(handles=legend_elements, fontsize=9, loc='lower right')
    save(fig3, 'stage_17_pathway_sensitivity.png')


# =============================================================================
# GAUGE DRIFT + LOO + SVD (analytical figures)
# =============================================================================

def gen_analytical():
    # Gauge drift ratios (hardcoded)
    fig, ax = plt.subplots(figsize=(12, 5.5))
    compounds_g = ['Canagliflozin', '17\u03b1E2', 'Aspirin', 'Glycine', 'Rapamycin', 'Acarbose']
    ratios = [0.151, 0.168, 0.177, 0.183, 0.191, 0.207]
    mean_val = np.mean(ratios)
    bar_colors = ['#4472c4' if r < mean_val else '#e74c3c' for r in ratios]

    bars = ax.bar(compounds_g, ratios, color=bar_colors, edgecolor='black', linewidth=0.5, width=0.6)
    for bar, v in zip(bars, ratios):
        ax.text(bar.get_x() + bar.get_width()/2, v + 0.002, f'{v:.3f}',
                ha='center', fontsize=10, fontweight='bold')
    ax.axhline(mean_val, color='black', linestyle='--', linewidth=1.5, label=f'Mean = {mean_val:.3f}')
    ax.set_ylabel('Drift / Baseline Ratio')
    ax.set_title(
        'Gauge Symmetry: Drift Ratios by Compound\n(+10% transcription_advantage perturbation)',
        fontsize=13, fontweight='bold'
    )
    ax.set_ylim(0.12, 0.24)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    # Gray arrow at acarbose
    ax.annotate('AMPK saturation\nbuffers acarbose', xy=(5, 0.207), xytext=(4.0, 0.232),
                fontsize=8, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))

    # Yellow text box (centered in x, 3/4 up in y, clear of all bars)
    ax.text(2.5, 0.225, '37% spread confirms\npartial (not full)\ngauge symmetry',
            fontsize=9, fontweight='bold', color='#8B6914', ha='center', va='bottom',
            bbox=dict(facecolor='#FFFFCC', edgecolor='#8B6914', boxstyle='round,pad=0.5'))

    save(fig, 'gauge_drift_ratios.png')

    # LOO Jacobian residuals (hardcoded)
    fig2, ax2 = plt.subplots(figsize=(12, 5.5))
    compounds_loo = ['Rapamycin', 'Acarbose', '17\u03b1E2', 'Canagliflozin', 'Glycine', 'Aspirin']
    capture_m = [95.4, 96.3, 99.9, 95.7, 79.2, 98.5]
    capture_f = [91.2, 89.5, 99.8, 85.1, 81.3, 94.2]
    x = np.arange(len(compounds_loo))
    w = 0.35
    bars_m = ax2.bar(x - w/2, capture_m, w, label='Male', color='#4472c4')
    bars_f = ax2.bar(x + w/2, capture_f, w, label='Female', color='#c0504d')
    for i in range(len(compounds_loo)):
        ax2.text(i - w/2, capture_m[i] + 0.5, f'{capture_m[i]:.1f}%', ha='center',
                 fontsize=9, fontweight='bold', color='#4472c4')
        ax2.text(i + w/2, capture_f[i] + 0.5, f'{capture_f[i]:.1f}%', ha='center',
                 fontsize=9, fontweight='bold', color='#c0504d')
    ax2.set_xticks(x)
    ax2.set_xticklabels(compounds_loo, rotation=30, ha='right')
    ax2.set_ylabel('Jacobian Capture (%)')
    ax2.set_title(
        'Leave-One-Out Vulnerability: Jacobian Analysis\nHigher = more implied by other compounds',
        fontsize=13, fontweight='bold'
    )
    ax2.set_ylim(70, 107)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')

    # "Most vulnerable" annotation at Glycine — in whitespace above bars
    ax2.annotate('Most vulnerable', xy=(4, 81.3), xytext=(4, 105),
                 fontsize=10, fontweight='bold', color='#8B6914', ha='center',
                 arrowprops=dict(arrowstyle='->', color='#8B6914', lw=2))
    ax2.set_ylim(70, 112)

    plt.tight_layout()
    save(fig2, 'loo_jacobian_residuals.png')

    # SVD spectrum (hardcoded)
    fig3, ax3 = plt.subplots(figsize=(12, 5.5))
    sigmas = [251, 68, 50, 30, 22, 15, 10, 6, 3, 1.2]
    labels_svd = [f'\u03c3{i+1}' for i in range(10)]
    subsystem_labels = [
        'Methylation\nclock', 'DNA damage\naccum.', 'mTORC1\ninhibition',
        'AMPK\nsaturation', 'Senescence\nclearance', 'NAD+\nhomeostasis',
        'Heteroplasmy\ndrift', 'SASP\ndynamics', 'Autophagy\nflux', 'Cancer\nsurveillance'
    ]
    bars = ax3.bar(labels_svd, sigmas, color='#4472c4', edgecolor='black', linewidth=0.5, width=0.6)
    ax3.set_yscale('log')
    ax3.axhline(1.2, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    ax3.text(9.5, 1.4, 'Observable threshold', fontsize=9, fontstyle='italic', color='red')

    # Subsystem labels
    for i, (bar, label) in enumerate(zip(bars, subsystem_labels)):
        ax3.text(bar.get_x() + bar.get_width()/2, sigmas[i] * 1.3, label,
                 ha='center', fontsize=7, rotation=0)

    ax3.set_ylim(top=1e3)
    ax3.set_xlabel('Singular Direction')
    ax3.set_ylabel('Singular Value')
    ax3.set_title(
        'SVD Spectrum: 10 Observable Degrees of Freedom\n'
        '(\u03ba = 14.6 full, \u03ba = 210 observable)',
        fontsize=13, fontweight='bold'
    )
    ax3.grid(True, alpha=0.3, axis='y')
    save(fig3, 'stage_A1_svd_spectrum.png')


# =============================================================================
# W_METH SWEEP
# =============================================================================

def gen_w_meth_sweep():
    orig_weights = {k: BIOAGE_PARAMS[k] for k in ['w_meth', 'w_damage', 'w_hetero', 'w_sen']}
    non_meth_total = orig_weights['w_damage'] + orig_weights['w_hetero'] + orig_weights['w_sen']

    w_values = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    mean_errors, max_errors, ra_preds = [], [], []

    for wm in w_values:
        remaining = 1.0 - wm
        scale = remaining / non_meth_total
        BIOAGE_PARAMS['w_meth'] = wm
        BIOAGE_PARAMS['w_damage'] = orig_weights['w_damage'] * scale
        BIOAGE_PARAMS['w_hetero'] = orig_weights['w_hetero'] * scale
        BIOAGE_PARAMS['w_sen'] = orig_weights['w_sen'] * scale

        errors = []
        for sex in ['M', 'F']:
            c = ctrl(sex)
            for cn, val in ITP_VALIDATION.items():
                r = treated(cn, sex)
                e = calculate_lifespan_extension(r, c)
                target = val.target_male if sex == 'M' else val.target_female
                errors.append(abs(e - target))
        mean_errors.append(np.mean(errors))
        max_errors.append(np.max(errors))
        ra_preds.append(run_combination_extension({'rapamycin': 1.0, 'acarbose': 1.0}, sex='M',
                                                    t_max=T_MAX, dt=DT, start_time=9.0/30.0))

    # Restore
    for k, v in orig_weights.items():
        BIOAGE_PARAMS[k] = v

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Left: calibration error
    ax1.plot(w_values, mean_errors, 'b-o', linewidth=2, markersize=8, label='Mean |error|')
    ax1.plot(w_values, max_errors, 'r--s', linewidth=2, markersize=8, label='Max |error|')
    ax1.axhline(y=1.0, color='gray', linestyle=':', alpha=0.7, label='\u00b11% threshold')
    acceptable = [w for w, mx in zip(w_values, max_errors) if mx <= 1.0]
    if acceptable:
        ax1.axvspan(min(acceptable) - 0.025, max(acceptable) + 0.025,
                    alpha=0.15, color='green', label='All 12 targets within \u00b11%')
    ax1.axvspan(0.05, 0.30, alpha=0.15, color='#7B2D8E', label='Clock decomposition range (0.05\u20130.30)')
    current_wm = orig_weights['w_meth']
    ax1.axvline(x=current_wm, color='#1a5276', linestyle='--', alpha=0.7, linewidth=1.5)
    ax1.annotate(f'current\n({current_wm:.2f})', xy=(current_wm, 1.3), fontsize=9, ha='center', color='#1a5276')
    ax1.set_xlabel('Methylation weight (w_meth)')
    ax1.set_ylabel('Calibration error (%)')
    ax1.set_title('BioAge Weight Sensitivity:\nCalibration Quality vs Methylation Weight', fontsize=12, fontweight='bold')
    ax1.set_xlim(0, 0.50)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Right: validation
    ax2.plot(w_values, ra_preds, 'g-D', linewidth=2, markersize=8)
    ax2.axhline(y=34.0, color='red', linestyle='--', linewidth=2, label='Observed: 34%')
    ax2.axhspan(33.0, 35.0, alpha=0.15, color='red', label='\u00b11% band')
    ax2.axvline(x=current_wm, color='#1a5276', linestyle='--', alpha=0.7, linewidth=1.5)
    ax2.annotate(f'current\n({current_wm:.2f})', xy=(current_wm, 29.5), fontsize=9, ha='center', color='#1a5276')
    ax2.set_xlim(0, 0.50)
    ax2.set_xlabel('Methylation weight (w_meth)')
    ax2.set_ylabel('Rapa+Acarbose prediction (%)')
    ax2.set_title('Validation Prediction:\nRapa+Acarbose Lifespan Extension', fontsize=12, fontweight='bold')
    ax2.legend(loc='lower right', fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    save(fig, 'w_meth_sweep.png')


# =============================================================================
# A6b: HETEROPLASMY & SENESCENCE WEIGHT SENSITIVITY
# =============================================================================

def _sweep_weight(weight_key, w_values):
    """Sweep a single BioAge weight, redistributing others proportionally.

    Returns dict with keys: w_values, max_errors, ra_preds, h_pct_at_death,
    mitophagy_marginal.
    """
    orig_weights = {k: BIOAGE_PARAMS[k] for k in ['w_meth', 'w_damage', 'w_hetero', 'w_sen']}
    non_target_keys = [k for k in orig_weights if k != weight_key]
    non_target_total = sum(orig_weights[k] for k in non_target_keys)

    max_errors = []
    ra_preds = []
    h_pct_at_death = []
    mitophagy_marginal = []

    T_MAX_SWEEP = 5.0

    # Build ITP base for whatif (all at 1x calibrated dose)
    base_compounds = [
        ('rapamycin', get_compound('rapamycin')),
        ('acarbose', get_compound('acarbose')),
        ('canagliflozin', get_compound('canagliflozin')),
        ('glycine', get_compound('glycine')),
    ]
    mito_inj = get_compound('urolithin_A')

    for wv in w_values:
        # Set weights
        remaining = 1.0 - wv
        scale = remaining / non_target_total
        BIOAGE_PARAMS[weight_key] = wv
        for k in non_target_keys:
            BIOAGE_PARAMS[k] = orig_weights[k] * scale

        # 1. Calibration error
        errors = []
        for sex in ['M', 'F']:
            c = ctrl(sex)
            for cn, val in ITP_VALIDATION.items():
                r = treated(cn, sex)
                e = calculate_lifespan_extension(r, c)
                target = val.target_male if sex == 'M' else val.target_female
                errors.append(abs(e - target))
        max_errors.append(np.max(errors))

        # 2. Rapa+Acarb validation
        ra = run_combination_extension({'rapamycin': 1.0, 'acarbose': 1.0}, sex='M',
                                       t_max=T_MAX, dt=DT, start_time=9.0/30.0)
        ra_preds.append(ra)

        # 3. H fraction of BioAge at death for ITP-base-treated male
        ctrl_m = run_control(sex='M', t_max=T_MAX_SWEEP, dt=DT_FINE)
        base_m_inj = _build_whatif_interventions(base_compounds, None, 'M')
        base_result = simulate(interventions=base_m_inj, sex='M', t_max=T_MAX_SWEEP, dt=DT_FINE)
        idx = min(int(round(base_result.t_death / DT_FINE)), len(base_result.t) - 1)
        meth_c = BIOAGE_PARAMS['w_meth'] * base_result.Methylation[idx] / BIOAGE_PARAMS['meth_norm']
        dmg_c = BIOAGE_PARAMS['w_damage'] * base_result.DNA_damage[idx] / BIOAGE_PARAMS['damage_norm']
        h_c = BIOAGE_PARAMS['w_hetero'] * base_result.Heteroplasmy[idx] / BIOAGE_PARAMS['H_norm']
        sen_c = BIOAGE_PARAMS['w_sen'] * base_result.SenCells[idx] / BIOAGE_PARAMS['sen_norm']
        total_ba = meth_c + dmg_c + h_c + sen_c
        h_pct_at_death.append(h_c / total_ba * 100 if total_ba > 0 else 0)

        # 4. Mitophagy marginal addition over ITP base
        base_ext = calculate_lifespan_extension(base_result, ctrl_m)
        mito_inj_m = _build_whatif_interventions(base_compounds, mito_inj, 'M')
        mito_result = simulate(interventions=mito_inj_m, sex='M', t_max=T_MAX_SWEEP, dt=DT_FINE)
        mito_ext = calculate_lifespan_extension(mito_result, ctrl_m)
        mitophagy_marginal.append(mito_ext - base_ext)

        print(f'  {weight_key}={wv:.2f}: max_err={max_errors[-1]:.2f}%, '
              f'RA={ra_preds[-1]:.1f}%, H%={h_pct_at_death[-1]:.0f}%, '
              f'mito_marginal=+{mitophagy_marginal[-1]:.1f}%')

    # Restore
    for k, v in orig_weights.items():
        BIOAGE_PARAMS[k] = v

    return {
        'w_values': w_values,
        'max_errors': max_errors,
        'ra_preds': ra_preds,
        'h_pct_at_death': h_pct_at_death,
        'mitophagy_marginal': mitophagy_marginal,
        'current': orig_weights[weight_key],
    }


def _nice_ceil(values, headroom=1.3):
    """Round up to a nice number ~30% above the max value."""
    peak = max(values) * headroom
    if peak <= 0:
        return 1.0
    mag = 10 ** int(np.floor(np.log10(peak)))
    return np.ceil(peak / mag) * mag


def gen_stage_A6b_weight_sensitivity():
    """Sensitivity of key conclusions to heteroplasmy and senescence weights."""

    # Sweep w_hetero
    print('  Sweeping w_hetero...')
    wh_values = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
    h_data = _sweep_weight('w_hetero', wh_values)

    # Sweep w_sen
    print('  Sweeping w_sen...')
    ws_values = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    s_data = _sweep_weight('w_sen', ws_values)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('BioAge Weight Sensitivity: Heteroplasmy & Senescence',
                 fontsize=14, fontweight='bold', y=0.98)

    # ── TOP ROW: w_hetero sweep ──
    ax = axes[0, 0]
    ax.plot(h_data['w_values'], h_data['max_errors'], 'r-s', linewidth=2, markersize=7)
    ax.axhline(y=1.0, color='gray', linestyle=':', alpha=0.7, label='±1% threshold')
    ax.axvline(x=h_data['current'], color='#1a5276', linestyle='--', alpha=0.7, linewidth=1.5)
    ax.annotate(f'current\n({h_data["current"]:.2f})',
                xy=(h_data['current'], 0.5), fontsize=8, ha='center', color='#1a5276')
    acceptable = [w for w, e in zip(h_data['w_values'], h_data['max_errors']) if e <= 1.0]
    if acceptable:
        ax.axvspan(min(acceptable) - 0.025, max(acceptable) + 0.025,
                   alpha=0.12, color='green', label='All 12 within ±1%')
    ax.set_xlabel('w_hetero')
    ax.set_ylabel('Max calibration error (%)')
    ax.set_title('Calibration vs w_hetero', fontweight='bold')
    ax.legend(fontsize=7, loc='upper left')
    ax.set_ylim(0, _nice_ceil(h_data['max_errors']))
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(h_data['w_values'], h_data['h_pct_at_death'], '#e74c3c', linewidth=2.5, marker='o', markersize=7)
    ax.axvline(x=h_data['current'], color='#1a5276', linestyle='--', alpha=0.7, linewidth=1.5)
    ax.set_xlabel('w_hetero')
    ax.set_ylabel('H as % of BioAge at death')
    ax.set_title('Bottleneck conclusion\n(ITP-base treated, male)', fontweight='bold')
    ax.set_ylim(0, _nice_ceil(h_data['h_pct_at_death']))
    ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    ax.plot(h_data['w_values'], h_data['mitophagy_marginal'], '#27ae60', linewidth=2.5, marker='D', markersize=7)
    ax.axvline(x=h_data['current'], color='#1a5276', linestyle='--', alpha=0.7, linewidth=1.5)
    ax.set_xlabel('w_hetero')
    ax.set_ylabel('Mitophagy marginal addition (%)')
    ax.set_title('Urolithin A value\n(marginal over ITP base, male)', fontweight='bold')
    ax.set_ylim(0, _nice_ceil(h_data['mitophagy_marginal']))
    ax.grid(True, alpha=0.3)

    # ── BOTTOM ROW: w_sen sweep ──
    ax = axes[1, 0]
    ax.plot(s_data['w_values'], s_data['max_errors'], 'r-s', linewidth=2, markersize=7)
    ax.axhline(y=1.0, color='gray', linestyle=':', alpha=0.7, label='±1% threshold')
    ax.axvline(x=s_data['current'], color='#1a5276', linestyle='--', alpha=0.7, linewidth=1.5)
    ax.annotate(f'current\n({s_data["current"]:.2f})',
                xy=(s_data['current'], 0.5), fontsize=8, ha='center', color='#1a5276')
    acceptable = [w for w, e in zip(s_data['w_values'], s_data['max_errors']) if e <= 1.0]
    if acceptable:
        ax.axvspan(min(acceptable) - 0.025, max(acceptable) + 0.025,
                   alpha=0.12, color='green', label='All 12 within ±1%')
    ax.set_xlabel('w_sen')
    ax.set_ylabel('Max calibration error (%)')
    ax.set_title('Calibration vs w_sen', fontweight='bold')
    ax.legend(fontsize=7, loc='upper left')
    ax.set_ylim(0, _nice_ceil(s_data['max_errors']))
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(s_data['w_values'], s_data['h_pct_at_death'], '#e74c3c', linewidth=2.5, marker='o', markersize=7)
    ax.axvline(x=s_data['current'], color='#1a5276', linestyle='--', alpha=0.7, linewidth=1.5)
    ax.set_xlabel('w_sen')
    ax.set_ylabel('H as % of BioAge at death')
    ax.set_title('Bottleneck conclusion\n(ITP-base treated, male)', fontweight='bold')
    ax.set_ylim(0, _nice_ceil(s_data['h_pct_at_death']))
    ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    ax.plot(s_data['w_values'], s_data['mitophagy_marginal'], '#27ae60', linewidth=2.5, marker='D', markersize=7)
    ax.axvline(x=s_data['current'], color='#1a5276', linestyle='--', alpha=0.7, linewidth=1.5)
    ax.set_xlabel('w_sen')
    ax.set_ylabel('Mitophagy marginal addition (%)')
    ax.set_title('Urolithin A value\n(marginal over ITP base, male)', fontweight='bold')
    ax.set_ylim(0, _nice_ceil(s_data['mitophagy_marginal']))
    ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, 'stage_A6b_weight_sensitivity.png')


# =============================================================================
# A10: PAIRWISE SYNERGY & NETWORK ANALYSIS
# =============================================================================

from matplotlib.colors import LinearSegmentedColormap
_SYNERGY_CMAP = LinearSegmentedColormap.from_list(
    'blue_white_green',
    ['#2166ac', '#92c5de', '#f7f7f7', '#78c679', '#1a9641'],
)

def gen_stage_A10():
    """Pairwise synergy matrices + network analysis (Taguchi 8 compounds)."""

    # Taguchi 8: 4 ITP unisex + 4 novel (all at 1x calibrated dose)
    TAGUCHI_COMPOUNDS = ['rapamycin', 'acarbose', 'canagliflozin', 'glycine',
                         'NMN', 'dasatinib_quercetin', 'urolithin_A']
    CD38I_KEY = 'cd38i'
    CD38I_INJ = {'cd38_inhibitor': 0.50}
    ALL_LABELS = ['Rapa', 'Acarb', 'Cana', 'Glycine',
                  'NMN', 'D+Q', 'UroA', 'CD38i']
    N = len(ALL_LABELS)  # 8
    NOVEL_KEYS = ['NMN', 'CD38i', 'D+Q', 'UroA']
    NOVEL_INJ = [
        get_compound('NMN'),
        CD38I_INJ,
        get_compound('dasatinib_quercetin'),
        get_compound('urolithin_A'),
    ]

    def _get_compound_interventions(compound, sex):
        """Get processed interventions for a named compound (dose-scaled, sex-modified)."""
        base = get_compound(compound)
        eff = effective_dose(compound, 1.0)
        scaled = {k: v * eff for k, v in base.items()}
        scaled = apply_sex_modifier(compound, scaled, sex)
        scaled['start_time'] = get_itp_start_time(compound)
        return scaled

    def _merge_interventions(*inj_dicts):
        """Merge multiple intervention dicts, summing shared keys. Earliest start_time wins."""
        merged = {}
        earliest = 1.0
        for d in inj_dicts:
            for k, v in d.items():
                if k == 'start_time':
                    if v < earliest:
                        earliest = v
                else:
                    merged[k] = merged.get(k, 0) + v
        merged['start_time'] = earliest
        return merged

    def _get_single_inj(idx, sex):
        """Get interventions for compound at index idx."""
        if idx < len(TAGUCHI_COMPOUNDS):
            return _get_compound_interventions(TAGUCHI_COMPOUNDS[idx], sex)
        else:
            # CD38i (index 7)
            novel = dict(CD38I_INJ)
            novel['start_time'] = 0.30
            return novel

    def _compute_ext(inj, sex, ctrl_cache):
        """Compute extension % for given interventions."""
        result = simulate(interventions=inj, sex=sex, t_max=T_MAX, dt=DT_FAST)
        return calculate_lifespan_extension(result, ctrl_cache[sex])

    # Pre-compute controls
    ctrl_cache = {
        'M': run_control(sex='M', t_max=T_MAX, dt=DT_FAST),
        'F': run_control(sex='F', t_max=T_MAX, dt=DT_FAST),
    }

    # --- Figure 1: Full 10x10 pairwise synergy matrix ---
    print('  Computing pairwise synergy matrix...')
    synergy = {}  # synergy[(sex, i, j)] = Δ_ij
    singles = {}  # singles[(sex, i)] = ext(i)

    for sex in ['M', 'F']:
        # Singles
        for i in range(N):
            inj_i = _get_single_inj(i, sex)
            singles[(sex, i)] = _compute_ext(inj_i, sex, ctrl_cache)
        # Pairs
        for i in range(N):
            for j in range(i + 1, N):
                inj_i = _get_single_inj(i, sex)
                inj_j = _get_single_inj(j, sex)
                inj_pair = _merge_interventions(inj_i, inj_j)
                ext_pair = _compute_ext(inj_pair, sex, ctrl_cache)
                synergy[(sex, i, j)] = ext_pair - singles[(sex, i)] - singles[(sex, j)]

    # Build matrices for plotting
    fig1, axes = plt.subplots(1, 2, figsize=(16, 7.5))
    for ax, sex, title_sex in zip(axes, ['M', 'F'], ['Male', 'Female']):
        mat = np.full((N, N), np.nan)
        # Diagonal: single-compound extensions
        for i in range(N):
            mat[i, i] = singles[(sex, i)]
        # Lower triangle: synergy values
        for i in range(N):
            for j in range(i + 1, N):
                mat[j, i] = synergy[(sex, i, j)]  # lower triangle: row > col

        # Plot heatmap
        vmax = np.nanmax(np.abs([v for k, v in synergy.items() if k[0] == sex]))
        vmax = max(vmax, 1.0)
        im = ax.imshow(mat, cmap=_SYNERGY_CMAP, vmin=-vmax, vmax=vmax, aspect='equal')

        # Annotate cells
        for i in range(N):
            for j in range(N):
                if not np.isnan(mat[i, j]):
                    val = mat[i, j]
                    color = 'white' if abs(val) > vmax * 0.6 else 'black'
                    fontsize = 7 if i == j else 8
                    fmt = f'{val:+.1f}' if i != j else f'{val:.1f}%'
                    ax.text(j, i, fmt, ha='center', va='center',
                            fontsize=fontsize, color=color, fontweight='bold')

        # Mask upper triangle
        for i in range(N):
            for j in range(i + 1, N):
                ax.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1,
                             facecolor='#f0f0f0', edgecolor='#ddd', linewidth=0.5))

        ax.set_xticks(range(N))
        ax.set_yticks(range(N))
        ax.set_xticklabels(ALL_LABELS, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(ALL_LABELS, fontsize=8)
        ax.set_title(f'{title_sex} Pairwise Synergy (\u0394 pp)', fontsize=12, fontweight='bold')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Synergy (pp)')

        # Divider line between ITP unisex and novel
        ax.axhline(3.5, color='black', linewidth=1.5, linestyle='--', alpha=0.5)
        ax.axvline(3.5, color='black', linewidth=1.5, linestyle='--', alpha=0.5)

    fig1.suptitle('Pairwise Synergy: \u0394\u1d62\u2c7c = ext(i+j) \u2212 ext(i) \u2212 ext(j)\n'
                  'Diagonal = single-compound extension (%). Off-diagonal = synergy (pp).',
                  fontsize=11, y=0.02)
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    save(fig1, 'stage_A10_synergy_full.png')

    # --- Figure 2: Conditional 4x4 novel synergy on ITP base ---
    print('  Computing conditional synergy (novel pairs on ITP base)...')
    T_MAX_WHATIF = 5.0
    ctrl_whatif = {
        'M': run_control(sex='M', t_max=T_MAX_WHATIF, dt=DT_FAST),
        'F': run_control(sex='F', t_max=T_MAX_WHATIF, dt=DT_FAST),
    }

    def _base_inj(sex):
        """Build ITP unisex base cocktail interventions (4 compounds at 1x)."""
        base_compounds = [
            ('rapamycin', get_compound('rapamycin')),
            ('acarbose', get_compound('acarbose')),
            ('canagliflozin', get_compound('canagliflozin')),
            ('glycine', get_compound('glycine')),
        ]
        return _build_whatif_interventions(base_compounds, None, sex)

    cond_synergy = {}  # cond_synergy[(sex, i, j)] for novel indices 0-3
    base_ext = {}
    base_plus_single = {}  # base_plus_single[(sex, i)]

    for sex in ['M', 'F']:
        bi = _base_inj(sex)
        base_result = simulate(interventions=bi, sex=sex, t_max=T_MAX_WHATIF, dt=DT_FAST)
        base_ext[sex] = calculate_lifespan_extension(base_result, ctrl_whatif[sex])

        # Base + each novel target
        for i in range(4):
            inj = _build_whatif_interventions(
                [('rapamycin', get_compound('rapamycin')),
                 ('acarbose', get_compound('acarbose')),
                 ('canagliflozin', get_compound('canagliflozin')),
                 ('glycine', get_compound('glycine'))],
                NOVEL_INJ[i], sex
            )
            r = simulate(interventions=inj, sex=sex, t_max=T_MAX_WHATIF, dt=DT_FAST)
            base_plus_single[(sex, i)] = calculate_lifespan_extension(r, ctrl_whatif[sex])

        # Base + novel pairs
        for i in range(4):
            for j in range(i + 1, 4):
                merged_novel = {}
                for k, v in NOVEL_INJ[i].items():
                    merged_novel[k] = merged_novel.get(k, 0) + v
                for k, v in NOVEL_INJ[j].items():
                    merged_novel[k] = merged_novel.get(k, 0) + v
                inj = _build_whatif_interventions(
                    [('rapamycin', get_compound('rapamycin')),
                     ('acarbose', get_compound('acarbose')),
                     ('canagliflozin', get_compound('canagliflozin')),
                     ('glycine', get_compound('glycine'))],
                    merged_novel, sex
                )
                r = simulate(interventions=inj, sex=sex, t_max=T_MAX_WHATIF, dt=DT_FAST)
                ext_pair = calculate_lifespan_extension(r, ctrl_whatif[sex])
                # Δ_ij|base = ext(base+i+j) - ext(base+i) - ext(base+j) + ext(base)
                cond_synergy[(sex, i, j)] = (ext_pair
                                             - base_plus_single[(sex, i)]
                                             - base_plus_single[(sex, j)]
                                             + base_ext[sex])

    # Plot conditional synergy
    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, sex, title_sex in zip(axes2, ['M', 'F'], ['Male', 'Female']):
        mat = np.full((4, 4), np.nan)
        # Diagonal: marginal addition of each novel target to base
        for i in range(4):
            mat[i, i] = base_plus_single[(sex, i)] - base_ext[sex]
        # Lower triangle: conditional synergy
        for i in range(4):
            for j in range(i + 1, 4):
                mat[j, i] = cond_synergy[(sex, i, j)]

        vmax = max(1.0, np.nanmax(np.abs([v for k, v in cond_synergy.items() if k[0] == sex])))
        im = ax.imshow(mat, cmap=_SYNERGY_CMAP, vmin=-vmax, vmax=vmax, aspect='equal')

        for i in range(4):
            for j in range(4):
                if not np.isnan(mat[i, j]):
                    val = mat[i, j]
                    color = 'white' if abs(val) > vmax * 0.6 else 'black'
                    fmt = f'{val:+.1f}' if i != j else f'+{val:.1f}%'
                    ax.text(j, i, fmt, ha='center', va='center',
                            fontsize=10, color=color, fontweight='bold')

        for i in range(4):
            for j in range(i + 1, 4):
                ax.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1,
                             facecolor='#f0f0f0', edgecolor='#ddd', linewidth=0.5))

        ax.set_xticks(range(4))
        ax.set_yticks(range(4))
        ax.set_xticklabels(NOVEL_KEYS, fontsize=9)
        ax.set_yticklabels(NOVEL_KEYS, fontsize=9)
        ax.set_title(f'{title_sex}: Conditional Synergy\n(given ITP base cocktail)',
                     fontsize=11, fontweight='bold')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Synergy (pp)')

    fig2.suptitle('\u0394\u1d62\u2c7c|base = ext(base+i+j) \u2212 ext(base+i) \u2212 ext(base+j) + ext(base)\n'
                  'Diagonal = marginal addition to ITP cocktail (pp)',
                  fontsize=10, y=0.02)
    plt.tight_layout(rect=[0, 0.06, 1, 1])
    save(fig2, 'stage_A10_synergy_conditional.png')

    # --- Figure 3: Network adjacency & centrality ---
    print('  Computing network adjacency and centrality...')

    # ODE state variables (nodes in the dependency graph)
    NODE_NAMES = ['AMPK', 'TSC2', 'mTORC1', 'NAD', 'SIRT1', 'FOXO3',
                  'Autophagy', 'DNA_Dmg', 'Meth', 'Hetero', 'UPRmt',
                  'SenCells', 'SASP', 'Immune', 'BioAge']
    n_nodes = len(NODE_NAMES)

    # Adjacency matrix: A[i,j] = influence strength of node j on node i
    # Derived from ODE coupling structure in model.py
    # Positive = activating, negative = inhibiting
    A = np.zeros((n_nodes, n_nodes))
    idx = {name: i for i, name in enumerate(NODE_NAMES)}

    # AMPK -> TSC2 (activation, Hill function)
    A[idx['TSC2'], idx['AMPK']] = 0.8
    # TSC2 -> mTORC1 (inhibition)
    A[idx['mTORC1'], idx['TSC2']] = -0.9
    # mTORC1 -> Autophagy (inhibition, so mTORC1 low = autophagy high)
    A[idx['Autophagy'], idx['mTORC1']] = -0.6
    # mTORC1 -> DNA_Dmg (inhibition of damage via proteostasis)
    A[idx['DNA_Dmg'], idx['mTORC1']] = -0.3
    # mTORC1 -> SenCells (inhibition)
    A[idx['SenCells'], idx['mTORC1']] = -0.3
    # NAD -> SIRT1 (activation)
    A[idx['SIRT1'], idx['NAD']] = 0.9
    # SIRT1 -> FOXO3 (activation)
    A[idx['FOXO3'], idx['SIRT1']] = 0.7
    # FOXO3 -> Autophagy (activation)
    A[idx['Autophagy'], idx['FOXO3']] = 0.5
    # FOXO3 -> Hetero (selection against mutants)
    A[idx['Hetero'], idx['FOXO3']] = -0.4
    # SIRT1 -> DNA_Dmg (repair boost)
    A[idx['DNA_Dmg'], idx['SIRT1']] = -0.5
    # SIRT1 -> Meth (EZH2 inhibition -> less methylation)
    A[idx['Meth'], idx['SIRT1']] = -0.4
    # SIRT1 -> SenCells (reduction)
    A[idx['SenCells'], idx['SIRT1']] = -0.3
    # DNA_Dmg -> Meth (damage accelerates epigenetic drift)
    A[idx['Meth'], idx['DNA_Dmg']] = 0.5
    # DNA_Dmg -> SenCells (damage-induced senescence)
    A[idx['SenCells'], idx['DNA_Dmg']] = 0.6
    # SenCells -> SASP (senescent cells produce SASP)
    A[idx['SASP'], idx['SenCells']] = 0.8
    # SASP -> NAD (CD38 upregulation consumes NAD)
    A[idx['NAD'], idx['SASP']] = -0.6
    # SASP -> DNA_Dmg (ROS production)
    A[idx['DNA_Dmg'], idx['SASP']] = 0.4
    # SASP -> SenCells (paracrine senescence)
    A[idx['SenCells'], idx['SASP']] = 0.5
    # SASP -> Meth (methylation acceleration)
    A[idx['Meth'], idx['SASP']] = 0.4
    # SASP -> Immune (suppression)
    A[idx['Immune'], idx['SASP']] = -0.5
    # Hetero -> NAD (OXPHOS efficiency loss)
    A[idx['NAD'], idx['Hetero']] = -0.5
    # Hetero -> DNA_Dmg (ROS from dysfunctional mitochondria)
    A[idx['DNA_Dmg'], idx['Hetero']] = 0.4
    # Hetero -> UPRmt (stress signal)
    A[idx['UPRmt'], idx['Hetero']] = 0.7
    # UPRmt -> Autophagy (mitophagy boost)
    A[idx['Autophagy'], idx['UPRmt']] = 0.3
    # Autophagy -> SenCells (clearance)
    A[idx['SenCells'], idx['Autophagy']] = -0.4
    # AMPK -> SenCells (reduction)
    A[idx['SenCells'], idx['AMPK']] = -0.3
    # AMPK -> NAD (NAMPT boost)
    A[idx['NAD'], idx['AMPK']] = 0.4
    # AMPK -> Meth (EZH2 inhibition)
    A[idx['Meth'], idx['AMPK']] = -0.3
    # BioAge composite inputs (from BIOAGE_PARAMS)
    A[idx['BioAge'], idx['Meth']] = 0.30
    A[idx['BioAge'], idx['DNA_Dmg']] = 0.30
    A[idx['BioAge'], idx['Hetero']] = 0.25
    A[idx['BioAge'], idx['SenCells']] = 0.15

    # Compute graph metrics
    abs_A = np.abs(A)
    in_degree = abs_A.sum(axis=1)   # sum of influences received
    out_degree = abs_A.sum(axis=0)  # sum of influences sent

    # Betweenness centrality via Floyd-Warshall on distance matrix
    # Distance = 1/|weight| for nonzero edges, inf otherwise
    dist = np.full((n_nodes, n_nodes), np.inf)
    np.fill_diagonal(dist, 0)
    nxt = np.full((n_nodes, n_nodes), -1, dtype=int)
    for i in range(n_nodes):
        for j in range(n_nodes):
            if abs_A[i, j] > 0:
                # A[i,j] = influence of j on i, so edge goes j → i
                dist[j, i] = 1.0 / abs_A[i, j]
                nxt[j, i] = i

    # Floyd-Warshall
    for k in range(n_nodes):
        for i in range(n_nodes):
            for j in range(n_nodes):
                if dist[i, k] + dist[k, j] < dist[i, j]:
                    dist[i, j] = dist[i, k] + dist[k, j]
                    nxt[i, j] = nxt[i, k]

    # Betweenness: count how often each node appears on shortest paths
    betweenness = np.zeros(n_nodes)
    for s in range(n_nodes):
        for t in range(n_nodes):
            if s == t or dist[s, t] == np.inf:
                continue
            # Trace path
            path = [s]
            cur = s
            while cur != t:
                cur = nxt[cur, t]
                if cur == -1:
                    break
                path.append(cur)
            # Count intermediaries
            for node in path[1:-1]:
                betweenness[node] += 1
    # Normalize
    if betweenness.max() > 0:
        betweenness = betweenness / betweenness.max()

    # Intervention entry points and shortest path to BioAge
    IV_POINTS = {
        'AMPK\n(Rapa,Acarb,Cana)': idx['AMPK'],
        'NAD\n(NMN,CD38i)': idx['NAD'],
        'SenCells\n(D+Q)': idx['SenCells'],
        'Hetero\n(Mitophagy)': idx['Hetero'],
        'Meth\n(17\u03b1E2)': idx['Meth'],
    }
    bioage_idx = idx['BioAge']
    path_lengths = {}
    for label, src in IV_POINTS.items():
        d = dist[src, bioage_idx]
        path_lengths[label] = d if d < np.inf else 999

    # --- Figure 3a: Adjacency heatmap ---
    fig3a, ax_adj = plt.subplots(figsize=(10, 8.5))
    vmax_adj = np.max(np.abs(A))
    im_adj = ax_adj.imshow(A, cmap=_SYNERGY_CMAP, vmin=-vmax_adj, vmax=vmax_adj, aspect='equal')
    ax_adj.set_xticks(range(n_nodes))
    ax_adj.set_yticks(range(n_nodes))
    ax_adj.set_xticklabels(NODE_NAMES, rotation=45, ha='right', fontsize=10)
    ax_adj.set_yticklabels(NODE_NAMES, fontsize=10)
    ax_adj.set_title('Weighted Adjacency Matrix\nA[i,j] = influence of column j on row i',
                     fontsize=13, fontweight='bold')
    plt.colorbar(im_adj, ax=ax_adj, fraction=0.046, pad=0.04,
                 label='Coupling strength (green = activating, blue = inhibiting)')
    for i in range(n_nodes):
        for j in range(n_nodes):
            if A[i, j] != 0:
                ax_adj.text(j, i, f'{A[i,j]:+.1f}', ha='center', va='center',
                           fontsize=8, color='white' if abs(A[i,j]) > 0.5 else 'black',
                           fontweight='bold')
    plt.tight_layout()
    save(fig3a, 'stage_A10_adjacency.png')

    # --- Figure 3b: In/Out degree bar chart ---
    fig3b, ax_deg = plt.subplots(figsize=(10, 7))
    x_deg = np.arange(n_nodes)
    w_deg = 0.35
    ax_deg.barh(x_deg - w_deg/2, in_degree, w_deg, color='#e74c3c',
                edgecolor='black', linewidth=0.5, label='In-degree (influence received)')
    ax_deg.barh(x_deg + w_deg/2, out_degree, w_deg, color='#3498db',
                edgecolor='black', linewidth=0.5, label='Out-degree (influence broadcast)')
    ax_deg.set_yticks(x_deg)
    ax_deg.set_yticklabels(NODE_NAMES, fontsize=11)
    ax_deg.set_xlabel('Weighted Degree', fontsize=11)
    ax_deg.set_title('Node Degree: Influence Received vs Broadcast',
                     fontsize=13, fontweight='bold')
    ax_deg.legend(fontsize=10, loc='lower right')
    ax_deg.invert_yaxis()
    ax_deg.grid(True, alpha=0.3, axis='x')
    for i in range(n_nodes):
        if in_degree[i] > 0.1:
            ax_deg.text(in_degree[i] + 0.03, i - w_deg/2, f'{in_degree[i]:.1f}',
                       va='center', fontsize=9, color='#e74c3c', fontweight='bold')
        if out_degree[i] > 0.1:
            ax_deg.text(out_degree[i] + 0.03, i + w_deg/2, f'{out_degree[i]:.1f}',
                       va='center', fontsize=9, color='#3498db', fontweight='bold')
    plt.tight_layout()
    save(fig3b, 'stage_A10_degree.png')

    # --- Figure 3c: Betweenness centrality ---
    fig3c, ax_btw = plt.subplots(figsize=(10, 7))
    sort_idx = np.argsort(betweenness)[::-1]
    bars_btw = ax_btw.barh(range(n_nodes), betweenness[sort_idx],
                           color='#2ecc71', edgecolor='black', linewidth=0.5)
    ax_btw.set_yticks(range(n_nodes))
    ax_btw.set_yticklabels([NODE_NAMES[i] for i in sort_idx], fontsize=11)
    ax_btw.set_xlabel('Normalized Betweenness Centrality', fontsize=11)
    ax_btw.set_title('Betweenness Centrality: Bottleneck Nodes',
                     fontsize=13, fontweight='bold')
    ax_btw.invert_yaxis()
    ax_btw.grid(True, alpha=0.3, axis='x')
    for i, si in enumerate(sort_idx):
        if betweenness[si] > 0.01:
            ax_btw.text(betweenness[si] + 0.01, i, f'{betweenness[si]:.2f}',
                       va='center', fontsize=10, fontweight='bold')
    plt.tight_layout()
    save(fig3c, 'stage_A10_betweenness.png')

    # --- Figure 3d: Shortest path from intervention points to BioAge ---
    fig3d, ax_path = plt.subplots(figsize=(10, 5.5))
    pl_labels = list(path_lengths.keys())
    pl_vals = list(path_lengths.values())
    colors_path = ['#e74c3c' if v < 3 else '#f39c12' if v < 5 else '#3498db' for v in pl_vals]
    ax_path.barh(range(len(pl_labels)), pl_vals, color=colors_path,
                 edgecolor='black', linewidth=0.5, height=0.6)
    ax_path.set_yticks(range(len(pl_labels)))
    ax_path.set_yticklabels(pl_labels, fontsize=12)
    ax_path.set_xlabel('Shortest Path Distance (1/|weight|)', fontsize=11)
    ax_path.set_title('Distance from Intervention Point \u2192 BioAge',
                      fontsize=13, fontweight='bold')
    ax_path.invert_yaxis()
    ax_path.grid(True, alpha=0.3, axis='x')
    for i, v in enumerate(pl_vals):
        ax_path.text(v + 0.1, i, f'{v:.1f}', va='center', fontsize=11, fontweight='bold')
    plt.tight_layout()
    save(fig3d, 'stage_A10_paths.png')


# =============================================================================
# A5 SENSITIVITY: RAPA+NMN PARAMETER TORNADO
# =============================================================================

def gen_stage_A5_sensitivity():
    """Parameter sensitivity tornado for Rapa+NMN prediction."""
    print('  Computing Rapa+NMN parameter sensitivity...')

    # Parameters to perturb: (dict, key, display_label)
    PARAMS_TO_SWEEP = [
        (NAD_PARAMS, 'k_nmn_synthesis', 'k_nmn_synthesis (NMN \u2192 NAD)'),
        (MTORC1_PARAMS, 'k_mtorc1_autophagy', 'k_mtorc1_autophagy (mTOR \u2192 autophagy)'),
        (SIRTUIN_PARAMS, 'Km_nad_sirt1', 'Km_nad_sirt1 (NAD affinity for SIRT1)'),
        (AUTOPHAGY_PARAMS, 'Km_autophagy', 'Km_autophagy (autophagy half-max)'),
        (TSC2_PARAMS, 'k_ampk_tsc2_max', 'k_ampk_tsc2_max (AMPK \u2192 TSC2)'),
        (DNA_REPAIR_PARAMS, 'k_ber_nad_dependence', 'k_ber_nad_dependence (NAD repair coupling)'),
        (BIOAGE_PARAMS, 'w_meth', 'w_meth (methylation weight)'),
        (BIOAGE_PARAMS, 'w_damage', 'w_damage (DNA damage weight)'),
        (NAD_PARAMS, 'k_cd38_base', 'k_cd38_base (CD38 NAD consumption)'),
        (MTORC1_PARAMS, 'k_mtorc1_senescence', 'k_mtorc1_senescence (mTOR \u2192 senescence)'),
    ]

    def _run_rapa_nmn(sex):
        """Run rapa alone, NMN alone, and rapa+NMN. Return (ext_combo, ext_rapa, ext_nmn)."""
        c = run_control(sex=sex, t_max=T_MAX, dt=DT_FAST)
        # Rapamycin
        rapa_inj = get_compound('rapamycin')
        rapa_inj = apply_sex_modifier('rapamycin', rapa_inj, sex)
        rapa_inj['start_time'] = get_itp_start_time('rapamycin')
        r_rapa = simulate(interventions=rapa_inj, sex=sex, t_max=T_MAX, dt=DT_FAST)
        ext_rapa = calculate_lifespan_extension(r_rapa, c)
        # NMN
        nmn_inj = dict(get_compound('NMN'))
        nmn_inj['start_time'] = get_itp_start_time('rapamycin')  # same start
        r_nmn = simulate(interventions=nmn_inj, sex=sex, t_max=T_MAX, dt=DT_FAST)
        ext_nmn = calculate_lifespan_extension(r_nmn, c)
        # Combo
        combo_inj = {}
        for k, v in rapa_inj.items():
            if k != 'start_time':
                combo_inj[k] = v
        for k, v in nmn_inj.items():
            if k != 'start_time':
                combo_inj[k] = combo_inj.get(k, 0) + v
        combo_inj['start_time'] = rapa_inj['start_time']
        r_combo = simulate(interventions=combo_inj, sex=sex, t_max=T_MAX, dt=DT_FAST)
        ext_combo = calculate_lifespan_extension(r_combo, c)
        return ext_combo, ext_rapa, ext_nmn

    # Baseline
    base_combo_m, base_rapa_m, base_nmn_m = _run_rapa_nmn('M')
    base_super_m = base_combo_m - base_rapa_m - base_nmn_m

    # Perturbation sweep
    PERTURB = 0.20  # ±20%
    results = []  # (label, combo_lo, combo_hi, super_lo, super_hi)

    for param_dict, key, label in PARAMS_TO_SWEEP:
        orig = param_dict[key]
        deltas_combo = []
        deltas_super = []
        for factor in [1 - PERTURB, 1 + PERTURB]:
            param_dict[key] = orig * factor
            combo, rapa, nmn = _run_rapa_nmn('M')
            superadd = combo - rapa - nmn
            deltas_combo.append(combo - base_combo_m)
            deltas_super.append(superadd - base_super_m)
            param_dict[key] = orig  # restore immediately

        results.append((label, deltas_combo[0], deltas_combo[1],
                        deltas_super[0], deltas_super[1]))

    # Sort by max absolute combo impact
    results.sort(key=lambda r: max(abs(r[1]), abs(r[2])))

    # Plot tornado with baseline reference bars
    # Each tornado rides on a dashed line at the baseline value,
    # so the reader sees perturbation size relative to the effect size.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 8))
    n_params = len(results)
    # Extra row at top for baseline reference bar
    y_ref = n_params  # top row index
    y_pos = np.arange(n_params)

    def _annotate_tornado(ax, baseline, i, lo, hi):
        """Place lo/hi annotations, staggering vertically if they'd overlap."""
        span = abs(hi - lo)
        tight = span < 0.6  # bars too close together for side-by-side labels
        for val, is_hi in [(hi, True), (lo, False)]:
            if abs(val) < 0.05:
                continue
            xpos = baseline + val
            if tight:
                # Stack: hi above center, lo below
                y_off = 0.18 if is_hi else -0.18
                # Both labels go to the right of the rightmost bar end
                x_anchor = baseline + max(abs(lo), abs(hi)) + 0.2
                ax.text(x_anchor, i + y_off, f'{val:+.1f}',
                        va='center', ha='left',
                        fontsize=8, fontweight='bold',
                        color='#2ecc71' if val > 0 else '#3498db')
            else:
                nudge = 0.15 if val > 0 else -0.15
                ax.text(xpos + nudge, i, f'{val:+.1f}',
                        va='center', ha='left' if val > 0 else 'right',
                        fontsize=8, fontweight='bold')

    # --- Upper tornado: total combo extension ---
    # Baseline reference bar (top row)
    ax1.barh(y_ref, base_combo_m, height=0.55, color='#2c3e50', alpha=0.85,
             edgecolor='black', linewidth=0.8, align='center')
    ax1.text(base_combo_m * 0.5, y_ref, f'Baseline: {base_combo_m:.1f}%',
             ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    # Dashed line dropping from baseline bar through all tornado rows
    ax1.axvline(base_combo_m, color='#2c3e50', linewidth=1.2, linestyle='--', alpha=0.6)

    # Tornado bars offset to ride at the baseline value
    for i, (label, lo, hi, _, _) in enumerate(results):
        ax1.barh(i, hi, left=base_combo_m, height=0.4,
                 color='#2ecc71' if hi > 0 else '#3498db',
                 edgecolor='black', linewidth=0.5, align='center')
        ax1.barh(i, lo, left=base_combo_m, height=0.4,
                 color='#2ecc71' if lo > 0 else '#3498db',
                 edgecolor='black', linewidth=0.5, align='center')

    ax1.set_yticks(list(y_pos) + [y_ref])
    ax1.set_yticklabels([r[0] for r in results] + [''], fontsize=9)
    ax1.set_xlabel('Rapa+NMN Extension (%)', fontsize=10)
    ax1.set_title('Total Combo Extension', fontsize=12, fontweight='bold')
    ax1.axvline(0, color='black', linewidth=0.5, alpha=0.3)
    ax1.set_xlim(0, None)
    ax1.grid(True, alpha=0.3, axis='x')

    # Annotate perturbation values
    for i, (_, lo, hi, _, _) in enumerate(results):
        _annotate_tornado(ax1, base_combo_m, i, lo, hi)

    # --- Lower tornado: superadditive component ---
    # Baseline reference bar (top row)
    ax2.barh(y_ref, base_super_m, height=0.55, color='#8e2323', alpha=0.85,
             edgecolor='black', linewidth=0.8, align='center')
    ax2.text(base_super_m * 0.5, y_ref,
             f'Baseline: {base_super_m:+.1f} pp',
             ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    # Dashed line
    ax2.axvline(base_super_m, color='#8e2323', linewidth=1.2, linestyle='--', alpha=0.6)

    for i, (label, _, _, lo, hi) in enumerate(results):
        ax2.barh(i, hi, left=base_super_m, height=0.4,
                 color='#e74c3c' if hi > 0 else '#3498db',
                 edgecolor='black', linewidth=0.5, align='center')
        ax2.barh(i, lo, left=base_super_m, height=0.4,
                 color='#e74c3c' if lo > 0 else '#3498db',
                 edgecolor='black', linewidth=0.5, align='center')

    ax2.set_yticks(list(y_pos) + [y_ref])
    ax2.set_yticklabels([r[0] for r in results] + [''], fontsize=9)
    ax2.set_xlabel('Superadditive Component (pp)', fontsize=10)
    ax2.set_title('Superadditive Component', fontsize=12, fontweight='bold')
    ax2.axvline(0, color='black', linewidth=0.5, alpha=0.3)
    ax2.set_xlim(0, None)
    ax2.grid(True, alpha=0.3, axis='x')

    for i, (_, _, _, lo, hi) in enumerate(results):
        _annotate_tornado(ax2, base_super_m, i, lo, hi)

    # Check if superadditivity survives all perturbations
    all_survive = all(
        (base_super_m + r[3] > 0) and (base_super_m + r[4] > 0) for r in results
    )
    verdict = ('Superadditivity survives all \u00b120% perturbations'
               if all_survive else
               'Superadditivity is fragile under some perturbations')
    fig.suptitle(f'Rapa+NMN Sensitivity: \u00b120% Parameter Perturbation (Male)\n{verdict}',
                 fontsize=13, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, 'stage_A5_sensitivity.png')


# =============================================================================
# STAGE A26: SERIAL RELIABILITY
# =============================================================================

def _norm_pdf(x, mu, sigma):
    """Gaussian PDF."""
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))

def _norm_cdf(x, mu, sigma):
    """Gaussian CDF via error function."""
    import math
    return 0.5 * (1 + np.vectorize(math.erf)((x - mu) / (sigma * np.sqrt(2))))

def _survival(x, mu, sigma):
    """Survival function for a single Gaussian failure mode."""
    return 1 - _norm_cdf(x, mu, sigma)

def gen_stage_A26_serial_reliability():
    """Serial reliability: why single-target testing is blind.

    Top panel: Three failure mode PDFs at different times, showing masking.
    Bottom panel: Competing-risks survival curves showing what you observe.
    """
    # Failure mode parameters (mouse lifespan scale, months)
    # B is wide enough that its left tail extends earlier than A's
    A_mu, A_sig = 30, 2.5     # earliest, narrow, tall — the dominant killer
    B_mu, B_sig = 37, 9       # middle, wider + earlier — left tail underruns A
    C_mu, C_sig = 52, 3.5     # latest, narrow — completely invisible

    COL_A = '#e74c3c'   # red
    COL_B = '#f39c12'   # orange
    COL_C = '#9b59b6'   # purple

    t = np.linspace(5, 70, 1200)

    # PDFs
    pdf_A = _norm_pdf(t, A_mu, A_sig)
    pdf_B = _norm_pdf(t, B_mu, B_sig)
    pdf_C = _norm_pdf(t, C_mu, C_sig)

    # Survival functions
    S_A = _survival(t, A_mu, A_sig)
    S_B = _survival(t, B_mu, B_sig)
    S_C = _survival(t, C_mu, C_sig)

    # Competing risks (series system): S = product of individual survivals
    S_all     = S_A * S_B * S_C           # baseline: all 3 active
    S_fix_B   = S_A * S_C                 # fix mode B only
    S_fix_A   = S_B * S_C                 # fix mode A only
    S_fix_AB  = S_C                       # fix A+B → mode C revealed

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'hspace': 0.32})
    fig.suptitle('Serial Reliability: Why Single-Intervention Testing Is Blind',
                 fontsize=14, fontweight='bold', y=0.98)

    # ── TOP PANEL: Raw Gaussians and the sliver ──
    # Mode A — the dominant killer (narrow, tall)
    ax1.fill_between(t, pdf_A, alpha=0.20, color=COL_A)
    ax1.plot(t, pdf_A, color=COL_A, linewidth=2.5,
             label='Mode A \u2014 earliest (e.g. cancer)')

    # Mode B — wider, left tail extends earlier than A's
    ax1.fill_between(t, pdf_B, alpha=0.12, color=COL_B)
    ax1.plot(t, pdf_B, color=COL_B, linewidth=2.5,
             label='Mode B \u2014 middle, wider (e.g. cardiovascular)')

    # THE SLIVER: on the early (left) side, B's wide tail sits above A's narrow tail.
    # These are deaths from B that happen before A would have killed the organism.
    # Fix B → these deaths are averted → tiny survival extension. That's the signal.
    # Find the left-side crossover (first point where A rises above B)
    cross_mask = (pdf_A[1:] > pdf_B[1:]) & (pdf_A[:-1] <= pdf_B[:-1])
    cross_indices = np.where(cross_mask)[0]
    left_cross_t = t[cross_indices[0]] if len(cross_indices) > 0 else A_mu

    ax1.fill_between(t, pdf_A, pdf_B,
                     where=(pdf_B > pdf_A) & (t < left_cross_t + 1),
                     interpolate=True,
                     alpha=0.65, color=COL_B, edgecolor='#b37400', linewidth=1.5,
                     label='The sliver \u2014 total signal from fixing B')

    # Mode C — completely invisible, no overlap
    ax1.plot(t, pdf_C, color=COL_C, linewidth=2, linestyle='--', alpha=0.7,
             label='Mode C \u2014 invisible (e.g. structural)')
    ax1.fill_between(t, pdf_C, alpha=0.08, color=COL_C)

    # Annotations — point arrow squarely into the shaded crescent
    # The sliver is thickest around t=22, between pdf_A≈0.001 and pdf_B≈0.011
    # Use a fixed target in the middle of the visible orange fill
    sliver_diff = np.where((pdf_B > pdf_A) & (t < left_cross_t + 1),
                           pdf_B - pdf_A, 0)
    sliver_best_idx = np.argmax(sliver_diff)
    arrow_t = t[sliver_best_idx]
    arrow_y = (pdf_A[sliver_best_idx] + pdf_B[sliver_best_idx]) / 2

    ax1.annotate('The sliver \u2014 fix B,\nthese deaths vanish.\nThat\u2019s the entire signal.',
                 xy=(arrow_t, arrow_y),
                 xytext=(7, 0.045),
                 fontsize=9, color='#b37400', fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color='#b37400', lw=2,
                                 shrinkB=0))

    ax1.annotate('No tail below A\n\u2014 undetectable until A\n(and preferably B) is removed',
                 xy=(49, max(pdf_C) * 0.5),
                 xytext=(55, max(pdf_A) * 0.55),
                 fontsize=9, color=COL_C, fontweight='bold', ha='left',
                 arrowprops=dict(arrowstyle='->', color=COL_C, lw=1.5))

    ax1.set_xlabel('Time (months)', fontsize=11)
    ax1.set_ylabel('Failure probability', fontsize=11)
    ax1.set_title('Three failure modes \u2014 the organism dies from whichever hits first',
                  fontsize=11, style='italic', pad=8)
    ax1.legend(loc='upper right', fontsize=8.5, framealpha=0.9)
    ax1.set_xlim(5, 70)
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)

    # ── BOTTOM PANEL: Survival curves ──
    ax2.plot(t, S_all * 100, color='black', linewidth=2.5, label='Baseline (all modes active)')
    ax2.plot(t, S_fix_B * 100, color=COL_A, linewidth=2, linestyle='--',
             label='Fix B only \u2192 almost no gain')

    # Shade the sliver in survival space
    ax2.fill_between(t, S_all * 100, S_fix_B * 100, alpha=0.4, color=COL_B)

    ax2.plot(t, S_fix_A * 100, color=COL_B, linewidth=2, linestyle='--',
             label='Fix A \u2192 mode B becomes bottleneck')
    ax2.plot(t, S_fix_AB * 100, color=COL_C, linewidth=3,
             label='Fix A + B \u2192 mode C finally revealed')

    # Annotations — point arrow into the widest part of the orange gap
    surv_gap = (S_fix_B - S_all) * 100
    gap_best_idx = np.argmax(surv_gap)
    gap_t = t[gap_best_idx]
    gap_mid_y = (S_all[gap_best_idx] + S_fix_B[gap_best_idx]) / 2 * 100
    ax2.annotate('Same sliver\nin survival space',
                 xy=(gap_t, gap_mid_y),
                 xytext=(12, 30),
                 fontsize=9, fontweight='bold', color='#b37400',
                 arrowprops=dict(arrowstyle='->', color='#b37400', lw=2,
                                 shrinkB=0))

    ax2.annotate('Combination testing\nreveals this',
                 xy=(49, S_fix_AB[np.searchsorted(t, 49)] * 100),
                 xytext=(55, 55),
                 fontsize=9, fontweight='bold', color=COL_C,
                 arrowprops=dict(arrowstyle='->', color=COL_C, lw=1.5))

    ax2.set_xlabel('Time (months)', fontsize=11)
    ax2.set_ylabel('Survival (%)', fontsize=11)
    ax2.set_title('What you actually observe \u2014 survival curves',
                  fontsize=11, style='italic', pad=8)
    ax2.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax2.set_xlim(5, 70)
    ax2.set_ylim(-2, 105)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, 'stage_A26_serial_reliability.png')


# =============================================================================
# AUDIT FIGURES (adapted from hill_audit.py)
# =============================================================================

# Taguchi 8 compound list for audit figures
_AUDIT_COMPOUNDS = ['rapamycin', 'acarbose', 'canagliflozin', 'glycine',
                    'NMN', 'dasatinib_quercetin', 'urolithin_A']
_AUDIT_CD38I = {'cd38_inhibitor': 0.50}

def _run_audit_fullstack(sex, t_max=5.0):
    """Run full Taguchi-8 stack at 1x dose with sex modifiers."""
    merged = {}
    for compound in _AUDIT_COMPOUNDS:
        base = get_compound(compound)
        if base is None:
            continue
        scaled = apply_sex_modifier(compound, base, sex)
        for k, v in scaled.items():
            merged[k] = merged.get(k, 0) + v
    for k, v in _AUDIT_CD38I.items():
        merged[k] = merged.get(k, 0) + v
    merged['start_time'] = 0.30
    return simulate(interventions=merged, sex=sex, t_max=t_max, dt=DT_FINE)


def gen_audit_trajectories():
    """8-panel trajectory comparison: Control vs Taguchi-8 full stack."""
    plot_vars = [
        ('NAD', 'NAD+ Level'), ('DNA_damage', 'DNA Damage'),
        ('Methylation', 'Methylation'), ('Heteroplasmy', 'Heteroplasmy'),
        ('SenCells', 'Senescent Cells'), ('SASP', 'SASP'),
        ('Autophagy', 'Autophagy'), ('BioAge', 'Biological Age'),
    ]

    for sex in ['M', 'F']:
        ctrl_r = run_control(sex=sex, t_max=5.0, dt=DT_FINE)
        full_r = _run_audit_fullstack(sex, t_max=5.0)
        ext_pct = calculate_lifespan_extension(full_r, ctrl_r)
        print(f'  Audit {sex}: +{ext_pct:.1f}%')

        t_plot_max = min(full_r.t_death + 0.3, 5.0)

        fig, axes = plt.subplots(4, 2, figsize=(14, 16))
        fig.suptitle(f'State Variable Trajectories \u2014 {sex}ales\n'
                     f'Control vs Taguchi-8 (1\u00d7 dose, +{ext_pct:.0f}%)',
                     fontsize=14, fontweight='bold', y=0.98)

        for idx, (var, label) in enumerate(plot_vars):
            ax = axes[idx // 2, idx % 2]
            traj_ctrl = getattr(ctrl_r, var)
            traj_full = getattr(full_r, var)

            ax.plot(ctrl_r.t, traj_ctrl, 'k-', linewidth=1.5, label='Control')
            ax.plot(full_r.t, traj_full, 'b-', linewidth=1.5, label='Taguchi-8')
            ax.axvline(ctrl_r.t_death, color='k', linestyle='--', alpha=0.4, linewidth=0.8)
            ax.axvline(full_r.t_death, color='b', linestyle='--', alpha=0.4, linewidth=0.8)

            ax.set_title(label, fontsize=11, fontweight='bold')
            ax.set_xlabel('Time (normalized)')
            ax.set_xlim(0, t_plot_max)
            ax.set_ylim(bottom=0)
            ax.legend(fontsize=8, loc='best')
            ax.grid(alpha=0.3)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        save(fig, f'audit_trajectories_{sex}.png')


def gen_audit_bioage_decomp():
    """Stacked area: BioAge components over time, Control vs Full-Stack."""
    components = [
        ('Methylation', BIOAGE_PARAMS['w_meth'], BIOAGE_PARAMS['meth_norm'], '#e74c3c'),
        ('DNA_damage', BIOAGE_PARAMS['w_damage'], BIOAGE_PARAMS['damage_norm'], '#3498db'),
        ('Heteroplasmy', BIOAGE_PARAMS['w_hetero'], BIOAGE_PARAMS['H_norm'], '#2ecc71'),
        ('SenCells', BIOAGE_PARAMS['w_sen'], BIOAGE_PARAMS['sen_norm'], '#9b59b6'),
    ]

    for sex in ['M', 'F']:
        ctrl_r = run_control(sex=sex, t_max=5.0, dt=DT_FINE)
        full_r = _run_audit_fullstack(sex, t_max=5.0)
        ext_pct = calculate_lifespan_extension(full_r, ctrl_r)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(f'BioAge Decomposition \u2014 {sex}ales (Taguchi-8 at 1\u00d7)',
                     fontsize=14, fontweight='bold')

        for panel_idx, (label_prefix, result) in enumerate([
            ('Control', ctrl_r), (f'Taguchi-8 (+{ext_pct:.0f}%)', full_r),
        ]):
            ax = axes[panel_idx]
            t_end = min(result.t_death + 0.2, 5.0)
            idx_end = int(t_end / DT_FINE)

            contribs = []
            labels = []
            colors = []
            for var, w, norm, color in components:
                if w == 0:
                    continue
                traj = getattr(result, var)[:idx_end]
                contrib = w * traj / norm
                contribs.append(contrib)
                labels.append(f'{var} (w={w})')
                colors.append(color)

            ax.stackplot(result.t[:idx_end], *contribs, labels=labels, colors=colors, alpha=0.7)
            ax.axhline(1.0, color='k', linestyle='--', linewidth=1.5, label='Death threshold')
            ax.axvline(result.t_death, color='k', linestyle=':', alpha=0.5)
            ax.set_title(f'{label_prefix} (t_death={result.t_death:.3f})', fontsize=11)
            ax.set_xlabel('Time (normalized)')
            ax.set_ylabel('BioAge contribution')
            ax.set_xlim(0, t_end)
            ax.set_ylim(0, 1.5)
            ax.legend(fontsize=8, loc='upper left')
            ax.grid(alpha=0.3)

        plt.tight_layout(rect=[0, 0, 1, 0.94])
        save(fig, f'audit_bioage_decomp_{sex}.png')


def gen_audit_pairwise():
    """8x8 pairwise interaction heatmap for the Taguchi 8 compounds."""
    all_names = _AUDIT_COMPOUNDS + ['cd38i']
    short_names = ['Rapa', 'Acarb', 'Cana', 'Gly', 'NMN', 'D+Q', 'UroA', 'CD38i']
    n_compounds = len(all_names)

    sex = 'M'
    ctrl_r = run_control(sex=sex, t_max=T_MAX, dt=DT_FINE)

    # Individual effects
    singles = {}
    for compound in _AUDIT_COMPOUNDS:
        result = simulate(compound=compound, sex=sex, t_max=T_MAX, dt=DT_FINE)
        singles[compound] = calculate_lifespan_extension(result, ctrl_r)
    cd38_result = simulate(
        interventions={**_AUDIT_CD38I, 'start_time': 0.30},
        sex=sex, t_max=T_MAX, dt=DT_FINE
    )
    singles['cd38i'] = calculate_lifespan_extension(cd38_result, ctrl_r)

    # Pairwise
    pair_matrix = np.zeros((n_compounds, n_compounds))
    ratio_matrix = np.ones((n_compounds, n_compounds))

    for i, c1 in enumerate(all_names):
        pair_matrix[i, i] = singles[c1]
        for j, c2 in enumerate(all_names):
            if j <= i:
                continue
            merged = {}
            for c in [c1, c2]:
                if c == 'cd38i':
                    for k, v in _AUDIT_CD38I.items():
                        merged[k] = merged.get(k, 0) + v
                else:
                    base = get_compound(c)
                    if base:
                        scaled = apply_sex_modifier(c, base, sex)
                        for k, v in scaled.items():
                            merged[k] = merged.get(k, 0) + v
            merged['start_time'] = 0.30
            result = simulate(interventions=merged, sex=sex, t_max=T_MAX, dt=DT_FINE)
            combo_ext = calculate_lifespan_extension(result, ctrl_r)
            pair_matrix[i, j] = combo_ext
            pair_matrix[j, i] = combo_ext

            naive_sum = singles[c1] + singles[c2]
            if naive_sum > 0:
                ratio = combo_ext / naive_sum
                ratio_matrix[i, j] = ratio
                ratio_matrix[j, i] = ratio

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle('Pairwise Interaction Audit (Males, 1\u00d7 dose)',
                 fontsize=14, fontweight='bold')

    im1 = axes[0].imshow(pair_matrix, cmap='YlOrRd', aspect='equal')
    axes[0].set_xticks(range(n_compounds))
    axes[0].set_xticklabels(short_names, rotation=45, ha='right', fontsize=9)
    axes[0].set_yticks(range(n_compounds))
    axes[0].set_yticklabels(short_names, fontsize=9)
    axes[0].set_title('Combination Extension (%)', fontsize=11)
    for i in range(n_compounds):
        for j in range(n_compounds):
            if i == j:
                axes[0].text(j, i, f'{pair_matrix[i,j]:.0f}',
                           ha='center', va='center', fontsize=8, fontweight='bold')
            elif j > i:
                axes[0].text(j, i, f'{pair_matrix[i,j]:.0f}',
                           ha='center', va='center', fontsize=8)
    plt.colorbar(im1, ax=axes[0], shrink=0.8)

    cmap = plt.cm.RdYlGn
    im2 = axes[1].imshow(ratio_matrix, cmap=cmap, vmin=0.5, vmax=1.5, aspect='equal')
    axes[1].set_xticks(range(n_compounds))
    axes[1].set_xticklabels(short_names, rotation=45, ha='right', fontsize=9)
    axes[1].set_yticks(range(n_compounds))
    axes[1].set_yticklabels(short_names, fontsize=9)
    axes[1].set_title('Ratio: combo / naive_sum\n(<1 = sub-additive, >1 = super-additive)',
                      fontsize=11)
    for i in range(n_compounds):
        for j in range(n_compounds):
            if j > i:
                r = ratio_matrix[i, j]
                color = 'white' if abs(r - 1.0) > 0.3 else 'black'
                axes[1].text(j, i, f'{r:.2f}',
                           ha='center', va='center', fontsize=8, color=color)
            elif i == j:
                axes[1].text(j, i, '\u2014', ha='center', va='center', fontsize=8, color='gray')
    plt.colorbar(im2, ax=axes[1], shrink=0.8)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, 'audit_pairwise_matrix.png')


# =============================================================================
# MAIN
# =============================================================================

GENERATORS = {
    'stage_00': gen_stage_00,
    'stage_01': gen_stage_01,
    'stage_02': gen_stage_02,
    'stage_03': gen_stage_03,
    'stage_04': gen_stage_04,
    'stage_05': gen_stage_05,
    'stage_06': gen_stage_06,
    'stage_07': gen_stage_07,
    'stage_08': gen_stage_08,
    'stage_09': gen_stage_09,
    'stage_10': gen_stage_10,
    'stage_11': gen_stage_11,
    'stage_14': gen_stage_14,
    'stage_15': gen_stage_15,
    'stage_16': gen_stage_16,
    'stage_17': gen_stage_17,
    'w_meth': gen_w_meth_sweep,
    'analytical': gen_analytical,
    'stage_A10': gen_stage_A10,
    'stage_A5': gen_stage_A5_sensitivity,
    'stage_A26': gen_stage_A26_serial_reliability,
    'stage_A6b': gen_stage_A6b_weight_sensitivity,
    'audit_trajectories': gen_audit_trajectories,
    'audit_bioage': gen_audit_bioage_decomp,
    'audit_pairwise': gen_audit_pairwise,
}

def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    if len(sys.argv) > 1:
        key = sys.argv[1]
        if key in GENERATORS:
            print(f'Generating {key} figures...')
            GENERATORS[key]()
        else:
            print(f'Unknown generator: {key}')
            print(f'Available: {", ".join(GENERATORS.keys())}')
            sys.exit(1)
    else:
        print(f'Generating all figures into {FIGURES_DIR}/')
        for key, gen in GENERATORS.items():
            print(f'\n[{key}]')
            gen()
        print('\nDone.')


if __name__ == '__main__':
    main()
