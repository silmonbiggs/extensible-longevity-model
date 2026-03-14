#!/usr/bin/env python3
"""
Generate sex-differences figure with error bars.

Model uncertainty (female predictions): RSS of OAT per-compound sensitivities,
scaled by each parameter's plausible range.

ITP target uncertainty: estimated SE of median extension from ITP multi-site
design (~150 mice/sex, 3 sites, bootstrap SE of median).

Usage:
    python scripts/generate_sex_diff_errorbars.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from elm.model import simulate, run_control, calculate_lifespan_extension
from elm.compounds import ITP_VALIDATION

FIGURES_DIR = Path(os.path.dirname(__file__)) / '..' / 'docs' / 'figures'
OAT_TSV = Path(os.path.dirname(__file__)) / '..' / 'docs' / 'oat_sensitivity.tsv'
DT = 0.002
T_MAX = 2.5
DPI = 300

# Compound display order (matches paper figure)
COMPOUNDS = ['aspirin', '17_alpha_estradiol', 'canagliflozin',
             'glycine', 'acarbose', 'rapamycin']
NAMES = {
    'aspirin': 'Aspirin', '17_alpha_estradiol': '17\u03b1-Estradiol',
    'canagliflozin': 'Canagliflozin', 'glycine': 'Glycine',
    'acarbose': 'Acarbose', 'rapamycin': 'Rapamycin',
}

# -------------------------------------------------------------------------
# Plausible half-ranges (as fraction of parameter value)
# From supplementary tables: "Set" parameters use stated ranges;
# literature parameters use +/- 20% (conservative default).
# -------------------------------------------------------------------------
PLAUSIBLE_HALF_RANGE = {
    # Heteroplasmy
    'heteroplasmy.transcription_advantage': 0.25,   # lit (Wallace 2010)
    'heteroplasmy.time_scale': 5.0 / 12.0,          # range 8-18, central 12
    'heteroplasmy.k_selection_foxo3': 0.075 / 0.10,  # range 0.05-0.20
    'heteroplasmy.H_0': 0.20,                        # lit (Stewart 2015)
    'heteroplasmy.Neff': 0.20,                       # lit (Shoubridge 2007)
    'heteroplasmy.k_selection_base': 0.50,           # range 0-0.01, wide
    # NAD+
    'nad.k_nampt_base': 0.20,                        # lit (Yoshino 2011)
    'nad.k_consumption_basal': 0.15 / 0.40,          # range 0.25-0.55
    'nad.k_cd38_base': 0.20,                         # lit (Chini 2020)
    'nad.k_cd38_age': 0.20,                          # lit (Camacho-Pereira 2016)
    'nad.k_cd38_sasp': 0.25 / 0.60,                  # range 0.3-0.8
    'nad.k_parp_consumption': 0.20,                  # lit (Massudi 2012)
    'nad.k_sirt_consumption': 0.25,                  # range 0.05-0.20
    'nad.k_nampt_age_decline': 0.20,                 # lit (Yoshino 2011)
    'nad.k_oxphos_nad_regen': 0.25,                  # estimated
    'nad.k_tryptophan': 0.25,                        # estimated
    'nad.Km_nad_consumption': 0.25,                  # estimated
    'nad.NAD_young': 0.10,                           # structural
    # mTORC1
    'mtorc1.k_mtorc1_damage_reduction': 0.25,       # estimated
    'mtorc1.k_mtorc1_senescence': 0.20,             # lit (Demidenko 2009)
    'mtorc1.k_mtorc1_autophagy': 0.20,              # lit (Kim 2011)
    # Sirtuin
    'sirtuin.k_sirt1_foxo3': 0.20,                  # lit (Brunet 2004)
    'sirtuin.Km_nad_sirt1': 0.15 / 0.40,            # range 0.3-0.6
    'sirtuin.n_hill_sirt1': 0.25,                    # range 1-3
    'sirtuin.foxo3_basal': 0.25,                     # estimated
    # TSC2
    'tsc2.k_tsc2_mtorc1_max': 0.20,                 # lit (Inoki 2003)
    'tsc2.Km_tsc2_mtorc1': 0.25,                    # estimated
    'tsc2.n_hill_tsc2_mtorc1': 0.25 / 1.2,          # range 1.0-2.0
    'tsc2.TSC2_basal': 0.20,                         # lit (Huang 2008)
    # DNA repair
    'dna_repair.k_damage_ros': 0.20 / 0.60,         # range 0.40-0.80
    'dna_repair.k_damage_replication': 0.25,         # range 0.05-0.25
    'dna_repair.k_damage_spontaneous': 0.25,         # estimated
    'dna_repair.k_ber_base': 0.20,                   # range 0.20-0.50
    'dna_repair.k_ber_nad_dependence': 0.20,         # range 0.4-1.0
    'dna_repair.k_antioxidant_ros_reduction': 0.20,  # range 0.35-0.75
    # Methylation
    'methylation.k_ezh2_base': 0.35 / 1.85,         # range 1.5-2.2
    'methylation.k_dnmt_base': 0.25,                 # range 0.15-0.40
    'methylation.k_sasp_methylation': 0.20 / 0.42,  # range 0.2-0.6
    'methylation.meth_young': 0.20,                  # lit (Horvath 2013)
    # Senescence
    'senescence.k_sen_from_damage': 0.25,            # range 0.05-0.30
    'senescence.k_sen_from_telomere': 0.25,          # range 0.05-0.20
    'senescence.k_sen_from_oncogene': 0.25,          # estimated
    'senescence.k_sen_natural_clear': 0.25,          # lit (Baker 2016)
    'senescence.k_sen_autophagy_clear': 0.20,        # range 0.40-1.00
    # Autophagy
    'autophagy.Km_autophagy': 0.20,                  # range 0.2-0.5
    'autophagy.w_foxo3_autophagy': 0.20,             # range 0.2-0.5
    'autophagy.Vmax_autophagy': 0.10,                # structural
    # UPRmt
    'uprmt.Vmax_uprmt': 0.10,                        # structural
    'uprmt.Km_uprmt': 0.20,                          # estimated
    'uprmt.k_uprmt_proteostasis': 0.25,              # estimated
    'uprmt.k_uprmt_nad': 0.25,                       # estimated
    'uprmt.k_uprmt_sen_reduction': 0.25,             # estimated
    'uprmt.k_uprmt_hetero': 0.25,                    # estimated
    'uprmt.k_uprmt_ros': 0.25,                       # estimated
    'uprmt.k_uprmt_autophagy': 0.25,                 # estimated
    # SASP
    'sasp.k_sasp_production': 0.20,                  # lit (Coppe 2008)
    'sasp.k_sasp_decay': 0.25,                       # range 0.30-1.00
}
DEFAULT_HALF_RANGE = 0.20  # for unlisted parameters

# -------------------------------------------------------------------------
# ITP target SE estimates (pp)
# Based on ~150 mice/sex across 3 ITP sites.
# SE(median extension) ~ 1.253 * SD / sqrt(n) for each group,
# SE(difference) ~ sqrt(2) * SE(single).
# With SD ~ 4 months, n ~ 150: SE ~ 2.0 pp.
# -------------------------------------------------------------------------
ITP_SE = {
    'rapamycin':         {'M': 2.5, 'F': 2.5},
    'acarbose':          {'M': 2.5, 'F': 2.5},
    'canagliflozin':     {'M': 2.5, 'F': 2.5},
    '17_alpha_estradiol': {'M': 2.5, 'F': 2.5},
    'aspirin':           {'M': 2.5, 'F': 2.5},
    'glycine':           {'M': 2.5, 'F': 2.5},
}


def compute_model_brackets():
    """Compute RSS model uncertainty brackets from OAT per-compound data."""
    # Read OAT TSV
    rows = []
    with open(OAT_TSV) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)

    # For each female prediction, compute RSS
    female_cols = {c: f'F_{c}' for c in COMPOUNDS}
    brackets = {}

    for compound in COMPOUNDS:
        col = female_cols[compound]
        ss = 0.0  # sum of squares
        for row in rows:
            param_name = row['parameter']
            oat_shift = float(row[col])
            if oat_shift == 0.0:
                continue
            # Get plausible half-range
            half_range = PLAUSIBLE_HALF_RANGE.get(param_name, DEFAULT_HALF_RANGE)
            # Scale: OAT was ±10%, actual range is ±half_range
            scaled = oat_shift * (half_range / 0.10)
            ss += scaled ** 2

        brackets[compound] = np.sqrt(ss)

    return brackets


def compute_predictions():
    """Compute model predictions for all compounds (current model = dual-action cana)."""
    preds = {}
    for sex in ['M', 'F']:
        ctrl = run_control(sex=sex, t_max=T_MAX, dt=DT)
        for c in COMPOUNDS:
            result = simulate(compound=c, sex=sex, t_max=T_MAX, dt=DT)
            ext = calculate_lifespan_extension(result, ctrl)
            preds[(c, sex)] = ext
    return preds


def compute_predictions_pure_ampk():
    """Compute predictions with canagliflozin as pure AMPK (baseline model)."""
    from scipy.optimize import brentq
    from elm.sex_mechanisms import apply_sex_modifier
    from elm.compounds import get_itp_start_time

    preds = compute_predictions()  # start with current model

    # Override canagliflozin with pure-AMPK
    cana_start = get_itp_start_time('canagliflozin')
    cana_target_m = ITP_VALIDATION['canagliflozin'].target_male

    ctrl_m = run_control(sex='M', t_max=T_MAX, dt=DT)
    ctrl_f = run_control(sex='F', t_max=T_MAX, dt=DT)

    # Find scale factor for pure-AMPK canagliflozin to hit male target
    def cana_error(s):
        inj = {'ampk': s}
        inj = apply_sex_modifier('canagliflozin', inj, 'M')
        inj['start_time'] = cana_start
        r = simulate(interventions=inj, sex='M', t_max=T_MAX, dt=DT)
        return calculate_lifespan_extension(r, ctrl_m) - cana_target_m

    s = brentq(cana_error, 0.01, 5.0, xtol=0.001)

    # Male prediction (should be ~14%)
    inj_m = apply_sex_modifier('canagliflozin', {'ampk': s}, 'M')
    inj_m['start_time'] = cana_start
    r_m = simulate(interventions=inj_m, sex='M', t_max=T_MAX, dt=DT)
    preds[('canagliflozin', 'M')] = calculate_lifespan_extension(r_m, ctrl_m)

    # Female prediction (pure-AMPK, attenuated)
    inj_f = apply_sex_modifier('canagliflozin', {'ampk': s}, 'F')
    inj_f['start_time'] = cana_start
    r_f = simulate(interventions=inj_f, sex='F', t_max=T_MAX, dt=DT)
    preds[('canagliflozin', 'F')] = calculate_lifespan_extension(r_f, ctrl_f)

    return preds


def plot_figure(preds, model_brackets, old_preds=None):
    """Generate journal-style sex differences figure with paired bars.

    If old_preds is provided, draw a ghost annotation on canagliflozin female
    showing the improvement from the old prediction to the new one.
    """
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
    })

    fig, ax = plt.subplots(figsize=(7.0, 4.0))

    n = len(COMPOUNDS)
    # Layout: 4 thin bars per compound, tightly packed
    # [M_ITP, M_ELM, F_ITP, F_ELM] — color distinguishes sex
    bw = 0.055         # thin bars
    gap = 0.005        # hairline gap between adjacent bars
    compound_spacing = 0.50  # generous space between compounds

    targets_m = [ITP_VALIDATION[c].target_male for c in COMPOUNDS]
    targets_f = [ITP_VALIDATION[c].target_female for c in COMPOUNDS]
    elm_m = [preds[(c, 'M')] for c in COMPOUNDS]
    elm_f = [preds[(c, 'F')] for c in COMPOUNDS]
    itp_se = [ITP_SE[c]['M'] for c in COMPOUNDS]  # uniform ±2.5
    model_se_f = [model_brackets[c] for c in COMPOUNDS]
    names = [NAMES[c] for c in COMPOUNDS]

    # 4 bars centered on x_center: positions -1.5, -0.5, +0.5, +1.5 units
    step = bw + gap
    x_centers = np.arange(n) * compound_spacing
    x_m_itp = x_centers - 1.5 * step
    x_m_elm = x_centers - 0.5 * step
    x_f_itp = x_centers + 0.5 * step
    x_f_elm = x_centers + 1.5 * step

    # Colors: ITP (data) gets saturated; ELM (model) gets muted
    c_itp_m = '#1565C0'   # strong blue — data
    c_elm_m = '#90CAF9'   # light blue — model
    c_itp_f = '#C62828'   # strong red — data
    c_elm_f = '#FFAB91'   # light salmon — model

    # Bars with error brackets
    ax.bar(x_m_itp, targets_m, bw, color=c_itp_m, edgecolor='#0D47A1',
           linewidth=0.5, yerr=itp_se, error_kw=dict(
               elinewidth=0.8, capsize=2, capthick=0.6, ecolor='#333333'),
           label='ITP male', zorder=2)
    ax.bar(x_m_elm, elm_m, bw, color=c_elm_m, edgecolor='#64B5F6',
           linewidth=0.5, label='ELM male', zorder=2)

    ax.bar(x_f_itp, targets_f, bw, color=c_itp_f, edgecolor='#B71C1C',
           linewidth=0.5, yerr=itp_se, error_kw=dict(
               elinewidth=0.8, capsize=2, capthick=0.6, ecolor='#333333'),
           label='ITP female', zorder=2)
    ax.bar(x_f_elm, elm_f, bw, color=c_elm_f, edgecolor='#FF8A65',
           linewidth=0.5, yerr=model_se_f, error_kw=dict(
               elinewidth=0.8, capsize=2, capthick=0.6, ecolor='#333333'),
           label='ELM female', zorder=2)

    # Ghost annotation: old canagliflozin prediction
    if old_preds is not None:
        cana_idx = COMPOUNDS.index('canagliflozin')
        old_f = old_preds[('canagliflozin', 'F')]
        new_f = elm_f[cana_idx]
        xbar = x_f_elm[cana_idx]  # x-center of the ELM female cana bar
        # Dashed ghost line at old bar height — extends past bar edge under arrow
        ghost_left = xbar - bw * 0.6
        ghost_right = xbar + bw * 1.4  # past the right edge, under the arrow
        ax.plot([ghost_left, ghost_right], [old_f, old_f],
                color='#888888', linestyle='--', linewidth=1.0, zorder=5)
        # Arrow from old to new, offset to the right of the bar
        x_arrow = xbar + bw * 0.9
        ax.annotate('', xy=(x_arrow, new_f), xytext=(x_arrow, old_f),
                    arrowprops=dict(arrowstyle='->', color='#555555',
                                    lw=1.2, shrinkA=1, shrinkB=1),
                    zorder=5)
        ax.text(x_arrow + bw * 0.4, (old_f + new_f) / 2,
                f'+{new_f - old_f:.1f}', fontsize=6.5, color='#555555',
                ha='left', va='center')

    ax.set_xticks(x_centers)
    ax.set_xticklabels(names, rotation=30, ha='right')
    ax.set_ylabel('Lifespan extension (%)')
    ax.set_ylim(-10, 32)
    ax.axhline(0, color='#cccccc', linewidth=0.5, zorder=1)

    ax.legend(fontsize=7, loc='upper left', frameon=False, ncol=2)

    # Mechanism brackets
    brackets = [
        ('COX-2/PGI2', '#cc00cc', 0, 0),
        ('Testosterone', '#e74c3c', 1, 1),
        ('AMPK sat.', '#333333', 2, 4),
        ('CYP3A PK', '#1565C0', 5, 5),
    ]
    margin = 2.0 * step + 0.02
    # Two tiers: AMPK on lower tier spanning ALL compounds; others on upper tier
    for label, color, i_start, i_end in brackets:
        if 'AMPK' in label:
            # Lower tier: spans ALL compounds
            y_ann = -7.0
            x0 = x_centers[0] - margin
            x1 = x_centers[-1] + margin
        else:
            # Upper tier
            y_ann = -3.2
            x0 = x_centers[i_start] - margin
            x1 = x_centers[i_end] + margin
        mid = (x0 + x1) / 2.0
        ax.plot([x0, x1], [y_ann, y_ann], color=color, lw=1.5,
                clip_on=False, solid_capstyle='butt')
        ax.plot([x0, x0], [y_ann + 0.7, y_ann], color=color, lw=1.5,
                clip_on=False, solid_capstyle='butt')
        ax.plot([x1, x1], [y_ann + 0.7, y_ann], color=color, lw=1.5,
                clip_on=False, solid_capstyle='butt')
        ax.text(mid, y_ann - 1.2, label, ha='center', fontsize=7,
                fontweight='bold', color=color)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out')

    fig.tight_layout()
    return fig


def print_table(label, preds, model_brackets):
    print(f'\n{label}:')
    print(f'  {"Compound":25s}  {"Pred F%":>8s}  {"Model SE":>9s}  {"ITP SE":>7s}')
    print(f'  {"-"*55}')
    for c in COMPOUNDS:
        pred_f = preds[(c, 'F')]
        target_f = ITP_VALIDATION[c].target_female
        m_se = model_brackets[c]
        i_se = ITP_SE[c]['F']
        combined = np.sqrt(m_se**2 + i_se**2)
        err = abs(pred_f - target_f)
        within = 'yes' if err <= combined else 'no'
        name = c.replace('17_alpha_estradiol', '17a-Estradiol')
        print(f'  {name:25s}  {pred_f:8.1f}  {m_se:9.2f}  {i_se:7.1f}'
              f'  | err={err:.1f}, combined={combined:.1f}, within={within}')


def main():
    print('Computing model uncertainty brackets from OAT...')
    model_brackets = compute_model_brackets()

    # Act 1: Pure-AMPK baseline (primary out-of-sample result)
    print('\nComputing pure-AMPK predictions (act 1)...')
    preds_baseline = compute_predictions_pure_ampk()
    print_table('Act 1: Pure-AMPK (out-of-sample)', preds_baseline, model_brackets)

    fig1 = plot_figure(preds_baseline, model_brackets)
    out1 = FIGURES_DIR / 'sex_diff_baseline.png'
    fig1.savefig(out1, dpi=DPI, bbox_inches='tight')
    plt.close(fig1)
    print(f'  Figure saved: {out1}')

    # Act 2: Dual-action canagliflozin (with FGF21 hypothesis)
    print('\nComputing dual-action predictions (act 2)...')
    preds_dual = compute_predictions()
    print_table('Act 2: Dual-action canagliflozin', preds_dual, model_brackets)

    fig2 = plot_figure(preds_dual, model_brackets, old_preds=preds_baseline)
    out2 = FIGURES_DIR / 'sex_diff_dual_action.png'
    fig2.savefig(out2, dpi=DPI, bbox_inches='tight')
    plt.close(fig2)
    print(f'  Figure saved: {out2}')


if __name__ == '__main__':
    main()
