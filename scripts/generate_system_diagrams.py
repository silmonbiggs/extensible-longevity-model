#!/usr/bin/env python3
"""
Generate static system diagram PNGs for Papers 1 and 2.

Produces:
  1. system_diagram_base.png       - The model architecture (Paper 1 Fig 1)
  2. system_diagram_itp.png        - With ITP compound entry points (Paper 1 Fig 1 alt)
  3. system_diagram_novel.png      - With ITP + novel compound entry points (Paper 2)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from pathlib import Path

FIGDIR = Path(__file__).parent.parent / 'docs' / 'figures'

# ===================================================================
# NODE POSITIONS (from deck.html, normalized to [0,1] range)
# ===================================================================

# Raw positions from deck
RAW = {
    'H': (120, 115), 'ROS': (120, 220), 'UPRmt': (120, 330),
    'NAMPT': (355, 115), 'NAD': (530, 180),
    'CD38': (355, 260), 'SIRT1': (530, 260),
    'Repair': (750, 180), 'DNA_Damage': (750, 440),
    'AMPK': (960, 115), 'TSC2': (960, 180), 'mTORC1': (960, 260),
    'Autophagy': (960, 330),
    'SenCells': (120, 550), 'SASP': (355, 550),
    'Methylation': (660, 550),
    'BioAge': (440, 660),
}

# Normalize to figure coords (flip y)
XMAX, YMAX = 1200, 720
def norm(pos):
    return {k: (v[0]/XMAX, 1 - v[1]/YMAX) for k, v in pos.items()}

POS = norm(RAW)

# ITP interventions
ITP_COMPOUNDS = {
    'Rapamycin': (1100/XMAX, 1 - (-10)/YMAX),
    'Acarbose': (960/XMAX, 1 - (-10)/YMAX),
    'Canagliflozin': (680/XMAX, 1 - (-10)/YMAX),
    '17aE2': (540/XMAX, 1 - (-10)/YMAX),
    'Aspirin': (400/XMAX, 1 - (-10)/YMAX),
    'Glycine': (820/XMAX, 1 - (-10)/YMAX),
}

# Novel compounds (positioned above existing nodes)
NOVEL_COMPOUNDS = {
    'NMN': (260/XMAX, 1 - (-10)/YMAX),
    'CD38i': (355/XMAX, 1 - (-60)/YMAX),
    'Urolithin A': (40/XMAX, 1 - (-10)/YMAX),
    'D+Q': (120/XMAX, 1 - (-60)/YMAX),
}

# ===================================================================
# REGION BOXES
# ===================================================================

REGIONS = {
    'MITOCHONDRIA': {'nodes': ['H', 'ROS', 'UPRmt'], 'color': '#FFCDD2', 'border': '#E57373'},
    'NAD+ HOMEOSTASIS': {'nodes': ['NAMPT', 'NAD', 'CD38', 'SIRT1'], 'color': '#C8E6C9', 'border': '#81C784'},
    'DNA DAMAGE': {'nodes': ['Repair', 'DNA_Damage'], 'color': '#FFE0B2', 'border': '#FFB74D'},
    'CORE SIGNALING': {'nodes': ['AMPK', 'TSC2', 'mTORC1', 'Autophagy'], 'color': '#BBDEFB', 'border': '#64B5F6'},
    'SENESCENCE': {'nodes': ['SenCells'], 'color': '#E1BEE7', 'border': '#BA68C8'},
    'INFLAMMATION': {'nodes': ['SASP'], 'color': '#FFCCBC', 'border': '#FF8A65'},
    'EPIGENETIC DRIFT': {'nodes': ['Methylation'], 'color': '#B2DFDB', 'border': '#4DB6AC'},
}

# ===================================================================
# EDGES
# ===================================================================

# Biological edges: (source, target, type)
# type: 'act' = activation (green), 'inh' = inhibition (red), 'dash' = coupling (gray dashed)
BIO_EDGES = [
    ('H', 'ROS', 'act'),
    ('ROS', 'DNA_Damage', 'act'),
    ('NAMPT', 'NAD', 'act'),
    ('NAD', 'SIRT1', 'act'),
    ('NAD', 'Repair', 'dash'),
    ('Repair', 'DNA_Damage', 'inh'),
    ('H', 'NAMPT', 'inh'),
    ('DNA_Damage', 'NAD', 'inh'),
    ('CD38', 'NAD', 'inh'),
    ('AMPK', 'TSC2', 'act'),
    ('TSC2', 'mTORC1', 'inh'),
    ('mTORC1', 'Autophagy', 'inh'),
    ('Autophagy', 'DNA_Damage', 'inh'),
    ('DNA_Damage', 'SenCells', 'act'),
    ('SenCells', 'SASP', 'act'),
    ('SASP', 'CD38', 'dash'),
    ('SASP', 'DNA_Damage', 'dash'),
    ('SASP', 'Methylation', 'dash'),
    ('SIRT1', 'Methylation', 'inh'),
    ('AMPK', 'NAMPT', 'dash'),
    # BioAge paths
    ('H', 'BioAge', 'bio'),
    ('DNA_Damage', 'BioAge', 'bio'),
    ('SenCells', 'BioAge', 'bio'),
    ('Methylation', 'BioAge', 'bio'),
]

# ITP intervention edges
ITP_EDGES = [
    ('Rapamycin', 'mTORC1', 'iv_inh'),
    ('Rapamycin', 'AMPK', 'iv_act'),
    ('Acarbose', 'AMPK', 'iv_act'),
    ('Acarbose', 'CD38', 'iv_inh'),
    ('Canagliflozin', 'AMPK', 'iv_act'),
    ('17aE2', 'AMPK', 'iv_act'),
    ('17aE2', 'DNA_Damage', 'iv_inh'),
    ('Aspirin', 'DNA_Damage', 'iv_inh'),
    ('Aspirin', 'SASP', 'iv_inh'),
    ('Glycine', 'AMPK', 'iv_act'),
    ('Glycine', 'DNA_Damage', 'iv_inh'),
]

# Novel compound edges
NOVEL_EDGES = [
    ('NMN', 'NAD', 'iv_act'),
    ('CD38i', 'CD38', 'iv_inh'),
    ('Urolithin A', 'H', 'iv_inh'),
    ('D+Q', 'SenCells', 'iv_inh'),
]


def draw_node(ax, name, pos, ntype='bio', highlight=False):
    x, y = pos
    if ntype == 'state':
        color = '#1565C0' if not highlight else '#E65100'
        ax.plot(x, y, 'o', markersize=18, color=color, zorder=5)
        ax.text(x, y, name.replace('_', '\n'), ha='center', va='center',
                fontsize=6, fontweight='bold', color='white', zorder=6)
    elif ntype == 'enzyme':
        ax.plot(x, y, 's', markersize=14, color='#546E7A', zorder=5)
        ax.text(x, y, name.replace('_', '\n'), ha='center', va='center',
                fontsize=5.5, color='white', zorder=6)
    elif ntype == 'output':
        ax.plot(x, y, 'D', markersize=20, color='#37474F', zorder=5)
        ax.text(x, y, 'BioAge', ha='center', va='center',
                fontsize=6, fontweight='bold', color='white', zorder=6)
    elif ntype == 'itp':
        ax.plot(x, y, 'o', markersize=14, color='#2E7D32', zorder=5)
        ax.text(x, y, name, ha='center', va='center',
                fontsize=4.5, fontweight='bold', color='white', zorder=6)
    elif ntype == 'novel':
        ax.plot(x, y, 'o', markersize=14, color='#C62828', zorder=5)
        ax.text(x, y, name, ha='center', va='center',
                fontsize=4.5, fontweight='bold', color='white', zorder=6)


def draw_edge(ax, p1, p2, etype='act'):
    x1, y1 = p1
    x2, y2 = p2
    styles = {
        'act': dict(color='#2E7D32', lw=1.2, ls='-', alpha=0.7),
        'inh': dict(color='#C62828', lw=1.2, ls='-', alpha=0.7),
        'dash': dict(color='#78909C', lw=0.8, ls='--', alpha=0.5),
        'bio': dict(color='#37474F', lw=1.0, ls=':', alpha=0.4),
        'iv_act': dict(color='#2E7D32', lw=1.5, ls='-', alpha=0.8),
        'iv_inh': dict(color='#C62828', lw=1.5, ls='-', alpha=0.8),
    }
    s = styles.get(etype, styles['act'])

    head = ']->' if 'inh' in etype else '->'
    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle='->', mutation_scale=10,
        color=s['color'], lw=s['lw'], linestyle=s['ls'], alpha=s['alpha'],
        zorder=3,
        connectionstyle='arc3,rad=0.05' if abs(x1 - x2) > 0.01 and abs(y1 - y2) > 0.01 else 'arc3,rad=0',
        shrinkA=10, shrinkB=10,
    )
    ax.add_patch(arrow)


def draw_region(ax, name, nodes, color, border):
    xs = [POS[n][0] for n in nodes if n in POS]
    ys = [POS[n][1] for n in nodes if n in POS]
    if not xs:
        return
    pad = 0.04
    x0, x1 = min(xs) - pad, max(xs) + pad
    y0, y1 = min(ys) - pad, max(ys) + pad
    rect = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                          boxstyle='round,pad=0.01',
                          facecolor=color, edgecolor=border,
                          alpha=0.25, lw=1.5, zorder=1)
    ax.add_patch(rect)
    ax.text(x0 + 0.005, y1 - 0.005, name, fontsize=5, color=border,
            fontweight='bold', va='top', zorder=2)


# Node types
STATE_VARS = ['H', 'NAD', 'DNA_Damage', 'SenCells', 'Methylation']
ENZYMES = ['ROS', 'UPRmt', 'NAMPT', 'CD38', 'SIRT1', 'Repair',
           'AMPK', 'TSC2', 'mTORC1', 'Autophagy', 'SASP']


def make_diagram(title, show_itp=False, show_novel=False, filename='diagram.png'):
    fig, ax = plt.subplots(1, 1, figsize=(10, 7.5))
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.02, 1.08)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)

    # Draw regions
    for rname, rdata in REGIONS.items():
        draw_region(ax, rname, rdata['nodes'], rdata['color'], rdata['border'])

    # Draw bio edges
    for src, tgt, etype in BIO_EDGES:
        if src in POS and tgt in POS:
            draw_edge(ax, POS[src], POS[tgt], etype)

    # Draw bio nodes
    for n in STATE_VARS:
        draw_node(ax, n, POS[n], 'state')
    for n in ENZYMES:
        draw_node(ax, n, POS[n], 'enzyme')
    draw_node(ax, 'BioAge', POS['BioAge'], 'output')

    # ITP compounds
    if show_itp:
        for name, pos in ITP_COMPOUNDS.items():
            draw_node(ax, name, pos, 'itp')
        for src, tgt, etype in ITP_EDGES:
            if src in ITP_COMPOUNDS and tgt in POS:
                draw_edge(ax, ITP_COMPOUNDS[src], POS[tgt], etype)

    # Novel compounds
    if show_novel:
        for name, pos in NOVEL_COMPOUNDS.items():
            draw_node(ax, name, pos, 'novel')
        for src, tgt, etype in NOVEL_EDGES:
            if src in NOVEL_COMPOUNDS and tgt in POS:
                draw_edge(ax, NOVEL_COMPOUNDS[src], POS[tgt], etype)

    # Legend
    legend_items = [
        mpatches.Patch(color='#1565C0', label='State variable'),
        mpatches.Patch(color='#546E7A', label='Signaling intermediate'),
        mpatches.Patch(color='#37474F', label='BioAge output'),
    ]
    if show_itp:
        legend_items.append(mpatches.Patch(color='#2E7D32', label='ITP compound'))
    if show_novel:
        legend_items.append(mpatches.Patch(color='#C62828', label='Novel target'))

    ax.legend(handles=legend_items, loc='lower right', fontsize=7,
              framealpha=0.9, edgecolor='#ccc')

    plt.tight_layout()
    outpath = FIGDIR / filename
    plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  {filename}')


if __name__ == '__main__':
    print('Generating system diagrams...')

    # Paper 1: Base model architecture
    make_diagram(
        'ELM System Architecture',
        show_itp=False, show_novel=False,
        filename='system_diagram_base.png'
    )

    # Paper 1 alt: With ITP interventions
    make_diagram(
        'ELM System Architecture with ITP Compound Entry Points',
        show_itp=True, show_novel=False,
        filename='system_diagram_itp.png'
    )

    # Paper 2: ITP + novel compounds
    make_diagram(
        'ELM System Architecture: ITP Compounds + Novel Targets',
        show_itp=True, show_novel=True,
        filename='system_diagram_all.png'
    )

    print('Done.')
