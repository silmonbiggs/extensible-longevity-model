"""
ELM Plotting Module

Publication-quality figure generation for ELM model results.
All plots use consistent styling suitable for journal publication.

Plot Types:
1. Trajectory plots - BioAge, NAD+, senescence over time
2. Survival curves - Mouse and human with uncertainty bands
3. Attribution charts - Intervention contributions
4. ITP calibration - Model vs experimental comparison
5. Sex difference visualization - Mechanism explanations
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path

from .model import SimulationResult, run_control, calculate_lifespan_extension
from .uncertainty import MonteCarloResult, HumanPrediction, AttributionResult
from .compounds import ITP_VALIDATION, HUMAN_BASELINE_LIFESPAN
from .pathways import HUMAN_SCALING_FACTOR


# =============================================================================
# PLOT STYLING
# =============================================================================

# Color palette for consistent styling
COLORS = {
    'control': '#2c3e50',      # Dark blue-gray
    'treated': '#27ae60',      # Green
    'ci_band': '#27ae60',      # Green (transparent)
    'male': '#3498db',         # Blue
    'female': '#e74c3c',       # Red
    'positive': '#27ae60',     # Green
    'negative': '#e74c3c',     # Red
    'neutral': '#95a5a6',      # Gray
    'threshold': '#e74c3c',    # Red
}

# Default figure settings
DEFAULT_FIGSIZE = (10, 6)
DEFAULT_DPI = 150
FONT_SIZE = 12


def _setup_style():
    """Apply consistent plot styling."""
    plt.rcParams.update({
        'font.size': FONT_SIZE,
        'axes.labelsize': FONT_SIZE,
        'axes.titlesize': FONT_SIZE + 2,
        'legend.fontsize': FONT_SIZE - 2,
        'xtick.labelsize': FONT_SIZE - 1,
        'ytick.labelsize': FONT_SIZE - 1,
    })


# =============================================================================
# TRAJECTORY PLOTS
# =============================================================================

def plot_bioage_trajectory(
    control: SimulationResult,
    treated: SimulationResult,
    mc_result: Optional[MonteCarloResult] = None,
    title: str = None,
    output_path: str = None,
    show_ci: bool = True,
    figsize: Tuple[int, int] = DEFAULT_FIGSIZE,
) -> plt.Figure:
    """
    Plot biological age trajectory comparing control and treated.

    Args:
        control: Control simulation result
        treated: Treated simulation result
        mc_result: Optional Monte Carlo result for uncertainty bands
        title: Plot title (auto-generated if None)
        output_path: Path to save figure (displays if None)
        show_ci: Whether to show confidence interval bands
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    _setup_style()
    fig, ax = plt.subplots(figsize=figsize)

    # Normalize time to control lifespan
    t_norm = control.t / control.t_death

    # Plot CI bands first (behind lines)
    if mc_result is not None and show_ci and mc_result.BioAge_ci_90 is not None:
        ax.fill_between(
            t_norm,
            mc_result.BioAge_ci_90[0],
            mc_result.BioAge_ci_90[1],
            alpha=0.2,
            color=COLORS['ci_band'],
            label='90% CI'
        )

    # Plot trajectories
    ax.plot(t_norm, control.BioAge, '--', color=COLORS['control'],
            linewidth=2, label='Control')
    ax.plot(t_norm, treated.BioAge, '-', color=COLORS['treated'],
            linewidth=2, label='Treated')

    # Death threshold
    ax.axhline(y=1.0, color=COLORS['threshold'], linestyle=':',
               alpha=0.7, label='Death threshold')

    # Markers for death times
    ax.axvline(x=1.0, color=COLORS['control'], linestyle=':', alpha=0.5)
    ax.axvline(x=treated.t_death / control.t_death, color=COLORS['treated'],
               linestyle=':', alpha=0.5)

    # Labels and formatting
    extension = calculate_lifespan_extension(treated, control)
    if title is None:
        title = f'Biological Age Trajectory\nLifespan Extension: +{extension:.1f}%'
    ax.set_title(title)
    ax.set_xlabel('Normalized Lifetime (fraction of control lifespan)')
    ax.set_ylabel('Biological Age (normalized)')
    ax.set_xlim(0, max(1.5, treated.t_death / control.t_death * 1.1))
    ax.set_ylim(0, 1.2)
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=DEFAULT_DPI, bbox_inches='tight')
        print(f"Saved: {output_path}")

    return fig


def plot_intermediates(
    control: SimulationResult,
    treated: SimulationResult,
    output_path: str = None,
    figsize: Tuple[int, int] = (15, 10),
) -> plt.Figure:
    """
    Plot intermediate model variables (NAD+, senescence, SASP, etc.).

    Args:
        control: Control simulation result
        treated: Treated simulation result
        output_path: Path to save figure
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    _setup_style()
    fig, axes = plt.subplots(2, 3, figsize=figsize)

    # Normalize time
    t_norm = control.t / control.t_death

    variables = [
        ('NAD', 'NAD+ Level', (0, 1.1)),
        ('SenCells', 'Senescent Cells', None),
        ('SASP', 'SASP Level', None),
        ('Methylation', 'Epigenetic Drift', None),
        ('DNA_damage', 'DNA Damage', None),
        ('BioAge', 'Biological Age', (0, 1.2)),
    ]

    for idx, (var, ylabel, ylim) in enumerate(variables):
        ax = axes.flatten()[idx]
        ax.plot(t_norm, getattr(control, var), '--', color=COLORS['control'],
                linewidth=2, label='Control')
        ax.plot(t_norm, getattr(treated, var), '-', color=COLORS['treated'],
                linewidth=2, label='Treated')

        if var == 'BioAge':
            ax.axhline(y=1.0, color=COLORS['threshold'], linestyle=':',
                       linewidth=2, label='Death threshold')

        ax.axvline(x=1.0, color='gray', linestyle=':', alpha=0.5)
        ax.set_xlabel('Normalized Lifetime')
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        if ylim:
            ax.set_ylim(ylim)
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)

    extension = calculate_lifespan_extension(treated, control)
    plt.suptitle(f'Model State Variables\nLifespan Extension: +{extension:.1f}%',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if output_path:
        plt.savefig(output_path, dpi=DEFAULT_DPI, bbox_inches='tight')
        print(f"Saved: {output_path}")

    return fig


# =============================================================================
# SURVIVAL CURVES
# =============================================================================

def plot_survival(
    control: SimulationResult,
    treated: SimulationResult,
    mc_result: Optional[MonteCarloResult] = None,
    species: str = 'mouse',
    sex: str = 'M',
    output_path: str = None,
    figsize: Tuple[int, int] = DEFAULT_FIGSIZE,
) -> plt.Figure:
    """
    Plot survival curves for mouse or human.

    Args:
        control: Control simulation result
        treated: Treated simulation result
        mc_result: Optional Monte Carlo result for uncertainty
        species: 'mouse' or 'human'
        sex: 'M' or 'F' (for human baseline)
        output_path: Path to save figure
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    _setup_style()
    fig, ax = plt.subplots(figsize=figsize)

    extension = calculate_lifespan_extension(treated, control)

    if species == 'mouse':
        # Mouse survival in years (baseline ~2.5 years)
        MOUSE_LIFESPAN = 2.5
        t_years = control.t * MOUSE_LIFESPAN / control.t_death
        control_death_years = MOUSE_LIFESPAN
        treated_death_years = treated.t_death * MOUSE_LIFESPAN / control.t_death

        ax.plot(t_years, control.Survival, '--', color=COLORS['control'],
                linewidth=2, label=f'Control ({control_death_years:.1f}y)')
        ax.plot(t_years, treated.Survival, '-', color=COLORS['treated'],
                linewidth=2, label=f'Treated ({treated_death_years:.1f}y, +{extension:.0f}%)')

        if mc_result is not None and mc_result.Survival_ci_90 is not None:
            ax.fill_between(t_years, mc_result.Survival_ci_90[0], mc_result.Survival_ci_90[1],
                            alpha=0.15, color=COLORS['ci_band'], label='90% CI')

        ax.set_xlim(0, max(5, treated_death_years * 1.3))
        ax.set_xlabel('Age (years)')
        ax.set_title('Mouse Survival Curves')

    else:  # human
        baseline = HUMAN_BASELINE_LIFESPAN[sex]
        human_ext = extension * HUMAN_SCALING_FACTOR
        human_lifespan = baseline * (1 + human_ext / 100)

        t_human = np.linspace(0, 150, 1000)
        sigma = baseline * 0.12

        control_surv = 1.0 / (1.0 + np.exp((t_human - baseline) / sigma))
        treated_surv = 1.0 / (1.0 + np.exp((t_human - human_lifespan) / (human_lifespan * 0.12)))

        ax.plot(t_human, control_surv, '--', color=COLORS['control'],
                linewidth=2, label=f'Control ({baseline:.0f}y)')
        ax.plot(t_human, treated_surv, '-', color=COLORS['treated'],
                linewidth=2, label=f'Treated ({human_lifespan:.1f}y, +{human_ext:.0f}%)')

        ax.axvline(x=baseline, color=COLORS['control'], linestyle=':', alpha=0.5)
        ax.axvline(x=human_lifespan, color=COLORS['treated'], linestyle=':', alpha=0.5)

        ax.set_xlim(0, 130)
        ax.set_xlabel('Age (years)')
        ax.set_title(f'Human Survival Curves ({sex})')

    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('Survival Probability')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='lower left')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=DEFAULT_DPI, bbox_inches='tight')
        print(f"Saved: {output_path}")

    return fig


# =============================================================================
# ATTRIBUTION PLOTS
# =============================================================================

def plot_attribution(
    attribution: Union[AttributionResult, Dict[str, float]],
    total_extension: float = None,
    output_path: str = None,
    figsize: Tuple[int, int] = DEFAULT_FIGSIZE,
) -> plt.Figure:
    """
    Plot intervention attribution (greedy subtraction decomposition).

    Args:
        attribution: AttributionResult or dict of name -> contribution
        total_extension: Total stack extension (used if dict provided)
        output_path: Path to save figure
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    _setup_style()
    fig, ax = plt.subplots(figsize=figsize)

    # Handle both AttributionResult and plain dict
    if isinstance(attribution, AttributionResult):
        contributions = attribution.contributions
        total = attribution.full_stack_extension
        interaction = attribution.interaction_effect
    else:
        contributions = attribution
        total = total_extension or sum(contributions.values())
        interaction = total - sum(contributions.values())

    # Filter and sort
    filtered = {k: v for k, v in contributions.items() if abs(v) >= 0.1}

    if not filtered:
        ax.text(0.5, 0.5, 'All contributions < 0.1%', ha='center', va='center',
                fontsize=12, transform=ax.transAxes)
        ax.axis('off')
    else:
        sorted_items = sorted(filtered.items(), key=lambda x: x[1], reverse=True)
        names = [x[0] for x in sorted_items]
        values = [x[1] for x in sorted_items]
        colors = [COLORS['positive'] if v > 0 else COLORS['negative'] for v in values]

        y_pos = np.arange(len(names))
        ax.barh(y_pos, values, color=colors, edgecolor='black', linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names)
        ax.set_xlabel('Lifespan Extension Contribution (%)')

        # Add value labels
        for i, v in enumerate(values):
            ax.text(v + 0.3 if v >= 0 else v - 0.3, i, f'{v:+.1f}%',
                    va='center', ha='left' if v >= 0 else 'right', fontsize=10)

        ax.axvline(x=0, color='black', linewidth=0.5)
        ax.grid(True, alpha=0.3, axis='x')

        human_ext = total * HUMAN_SCALING_FACTOR
        ax.set_title(f'Intervention Attribution\n'
                     f'Mouse: +{total:.1f}% | Human (est): +{human_ext:.1f}%')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=DEFAULT_DPI, bbox_inches='tight')
        print(f"Saved: {output_path}")

    return fig


# =============================================================================
# ITP CALIBRATION PLOTS
# =============================================================================

def plot_itp_calibration(
    model_results: Dict[str, Dict[str, float]],
    output_path: str = None,
    figsize: Tuple[int, int] = (12, 6),
) -> plt.Figure:
    """
    Plot ITP calibration: model predictions vs experimental targets.

    Args:
        model_results: Dict of compound -> {'M': ext_m, 'F': ext_f}
        output_path: Path to save figure
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    _setup_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    compounds = list(model_results.keys())
    n = len(compounds)
    x = np.arange(n)
    width = 0.35

    # Get targets from ITP_VALIDATION
    targets_m = []
    targets_f = []
    model_m = []
    model_f = []

    for compound in compounds:
        key = compound.lower().replace('-', '_').replace(' ', '_')
        if key in ITP_VALIDATION:
            targets_m.append(ITP_VALIDATION[key].target_male)
            targets_f.append(ITP_VALIDATION[key].target_female)
        else:
            targets_m.append(0)
            targets_f.append(0)
        model_m.append(model_results[compound]['M'])
        model_f.append(model_results[compound]['F'])

    # Bar chart comparison
    ax1.bar(x - width/2, targets_m, width, label='ITP Target', color=COLORS['control'], alpha=0.7)
    ax1.bar(x + width/2, model_m, width, label='ELM Model', color=COLORS['treated'], alpha=0.7)
    ax1.set_ylabel('Lifespan Extension (%)')
    ax1.set_title('Male')
    ax1.set_xticks(x)
    ax1.set_xticklabels(compounds, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

    ax2.bar(x - width/2, targets_f, width, label='ITP Target', color=COLORS['control'], alpha=0.7)
    ax2.bar(x + width/2, model_f, width, label='ELM Model', color=COLORS['treated'], alpha=0.7)
    ax2.set_ylabel('Lifespan Extension (%)')
    ax2.set_title('Female')
    ax2.set_xticks(x)
    ax2.set_xticklabels(compounds, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')

    # Calculate MAE
    all_targets = targets_m + targets_f
    all_model = model_m + model_f
    mae = np.mean(np.abs(np.array(all_targets) - np.array(all_model)))

    plt.suptitle(f'ITP Calibration: ELM vs Experimental\nMAE: {mae:.2f}%',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    if output_path:
        plt.savefig(output_path, dpi=DEFAULT_DPI, bbox_inches='tight')
        print(f"Saved: {output_path}")

    return fig


def plot_itp_scatter(
    model_results: Dict[str, Dict[str, float]],
    output_path: str = None,
    figsize: Tuple[int, int] = (8, 8),
) -> plt.Figure:
    """
    Plot ITP calibration as scatter (model vs target).

    Args:
        model_results: Dict of compound -> {'M': ext_m, 'F': ext_f}
        output_path: Path to save figure
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    _setup_style()
    fig, ax = plt.subplots(figsize=figsize)

    targets = []
    predictions = []
    labels = []
    colors_list = []

    for compound, results in model_results.items():
        key = compound.lower().replace('-', '_').replace(' ', '_')
        if key in ITP_VALIDATION:
            # Male
            targets.append(ITP_VALIDATION[key].target_male)
            predictions.append(results['M'])
            labels.append(f'{compound} M')
            colors_list.append(COLORS['male'])
            # Female
            targets.append(ITP_VALIDATION[key].target_female)
            predictions.append(results['F'])
            labels.append(f'{compound} F')
            colors_list.append(COLORS['female'])

    targets = np.array(targets)
    predictions = np.array(predictions)

    # Scatter plot
    ax.scatter(targets, predictions, c=colors_list, s=100, alpha=0.7, edgecolors='black')

    # Identity line
    max_val = max(max(targets), max(predictions)) * 1.1
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='Perfect calibration')

    # Add labels
    for i, label in enumerate(labels):
        ax.annotate(label, (targets[i], predictions[i]), fontsize=8,
                    xytext=(5, 5), textcoords='offset points')

    # Statistics
    mae = np.mean(np.abs(targets - predictions))
    r2 = 1 - np.sum((targets - predictions)**2) / np.sum((targets - np.mean(targets))**2)

    ax.set_xlabel('ITP Target (%)')
    ax.set_ylabel('ELM Model (%)')
    ax.set_title(f'ITP Calibration\nMAE: {mae:.2f}%, R²: {r2:.3f}')
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # Legend for sex
    legend_elements = [
        Patch(facecolor=COLORS['male'], label='Male'),
        Patch(facecolor=COLORS['female'], label='Female'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=DEFAULT_DPI, bbox_inches='tight')
        print(f"Saved: {output_path}")

    return fig


# =============================================================================
# SEX DIFFERENCE PLOTS
# =============================================================================

def plot_sex_comparison(
    compounds: List[str],
    results: Dict[str, Dict[str, float]],
    output_path: str = None,
    figsize: Tuple[int, int] = (10, 6),
) -> plt.Figure:
    """
    Plot sex-specific effects for multiple compounds.

    Args:
        compounds: List of compound names
        results: Dict of compound -> {'M': ext_m, 'F': ext_f}
        output_path: Path to save figure
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    _setup_style()
    fig, ax = plt.subplots(figsize=figsize)

    n = len(compounds)
    x = np.arange(n)
    width = 0.35

    males = [results[c]['M'] for c in compounds]
    females = [results[c]['F'] for c in compounds]

    ax.bar(x - width/2, males, width, label='Male', color=COLORS['male'], alpha=0.8)
    ax.bar(x + width/2, females, width, label='Female', color=COLORS['female'], alpha=0.8)

    ax.set_ylabel('Lifespan Extension (%)')
    ax.set_title('Sex-Specific Effects')
    ax.set_xticks(x)
    ax.set_xticklabels(compounds, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    ax.axhline(y=0, color='black', linewidth=0.5)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=DEFAULT_DPI, bbox_inches='tight')
        print(f"Saved: {output_path}")

    return fig


# =============================================================================
# HUMAN PREDICTION PLOTS
# =============================================================================

def plot_human_prediction(
    prediction: HumanPrediction,
    output_path: str = None,
    figsize: Tuple[int, int] = (10, 5),
) -> plt.Figure:
    """
    Plot human prediction uncertainty distribution.

    Args:
        prediction: HumanPrediction from predict_human_extension
        output_path: Path to save figure
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    _setup_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Histogram of human extension samples
    ax1.hist(prediction.samples, bins=50, density=True, alpha=0.7, color='steelblue')
    ax1.axvline(prediction.human_extension_median, color='red', linewidth=2,
                label=f'Median: {prediction.human_extension_median:.1f}%')
    ax1.axvline(prediction.human_extension_ci_90[0], color='orange', linewidth=1.5,
                linestyle='--', label=f'90% CI: [{prediction.human_extension_ci_90[0]:.1f}, '
                                      f'{prediction.human_extension_ci_90[1]:.1f}]%')
    ax1.axvline(prediction.human_extension_ci_90[1], color='orange', linewidth=1.5,
                linestyle='--')
    ax1.set_xlabel('Human Lifespan Extension (%)')
    ax1.set_ylabel('Probability Density')
    ax1.set_title('Human Prediction Uncertainty')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Box plot comparing mouse vs human
    data = [
        [prediction.mouse_extension],
        prediction.samples
    ]
    bp = ax2.boxplot(data, labels=['Mouse\n(observed)', 'Human\n(predicted)'],
                     patch_artist=True)
    bp['boxes'][0].set_facecolor(COLORS['treated'])
    bp['boxes'][1].set_facecolor('steelblue')
    ax2.set_ylabel('Lifespan Extension (%)')
    ax2.set_title(f'Cross-Species Translation\n(dampening: {prediction.dampening_factor:.2f})')
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=DEFAULT_DPI, bbox_inches='tight')
        print(f"Saved: {output_path}")

    return fig


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def save_all_paper_figures(
    output_dir: str,
    control: SimulationResult,
    treated: SimulationResult,
    mc_result: Optional[MonteCarloResult] = None,
    attribution: Optional[AttributionResult] = None,
    itp_results: Optional[Dict[str, Dict[str, float]]] = None,
) -> None:
    """
    Generate all standard paper figures.

    Args:
        output_dir: Directory to save figures
        control: Control simulation
        treated: Treated simulation
        mc_result: Optional Monte Carlo results
        attribution: Optional attribution results
        itp_results: Optional ITP calibration results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Generating paper figures...")

    # Trajectory
    plot_bioage_trajectory(control, treated, mc_result,
                           output_path=output_dir / 'bioage_trajectory.png')

    # Intermediates
    plot_intermediates(control, treated,
                       output_path=output_dir / 'model_intermediates.png')

    # Survival
    plot_survival(control, treated, mc_result, species='mouse',
                  output_path=output_dir / 'survival_mouse.png')
    plot_survival(control, treated, species='human', sex='M',
                  output_path=output_dir / 'survival_human.png')

    # Attribution
    if attribution:
        plot_attribution(attribution,
                         output_path=output_dir / 'attribution.png')

    # ITP calibration
    if itp_results:
        plot_itp_calibration(itp_results,
                             output_path=output_dir / 'itp_calibration.png')
        plot_itp_scatter(itp_results,
                         output_path=output_dir / 'itp_scatter.png')

    print(f"Saved all figures to: {output_dir}")
