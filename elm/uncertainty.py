"""
ELM Uncertainty Quantification

This module provides Monte Carlo methods for quantifying uncertainty in
ELM model predictions. Key features:

1. Trajectory uncertainty - Confidence bands for biological age trajectories
2. Lifespan uncertainty - Confidence intervals for lifespan extension
3. Human translation - Cross-species prediction with uncertainty
4. Intervention attribution - Greedy subtraction decomposition

Uncertainty Sources:
- Individual variation (CV ~12% for lifespan)
- Trajectory stochasticity (CV ~15% for biomarkers)
- Cross-species translation (CV ~40% for dampening factors)
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

from .model import simulate, run_control, calculate_lifespan_extension
from .pathways import CROSS_SPECIES_DAMPENING, INTERVENTION_DAMPENING, HUMAN_SCALING_FACTOR


# =============================================================================
# RESULT DATACLASSES
# =============================================================================

@dataclass
class MonteCarloResult:
    """Container for Monte Carlo simulation results."""
    # Lifespan statistics
    t_death_median: float
    t_death_ci_90: Tuple[float, float]  # (5th, 95th percentile)
    t_death_ci_50: Tuple[float, float]  # (25th, 75th percentile)
    t_death_samples: np.ndarray

    # Extension statistics (if control provided)
    extension_median: Optional[float] = None
    extension_ci_90: Optional[Tuple[float, float]] = None

    # Trajectory statistics (time series)
    t: Optional[np.ndarray] = None
    BioAge_median: Optional[np.ndarray] = None
    BioAge_ci_90: Optional[Tuple[np.ndarray, np.ndarray]] = None
    NAD_median: Optional[np.ndarray] = None
    NAD_ci_90: Optional[Tuple[np.ndarray, np.ndarray]] = None
    Survival_median: Optional[np.ndarray] = None
    Survival_ci_90: Optional[Tuple[np.ndarray, np.ndarray]] = None


@dataclass
class HumanPrediction:
    """Container for human lifespan prediction."""
    mouse_extension: float
    human_extension_median: float
    human_extension_ci_90: Tuple[float, float]
    human_extension_ci_50: Tuple[float, float]
    dampening_factor: float
    intervention_type: str
    samples: np.ndarray


@dataclass
class AttributionResult:
    """Container for intervention attribution results."""
    contributions: Dict[str, float]
    full_stack_extension: float
    interaction_effect: float


# =============================================================================
# MONTE CARLO PARAMETERS
# =============================================================================

# Coefficients of variation for uncertainty sources
LIFESPAN_CV = 0.12      # Individual lifespan variation
TRAJECTORY_CV = 0.15    # Biomarker trajectory variation
DAMPENING_CV = 0.40     # Cross-species translation uncertainty


# =============================================================================
# CORE MONTE CARLO FUNCTIONS
# =============================================================================

def run_monte_carlo(
    interventions: Dict[str, float] = None,
    compound: str = None,
    sex: str = 'M',
    n_samples: int = 500,
    seed: int = 42,
    t_max: float = 1.5,
    dt: float = 0.001,
    include_trajectories: bool = True,
    verbose: bool = False
) -> MonteCarloResult:
    """
    Run Monte Carlo simulation for uncertainty quantification.

    Uses a fast resampling approach: run one deterministic simulation,
    then resample with time-scaling and trajectory perturbations.

    Args:
        interventions: Dict of pathway -> injection value
        compound: Compound name (alternative to interventions)
        sex: 'M' or 'F'
        n_samples: Number of Monte Carlo samples
        seed: Random seed for reproducibility
        t_max: Maximum simulation time
        dt: Time step
        include_trajectories: Whether to compute trajectory statistics
        verbose: Print progress

    Returns:
        MonteCarloResult with lifespan and trajectory statistics

    Example:
        >>> mc = run_monte_carlo(compound='rapamycin', sex='M', n_samples=1000)
        >>> print(f"Extension: {mc.extension_median:.1f}% "
        ...       f"(90% CI: {mc.extension_ci_90[0]:.1f}-{mc.extension_ci_90[1]:.1f}%)")
    """
    rng = np.random.default_rng(seed)

    # Run base simulation
    if compound is not None:
        base_result = simulate(compound=compound, sex=sex, t_max=t_max, dt=dt)
    else:
        base_result = simulate(interventions=interventions, sex=sex, t_max=t_max, dt=dt)

    # Run control for extension calculation
    control = run_control(sex=sex, t_max=t_max, dt=dt)

    t_array = base_result.t
    n = len(t_array)
    base_t_death = base_result.t_death

    # Sample arrays
    t_death_samples = np.zeros(n_samples)

    if include_trajectories:
        bioage_samples = np.zeros((n_samples, n))
        nad_samples = np.zeros((n_samples, n))
        survival_samples = np.zeros((n_samples, n))

    # Generate samples via time-scaling and perturbation
    for i in range(n_samples):
        if verbose and (i + 1) % 100 == 0:
            print(f"  MC sample {i+1}/{n_samples}")

        # Lifespan variation (log-normal)
        lifespan_mult = np.exp(rng.normal(0, LIFESPAN_CV))
        t_death_samples[i] = base_t_death * lifespan_mult

        if include_trajectories:
            # Trajectory variation
            trajectory_mult = np.exp(rng.normal(0, TRAJECTORY_CV))

            # Time-scale trajectories
            time_scale = 1.0 / lifespan_mult
            t_scaled = np.clip(t_array * time_scale, 0, t_array[-1])

            bioage_samples[i] = np.interp(t_scaled, t_array, base_result.BioAge) * trajectory_mult
            nad_samples[i] = np.interp(t_scaled, t_array, base_result.NAD)

            # Recompute survival for this sample's death time
            sigma = 0.12
            survival_samples[i] = 1.0 / (1.0 + np.exp((t_array - t_death_samples[i]) / sigma))

    # Compute extension samples
    extension_samples = (t_death_samples - control.t_death) / control.t_death * 100

    # Build result
    result = MonteCarloResult(
        t_death_median=np.median(t_death_samples),
        t_death_ci_90=(np.percentile(t_death_samples, 5), np.percentile(t_death_samples, 95)),
        t_death_ci_50=(np.percentile(t_death_samples, 25), np.percentile(t_death_samples, 75)),
        t_death_samples=t_death_samples,
        extension_median=np.median(extension_samples),
        extension_ci_90=(np.percentile(extension_samples, 5), np.percentile(extension_samples, 95)),
        t=t_array if include_trajectories else None,
    )

    if include_trajectories:
        result.BioAge_median = np.median(bioage_samples, axis=0)
        result.BioAge_ci_90 = (
            np.percentile(bioage_samples, 5, axis=0),
            np.percentile(bioage_samples, 95, axis=0)
        )
        result.NAD_median = np.median(nad_samples, axis=0)
        result.NAD_ci_90 = (
            np.percentile(nad_samples, 5, axis=0),
            np.percentile(nad_samples, 95, axis=0)
        )
        result.Survival_median = np.median(survival_samples, axis=0)
        result.Survival_ci_90 = (
            np.percentile(survival_samples, 5, axis=0),
            np.percentile(survival_samples, 95, axis=0)
        )

    return result


# =============================================================================
# HUMAN TRANSLATION
# =============================================================================

def predict_human_extension(
    mouse_extension: float,
    intervention_type: str = 'metabolic',
    n_samples: int = 1000,
    seed: int = 42
) -> HumanPrediction:
    """
    Predict human lifespan extension from mouse data with uncertainty.

    Uses pathway-specific dampening factors derived from comparative
    biology and allometric scaling considerations.

    Args:
        mouse_extension: Observed mouse lifespan extension (%)
        intervention_type: Type of intervention for dampening lookup.
            Options: 'mtor_pathway', 'ampk_pathway', 'nad_pathway',
                     'senolytic', 'antioxidant', 'antiinflam', 'metabolic'
        n_samples: Number of samples for uncertainty
        seed: Random seed

    Returns:
        HumanPrediction with median and confidence intervals

    Example:
        >>> pred = predict_human_extension(23.0, 'mtor_pathway')
        >>> print(f"Human prediction: {pred.human_extension_median:.1f}% "
        ...       f"(90% CI: {pred.human_extension_ci_90[0]:.1f}-"
        ...       f"{pred.human_extension_ci_90[1]:.1f}%)")
    """
    rng = np.random.default_rng(seed)

    # Get base dampening factor
    base_dampening = INTERVENTION_DAMPENING.get(intervention_type, HUMAN_SCALING_FACTOR)

    # Sample dampening factors (log-normal distribution)
    dampening_samples = base_dampening * np.exp(rng.normal(0, DAMPENING_CV, n_samples))

    # Compute human extension samples
    human_samples = mouse_extension * dampening_samples

    return HumanPrediction(
        mouse_extension=mouse_extension,
        human_extension_median=np.median(human_samples),
        human_extension_ci_90=(np.percentile(human_samples, 5), np.percentile(human_samples, 95)),
        human_extension_ci_50=(np.percentile(human_samples, 25), np.percentile(human_samples, 75)),
        dampening_factor=base_dampening,
        intervention_type=intervention_type,
        samples=human_samples,
    )


def predict_human_lifespan(
    mouse_extension: float,
    intervention_type: str = 'metabolic',
    sex: str = 'M',
    n_samples: int = 1000,
    seed: int = 42
) -> Dict[str, float]:
    """
    Predict absolute human lifespan with intervention.

    Args:
        mouse_extension: Mouse lifespan extension (%)
        intervention_type: Intervention type for dampening
        sex: 'M' or 'F' for baseline lifespan
        n_samples: Number of samples
        seed: Random seed

    Returns:
        Dict with baseline, predicted lifespan, and extension statistics
    """
    from .compounds import HUMAN_BASELINE_LIFESPAN

    baseline = HUMAN_BASELINE_LIFESPAN[sex]
    pred = predict_human_extension(mouse_extension, intervention_type, n_samples, seed)

    return {
        'baseline_years': baseline,
        'extension_percent_median': pred.human_extension_median,
        'extension_percent_ci_90': pred.human_extension_ci_90,
        'predicted_lifespan_median': baseline * (1 + pred.human_extension_median / 100),
        'predicted_lifespan_ci_90': (
            baseline * (1 + pred.human_extension_ci_90[0] / 100),
            baseline * (1 + pred.human_extension_ci_90[1] / 100)
        ),
        'years_gained_median': baseline * pred.human_extension_median / 100,
    }


# =============================================================================
# INTERVENTION ATTRIBUTION
# =============================================================================

def compute_attribution(
    intervention_components: Dict[str, Dict[str, float]],
    sex: str = 'M',
    t_max: float = 1.5,
    dt: float = 0.001,
    start_time: float = 0.56,
    verbose: bool = True
) -> AttributionResult:
    """
    Compute intervention contributions via greedy subtraction.

    For a stack of interventions, determines how much each component
    contributes by measuring the reduction when that component is removed.

    Args:
        intervention_components: Dict mapping names to pathway injections.
            Example: {'Rapamycin': {'mtorc1_inhibition': 0.85, 'ampk': 0.4},
                      'NMN': {'nmn': 0.5}}
        sex: 'M' or 'F'
        t_max: Maximum simulation time
        dt: Time step
        verbose: Print decomposition results

    Returns:
        AttributionResult with per-intervention contributions

    Example:
        >>> components = {
        ...     'Rapamycin': {'mtorc1_inhibition': 0.85, 'ampk': 0.4},
        ...     'NMN': {'nmn': 0.5},
        ...     'Exercise': {'ampk': 0.16, 'senolytic': 0.06}
        ... }
        >>> attr = compute_attribution(components, sex='M')
        >>> for name, contrib in attr.contributions.items():
        ...     print(f"{name}: {contrib:+.1f}%")
    """
    intervention_names = list(intervention_components.keys())
    n_int = len(intervention_names)

    if verbose:
        print(f"\nComputing greedy subtraction for {n_int} interventions...")
        print(f"  Requires {n_int + 1} simulations")

    # Run control
    control = run_control(sex=sex, t_max=t_max, dt=dt)

    # Build full stack
    full_stack = {'start_time': start_time}
    for component in intervention_components.values():
        for key, value in component.items():
            if key == 'start_time':
                continue
            if key in full_stack:
                full_stack[key] = max(full_stack[key], value)
            else:
                full_stack[key] = value

    # Run full stack
    full_result = simulate(interventions=full_stack, sex=sex, t_max=t_max, dt=dt)
    full_stack_ext = calculate_lifespan_extension(full_result, control)

    if verbose:
        print(f"  Full Stack: +{full_stack_ext:.1f}%")

    # Compute contribution of each intervention
    contributions = {}
    for intervention_name in intervention_names:
        # Build reduced stack (without this intervention)
        reduced_stack = {'start_time': start_time}
        for name, component in intervention_components.items():
            if name == intervention_name:
                continue
            for key, value in component.items():
                if key == 'start_time':
                    continue
                if key in reduced_stack:
                    reduced_stack[key] = max(reduced_stack[key], value)
                else:
                    reduced_stack[key] = value

        # Run reduced stack
        result = simulate(interventions=reduced_stack, sex=sex, t_max=t_max, dt=dt)
        reduced_ext = calculate_lifespan_extension(result, control)

        # Contribution = full - reduced
        contribution = full_stack_ext - reduced_ext
        contributions[intervention_name] = contribution

    # Compute interaction effect
    total_contributions = sum(contributions.values())
    interaction_effect = full_stack_ext - total_contributions

    if verbose:
        print("\nGreedy Subtraction Decomposition:")
        print("-" * 50)
        for name, value in sorted(contributions.items(), key=lambda x: -x[1]):
            print(f"  {name:20s}: {value:+6.1f}%")
        print("-" * 50)
        print(f"  {'SUM':20s}: {total_contributions:+6.1f}%")
        print(f"  {'Full Stack (actual)':20s}: {full_stack_ext:+6.1f}%")
        print(f"  {'Interaction/Overlap':20s}: {interaction_effect:+6.1f}%")

    return AttributionResult(
        contributions=contributions,
        full_stack_extension=full_stack_ext,
        interaction_effect=interaction_effect,
    )


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def compute_confidence_interval(
    samples: np.ndarray,
    confidence: float = 0.90
) -> Tuple[float, float]:
    """
    Compute confidence interval from samples.

    Args:
        samples: Array of samples
        confidence: Confidence level (0-1)

    Returns:
        (lower, upper) bounds
    """
    alpha = 1 - confidence
    lower = np.percentile(samples, 100 * alpha / 2)
    upper = np.percentile(samples, 100 * (1 - alpha / 2))
    return (lower, upper)


def summarize_uncertainty(mc_result: MonteCarloResult) -> str:
    """
    Generate human-readable summary of Monte Carlo results.

    Args:
        mc_result: MonteCarloResult from run_monte_carlo

    Returns:
        Formatted string summary
    """
    lines = [
        "Monte Carlo Uncertainty Summary",
        "=" * 40,
        f"Lifespan (normalized):",
        f"  Median: {mc_result.t_death_median:.3f}",
        f"  90% CI: {mc_result.t_death_ci_90[0]:.3f} - {mc_result.t_death_ci_90[1]:.3f}",
        f"  IQR:    {mc_result.t_death_ci_50[0]:.3f} - {mc_result.t_death_ci_50[1]:.3f}",
    ]

    if mc_result.extension_median is not None:
        lines.extend([
            f"\nLifespan Extension:",
            f"  Median: {mc_result.extension_median:.1f}%",
            f"  90% CI: {mc_result.extension_ci_90[0]:.1f}% - {mc_result.extension_ci_90[1]:.1f}%",
        ])

    return "\n".join(lines)
