"""
ELM: Extensible Longevity Model

A mechanistic model of mammalian longevity calibrated to ITP intervention data.
Explains sex-specific effects through three mechanisms: AMPK saturation,
pharmacokinetics, and testosterone dependence.

Basic usage:
    from elm import simulate, run_itp_compound

    # Simulate an ITP compound
    result, extension = run_itp_compound('rapamycin', sex='M')
    print(f"Lifespan extension: {extension:.1f}%")

    # Custom intervention
    result = simulate(interventions={'ampk': 0.5, 'mtorc1_inhibition': 0.3})

For more examples, see the examples/ directory.
"""

__version__ = "1.0.0"
__author__ = "ELM Development Team"

# Core simulation
from .model import (
    simulate,
    run_control,
    run_itp_compound,
    calculate_lifespan_extension,
    SimulationResult,
)

# Compound database
from .compounds import (
    COMPOUNDS,
    ITP_VALIDATION,
    get_compound,
    get_itp_targets,
    list_compounds,
)

# Sex mechanisms
from .sex_mechanisms import (
    MechanismType,
    SEX_MECHANISMS,
    apply_sex_modifier,
    get_sex_effect,
    get_mechanism_info,
)

# Pathway parameters
from .pathways import (
    ALL_PATHWAY_PARAMS,
    NAD_PARAMS,
    AMPK_PARAMS,
    MTORC1_PARAMS,
    SIRTUIN_PARAMS,
    BIOAGE_PARAMS,
)

# Uncertainty quantification
from .uncertainty import (
    run_monte_carlo,
    predict_human_extension,
    predict_human_lifespan,
    compute_attribution,
    MonteCarloResult,
    HumanPrediction,
    AttributionResult,
)

# Plotting
from .plotting import (
    plot_bioage_trajectory,
    plot_intermediates,
    plot_survival,
    plot_attribution,
    plot_itp_calibration,
    plot_itp_scatter,
    plot_sex_comparison,
    plot_human_prediction,
    save_all_paper_figures,
)

__all__ = [
    # Core
    'simulate',
    'run_control',
    'run_itp_compound',
    'calculate_lifespan_extension',
    'SimulationResult',
    # Compounds
    'COMPOUNDS',
    'ITP_VALIDATION',
    'get_compound',
    'get_itp_targets',
    'list_compounds',
    # Sex mechanisms
    'MechanismType',
    'SEX_MECHANISMS',
    'apply_sex_modifier',
    'get_sex_effect',
    'get_mechanism_info',
    # Parameters
    'ALL_PATHWAY_PARAMS',
    # Uncertainty
    'run_monte_carlo',
    'predict_human_extension',
    'predict_human_lifespan',
    'compute_attribution',
    'MonteCarloResult',
    'HumanPrediction',
    'AttributionResult',
    # Plotting
    'plot_bioage_trajectory',
    'plot_intermediates',
    'plot_survival',
    'plot_attribution',
    'plot_itp_calibration',
    'plot_itp_scatter',
    'plot_sex_comparison',
    'plot_human_prediction',
    'save_all_paper_figures',
]
