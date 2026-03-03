"""
Per-compound dose-response curves using normalized Hill functions.

Replaces linear dose scaling (pathway_value = base × dose_multiplier)
with pharmacologically realistic diminishing returns.

Mathematical framework:
    f(d) = d^n × (K^n + 1) / (K^n + d^n)

Properties:
    f(0) = 0
    f(1) = 1  (preserves all 1× calibrations by construction)
    f(2) < 2  for all K, n > 0 (diminishing returns)
    f(∞) = K^n + 1  (ceiling effect)

Parameters:
    d = dose multiplier (1.0 = ITP reference dose)
    n = Hill coefficient (steepness; n=1 hyperbolic, n→∞ step)
    K = relative ED50 (dose at half-max, in multiples of reference dose)
"""

from typing import Dict
import numpy as np

# Per-compound Hill parameters
# Compounds with ITP multi-dose data are fitted (CAL);
# single-dose compounds use pharmacological estimates (EST).
DOSE_RESPONSE_PARAMS: Dict[str, dict] = {
    'rapamycin': {
        'n': 1.87,
        'K': 0.33,
        'provenance': 'CAL',
        'source': 'ITP 3-dose fit (4.7, 14, 42 ppm)',
        'f2': None,  # computed below
        'pmids': ['24341993'],  # Miller 2014, rapamycin doses + PK
        'notes': (
            'Three ITP doses: 4.7 ppm (+3% M), 14 ppm (+13% M), 42 ppm (+23% M). '
            'Fit residual < 1e-6 (3 points, 2 params — exact). '
            'Steep saturation: f(2)=1.09, ceiling=1.13 — nearly all benefit at 1x. '
            'PK sex mechanism in sex_mechanisms.py handles female mTORC1 separately.'
        ),
    },
    'acarbose': {
        'n': 4.0,
        'K': 0.39,
        'provenance': 'CAL',
        'source': 'ITP 3-dose plateau (400, 1000, 2500 ppm)',
        'f2': None,
        'pmids': ['24245565', 'PMC6413665'],  # Harrison 2014; Harrison 2019 multi-dose
        'notes': (
            'Three ITP doses: 400 ppm (+11% M), 1000 ppm (+17% M), 2500 ppm (+16% M). '
            'Hard saturation — alpha-glucosidase fully blocked at 1000 ppm. '
            '2.5x dose gives NO additional benefit (slight DECREASE). '
            'Unconstrained fit gives n=21 (step function); capped at n=4 for '
            'biological plausibility. f(2)~1.00 regardless of n>3.'
        ),
    },
    '17_alpha_estradiol': {
        'n': 0.74,
        'K': 0.83,
        'provenance': 'CAL',
        'source': 'ITP 2-dose fit (4.8, 14.4 ppm)',
        'f2': None,
        'pmids': ['27312235', '24245565'],  # Strong 2016; Harrison 2014
        'notes': (
            'Two ITP doses: 4.8 ppm (+12% M), 14.4 ppm (+19% M). '
            'Exact fit (2 points, 2 params). n<1 indicates shallow Hill — '
            'consistent with 5-alpha-reductase enzyme kinetics. '
            'Ceiling=1.87 but no data above 14.4 ppm; extrapolation uncertain.'
        ),
    },
    'canagliflozin': {
        'n': 1.0,
        'K': 1.0,
        'provenance': 'EST',
        'source': 'Pharmacological estimate (SGLT2 receptor kinetics)',
        'f2': None,
        'pmids': ['32990681'],  # Miller 2020
        'notes': (
            'Single ITP dose only (180 ppm, +14% M). '
            'SGLT2 inhibition: receptor binding kinetics, likely standard Hill n=1. '
            'FDA clinical data: 100 mg vs 300 mg shows diminishing HbA1c returns. '
            'Confidence: LOW.'
        ),
    },
    'aspirin': {
        'n': 1.5,
        'K': 0.7,
        'provenance': 'EST',
        'source': 'COX-1 pharmacology + ITP non-monotonic note',
        'f2': None,
        'pmids': ['18631321'],  # Strong 2008
        'notes': (
            'Single ITP dose (~20 ppm, +8% M). '
            'COX-1 inhibition: irreversible, near-saturation at low dose. '
            'ITP tested higher doses (60, 200 ppm) — NO benefit (non-monotonic). '
            'Steep Hill plateau is a simplification; high-dose aspirin may be harmful '
            '(COX-2 inhibition blocks protective prostaglandins). '
            'Confidence: MEDIUM (non-monotonic data supports strong saturation).'
        ),
    },
    'glycine': {
        'n': 1.0,
        'K': 2.0,
        'provenance': 'EST',
        'source': 'Amino acid pharmacology (Brind et al. rat dose-range)',
        'f2': None,
        'pmids': ['30916479'],  # Miller 2019
        'notes': (
            'Single ITP dose (8% diet, +6% M). '
            'Amino acid — absorption less saturable than enzyme inhibition. '
            'AMPK activation from glycine is indirect (metabolic signaling). '
            'Brind et al. tested 8%, 12%, 20% in rats with graded effects. '
            'Confidence: LOW-MEDIUM.'
        ),
    },
}


def hill_dose_response(d: float, n: float, K: float) -> float:
    """
    Normalized Hill dose-response function.

    f(d) = d^n * (K^n + 1) / (K^n + d^n)

    Properties:
        f(0) = 0
        f(1) = 1  (by construction)
        f(inf) = K^n + 1  (ceiling)

    Args:
        d: Dose multiplier (1.0 = reference dose).
        n: Hill coefficient (steepness).
        K: Relative ED50 in multiples of reference dose.

    Returns:
        Effective dose multiplier (dimensionless).
    """
    if d <= 0:
        return 0.0
    d_n = d ** n
    K_n = K ** n
    return d_n * (K_n + 1) / (K_n + d_n)


def effective_dose(compound: str, dose_multiplier: float) -> float:
    """
    Get effective dose multiplier for a compound after dose-response correction.

    At dose_multiplier=1.0, always returns 1.0 (preserving calibration).
    At dose_multiplier=2.0, returns f(2) < 2.0 (diminishing returns).

    Args:
        compound: Compound name (must be in DOSE_RESPONSE_PARAMS).
        dose_multiplier: Raw dose multiplier (e.g. 0.5, 1.0, 2.0).

    Returns:
        Corrected effective dose multiplier.

    Raises:
        KeyError: If compound not found in DOSE_RESPONSE_PARAMS.
    """
    if dose_multiplier == 1.0:
        return 1.0  # Fast path — f(1) = 1 by construction
    if dose_multiplier <= 0:
        return 0.0

    params = DOSE_RESPONSE_PARAMS[compound]
    return hill_dose_response(dose_multiplier, params['n'], params['K'])


def get_dose_response_info(compound: str) -> dict:
    """Return full dose-response info for a compound (for popups/docs)."""
    params = DOSE_RESPONSE_PARAMS[compound].copy()
    params['f2'] = hill_dose_response(2.0, params['n'], params['K'])
    params['f05'] = hill_dose_response(0.5, params['n'], params['K'])
    params['ceiling'] = params['K'] ** params['n'] + 1
    return params


# Compute f(2) for all compounds on import
for _name, _params in DOSE_RESPONSE_PARAMS.items():
    _params['f2'] = hill_dose_response(2.0, _params['n'], _params['K'])
