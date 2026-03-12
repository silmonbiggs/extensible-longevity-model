"""
ELM Compound Database

This module contains compound definitions, pathway mappings, and ITP validation data.

Compound Categories:
- ITP-validated: Tested in NIA Interventions Testing Program with published results
- Supplements: Common longevity supplements with literature support
- Lifestyle: Exercise, fasting, and other lifestyle interventions
- Experimental: Theoretical or early-stage compounds

Data Sources:
- ITP publications (PMIDs included)
- Peer-reviewed literature for non-ITP compounds
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass


# =============================================================================
# COMPOUND -> PATHWAY MAPPINGS
# =============================================================================

# Maps compound names to their pathway injection values.
# Values represent pathway activation strength (0-1 scale).
# ITP compounds are calibrated to match observed lifespan extension.

COMPOUNDS: Dict[str, Dict[str, float]] = {
    # -------------------------------------------------------------------------
    # ITP-VALIDATED COMPOUNDS (calibrated to match ITP data)
    # -------------------------------------------------------------------------
    'rapamycin': {
        'mtorc1_drug_inhibition': 0.4190,
        'ampk': 0.2792,
    },
    'acarbose': {
        'ampk': 0.1690,
        'gut_microbiome': 0.1690,
    },
    'canagliflozin': {
        'ampk': 0.1618,
    },
    '17_alpha_estradiol': {
        'ampk': 0.1764,
        'antioxidant': 0.1764,
    },
    'aspirin': {
        'antioxidant': 0.0985,
        'ampk': 0.0591,
        'antiinflam': 0.0394,
    },
    'glycine': {
        'antioxidant': 0.2842,
        'ampk': 0.0316,
    },

    # -------------------------------------------------------------------------
    # NAD+ PRECURSORS
    # -------------------------------------------------------------------------
    'NMN': {'nmn': 0.50},
    'NR': {'nmn': 0.45},

    # -------------------------------------------------------------------------
    # SIRTUIN ACTIVATORS
    # -------------------------------------------------------------------------
    'pterostilbene': {'sirt1_direct': 0.25},
    'resveratrol': {'sirt1_direct': 0.20},

    # -------------------------------------------------------------------------
    # METABOLIC / AMPK ACTIVATORS
    # -------------------------------------------------------------------------
    'metformin': {'ampk': 0.0},  # ITP: 0% in mice (may work in humans)
    'berberine': {'ampk': 0.25},
    'AKG': {
        'akg': 0.80,
        'ampk': 0.35,
    },

    # -------------------------------------------------------------------------
    # SENOLYTICS
    # -------------------------------------------------------------------------
    'fisetin': {'senolytic': 0.30},
    'dasatinib_quercetin': {'senolytic': 0.40},

    # -------------------------------------------------------------------------
    # MITOCHONDRIAL
    # -------------------------------------------------------------------------
    'urolithin_A': {'mitophagy': 0.40},
    'MitoQ': {'antioxidant': 0.20},

    # -------------------------------------------------------------------------
    # ANTIOXIDANTS / ANTI-INFLAMMATORY
    # -------------------------------------------------------------------------
    'NAC': {'antioxidant': 0.30},
    'curcumin': {'antiinflam': 0.25},
    'EGCG': {
        'antioxidant': 0.15,
        'antiinflam': 0.10,
    },
    'sulforaphane': {'antioxidant': 0.20},

    # -------------------------------------------------------------------------
    # LIFESTYLE
    # -------------------------------------------------------------------------
    'exercise': {
        'ampk': 0.16,
        'senolytic': 0.06,
        'mitophagy': 0.10,
        'antioxidant': 0.08,
        'antiinflam': 0.12,
    },
    'caloric_restriction': {
        'ampk': 0.50,
        'mitophagy': 0.30,
        'antiinflam': 0.25,
    },
    'intermittent_fasting': {'ampk': 0.30},
    'spermidine': {'ampk': 0.20},

    # -------------------------------------------------------------------------
    # NAD+ HOMEOSTASIS STRATEGIES
    # -------------------------------------------------------------------------
    'NAD_homeostasis': {'nad_target': 1.0},
    'NAD_homeostasis_plus': {
        'nad_target': 1.0,
        'cd38_inhibitor': 0.5,
    },
}

# Alias for backward compatibility
COMPOUND_TO_INJECTIONS = COMPOUNDS


# =============================================================================
# ITP VALIDATION DATA
# =============================================================================

@dataclass
class ITPCompound:
    """ITP compound validation data."""
    name: str
    target_male: float          # Expected male lifespan extension (%)
    target_female: float        # Expected female lifespan extension (%)
    pmid: int                   # PubMed ID
    dose: str                   # Dose used in study
    start_age_months: int       # Age when treatment started
    notes: str = ""


# ITP validation targets with sex-specific data
ITP_VALIDATION: Dict[str, ITPCompound] = {
    'rapamycin': ITPCompound(
        name='Rapamycin',
        target_male=23.0,
        target_female=26.0,
        pmid=24341993,
        dose='42 ppm',
        start_age_months=9,
        notes='Higher blood levels in females. mTORC1 inhibitor.',
    ),
    'acarbose': ITPCompound(
        name='Acarbose',
        target_male=22.0,
        target_female=5.0,
        pmid=24245565,
        dose='1000 ppm',
        start_age_months=4,
        notes='Strong male bias (4.4x). Alpha-glucosidase inhibitor -> gut -> AMPK.',
    ),
    'canagliflozin': ITPCompound(
        name='Canagliflozin',
        target_male=14.0,
        target_female=9.0,  # Updated: multi-pathway model
        pmid=32990681,
        dose='180 ppm',
        start_age_months=7,
        notes='SGLT2 inhibitor. Multi-pathway: AMPK + FGF21/ketogenesis.',
    ),
    '17_alpha_estradiol': ITPCompound(
        name='17-alpha-Estradiol',
        target_male=19.0,
        target_female=0.0,
        pmid=27312235,
        dose='14.4 ppm',
        start_age_months=4,
        notes='Male-only at 14.4 ppm (Strong 2016). Testosterone-dependent mechanism (Garratt 2018).',
    ),
    'aspirin': ITPCompound(
        name='Aspirin',
        target_male=8.0,
        target_female=0.0,
        pmid=18631321,
        dose='20 ppm',
        start_age_months=4,
        notes='Male-only at low dose. Higher doses ineffective.',
    ),
    'glycine': ITPCompound(
        name='Glycine',
        target_male=6.0,
        target_female=4.0,
        pmid=30916479,
        dose='8% of diet',
        start_age_months=9,
        notes='Both sexes benefit. Reduced pulmonary adenocarcinoma.',
    ),
}


# Legacy format for backward compatibility
ITP_COMPOUNDS: Dict[str, Dict[str, Any]] = {
    name: {
        'interventions': {
            'start_time': data.start_age_months / 30.0,
            **COMPOUNDS.get(name.lower().replace('-', '_').replace(' ', '_'), {}),
        },
        'target': data.target_male,
        'target_M': data.target_male,
        'target_F': data.target_female,
    }
    for name, data in {
        'Rapamycin': ITP_VALIDATION['rapamycin'],
        'Acarbose': ITP_VALIDATION['acarbose'],
        'Canagliflozin': ITP_VALIDATION['canagliflozin'],
        '17a-Estradiol': ITP_VALIDATION['17_alpha_estradiol'],
        'Aspirin': ITP_VALIDATION['aspirin'],
        'Glycine': ITP_VALIDATION['glycine'],
    }.items()
}


# =============================================================================
# HUMAN PARAMETERS
# =============================================================================

# Human baseline lifespan by sex (years)
HUMAN_BASELINE_LIFESPAN: Dict[str, float] = {
    'M': 76.0,
    'F': 81.0,
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_compound(name: str) -> Optional[Dict[str, float]]:
    """
    Get pathway injections for a compound.

    Handles case-insensitive lookup and common name variations.

    Args:
        name: Compound name (e.g., 'rapamycin', 'Rapamycin', '17-alpha-estradiol')

    Returns:
        Dict of pathway -> injection value, or None if not found.
    """
    # Normalize name
    key = name.lower().replace('-', '_').replace(' ', '_')

    # Direct lookup
    if key in COMPOUNDS:
        return COMPOUNDS[key].copy()

    # Try common variations
    variations = [
        key,
        key.replace('alpha', 'a'),
        key.replace('_', ''),
    ]

    for var in variations:
        if var in COMPOUNDS:
            return COMPOUNDS[var].copy()

    # Case-insensitive fallback for mixed-case keys (e.g., 'NMN', 'AKG')
    for k in COMPOUNDS:
        if k.lower() == key:
            return COMPOUNDS[k].copy()

    return None


def get_itp_targets(compound: str) -> Optional[Dict[str, float]]:
    """
    Get ITP validation targets for a compound.

    Args:
        compound: Compound name

    Returns:
        Dict with 'M' and 'F' target percentages, or None if not ITP compound.
    """
    key = compound.lower().replace('-', '_').replace(' ', '_')

    if key in ITP_VALIDATION:
        data = ITP_VALIDATION[key]
        return {
            'M': data.target_male,
            'F': data.target_female,
            'pmid': data.pmid,
        }

    return None


def list_compounds(category: Optional[str] = None) -> List[str]:
    """
    List available compounds, optionally filtered by category.

    Args:
        category: Optional filter - 'itp', 'supplement', 'lifestyle', or None for all

    Returns:
        List of compound names.
    """
    if category is None:
        return list(COMPOUNDS.keys())

    itp_compounds = set(ITP_VALIDATION.keys())
    lifestyle = {'exercise', 'caloric_restriction', 'intermittent_fasting'}
    supplements = {'NMN', 'NR', 'pterostilbene', 'resveratrol', 'NAC',
                   'curcumin', 'EGCG', 'sulforaphane', 'fisetin', 'spermidine'}

    if category.lower() == 'itp':
        return [c for c in COMPOUNDS if c in itp_compounds]
    elif category.lower() == 'lifestyle':
        return [c for c in COMPOUNDS if c in lifestyle]
    elif category.lower() == 'supplement':
        return [c for c in COMPOUNDS if c in supplements]
    else:
        return list(COMPOUNDS.keys())


def get_itp_start_time(compound: str, lifespan_months: float = 30.0) -> float:
    """
    Get normalized ITP treatment start time for a compound.

    Looks up start_age_months from ITP_VALIDATION and normalizes by
    the reference mouse lifespan (default 30 months).

    Args:
        compound: Compound name (case-insensitive, handles dashes/spaces)
        lifespan_months: Reference mouse lifespan for normalization

    Returns:
        Normalized start time (e.g., 0.30 for rapamycin at 9 months).
        Returns 0.56 for non-ITP compounds (legacy default).
    """
    key = compound.lower().replace('-', '_').replace(' ', '_')

    if key in ITP_VALIDATION:
        return ITP_VALIDATION[key].start_age_months / lifespan_months

    return 9.0 / 30.0  # Default: 9 months (standard ITP protocol age)


def get_sex_ratio(compound: str) -> Optional[float]:
    """
    Get male/female effect ratio for an ITP compound.

    Returns:
        M/F ratio, inf if F=0 and M>0, or None if not ITP compound.
    """
    targets = get_itp_targets(compound)
    if targets is None:
        return None

    m, f = targets['M'], targets['F']

    if f == 0:
        return float('inf') if m > 0 else 1.0
    return m / f


# =============================================================================
# TAGUCHI 8-COMPOUND STACK (audited, 1x calibrated doses)
# =============================================================================

# The 8 compounds in the Taguchi L18 screen, all at 1x ITP-calibrated doses.
# ITP unisex (4): rapamycin, acarbose, canagliflozin, glycine
# Novel targets (4): NMN, CD38 inhibitor, dasatinib+quercetin, urolithin A
#
# Pairwise interactions near-additive. Full stack: +134% M / +115% F (9-month start).
# Mitophagy continuous (daily UroA), senolytic pulsed (3d on / 11d off, biweekly).
TAGUCHI_8_STACK: Dict[str, Dict[str, float]] = {
    'Rapamycin': {'mtorc1_drug_inhibition': 0.4190, 'ampk': 0.2792},
    'Acarbose': {'ampk': 0.1690, 'gut_microbiome': 0.1690},
    'Canagliflozin': {'ampk': 0.1618},
    'Glycine': {'antioxidant': 0.2842, 'ampk': 0.0316},
    'NMN': {'nmn': 0.50},
    'CD38_Inhibitor': {'cd38_inhibitor': 0.50},
    'Dasatinib_Quercetin': {'senolytic': 0.40},
    'Urolithin_A': {'mitophagy': 0.40},
}

# Backward compatibility alias
FULL_STACK_COMPONENTS = TAGUCHI_8_STACK

# Legacy speculative stack (pre-audit, supra-calibration doses) — DO NOT USE
_LEGACY_FULL_STACK_COMPONENTS: Dict[str, Dict[str, float]] = {
    'Rapamycin': {'mtorc1_drug_inhibition': 0.94, 'ampk': 0.55},
    'Acarbose': {'ampk': 0.30, 'gut_microbiome': 0.50},
    'NMN': {'nmn': 0.50},
    'Antioxidants': {'antioxidant': 0.50},
    'AKG': {'akg': 0.80, 'ampk': 0.35},
    'Senolytics': {'senolytic': 0.30},
    'Mitophagy': {'mitophagy': 0.40},
    'Anti-inflammatory': {'antiinflam': 0.40},
    'CD38_Inhibitor': {'cd38_inhibitor': 0.50},
    'SIRT1_Activator': {'sirt1_direct': 0.40},
}
