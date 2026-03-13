"""
ELM Sex-Specific Mechanisms

This module implements mechanistic explanations for sex differences in
longevity intervention effects, as observed in ITP data.

Four Mechanism Types:
1. AMPK Saturation - Females have higher baseline AMPK (LKB1), causing
   diminishing returns for AMPK activators (acarbose, canagliflozin, glycine)

2. Pharmacokinetic - Sex differences in drug metabolism lead to different
   blood levels (rapamycin: females have 3.48x higher levels)

3. Testosterone Dependent - Some interventions require testosterone as a
   substrate or cofactor (17-alpha-estradiol works via 5α-reductase)

4. Estrogen Interference - Aspirin suppresses the COX-2→PGI2 pathway that
   estrogen uses for cardioprotection in females, cancelling its benefit

Literature Support:
- AMPK/LKB1: McInnes 2011 (PMID 21876548) - estrogen regulation of LKB1
- Rapamycin PK: Miller 2014 (PMID 24341993) - sex-specific blood levels
- 17α-E2 mechanism: Garratt 2018, 2019 (PMIDs 29806096, 30740872)
- Aspirin/estrogen: Egan 2004 (PMID 15550624) - COX-2→PGI2 atheroprotection
  Grosser 2006 (PMID 17518513) - NSAIDs reverse HRT cardioprotection
"""

from typing import Dict, Optional, Any
from enum import Enum
from dataclasses import dataclass


# =============================================================================
# MECHANISM TYPES
# =============================================================================

class MechanismType(str, Enum):
    """Types of sex-specific mechanisms."""
    NONE = 'none'
    PHARMACOKINETIC = 'pharmacokinetic'
    AMPK_SATURATION = 'ampk_saturation'
    TESTOSTERONE_DEPENDENT = 'testosterone_dependent'
    ESTROGEN_INTERFERENCE = 'estrogen_interference'
    UNKNOWN = 'unknown'


# =============================================================================
# SEX MECHANISM DATA
# =============================================================================

@dataclass
class SexMechanism:
    """Sex-specific mechanism for a compound."""
    mechanism_type: MechanismType
    male_effect: float
    female_effect: float
    explanation: str
    pmids: list
    confidence: str  # 'high', 'medium', 'low'

    # Type-specific parameters
    pk_multiplier: Optional[float] = None      # For PHARMACOKINETIC (blood level ratio)
    pk_ampk_boost: Optional[float] = None      # For PHARMACOKINETIC (AMPK multiplier from higher exposure)
    ampk_baseline: Optional[float] = None      # For AMPK_SATURATION
    ampk_k_diminishing: Optional[float] = None # For AMPK_SATURATION
    requires_testosterone: bool = False         # For TESTOSTERONE_DEPENDENT


# Mechanism data for ITP compounds
SEX_MECHANISMS: Dict[str, SexMechanism] = {
    'rapamycin': SexMechanism(
        mechanism_type=MechanismType.PHARMACOKINETIC,
        male_effect=1.0,
        female_effect=1.0,  # Base effect same; PK handles difference
        pk_multiplier=3.48,
        pk_ampk_boost=1.06,  # Higher blood levels also boost AMPK activation
        explanation=(
            'Females have 3.48x higher blood levels due to CYP3A metabolism '
            'differences. Higher exposure leads to greater mTORC1 inhibition '
            'and stronger AMPK activation.'
        ),
        pmids=[24341993],
        confidence='high',
    ),

    'acarbose': SexMechanism(
        mechanism_type=MechanismType.AMPK_SATURATION,
        male_effect=1.0,
        female_effect=0.29,
        ampk_baseline=1.3,
        ampk_k_diminishing=8.0,
        explanation=(
            'Estrogen upregulates LKB1 (AMPK activator). Females have ~1.3x '
            'higher baseline AMPK, causing diminishing returns from activation. '
            'Global AMPK saturation parameters (baseline=1.3, k=8) from '
            'LKB1 literature, not fitted to ITP sex data.'
        ),
        pmids=[21876548, 24245565],
        confidence='medium',
    ),

    'canagliflozin': SexMechanism(
        mechanism_type=MechanismType.AMPK_SATURATION,
        male_effect=1.0,
        female_effect=0.29,
        ampk_baseline=1.3,
        ampk_k_diminishing=8.0,
        explanation=(
            'Dual-action: 88% antioxidant/FGF21 (sex-independent) + 12% AMPK '
            '(subject to female saturation). The antioxidant fraction represents '
            'FGF21/ketogenesis signaling (Ferrannini 2016, Osataphan 2019, '
            'Packer 2020), which passes through the sex modifier at full '
            'strength. Sweep-calibrated f=0.88 minimizes female prediction '
            'error (0.2 pp vs 4.1 pp with pure AMPK).'
        ),
        pmids=[32990681, 30843877, 37406767, 21876548],
        confidence='medium',
    ),

    '17_alpha_estradiol': SexMechanism(
        mechanism_type=MechanismType.TESTOSTERONE_DEPENDENT,
        male_effect=1.0,
        female_effect=0.0,
        requires_testosterone=True,
        explanation=(
            'Works via 5α-reductase inhibition (testosterone→DHT conversion). '
            'Requires testosterone substrate. Castrated males and females '
            'show zero effect.'
        ),
        pmids=[29806096, 30740872],
        confidence='high',
    ),

    'aspirin': SexMechanism(
        mechanism_type=MechanismType.ESTROGEN_INTERFERENCE,
        male_effect=1.0,
        female_effect=0.0,
        explanation=(
            'Male-only effect at low dose (20 ppm). Aspirin is redundant with '
            'or interferes with estrogen\'s natural cardioprotective pathway. '
            'Estrogen activates COX-2→PGI2 for atheroprotection in females '
            '(Egan 2004); aspirin\'s COX inhibition suppresses this. Deleting '
            'the PGI2 receptor removes estrogen\'s protective effect. In '
            'postmenopausal women on HRT, concurrent NSAID use reversed '
            'cardioprotection (OR 0.66→1.50, Grosser 2006). In males (low '
            'estrogen, no COX-2→PGI2 protection to lose), aspirin provides '
            'net anti-inflammatory benefit.'
        ),
        pmids=[18631321, 15550624, 17518513],
        confidence='medium',
    ),

    'glycine': SexMechanism(
        mechanism_type=MechanismType.AMPK_SATURATION,
        male_effect=1.0,
        female_effect=0.29,
        ampk_baseline=1.3,
        ampk_k_diminishing=8.0,
        explanation=(
            'Mild AMPK activator with antioxidant component. Same global '
            'AMPK saturation parameters (baseline=1.3, k=8). The antioxidant '
            'pathway is not subject to AMPK saturation, preserving most of '
            'the female effect.'
        ),
        pmids=[30916479, 21876548],
        confidence='medium',
    ),
}


# =============================================================================
# SEX MODIFIER FUNCTIONS
# =============================================================================

def apply_sex_modifier(
    compound: str,
    injections: Dict[str, float],
    sex: str
) -> Dict[str, float]:
    """
    Apply sex-specific modifier to pathway injections.

    Uses mechanistic models instead of simple multipliers:
    - AMPK_SATURATION: Reduces AMPK-related pathways based on female baseline
    - PHARMACOKINETIC: Adjusts mTORC1 inhibition based on blood level ratio
    - TESTOSTERONE_DEPENDENT: Zero effect for females
    - UNKNOWN: Uses simple effect multiplier

    Args:
        compound: Compound name
        injections: Dict of pathway -> injection value
        sex: 'M' or 'F'

    Returns:
        Modified injections dict
    """
    # Normalize compound name
    key = compound.lower().replace('-', '_').replace(' ', '_')

    if key not in SEX_MECHANISMS:
        return injections.copy()

    mech = SEX_MECHANISMS[key]

    # Males always get full effect
    if sex == 'M':
        return injections.copy()

    # Handle each mechanism type
    if mech.mechanism_type == MechanismType.TESTOSTERONE_DEPENDENT:
        # Categorical zero for females
        return {k: 0.0 for k in injections}

    elif mech.mechanism_type == MechanismType.PHARMACOKINETIC:
        # Higher blood levels boost all pathways proportionally
        # mTORC1: Hill saturation from higher exposure
        # AMPK: linear proportional boost (more drug -> more AMPK activation)
        modified = injections.copy()
        mtorc1_key = ('mtorc1_drug_inhibition' if 'mtorc1_drug_inhibition' in injections
                       else 'mtorc1_inhibition' if 'mtorc1_inhibition' in injections
                       else None)
        if mtorc1_key and mech.pk_multiplier:
            base_inhib = injections[mtorc1_key]
            ec50 = (1.0 - base_inhib) / base_inhib if base_inhib > 0 else 0.176
            new_inhib = mech.pk_multiplier / (ec50 + mech.pk_multiplier)
            modified[mtorc1_key] = min(0.98, new_inhib)
        if mech.pk_ampk_boost and 'ampk' in injections:
            modified['ampk'] = min(1.0, injections['ampk'] * mech.pk_ampk_boost)
        return modified

    elif mech.mechanism_type == MechanismType.AMPK_SATURATION:
        # Apply AMPK diminishing returns
        if mech.ampk_baseline and mech.ampk_k_diminishing:
            dim_factor = 1.0 / (1.0 + mech.ampk_k_diminishing * (mech.ampk_baseline - 1.0))
            modified = {}
            for k, v in injections.items():
                if k in ['ampk', 'gut_microbiome']:
                    modified[k] = v * dim_factor
                else:
                    modified[k] = v
            return modified
        return injections.copy()

    else:
        # Fallback: simple multiplier
        modifier = mech.female_effect
        return {k: v * modifier for k, v in injections.items()}


def get_sex_effect(compound: str, sex: str) -> float:
    """
    Get the sex-specific effect multiplier for a compound.

    Args:
        compound: Compound name
        sex: 'M' or 'F'

    Returns:
        Effect multiplier (1.0 for male, varies for female)
    """
    key = compound.lower().replace('-', '_').replace(' ', '_')

    if key not in SEX_MECHANISMS:
        return 1.0

    mech = SEX_MECHANISMS[key]
    return mech.male_effect if sex == 'M' else mech.female_effect


def get_mechanism_info(compound: str) -> Optional[Dict[str, Any]]:
    """
    Get detailed mechanism information for a compound.

    Args:
        compound: Compound name

    Returns:
        Dict with mechanism type, effects, explanation, PMIDs, or None.
    """
    key = compound.lower().replace('-', '_').replace(' ', '_')

    if key not in SEX_MECHANISMS:
        return None

    mech = SEX_MECHANISMS[key]
    return {
        'mechanism_type': mech.mechanism_type.value,
        'male_effect': mech.male_effect,
        'female_effect': mech.female_effect,
        'explanation': mech.explanation,
        'pmids': mech.pmids,
        'confidence': mech.confidence,
    }


def compute_ampk_saturation_factor(
    female_baseline: float = 1.3,
    k_diminishing: float = 8.0
) -> float:
    """
    Compute the AMPK saturation diminishing factor for females.

    The model: if females have higher baseline AMPK, additional AMPK
    activation has diminishing returns.

    Factor = 1 / (1 + k * (baseline - 1))

    With default parameters (baseline=1.3, k=8):
    Factor = 1 / (1 + 8 * 0.3) = 1/3.4 ≈ 0.294

    Args:
        female_baseline: Female AMPK baseline relative to male (default 1.3)
        k_diminishing: Diminishing returns strength (default 8.0)

    Returns:
        Diminishing factor (0-1)
    """
    return 1.0 / (1.0 + k_diminishing * (female_baseline - 1.0))


# =============================================================================
# SUMMARY FUNCTIONS
# =============================================================================

def print_mechanism_summary():
    """Print a summary of all sex mechanisms."""
    print("=" * 80)
    print("ELM Sex-Specific Mechanisms")
    print("=" * 80)
    print()

    for compound, mech in SEX_MECHANISMS.items():
        print(f"{compound}")
        print(f"  Type: {mech.mechanism_type.value}")
        print(f"  M effect: {mech.male_effect:.2f}, F effect: {mech.female_effect:.2f}")
        print(f"  Confidence: {mech.confidence}")
        print(f"  PMIDs: {mech.pmids}")
        print()


if __name__ == '__main__':
    print_mechanism_summary()
