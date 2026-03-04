"""
ELM Core Simulation Model

This module contains the core ODE-based simulation for the ELM longevity model.
The model integrates 8 biological subsystems:

1. NAD+ Homeostasis - Synthesis, consumption, and supplementation
2. Sirtuins - NAD+-dependent deacetylases (SIRT1, FOXO3)
3. mTORC1 - Growth signaling and autophagy regulation
4. DNA Damage & Repair - ROS damage and NAD+-dependent repair
5. Epigenetics - Methylation drift and demethylation
6. Mitochondria - Heteroplasmy dynamics
7. Senescence - Cell entry, SASP, and clearance
8. Cancer - Mutation accumulation and immune surveillance

Simulation outputs biological age trajectory and lifespan.
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, Any

# Import pathway parameters
from .pathways import (
    NAD_PARAMS,
    NMN_SATURATION_PARAMS,
    CD38_SATURATION_PARAMS,
    AMPK_PARAMS,
    MTORC1_PARAMS,
    SIRTUIN_PARAMS,
    AUTOPHAGY_PARAMS,
    DNA_REPAIR_PARAMS,
    METHYLATION_PARAMS,
    HETEROPLASMY_PARAMS,
    UPRMT_PARAMS,
    SENESCENCE_PARAMS,
    SASP_PARAMS,
    CANCER_PARAMS,
    BIOAGE_PARAMS,
    TSC2_PARAMS,
)

# Import compound and sex mechanism data
from .compounds import COMPOUNDS, get_compound, get_itp_start_time
from .sex_mechanisms import apply_sex_modifier, get_sex_effect


# =============================================================================
# RESULT DATACLASS
# =============================================================================

@dataclass
class SimulationResult:
    """Container for simulation outputs."""
    # Time and lifespan
    t: np.ndarray
    t_death: float
    extension_percent: float

    # Core trajectories
    NAD: np.ndarray
    BioAge: np.ndarray
    Survival: np.ndarray

    # State variables
    SIRT1: np.ndarray
    FOXO3: np.ndarray
    Autophagy: np.ndarray
    DNA_damage: np.ndarray
    Methylation: np.ndarray
    Heteroplasmy: np.ndarray
    SenCells: np.ndarray
    SASP: np.ndarray
    Mutations: np.ndarray

    # Algebraic intermediates
    TSC2: np.ndarray
    mTORC1_activity: np.ndarray

    # Additional outputs
    Cancer_prob: np.ndarray
    cancer_occurred: bool
    cancer_time: Optional[float]

    # Intervention tracking
    interventions: Dict[str, float]

    # Monte Carlo (optional)
    ci_90: Optional[Tuple[float, float]] = None


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def hill_function(x: float, Km: float, n: float, Vmax: float = 1.0) -> float:
    """Hill equation for cooperative binding."""
    if x <= 0:
        return 0.0
    return Vmax * (x ** n) / (Km ** n + x ** n)


def michaelis_menten(substrate: float, Km: float, Vmax: float = 1.0) -> float:
    """Michaelis-Menten saturation kinetics."""
    if substrate <= 0:
        return 0.0
    return Vmax * substrate / (Km + substrate)


# =============================================================================
# AMPK FUNCTIONS
# =============================================================================

def calculate_effective_ampk(
    ampk_direct: float,
    gut_contribution: float,
    params: Dict = None
) -> float:
    """
    Calculate effective AMPK activity with saturation.

    AMPK activation follows Michaelis-Menten kinetics to model
    the saturation observed in ITP compounds (especially in females).
    """
    if params is None:
        params = AMPK_PARAMS
    total_input = ampk_direct + params['k_gut_ampk'] * gut_contribution
    return michaelis_menten(total_input, params['Km_ampk'], params['Vmax_ampk'])


# =============================================================================
# TSC2 / mTORC1 FUNCTIONS
# =============================================================================

def calculate_tsc2(ampk: float, params: Dict = None) -> float:
    """
    Calculate TSC2 activity as an algebraic intermediate.

    TSC2 has a basal GAP activity that is enhanced by AMPK phosphorylation
    at Thr1271 and Ser1387. Uses a Hill function for cooperative activation
    (AMPK priming enables GSK3-beta sequential phosphorylation).

    References:
        Inoki, Zhu & Guan 2003, Cell (PMID 14651849)
        Inoki et al. 2006, Cell (PMID 16959574)
        Huang & Manning 2008, Biochem J (PMID 18466115)
    """
    if params is None:
        params = TSC2_PARAMS

    ampk_activation = hill_function(
        ampk,
        Km=params['Km_ampk_tsc2'],
        n=params['n_hill_ampk_tsc2'],
        Vmax=params['k_ampk_tsc2_max']
    )

    tsc2 = params['TSC2_basal'] + ampk_activation
    return np.clip(tsc2, params['TSC2_min'], 1.0)


def calculate_mtorc1_activity(
    tsc2_activity: float,
    drug_inhibition: float,
    params: Dict = None
) -> float:
    """
    Calculate effective mTORC1 activity.

    mTORC1 is inhibited by two independent mechanisms that stack multiplicatively:
    1. TSC2 GAP activity on Rheb (AMPK -> TSC2 -> Rheb-GDP -> mTORC1 off)
    2. Direct drug inhibition (rapamycin-FKBP12 -> mTOR FRB domain)

    References:
        TSC2->Rheb: Inoki et al. 2003 (PMID 14651849)
        Rapamycin: Mukhopadhyay et al. 2016 (PMID 26839163)
        Dual requirement: Kazyken et al. 2020 (PMID 32912901)
    """
    if params is None:
        params = TSC2_PARAMS

    # Arm 1: TSC2-mediated inhibition (via Rheb inactivation)
    tsc2_inhibition = hill_function(
        tsc2_activity,
        Km=params['Km_tsc2_mtorc1'],
        n=params['n_hill_tsc2_mtorc1'],
        Vmax=params['k_tsc2_mtorc1_max']
    )

    # Arm 2: Direct drug inhibition (rapamycin)
    # Multiplicative stacking: independent mechanisms
    mtorc1 = (1.0 - tsc2_inhibition) * (1.0 - drug_inhibition)

    return np.clip(mtorc1, 0.05, 1.0)


def calculate_mtorc1_effects(mtorc1_inhibition: float, params: Dict = None) -> Dict:
    """Calculate downstream effects of mTORC1 inhibition."""
    if params is None:
        params = MTORC1_PARAMS
    return {
        'autophagy_boost': params['k_mtorc1_autophagy'] * mtorc1_inhibition,
        'proteostasis_boost': params['k_mtorc1_proteostasis'] * mtorc1_inhibition,
        'tfeb_activation': params['k_mtorc1_tfeb'] * mtorc1_inhibition,
        'senescence_reduction': params['k_mtorc1_senescence'] * mtorc1_inhibition,
        'damage_reduction': params['k_mtorc1_damage_reduction'] * mtorc1_inhibition,
    }


# =============================================================================
# SASP FUNCTIONS
# =============================================================================

def calculate_sasp_production(
    sen_cells: float,
    sen_duration: float,
    params: Dict = None
) -> float:
    """Calculate SASP production from senescent cells."""
    if params is None:
        params = SASP_PARAMS
    maturation_factor = min(1.0, sen_duration / params['k_sasp_delay']) if params['k_sasp_delay'] > 0 else 1.0
    return params['k_sasp_production'] * sen_cells * maturation_factor


def calculate_sasp_clearance(
    sasp: float,
    liver_function: float,
    antiinflam: float,
    params: Dict = None
) -> float:
    """Calculate SASP clearance rate."""
    if params is None:
        params = SASP_PARAMS
    decay = params['k_sasp_decay'] * sasp
    liver = params['k_sasp_liver_clear'] * liver_function * sasp
    antiinflam_clear = params['k_antiinflam_sasp_reduction'] * antiinflam * sasp
    return decay + liver + antiinflam_clear


def get_sasp_effects(sasp: float, params: Dict = None) -> Dict:
    """Get downstream effects of SASP."""
    if params is None:
        params = SASP_PARAMS
    return {
        'cd38_upregulation': params['k_sasp_cd38'] * sasp,
        'paracrine_senescence': params['k_sasp_paracrine_sen'] * sasp,
        'ros_production': params['k_sasp_ros'] * sasp,
        'immune_suppression': params['k_sasp_immune_suppress'] * sasp,
        'meth_acceleration': params['k_sasp_meth_drift'] * sasp,
    }


# =============================================================================
# UPRmt FUNCTIONS
# =============================================================================

def calculate_uprmt(
    heteroplasmy: float,
    ros_level: float,
    nad_ratio: float,
    params: Dict = None
) -> float:
    """Calculate mitochondrial unfolded protein response."""
    if params is None:
        params = UPRMT_PARAMS
    hetero_stress = params['k_uprmt_hetero'] * heteroplasmy
    ros_stress = params['k_uprmt_ros'] * max(0, ros_level - 1.0)
    nad_stress = params['k_uprmt_nad'] * max(0, 1 - nad_ratio)
    total_stress = hetero_stress + ros_stress + nad_stress
    return michaelis_menten(total_stress, params['Km_uprmt'], params['Vmax_uprmt'])


def get_uprmt_effects(uprmt: float, params: Dict = None) -> Dict:
    """Get downstream effects of UPRmt."""
    if params is None:
        params = UPRMT_PARAMS
    return {
        'proteostasis_boost': params['k_uprmt_proteostasis'] * uprmt,
        'sen_reduction': params['k_uprmt_sen_reduction'] * uprmt,
        'autophagy_boost': params['k_uprmt_autophagy'] * uprmt,
    }


# =============================================================================
# NAD+ FUNCTIONS
# =============================================================================

def compute_nmn_boost_saturated(
    nmn_supplement: float,
    age_factor: float,
    params: Dict = None
) -> Tuple[float, float]:
    """
    Compute NMN boost with Michaelis-Menten saturation.

    Returns (actual_boost, effective_dose).
    """
    if params is None:
        params = NMN_SATURATION_PARAMS
    Vmax = params['Vmax_nmn']
    Km = params['Km_nmn']
    k_decline = params['k_nmn_absorption_decline']

    absorption_efficiency = 1.0 - k_decline * age_factor
    effective_dose = nmn_supplement
    absorbed_nmn = effective_dose * absorption_efficiency

    actual_boost = Vmax * absorbed_nmn / (Km + absorbed_nmn) if absorbed_nmn > 0 else 0
    return actual_boost, effective_dose


def compute_cd38_reduction_saturated(
    cd38_inhibitor: float,
    gut_antiinflam: float = 0,
    params_cd38: Dict = None,
    params_nad: Dict = None
) -> float:
    """Compute CD38 reduction with Hill equation and ceiling."""
    if params_cd38 is None:
        params_cd38 = CD38_SATURATION_PARAMS
    if params_nad is None:
        params_nad = NAD_PARAMS

    Max_inh = params_cd38['Max_cd38_inhibition']
    IC50 = params_cd38['IC50_cd38_inhibitor']
    n = params_cd38['Hill_n_cd38']

    gut_reduction = gut_antiinflam * params_nad['k_gut_antiinflam']

    if cd38_inhibitor > 0:
        inhibitor_reduction = Max_inh * (cd38_inhibitor ** n) / (IC50 ** n + cd38_inhibitor ** n)
    else:
        inhibitor_reduction = 0

    total_reduction = min(Max_inh, inhibitor_reduction + gut_reduction * 0.5)
    return total_reduction


def compute_nmn_dose_for_nad_target(
    current_nad: float,
    target_nad: float,
    age_factor: float,
    base_synthesis: float,
    consumption: float,
    params: Dict = None
) -> Tuple[float, bool]:
    """
    Calculate NMN dose required to maintain target NAD+ level.

    Uses bang-bang control: dose at maximum when below target,
    steady-state dose when at target.

    Returns (required_dose, is_achievable).
    """
    if params is None:
        params = NMN_SATURATION_PARAMS

    Vmax = params['Vmax_nmn']
    Km = params['Km_nmn']
    k_decline = params['k_nmn_absorption_decline']

    # Steady-state boost needed
    steady_state_boost = consumption - base_synthesis

    # Check achievability
    if steady_state_boost >= Vmax * 0.99:
        return 100.0, False

    # Bang-bang control
    error = target_nad - current_nad

    if error > 0.01:
        required_boost = Vmax * 0.95
    elif error < -0.01:
        required_boost = max(0, steady_state_boost * 0.5)
    else:
        required_boost = steady_state_boost

    if required_boost <= 0:
        return 0.0, True

    # Inverse Michaelis-Menten
    required_absorbed = required_boost * Km / (Vmax - required_boost)
    absorption_eff = 1.0 - k_decline * age_factor
    required_dose = required_absorbed / max(0.3, absorption_eff)

    return required_dose, True


def calculate_nad_synthesis(
    ampk: float,
    nmn_supplement: float,
    age_factor: float,
    heteroplasmy: float,
    params: Dict = None
) -> float:
    """Calculate NAD+ synthesis rate."""
    if params is None:
        params = NAD_PARAMS

    nampt_activity = (params['k_nampt_base']
                      * np.exp(-params['k_nampt_age_decline'] * age_factor)
                      * (1 + params['k_nampt_ampk'] * ampk))
    nampt_activity = max(0.1, nampt_activity)

    oxphos_efficiency = 1 - params['k_hetero_oxphos'] * heteroplasmy
    oxphos_efficiency = max(0.30, oxphos_efficiency)
    oxphos_nad_regen = params['k_oxphos_nad_regen'] * oxphos_efficiency

    de_novo = params['k_tryptophan']
    nmn_boost = params['k_nmn_synthesis'] * nmn_supplement

    return nampt_activity * (0.7 + 0.3 * oxphos_efficiency) + de_novo + nmn_boost + oxphos_nad_regen


def calculate_nad_consumption(
    nad: float,
    dna_damage: float,
    sirt1: float,
    sasp: float,
    age_factor: float,
    cd38_inhibitor: float,
    gut_antiinflam: float,
    params: Dict = None
) -> float:
    """Calculate NAD+ consumption rate.

    Consumption is substrate-dependent: rate × NAD/(NAD + Km).
    This prevents unphysical NAD depletion and creates a natural
    equilibrium at late ages (enzymes can't consume absent substrate).
    """
    if params is None:
        params = NAD_PARAMS

    cd38 = (params['k_cd38_base']
            + params['k_cd38_age'] * age_factor
            + params['k_cd38_sasp'] * sasp)
    cd38_reduction = 1 - cd38_inhibitor * params['k_cd38_inhibitor'] - gut_antiinflam * params['k_gut_antiinflam']
    cd38 *= max(0.2, cd38_reduction)
    cd38 = max(0, cd38)

    parp = params['k_parp_consumption'] * dna_damage
    sirt = params['k_sirt_consumption'] * sirt1
    basal = params['k_consumption_basal']

    rate = basal + cd38 + parp + sirt

    # Michaelis-Menten substrate dependence: consumption → 0 as NAD → 0
    Km_nad_consumption = params.get('Km_nad_consumption', 0.05)
    substrate_factor = nad / (nad + Km_nad_consumption) if nad > 0 else 0.0

    return rate * substrate_factor


# =============================================================================
# DNA DAMAGE & REPAIR
# =============================================================================

def calculate_damage_rate(
    ros_level: float,
    proliferation: float,
    antioxidant: float,
    ampk: float,
    uprmt_proteostasis: float,
    mtorc1_damage_reduction: float,
    sasp_ros: float,
    params: Dict = None
) -> float:
    """Calculate DNA damage accumulation rate."""
    if params is None:
        params = DNA_REPAIR_PARAMS

    ros_protection = 1 - params['k_antioxidant_ros_reduction'] * antioxidant
    ros_protection = max(0.2, ros_protection)

    total_ros = ros_level + sasp_ros
    ros_damage = params['k_damage_ros'] * total_ros * ros_protection
    replication_damage = params['k_damage_replication'] * (1 + proliferation)
    spontaneous = params['k_damage_spontaneous']

    ampk_reduction = max(0.3, 1 - params['k_ampk_damage_reduction'] * ampk)
    uprmt_protection = max(0.7, 1 - uprmt_proteostasis)
    mtorc1_protection = max(0.7, 1 - mtorc1_damage_reduction)

    return (ros_damage + replication_damage + spontaneous) * ampk_reduction * uprmt_protection * mtorc1_protection


def calculate_repair_rate(
    nad: float,
    sirt1: float,
    dna_damage: float,
    params: Dict = None
) -> float:
    """Calculate DNA repair capacity."""
    if params is None:
        params = DNA_REPAIR_PARAMS

    nad_ratio = nad / NAD_PARAMS['NAD_young']
    ber = params['k_ber_base'] * (nad_ratio ** params['k_ber_nad_dependence'])
    ner = params['k_ner_base'] * (nad_ratio ** params['k_ner_nad_dependence'])
    hr = params['k_hr_base'] * (1 + params['k_hr_sirt1_boost'] * sirt1)

    return (ber + ner + hr) * dna_damage


# =============================================================================
# METHYLATION
# =============================================================================

def calculate_methylation_rate(
    ezh2: float,
    dnmt: float,
    damage: float,
    sasp: float,
    sasp_meth_accel: float,
    gut: float,
    params: Dict = None
) -> float:
    """Calculate methylation accumulation rate."""
    if params is None:
        params = METHYLATION_PARAMS

    ezh2_meth = params['k_ezh2_base'] * ezh2 * (0.1 + damage)

    gut_dnmt_factor = max(0.3, 1 - gut * params['k_gut_dnmt_inhibition'])
    dnmt_meth = params['k_dnmt_base'] * dnmt * gut_dnmt_factor

    sasp_meth = params['k_sasp_methylation'] * sasp * (1 + sasp_meth_accel)

    return ezh2_meth + dnmt_meth + sasp_meth


def calculate_demethylation_rate(
    meth: float,
    tet_activity: float,
    osk_active: bool,
    params: Dict = None
) -> float:
    """Calculate methylation removal rate."""
    if params is None:
        params = METHYLATION_PARAMS

    tet_demeth = params['k_tet_base'] * tet_activity * meth

    if osk_active:
        osk_demeth = params['k_osk_demethylation'] * params['k_osk_efficiency'] * meth
    else:
        osk_demeth = 0.0

    return tet_demeth + osk_demeth


# =============================================================================
# HETEROPLASMY
# =============================================================================

def calculate_heteroplasmy_dynamics(
    H: float,
    foxo3: float,
    mitophagy_active: bool,
    mitophagy_strength: float,
    params: Dict = None
) -> float:
    """Kowald-Kirkwood heteroplasmy dynamics."""
    if params is None:
        params = HETEROPLASMY_PARAMS

    time_scale = params['time_scale']
    Neff = params['Neff']

    drift = H * (1 - H) / Neff
    mutation_pressure = params['transcription_advantage'] * H * (1 - H)

    selection = params['k_selection_base'] + params['k_selection_foxo3'] * foxo3
    if mitophagy_active:
        selection += params['k_selection_mitophagy'] * mitophagy_strength
    selection *= H

    return (drift + mutation_pressure - selection) * time_scale


# =============================================================================
# SENESCENCE
# =============================================================================

def calculate_senescence_rate(
    damage: float,
    sirt1: float,
    ampk: float,
    gut: float,
    uprmt_sen_reduction: float,
    mtorc1_sen_reduction: float,
    sasp_paracrine: float,
    params: Dict = None
) -> float:
    """Calculate senescence entry rate."""
    if params is None:
        params = SENESCENCE_PARAMS

    threshold = params.get('damage_senescence_threshold', 0.3)
    width = params.get('damage_senescence_width', 0.05)
    sigmoid_gate = 1.0 / (1.0 + np.exp(-(damage - threshold) / width))
    damage_sen = params['k_sen_from_damage'] * max(0, damage - threshold) * sigmoid_gate
    replicative = params['k_sen_from_telomere']
    oncogene = params['k_sen_from_oncogene']
    paracrine_sen = params['k_sen_from_sasp'] * sasp_paracrine

    sirt1_factor = 1 + 0.5 * (1 - sirt1)
    ampk_reduction = max(0.3, 1 - params['k_ampk_sen_reduction'] * ampk)
    gut_reduction = max(0.3, 1 - params['k_gut_sen_reduction'] * gut)
    uprmt_reduction = max(0.7, 1 - uprmt_sen_reduction)
    mtorc1_reduction = max(0.7, 1 - mtorc1_sen_reduction)

    base_sen = (damage_sen + replicative + oncogene) * sirt1_factor * ampk_reduction * gut_reduction * uprmt_reduction * mtorc1_reduction
    return base_sen + paracrine_sen


def calculate_senescence_clearance(
    sen: float,
    autophagy: float,
    senolytic_active: bool,
    params: Dict = None
) -> float:
    """Calculate senescent cell clearance rate."""
    if params is None:
        params = SENESCENCE_PARAMS

    natural = params['k_sen_natural_clear'] * sen
    auto_clear = params['k_sen_autophagy_clear'] * autophagy * sen

    senolytic = params['k_senolytic_clear'] * sen if senolytic_active else 0.0

    return natural + auto_clear + senolytic


# =============================================================================
# SIRTUINS
# =============================================================================

def calculate_sirt1(nad: float, sirt1_direct: float, params: Dict = None) -> float:
    """Calculate SIRT1 activity."""
    if params is None:
        params = SIRTUIN_PARAMS

    sirt1_nad = hill_function(nad, params['Km_nad_sirt1'], params['n_hill_sirt1'])
    sirt1_boost = sirt1_direct * params['k_sirt1_direct_boost'] * (1 - sirt1_nad)

    return min(1.0, sirt1_nad + sirt1_boost)


def calculate_foxo3(sirt1: float, params: Dict = None) -> float:
    """Calculate FOXO3 activity."""
    if params is None:
        params = SIRTUIN_PARAMS

    return min(1.0, params['foxo3_basal'] + params['k_sirt1_foxo3'] * sirt1)


# =============================================================================
# AUTOPHAGY
# =============================================================================

def calculate_autophagy(
    ampk: float,
    foxo3: float,
    gut: float,
    uprmt_autophagy: float,
    mtorc1_autophagy: float,
    params: Dict = None
) -> float:
    """Calculate autophagy activity."""
    if params is None:
        params = AUTOPHAGY_PARAMS

    autophagy_input = (ampk * params['w_ampk_autophagy'] +
                       foxo3 * params['w_foxo3_autophagy'] +
                       gut * params['w_gut_autophagy'] +
                       uprmt_autophagy + mtorc1_autophagy)

    return michaelis_menten(autophagy_input, params['Km_autophagy'], params['Vmax_autophagy'])


# =============================================================================
# CANCER
# =============================================================================

def calculate_mutation_rate(
    unrepaired_damage: float,
    proliferation: float,
    ros: float,
    params: Dict = None
) -> float:
    """Calculate mutation accumulation rate."""
    if params is None:
        params = CANCER_PARAMS

    damage_mutations = params['k_mut_from_damage'] * unrepaired_damage
    replication_mutations = params['k_mut_from_replication'] * (1 + params['k_mut_proliferation'] * proliferation)
    ros_mutations = params['k_mut_from_ros'] * ros

    return damage_mutations + replication_mutations + ros_mutations


def calculate_mutation_clearance(
    mutations: float,
    immune_function: float,
    p53_activity: float,
    senescent_cells: float,
    params: Dict = None
) -> float:
    """Calculate mutation/pre-cancer clearance rate."""
    if params is None:
        params = CANCER_PARAMS

    immune_clear = params['k_immune_surveillance'] * immune_function * mutations
    apoptosis = params['k_p53_apoptosis'] * p53_activity * mutations

    return immune_clear + apoptosis


def calculate_cancer_probability(mutations: float, params: Dict = None) -> float:
    """Calculate instantaneous cancer probability."""
    if params is None:
        params = CANCER_PARAMS

    if mutations <= 0:
        return 0.0

    ratio = mutations / params['cancer_threshold']
    return 1 - np.exp(-params['k_cancer_rate'] * (ratio ** params['cancer_hits']))


# =============================================================================
# BIOLOGICAL AGE
# =============================================================================

def calculate_bioage(
    meth: float,
    damage: float,
    H: float,
    sen: float,
    sasp: float,
    params: Dict = None
) -> float:
    """Calculate biological age from component values."""
    if params is None:
        params = BIOAGE_PARAMS

    meth_contrib = params['w_meth'] * meth / params['meth_norm']
    damage_contrib = params['w_damage'] * damage / params['damage_norm']
    h_contrib = params['w_hetero'] * H / params['H_norm']
    sen_contrib = params['w_sen'] * sen / params['sen_norm']
    sasp_contrib = params['w_sasp'] * sasp / params['SASP_norm']

    return meth_contrib + damage_contrib + h_contrib + sen_contrib + sasp_contrib


def derive_normalization_constants(dt: float = 0.001) -> Dict[str, float]:
    """
    Derive BioAge normalization constants from control simulation.

    Runs an untreated male control, reads each BioAge component at t=1.0,
    and returns the values as normalization constants. Setting each norm to
    its component's value at t=1.0 guarantees BioAge(t=1.0) = 1.0 by
    construction. No circularity: ODE dynamics are independent of norms.

    Returns:
        Dict with keys meth_norm, damage_norm, H_norm, sen_norm, SASP_norm.
    """
    ctrl = simulate(interventions={'start_time': 0.56}, sex='M', t_max=2.5, dt=dt)
    idx = int(round(1.0 / dt))
    return {
        'meth_norm': float(ctrl.Methylation[idx]),
        'damage_norm': float(ctrl.DNA_damage[idx]),
        'H_norm': float(ctrl.Heteroplasmy[idx]),
        'sen_norm': float(ctrl.SenCells[idx]),
        'SASP_norm': float(ctrl.SASP[idx]),
    }


# =============================================================================
# CORE SIMULATION
# =============================================================================

def simulate(
    interventions: Dict[str, float] = None,
    compound: str = None,
    sex: str = 'M',
    t_max: float = 2.5,
    dt: float = 0.001,
    verbose: bool = False,
    seed: int = 42
) -> SimulationResult:
    """
    Run ELM longevity simulation.

    Args:
        interventions: Dict of pathway -> injection value (0-1 scale)
        compound: Name of compound from COMPOUNDS database (alternative to interventions)
        sex: 'M' or 'F' for sex-specific effects
        t_max: Maximum simulation time (normalized; 1.0 = baseline lifespan)
        dt: Time step
        verbose: Print progress
        seed: Random seed for cancer stochastic events (default 42)

    Returns:
        SimulationResult with trajectories and lifespan extension

    Example:
        >>> result = simulate(compound='rapamycin', sex='M')
        >>> print(f"Extension: {result.extension_percent:.1f}%")
    """
    rng = np.random.default_rng(seed)
    # Build time array
    t_array = np.arange(0, t_max, dt)
    n = len(t_array)

    # Get interventions from compound if specified
    if compound is not None:
        interventions = get_compound(compound)
        if interventions is None:
            raise ValueError(f"Unknown compound: {compound}")
        # Apply sex-specific modifier
        interventions = apply_sex_modifier(compound, interventions, sex)
        interventions['start_time'] = get_itp_start_time(compound)
    elif interventions is None:
        interventions = {}

    # Extract intervention parameters
    start_time = interventions.get('start_time', 0.56)
    ampk_boost = interventions.get('ampk', 0)
    mtorc1_inhibition = interventions.get('mtorc1_drug_inhibition',
                                          interventions.get('mtorc1_inhibition', 0))
    nmn_supplement = interventions.get('nmn', 0)
    cd38_inhibitor = interventions.get('cd38_inhibitor', 0)
    antioxidant = interventions.get('antioxidant', 0)
    sirt1_direct = interventions.get('sirt1_direct', 0)
    akg_supplement = interventions.get('akg', 0)
    gut_microbiome = interventions.get('gut_microbiome', 0)
    antiinflam = interventions.get('antiinflam', 0)
    senolytic = interventions.get('senolytic', 0)
    mitophagy = interventions.get('mitophagy', 0)
    osk_reprogramming = interventions.get('osk', 0)
    nad_target = interventions.get('nad_target', 0)

    # State arrays
    NAD = np.zeros(n)
    SIRT1 = np.zeros(n)
    FOXO3 = np.zeros(n)
    EZH2 = np.zeros(n)
    Autophagy = np.zeros(n)
    DNA_damage = np.zeros(n)
    Methylation = np.zeros(n)
    Heteroplasmy = np.zeros(n)
    SenCells = np.zeros(n)
    SASP = np.zeros(n)
    Mutations = np.zeros(n)
    BioAge = np.zeros(n)
    Cancer_prob = np.zeros(n)
    Immune_function = np.zeros(n)
    UPRmt = np.zeros(n)
    TSC2 = np.zeros(n)
    mTORC1_activity = np.zeros(n)

    # NAD+ targeting tracking
    NAD_target_dose = np.zeros(n)
    NAD_target_achievable = np.ones(n, dtype=bool)

    # Initial conditions
    NAD[0] = 1.0
    SIRT1[0] = calculate_sirt1(NAD[0], 0, SIRTUIN_PARAMS)
    FOXO3[0] = calculate_foxo3(SIRT1[0], SIRTUIN_PARAMS)
    EZH2[0] = 1.0 - 0.5 * SIRT1[0]
    DNA_damage[0] = 0.05
    Methylation[0] = METHYLATION_PARAMS['meth_young']
    Heteroplasmy[0] = HETEROPLASMY_PARAMS['H_0']
    SenCells[0] = 0.01
    SASP[0] = SASP_PARAMS['SASP_young']
    Mutations[0] = 5.0
    Immune_function[0] = 1.0
    UPRmt[0] = 0.05
    TSC2[0] = calculate_tsc2(0, TSC2_PARAMS)  # Basal TSC2 at zero AMPK
    mTORC1_activity[0] = calculate_mtorc1_activity(TSC2[0], 0, TSC2_PARAMS)
    Autophagy[0] = calculate_autophagy(0, FOXO3[0], 0, 0, 0, AUTOPHAGY_PARAMS)

    # Pulse tracking
    last_senolytic_pulse = -1.0
    last_mitophagy_pulse = -1.0
    last_osk_pulse = -1.0
    sen_duration = 0.0

    # Cancer tracking
    cancer_occurred = False
    cancer_time = None

    # Main simulation loop
    for i in range(1, n):
        t = t_array[i]
        age_factor = t
        active = 1.0 if t >= start_time else 0.0

        # Effective AMPK
        ampk_eff = calculate_effective_ampk(
            ampk_direct=ampk_boost * active,
            gut_contribution=gut_microbiome * active,
            params=AMPK_PARAMS
        )

        # TSC2 and mTORC1 (canonical AMPK -> TSC2 -> Rheb -> mTORC1 pathway)
        TSC2[i] = calculate_tsc2(ampk_eff, TSC2_PARAMS)
        drug_inh = mtorc1_inhibition * active
        mTORC1_activity[i] = calculate_mtorc1_activity(TSC2[i], drug_inh, TSC2_PARAMS)
        mtorc1_inh = 1.0 - mTORC1_activity[i]
        mtorc1_effects = calculate_mtorc1_effects(mtorc1_inh, MTORC1_PARAMS)

        # SASP effects
        sasp_effects = get_sasp_effects(SASP[i-1], SASP_PARAMS)

        # OXPHOS efficiency
        oxphos_eff = max(0.30, 1 - NAD_PARAMS['k_hetero_oxphos'] * Heteroplasmy[i-1])

        # UPRmt
        ros_level = 1.0 + 0.5 * Heteroplasmy[i-1]
        nad_ratio = NAD[i-1] / NAD_PARAMS['NAD_young']
        UPRmt[i] = calculate_uprmt(Heteroplasmy[i-1], ros_level, nad_ratio, UPRMT_PARAMS)
        uprmt_effects = get_uprmt_effects(UPRmt[i], UPRMT_PARAMS)

        # NAD+ Homeostasis
        consumption = calculate_nad_consumption(
            nad=NAD[i-1],
            dna_damage=DNA_damage[i-1],
            sirt1=SIRT1[i-1],
            sasp=SASP[i-1],
            age_factor=age_factor,
            cd38_inhibitor=cd38_inhibitor * active,
            gut_antiinflam=gut_microbiome * active,
            params=NAD_PARAMS
        )

        # NAD+ targeting mode (closed-loop control)
        if nad_target > 0 and active > 0:
            base_synthesis = calculate_nad_synthesis(
                ampk=ampk_eff,
                nmn_supplement=0,
                age_factor=age_factor,
                heteroplasmy=Heteroplasmy[i-1],
                params=NAD_PARAMS
            )

            required_dose, achievable = compute_nmn_dose_for_nad_target(
                current_nad=NAD[i-1],
                target_nad=nad_target,
                age_factor=age_factor,
                base_synthesis=base_synthesis,
                consumption=consumption,
                params=NMN_SATURATION_PARAMS
            )

            NAD_target_dose[i] = required_dose
            NAD_target_achievable[i] = achievable

            if achievable:
                NAD[i] = nad_target
            else:
                nmn_boost, _ = compute_nmn_boost_saturated(100.0, age_factor, NMN_SATURATION_PARAMS)
                synthesis = base_synthesis + nmn_boost
                dNAD = synthesis - consumption
                NAD[i] = np.clip(NAD[i-1] + dNAD * dt, NAD_PARAMS['NAD_min'], nad_target)
        else:
            synthesis = calculate_nad_synthesis(
                ampk=ampk_eff,
                nmn_supplement=nmn_supplement * active,
                age_factor=age_factor,
                heteroplasmy=Heteroplasmy[i-1],
                params=NAD_PARAMS
            )
            dNAD = synthesis - consumption
            # No artificial floor needed: substrate-dependent consumption
            # creates a natural equilibrium as NAD → 0
            NAD[i] = min(NAD[i-1] + dNAD * dt, NAD_PARAMS['NAD_max'])
            NAD[i] = max(NAD[i], 0.001)  # safety net only (should never bind)

        # SIRT1 and FOXO3
        SIRT1[i] = calculate_sirt1(NAD[i], sirt1_direct * active, SIRTUIN_PARAMS)
        FOXO3[i] = calculate_foxo3(SIRT1[i], SIRTUIN_PARAMS)
        ampk_ezh2_inhib = METHYLATION_PARAMS['k_ezh2_ampk_inhibition'] * ampk_eff * active
        EZH2[i] = np.clip(1.0 - 0.5 * SIRT1[i] - ampk_ezh2_inhib, 0.1, 1.0)

        # Autophagy
        Autophagy[i] = calculate_autophagy(
            ampk=ampk_eff,
            foxo3=FOXO3[i],
            gut=gut_microbiome * active,
            uprmt_autophagy=uprmt_effects['autophagy_boost'],
            mtorc1_autophagy=mtorc1_effects['autophagy_boost'],
            params=AUTOPHAGY_PARAMS
        )

        # DNA Damage
        osk_pulse_active = False
        if osk_reprogramming > 0 and active > 0:
            if t - last_osk_pulse >= METHYLATION_PARAMS['osk_pulse_interval']:
                last_osk_pulse = t
                osk_pulse_active = True
            elif t - last_osk_pulse < METHYLATION_PARAMS['osk_pulse_duration']:
                osk_pulse_active = True

        proliferation = CANCER_PARAMS['k_osk_proliferation'] if osk_pulse_active else 0.0

        damage_rate = calculate_damage_rate(
            ros_level=ros_level,
            proliferation=proliferation,
            antioxidant=antioxidant * active,
            ampk=ampk_eff * active,
            uprmt_proteostasis=uprmt_effects['proteostasis_boost'],
            mtorc1_damage_reduction=mtorc1_effects['damage_reduction'],
            sasp_ros=sasp_effects['ros_production'],
            params=DNA_REPAIR_PARAMS
        )

        repair_rate = calculate_repair_rate(
            nad=NAD[i],
            sirt1=SIRT1[i],
            dna_damage=DNA_damage[i-1],
            params=DNA_REPAIR_PARAMS
        )

        dDamage = damage_rate - repair_rate
        DNA_damage[i] = np.clip(DNA_damage[i-1] + dDamage * dt, 0, DNA_REPAIR_PARAMS['damage_max'])
        unrepaired = max(0, damage_rate - repair_rate)

        # Methylation
        meth_rate = calculate_methylation_rate(
            ezh2=EZH2[i],
            dnmt=1.0,
            damage=DNA_damage[i],
            sasp=SASP[i-1],
            sasp_meth_accel=sasp_effects['meth_acceleration'],
            gut=gut_microbiome * active,
            params=METHYLATION_PARAMS
        )
        tet_activity = 1.0 + METHYLATION_PARAMS['k_tet_akg_boost'] * akg_supplement * active
        demeth_rate = calculate_demethylation_rate(
            meth=Methylation[i-1],
            tet_activity=tet_activity,
            osk_active=osk_pulse_active,
            params=METHYLATION_PARAMS
        )
        dMeth = meth_rate - demeth_rate
        Methylation[i] = max(METHYLATION_PARAMS['meth_young'] * 0.5, Methylation[i-1] + dMeth * dt)

        # Heteroplasmy
        mitophagy_pulse_active = False
        if mitophagy > 0 and active > 0:
            if t - last_mitophagy_pulse >= HETEROPLASMY_PARAMS['mitophagy_pulse_interval']:
                last_mitophagy_pulse = t
                mitophagy_pulse_active = True
            elif t - last_mitophagy_pulse < HETEROPLASMY_PARAMS['mitophagy_pulse_duration']:
                mitophagy_pulse_active = True

        dH = calculate_heteroplasmy_dynamics(
            H=Heteroplasmy[i-1],
            foxo3=FOXO3[i],
            mitophagy_active=mitophagy_pulse_active,
            mitophagy_strength=mitophagy,
            params=HETEROPLASMY_PARAMS
        )
        Heteroplasmy[i] = np.clip(Heteroplasmy[i-1] + dH * dt, HETEROPLASMY_PARAMS['H_floor'], 1.0)

        # Senescent Cells
        sen_entry = calculate_senescence_rate(
            damage=DNA_damage[i],
            sirt1=SIRT1[i],
            ampk=ampk_eff * active,
            gut=gut_microbiome * active,
            uprmt_sen_reduction=uprmt_effects['sen_reduction'],
            mtorc1_sen_reduction=mtorc1_effects['senescence_reduction'],
            sasp_paracrine=sasp_effects['paracrine_senescence'],
            params=SENESCENCE_PARAMS
        )

        senolytic_pulse_active = False
        senolytic_pulse_just_started = False
        if senolytic > 0 and active > 0:
            if t - last_senolytic_pulse >= SENESCENCE_PARAMS['senolytic_pulse_interval']:
                last_senolytic_pulse = t
                senolytic_pulse_active = True
                senolytic_pulse_just_started = True
            elif t - last_senolytic_pulse < SENESCENCE_PARAMS['senolytic_pulse_duration']:
                senolytic_pulse_active = True

        if senolytic_pulse_just_started:
            kill_fraction = SENESCENCE_PARAMS['k_senolytic_kill_fraction'] * senolytic
            SenCells[i-1] *= (1 - kill_fraction)
            SenCells[i-1] = max(SenCells[i-1], SENESCENCE_PARAMS['sen_min_residual'])

        sen_clear = calculate_senescence_clearance(
            sen=SenCells[i-1],
            autophagy=Autophagy[i],
            senolytic_active=senolytic_pulse_active,
            params=SENESCENCE_PARAMS
        )

        dSen = sen_entry - sen_clear
        SenCells[i] = np.clip(SenCells[i-1] + dSen * dt, SENESCENCE_PARAMS['sen_min_residual'], 1.0)

        if SenCells[i] > SenCells[i-1]:
            sen_duration += dt
        else:
            sen_duration = max(0, sen_duration - dt * 0.5)

        # SASP dynamics
        sasp_prod = calculate_sasp_production(
            sen_cells=SenCells[i],
            sen_duration=sen_duration,
            params=SASP_PARAMS
        )
        liver_function = max(0.5, 1.0 - 0.3 * age_factor)
        sasp_clear = calculate_sasp_clearance(
            sasp=SASP[i-1],
            liver_function=liver_function,
            antiinflam=antiinflam * active,
            params=SASP_PARAMS
        )
        if senolytic_pulse_active:
            sasp_clear += SASP_PARAMS['k_senolytic_sasp_reduction'] * SASP[i-1]

        dSASP = sasp_prod - sasp_clear
        SASP[i] = np.clip(SASP[i-1] + dSASP * dt, SASP_PARAMS['SASP_young'] * 0.5, 2.0)

        # Immune Function
        Immune_function[i] = max(0.2, 1.0 - CANCER_PARAMS['k_immune_age_decline'] * age_factor - sasp_effects['immune_suppression'])

        # Mutations
        mut_rate = calculate_mutation_rate(
            unrepaired_damage=unrepaired,
            proliferation=proliferation,
            ros=ros_level + sasp_effects['ros_production'],
            params=CANCER_PARAMS
        )
        p53_activity = 0.5 + 0.5 * SIRT1[i]
        mut_clear = calculate_mutation_clearance(
            mutations=Mutations[i-1],
            immune_function=Immune_function[i],
            p53_activity=p53_activity,
            senescent_cells=SenCells[i],
            params=CANCER_PARAMS
        )
        dMut = mut_rate - mut_clear
        Mutations[i] = max(0, Mutations[i-1] + dMut * dt)
        Cancer_prob[i] = calculate_cancer_probability(Mutations[i], CANCER_PARAMS)

        if not cancer_occurred and rng.random() < Cancer_prob[i] * dt:
            cancer_occurred = True
            cancer_time = t

        # Biological Age
        BioAge[i] = calculate_bioage(
            meth=Methylation[i],
            damage=DNA_damage[i],
            H=Heteroplasmy[i],
            sen=SenCells[i],
            sasp=SASP[i],
            params=BIOAGE_PARAMS
        )

    # Find death time
    death_idx = np.argmax(BioAge >= BIOAGE_PARAMS['death_threshold'])
    if death_idx == 0 and BioAge[0] < BIOAGE_PARAMS['death_threshold']:
        death_idx = n - 1
    t_death = t_array[death_idx]

    # Survival curve
    sigma = 0.12
    Survival = 1.0 / (1.0 + np.exp((t_array - t_death) / sigma))

    # Extension will be calculated by caller (avoids recursion)
    extension_percent = 0.0

    return SimulationResult(
        t=t_array,
        t_death=t_death,
        extension_percent=extension_percent,
        NAD=NAD,
        BioAge=BioAge,
        Survival=Survival,
        SIRT1=SIRT1,
        FOXO3=FOXO3,
        Autophagy=Autophagy,
        DNA_damage=DNA_damage,
        Methylation=Methylation,
        Heteroplasmy=Heteroplasmy,
        SenCells=SenCells,
        SASP=SASP,
        TSC2=TSC2,
        mTORC1_activity=mTORC1_activity,
        Mutations=Mutations,
        Cancer_prob=Cancer_prob,
        cancer_occurred=cancer_occurred,
        cancer_time=cancer_time,
        interventions=interventions or {},
    )


def calculate_lifespan_extension(treated: SimulationResult, control: SimulationResult) -> float:
    """Calculate percent lifespan extension."""
    if control.t_death == 0:
        return 0.0
    return (treated.t_death - control.t_death) / control.t_death * 100


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def run_control(sex: str = 'M', t_max: float = 2.5, dt: float = 0.001) -> SimulationResult:
    """Run untreated control simulation."""
    return simulate(interventions={'start_time': 0.56}, sex=sex, t_max=t_max, dt=dt)


def run_itp_compound(
    compound: str,
    sex: str = 'M',
    t_max: float = 2.5,
    dt: float = 0.001
) -> Tuple[SimulationResult, float]:
    """
    Run simulation for an ITP compound and return extension.

    Returns (result, extension_percent).
    """
    control = run_control(sex=sex, t_max=t_max, dt=dt)
    treated = simulate(compound=compound, sex=sex, t_max=t_max, dt=dt)
    extension = calculate_lifespan_extension(treated, control)
    # Update the result with calculated extension
    treated.extension_percent = extension
    return treated, extension


def run_combination(
    compound_doses: Dict[str, float],
    sex: str = 'M',
    t_max: float = 2.5,
    dt: float = 0.001,
    start_time: Optional[float] = None
) -> SimulationResult:
    """
    Run simulation for a combination of compounds at specified dose levels.

    Each compound's pathway injections are scaled by its dose multiplier,
    sex-modified independently, then merged by summing shared pathway keys.

    Args:
        compound_doses: Dict of compound_name -> dose_multiplier
            e.g. {'rapamycin': 1.0, 'acarbose': 0.5}
        sex: 'M' or 'F'
        t_max: Maximum simulation time
        dt: Time step
        start_time: Normalized treatment start time. If provided, overrides
            per-compound ITP start times (use for combo studies where all
            drugs begin at the same age). If None, uses the earliest
            individual ITP start time across compounds.

    Returns:
        SimulationResult with extension_percent set vs control.
    """
    merged = {}
    earliest_start = 1.0  # Track earliest start time across compounds
    for compound, dose in compound_doses.items():
        if dose <= 0:
            continue
        base = get_compound(compound)
        if base is None:
            raise ValueError(f"Unknown compound: {compound}")
        # Track earliest ITP start time (used as fallback)
        start = get_itp_start_time(compound)
        if start < earliest_start:
            earliest_start = start
        # Apply dose-response curve (diminishing returns at high doses)
        from .dose_response import effective_dose
        eff_dose = effective_dose(compound, dose)
        scaled = {k: v * eff_dose for k, v in base.items()}
        # Apply sex modifier per compound (sex acts on each drug independently)
        scaled = apply_sex_modifier(compound, scaled, sex)
        # Sum into merged dict
        for k, v in scaled.items():
            merged[k] = merged.get(k, 0) + v

    merged['start_time'] = start_time if start_time is not None else earliest_start

    control = run_control(sex=sex, t_max=t_max, dt=dt)
    treated = simulate(interventions=merged, sex=sex, t_max=t_max, dt=dt)
    extension = calculate_lifespan_extension(treated, control)
    treated.extension_percent = extension
    return treated


def run_combination_extension(
    compound_doses: Dict[str, float],
    sex: str = 'M',
    t_max: float = 2.5,
    dt: float = 0.001,
    start_time: Optional[float] = None
) -> float:
    """Convenience wrapper: returns just the extension_percent float."""
    result = run_combination(compound_doses, sex=sex, t_max=t_max, dt=dt,
                             start_time=start_time)
    return result.extension_percent
