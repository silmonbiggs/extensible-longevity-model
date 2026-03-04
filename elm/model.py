"""
ELM Core Simulation Model

This module contains the core ODE-based simulation for the ELM longevity model.
The model tracks 8 ODE state variables as interconnected subsystems:

1. Heteroplasmy - Mitochondrial DNA damage fraction
2. NAD+ - Nicotinamide adenine dinucleotide level
3. DNA damage - Genomic lesion accumulation
4. AMPK/mTORC1 signaling - Nutrient-sensing axis
5. Autophagy - Cellular recycling flux
6. Senescent cell burden
7. SASP - Senescence-associated secretory phenotype
8. Methylation drift - Epigenetic age

Additional algebraic variables (SIRT1, FOXO3, EZH2, UPRmt, immune
function) are computed from the state at each time step. The cancer
module (mutation accumulation) is present but INACTIVE in all current
predictions -- it does not influence BioAge or lifespan.

Simulation outputs biological age trajectory and lifespan.

Integration uses scipy.integrate.solve_ivp (RK45, adaptive stepping).

Sensitivity analysis: see scripts/generate_figures.py (tornado charts,
weight sweeps, pairwise synergy matrix) and tests/test_pathways.py
(numerical convergence, null model comparison).
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, Any
from scipy.integrate import solve_ivp

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

    # -----------------------------------------------------------------
    # State vector: 8 ODE variables
    #   [0] NAD, [1] DNA_damage, [2] Methylation, [3] Heteroplasmy,
    #   [4] SenCells, [5] SASP, [6] Mutations, [7] sen_duration
    # -----------------------------------------------------------------
    I_NAD, I_DMG, I_METH, I_H, I_SEN, I_SASP, I_MUT, I_SENDUR = range(8)

    y0 = np.array([
        1.0,                                # NAD
        0.05,                               # DNA_damage
        METHYLATION_PARAMS['meth_young'],   # Methylation
        HETEROPLASMY_PARAMS['H_0'],         # Heteroplasmy
        0.01,                               # SenCells
        SASP_PARAMS['SASP_young'],          # SASP
        5.0,                                # Mutations
        0.0,                                # sen_duration
    ])

    # Pulse schedule helpers
    def _is_pulse_active(t, interval, duration):
        if t < start_time:
            return False
        phase = (t - start_time) % interval
        return phase < duration

    # Convert senolytic instantaneous kill to continuous fast clearance
    # rate over the pulse duration: -ln(1-f)/d gives equivalent total kill.
    sen_pulse_interval = SENESCENCE_PARAMS.get('senolytic_pulse_interval', 0.1)
    sen_pulse_duration = SENESCENCE_PARAMS.get('senolytic_pulse_duration', 0.01)
    if senolytic > 0 and sen_pulse_duration > 0:
        kill_frac = SENESCENCE_PARAMS['k_senolytic_kill_fraction'] * senolytic
        senolytic_fast_rate = -np.log(max(1e-12, 1 - kill_frac)) / sen_pulse_duration
    else:
        senolytic_fast_rate = 0.0

    mito_interval = HETEROPLASMY_PARAMS.get('mitophagy_pulse_interval', 0.1)
    mito_duration = HETEROPLASMY_PARAMS.get('mitophagy_pulse_duration', 0.03)
    osk_interval = METHYLATION_PARAMS.get('osk_pulse_interval', 0.1)
    osk_duration = METHYLATION_PARAMS.get('osk_pulse_duration', 0.01)

    def rhs(t, y):
        """ODE right-hand side for the 8-variable state vector."""
        nad     = max(y[I_NAD], 0.001)
        damage  = max(y[I_DMG], 0.0)
        meth    = max(y[I_METH], METHYLATION_PARAMS['meth_young'] * 0.5)
        hetero  = np.clip(y[I_H], HETEROPLASMY_PARAMS['H_floor'], 1.0)
        sen     = max(y[I_SEN], SENESCENCE_PARAMS['sen_min_residual'])
        sasp    = max(y[I_SASP], SASP_PARAMS['SASP_young'] * 0.5)
        mut     = max(y[I_MUT], 0.0)
        sendur  = max(y[I_SENDUR], 0.0)

        age_factor = t
        active = 1.0 if t >= start_time else 0.0

        # AMPK
        ampk_eff = calculate_effective_ampk(
            ampk_direct=ampk_boost * active,
            gut_contribution=gut_microbiome * active,
            params=AMPK_PARAMS
        )

        # TSC2 / mTORC1
        tsc2 = calculate_tsc2(ampk_eff, TSC2_PARAMS)
        drug_inh = mtorc1_inhibition * active
        mtorc1_act = calculate_mtorc1_activity(tsc2, drug_inh, TSC2_PARAMS)
        mtorc1_effects = calculate_mtorc1_effects(1.0 - mtorc1_act, MTORC1_PARAMS)

        # SASP effects
        sasp_effects = get_sasp_effects(sasp, SASP_PARAMS)

        # ROS and UPRmt
        ros_level = 1.0 + 0.5 * hetero
        nad_ratio = nad / NAD_PARAMS['NAD_young']
        uprmt = calculate_uprmt(hetero, ros_level, nad_ratio, UPRMT_PARAMS)
        uprmt_effects = get_uprmt_effects(uprmt, UPRMT_PARAMS)

        # SIRT1 / FOXO3 / EZH2 (algebraic)
        sirt1 = calculate_sirt1(nad, sirt1_direct * active, SIRTUIN_PARAMS)
        foxo3 = calculate_foxo3(sirt1, SIRTUIN_PARAMS)
        ampk_ezh2_inhib = METHYLATION_PARAMS['k_ezh2_ampk_inhibition'] * ampk_eff * active
        ezh2 = np.clip(1.0 - 0.5 * sirt1 - ampk_ezh2_inhib, 0.1, 1.0)

        # Autophagy (algebraic)
        autophagy = calculate_autophagy(
            ampk=ampk_eff, foxo3=foxo3,
            gut=gut_microbiome * active,
            uprmt_autophagy=uprmt_effects['autophagy_boost'],
            mtorc1_autophagy=mtorc1_effects['autophagy_boost'],
            params=AUTOPHAGY_PARAMS
        )

        # --- dNAD ---
        consumption = calculate_nad_consumption(
            nad=nad, dna_damage=damage, sirt1=sirt1, sasp=sasp,
            age_factor=age_factor,
            cd38_inhibitor=cd38_inhibitor * active,
            gut_antiinflam=gut_microbiome * active,
            params=NAD_PARAMS
        )
        if nad_target > 0 and active > 0:
            base_syn = calculate_nad_synthesis(
                ampk=ampk_eff, nmn_supplement=0,
                age_factor=age_factor, heteroplasmy=hetero, params=NAD_PARAMS
            )
            _, achievable = compute_nmn_dose_for_nad_target(
                current_nad=nad, target_nad=nad_target,
                age_factor=age_factor, base_synthesis=base_syn,
                consumption=consumption, params=NMN_SATURATION_PARAMS
            )
            if achievable:
                dNAD = (nad_target - nad) * 100.0  # fast pull toward target
            else:
                nmn_boost, _ = compute_nmn_boost_saturated(100.0, age_factor, NMN_SATURATION_PARAMS)
                dNAD = base_syn + nmn_boost - consumption
        else:
            synthesis = calculate_nad_synthesis(
                ampk=ampk_eff, nmn_supplement=nmn_supplement * active,
                age_factor=age_factor, heteroplasmy=hetero, params=NAD_PARAMS
            )
            dNAD = synthesis - consumption
        # Soft bounds
        if nad >= NAD_PARAMS['NAD_max'] and dNAD > 0:
            dNAD = 0.0
        if nad <= 0.001 and dNAD < 0:
            dNAD = 0.0

        # --- dDamage ---
        osk_pulse = _is_pulse_active(t, osk_interval, osk_duration) if (osk_reprogramming > 0 and active > 0) else False
        proliferation = CANCER_PARAMS['k_osk_proliferation'] if osk_pulse else 0.0
        damage_rate = calculate_damage_rate(
            ros_level=ros_level, proliferation=proliferation,
            antioxidant=antioxidant * active, ampk=ampk_eff * active,
            uprmt_proteostasis=uprmt_effects['proteostasis_boost'],
            mtorc1_damage_reduction=mtorc1_effects['damage_reduction'],
            sasp_ros=sasp_effects['ros_production'],
            params=DNA_REPAIR_PARAMS
        )
        repair_rate = calculate_repair_rate(
            nad=nad, sirt1=sirt1, dna_damage=damage, params=DNA_REPAIR_PARAMS
        )
        dDamage = damage_rate - repair_rate
        if damage >= DNA_REPAIR_PARAMS['damage_max'] and dDamage > 0:
            dDamage = 0.0
        if damage <= 0 and dDamage < 0:
            dDamage = 0.0

        # --- dMeth ---
        meth_rate = calculate_methylation_rate(
            ezh2=ezh2, dnmt=1.0, damage=damage, sasp=sasp,
            sasp_meth_accel=sasp_effects['meth_acceleration'],
            gut=gut_microbiome * active, params=METHYLATION_PARAMS
        )
        tet_act = 1.0 + METHYLATION_PARAMS['k_tet_akg_boost'] * akg_supplement * active
        demeth_rate = calculate_demethylation_rate(
            meth=meth, tet_activity=tet_act, osk_active=osk_pulse,
            params=METHYLATION_PARAMS
        )
        dMeth = meth_rate - demeth_rate
        if meth <= METHYLATION_PARAMS['meth_young'] * 0.5 and dMeth < 0:
            dMeth = 0.0

        # --- dH ---
        mito_pulse = _is_pulse_active(t, mito_interval, mito_duration) if (mitophagy > 0 and active > 0) else False
        dH = calculate_heteroplasmy_dynamics(
            H=hetero, foxo3=foxo3,
            mitophagy_active=mito_pulse, mitophagy_strength=mitophagy,
            params=HETEROPLASMY_PARAMS
        )
        if hetero >= 1.0 and dH > 0:
            dH = 0.0
        if hetero <= HETEROPLASMY_PARAMS['H_floor'] and dH < 0:
            dH = 0.0

        # --- dSen ---
        sen_entry = calculate_senescence_rate(
            damage=damage, sirt1=sirt1, ampk=ampk_eff * active,
            gut=gut_microbiome * active,
            uprmt_sen_reduction=uprmt_effects['sen_reduction'],
            mtorc1_sen_reduction=mtorc1_effects['senescence_reduction'],
            sasp_paracrine=sasp_effects['paracrine_senescence'],
            params=SENESCENCE_PARAMS
        )
        sen_pulse = _is_pulse_active(t, sen_pulse_interval, sen_pulse_duration) if (senolytic > 0 and active > 0) else False
        sen_clear = calculate_senescence_clearance(
            sen=sen, autophagy=autophagy, senolytic_active=sen_pulse,
            params=SENESCENCE_PARAMS
        )
        dSen = sen_entry - sen_clear
        # Senolytic fast kill (continuous equivalent of instantaneous fraction kill)
        if sen_pulse and senolytic_fast_rate > 0:
            dSen -= senolytic_fast_rate * sen
        if sen <= SENESCENCE_PARAMS['sen_min_residual'] and dSen < 0:
            dSen = 0.0

        # --- dSASP ---
        sasp_prod = calculate_sasp_production(
            sen_cells=sen, sen_duration=sendur, params=SASP_PARAMS
        )
        liver_fn = max(0.5, 1.0 - 0.3 * age_factor)
        sasp_clr = calculate_sasp_clearance(
            sasp=sasp, liver_function=liver_fn,
            antiinflam=antiinflam * active, params=SASP_PARAMS
        )
        if sen_pulse:
            sasp_clr += SASP_PARAMS['k_senolytic_sasp_reduction'] * sasp
        dSASP = sasp_prod - sasp_clr
        if sasp >= 2.0 and dSASP > 0:
            dSASP = 0.0
        if sasp <= SASP_PARAMS['SASP_young'] * 0.5 and dSASP < 0:
            dSASP = 0.0

        # --- dMut ---
        unrepaired = max(0, damage_rate - repair_rate)
        immune_fn = max(0.2, 1.0 - CANCER_PARAMS['k_immune_age_decline'] * age_factor
                        - sasp_effects['immune_suppression'])
        mut_rate = calculate_mutation_rate(
            unrepaired_damage=unrepaired, proliferation=proliferation,
            ros=ros_level + sasp_effects['ros_production'],
            params=CANCER_PARAMS
        )
        p53_act = 0.5 + 0.5 * sirt1
        mut_clear = calculate_mutation_clearance(
            mutations=mut, immune_function=immune_fn,
            p53_activity=p53_act, senescent_cells=sen,
            params=CANCER_PARAMS
        )
        dMut = mut_rate - mut_clear
        if mut <= 0 and dMut < 0:
            dMut = 0.0

        # --- dSenDur (tracks how long SenCells has been rising) ---
        dSenDur = 1.0 if dSen > 0 else -0.5
        if sendur <= 0 and dSenDur < 0:
            dSenDur = 0.0

        return [dNAD, dDamage, dMeth, dH, dSen, dSASP, dMut, dSenDur]

    # -----------------------------------------------------------------
    # Integrate with scipy RK45
    # -----------------------------------------------------------------
    sol = solve_ivp(
        rhs, [t_array[0], t_array[-1]], y0,
        method='RK45', t_eval=t_array,
        rtol=1e-8, atol=1e-10, max_step=0.01
    )

    if sol.status != 0:
        raise RuntimeError(f"ODE integration failed: {sol.message}")

    # Extract state trajectories
    NAD          = sol.y[I_NAD]
    DNA_damage   = sol.y[I_DMG]
    Methylation  = sol.y[I_METH]
    Heteroplasmy = sol.y[I_H]
    SenCells     = sol.y[I_SEN]
    SASP_arr     = sol.y[I_SASP]
    Mutations    = sol.y[I_MUT]

    # -----------------------------------------------------------------
    # Post-process: compute algebraic quantities on output grid
    # -----------------------------------------------------------------
    SIRT1          = np.zeros(n)
    FOXO3          = np.zeros(n)
    Autophagy      = np.zeros(n)
    TSC2_arr       = np.zeros(n)
    mTORC1_act_arr = np.zeros(n)
    BioAge         = np.zeros(n)
    Cancer_prob    = np.zeros(n)

    cancer_occurred = False
    cancer_time = None

    for i in range(n):
        t = t_array[i]
        active = 1.0 if t >= start_time else 0.0
        ampk_eff = calculate_effective_ampk(
            ampk_direct=ampk_boost * active,
            gut_contribution=gut_microbiome * active,
            params=AMPK_PARAMS
        )
        TSC2_arr[i] = calculate_tsc2(ampk_eff, TSC2_PARAMS)
        mTORC1_act_arr[i] = calculate_mtorc1_activity(
            TSC2_arr[i], mtorc1_inhibition * active, TSC2_PARAMS
        )
        SIRT1[i] = calculate_sirt1(NAD[i], sirt1_direct * active, SIRTUIN_PARAMS)
        FOXO3[i] = calculate_foxo3(SIRT1[i], SIRTUIN_PARAMS)

        mtorc1_eff = calculate_mtorc1_effects(1.0 - mTORC1_act_arr[i], MTORC1_PARAMS)
        uprmt_val = calculate_uprmt(
            Heteroplasmy[i], 1.0 + 0.5 * Heteroplasmy[i],
            NAD[i] / NAD_PARAMS['NAD_young'], UPRMT_PARAMS
        )
        uprmt_eff = get_uprmt_effects(uprmt_val, UPRMT_PARAMS)
        Autophagy[i] = calculate_autophagy(
            ampk=ampk_eff, foxo3=FOXO3[i],
            gut=gut_microbiome * active,
            uprmt_autophagy=uprmt_eff['autophagy_boost'],
            mtorc1_autophagy=mtorc1_eff['autophagy_boost'],
            params=AUTOPHAGY_PARAMS
        )

        BioAge[i] = calculate_bioage(
            meth=Methylation[i], damage=DNA_damage[i],
            H=Heteroplasmy[i], sen=SenCells[i], sasp=SASP_arr[i],
            params=BIOAGE_PARAMS
        )
        Cancer_prob[i] = calculate_cancer_probability(Mutations[i], CANCER_PARAMS)

        if not cancer_occurred and rng.random() < Cancer_prob[i] * dt:
            cancer_occurred = True
            cancer_time = t

    # Find death time
    death_idx = np.argmax(BioAge >= BIOAGE_PARAMS['death_threshold'])
    if death_idx == 0 and BioAge[0] < BIOAGE_PARAMS['death_threshold']:
        death_idx = n - 1
    t_death = t_array[death_idx]

    # Survival curve
    sigma = 0.12
    Survival = 1.0 / (1.0 + np.exp((t_array - t_death) / sigma))

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
        SASP=SASP_arr,
        TSC2=TSC2_arr,
        mTORC1_activity=mTORC1_act_arr,
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
