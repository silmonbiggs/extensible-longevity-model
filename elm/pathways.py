"""
ELM Pathway Parameters

This module contains all parameter dictionaries for the ELM model pathways.
Parameters are organized by biological system and include literature references
where available.

Parameter Provenance:
- ITP: Calibrated to match ITP intervention data
- Literature: Derived from published experimental data
- Estimated: Model estimates with physiological constraints

All parameters are dimensionless ratios or rate constants normalized to
a baseline mouse lifespan of ~2.5 years.
"""

from typing import Dict

# =============================================================================
# NAD+ HOMEOSTASIS
# =============================================================================

NAD_PARAMS: Dict[str, float] = {
    # Synthesis
    'k_nampt_base': 0.55,          # Basal NAMPT activity
    'k_nampt_ampk': 0.50,          # AMPK boost to NAMPT
    'k_nampt_age_decline': 0.431,  # Age-related NAMPT decline (exponential: exp(-0.431)=0.65 at t=1)
    'k_tryptophan': 0.05,          # De novo synthesis contribution

    # Consumption
    'k_consumption_basal': 0.40,   # Baseline NAD+ turnover
    'k_cd38_base': 0.35,           # Basal CD38 consumption
    'k_cd38_age': 0.50,            # Age-related CD38 increase
    'k_cd38_sasp': 0.60,           # SASP-induced CD38
    'k_parp_consumption': 0.70,    # PARP consumption (damage-dependent)
    'k_sirt_consumption': 0.10,    # Sirtuin consumption

    # Interventions
    'k_nmn_synthesis': 0.80,       # NMN supplementation effect
    'k_cd38_inhibitor': 0.70,      # CD38 inhibitor effect
    'k_gut_antiinflam': 0.75,      # Gut microbiome anti-inflammatory

    # Mitochondrial coupling
    'k_hetero_oxphos': 0.35,       # Heteroplasmy impact on OXPHOS
    'k_oxphos_nad_regen': 0.15,    # OXPHOS contribution to NAD+ regen

    # Bounds
    'NAD_min': 0.15,               # Minimum viable NAD+ level
    'NAD_max': 1.2,                # Maximum NAD+ level
    'NAD_young': 1.0,              # Youthful NAD+ reference (normalized)
    'Km_nad_consumption': 0.05,    # Substrate Km for NAD+ consumption (Michaelis-Menten)
}

# NMN dose-response (Michaelis-Menten saturation)
NMN_SATURATION_PARAMS: Dict[str, float] = {
    'Vmax_nmn': 0.80,              # Maximum NMN boost
    'Km_nmn': 0.60,                # Half-max dose
    'k_nmn_absorption_decline': 0.20,  # Age-related absorption decline
}

# CD38 inhibitor dose-response (Hill function with ceiling)
CD38_SATURATION_PARAMS: Dict[str, float] = {
    'Max_cd38_inhibition': 0.75,   # Maximum achievable inhibition
    'IC50_cd38_inhibitor': 0.30,   # Half-max dose
    'Hill_n_cd38': 1.5,            # Hill coefficient
}

# =============================================================================
# AMPK / ENERGY SENSING
# =============================================================================

AMPK_PARAMS: Dict[str, float] = {
    'Km_ampk': 0.40,               # Saturation constant
    'Vmax_ampk': 1.0,              # Maximum AMPK activity
    'k_gut_ampk': 0.20,            # Gut microbiome contribution
}

# =============================================================================
# mTORC1 / GROWTH SIGNALING
# =============================================================================

MTORC1_PARAMS: Dict[str, float] = {
    'k_mtorc1_autophagy': 0.70,    # mTORC1 inhibition -> autophagy boost
    'k_mtorc1_proteostasis': 0.45, # mTORC1 inhibition -> proteostasis
    'k_mtorc1_tfeb': 0.35,         # mTORC1 inhibition -> TFEB activation
    'k_mtorc1_senescence': 0.35,   # mTORC1 inhibition -> senescence reduction
    'k_mtorc1_damage_reduction': 0.25,  # mTORC1 inhibition -> damage reduction
}

# =============================================================================
# SIRTUINS
# =============================================================================

SIRTUIN_PARAMS: Dict[str, float] = {
    'Km_nad_sirt1': 0.4,           # NAD+ affinity for SIRT1
    'n_hill_sirt1': 2.0,           # Hill coefficient (cooperativity)
    'k_sirt1_direct_boost': 0.5,   # Direct SIRT1 activator effect
    'foxo3_basal': 0.2,            # Basal FOXO3 activity
    'k_sirt1_foxo3': 0.8,          # SIRT1 -> FOXO3 activation
}

# =============================================================================
# AUTOPHAGY
# =============================================================================

AUTOPHAGY_PARAMS: Dict[str, float] = {
    'Vmax_autophagy': 1.0,         # Maximum autophagy rate
    'Km_autophagy': 0.35,          # Half-max activation
    'w_ampk_autophagy': 0.45,      # AMPK contribution weight
    'w_foxo3_autophagy': 0.35,     # FOXO3 contribution weight
    'w_gut_autophagy': 0.40,       # Gut microbiome contribution
    'k_autophagy_damage_clear': 0.30,  # Autophagy damage clearance
}

# =============================================================================
# DNA DAMAGE & REPAIR
# =============================================================================

DNA_REPAIR_PARAMS: Dict[str, float] = {
    # Damage sources
    'k_damage_ros': 0.60,          # ROS-induced damage
    'k_damage_replication': 0.12,  # Replication errors
    'k_damage_spontaneous': 0.18,  # Spontaneous damage
    'k_damage_sasp': 0.25,         # SASP-induced damage

    # Repair pathways
    'k_ber_base': 0.35,            # Base excision repair
    'k_ber_nad_dependence': 0.7,   # BER NAD+ dependence
    'k_ner_base': 0.18,            # Nucleotide excision repair
    'k_ner_nad_dependence': 0.4,   # NER NAD+ dependence
    'k_hr_base': 0.10,             # Homologous recombination
    'k_hr_sirt1_boost': 0.5,       # HR SIRT1 dependence

    # Downstream effects
    'k_damage_to_mutation': 0.02,  # Damage -> mutation conversion
    'k_damage_to_senescence': 0.15,  # Damage -> senescence
    'damage_senescence_threshold': 0.3,  # Midpoint for senescence sigmoid
    'damage_senescence_width': 0.05,     # Sigmoid transition width (smooth over D±2w)

    # Interventions
    'k_antioxidant_ros_reduction': 0.55,  # Antioxidant effect
    'k_ampk_damage_reduction': 0.50,      # AMPK protective effect

    # Bounds
    'damage_max': 2.0,             # Maximum damage level
}

# =============================================================================
# EPIGENETICS / METHYLATION
# =============================================================================

METHYLATION_PARAMS: Dict[str, float] = {
    # Writers
    'k_ezh2_base': 1.85,           # Basal EZH2 activity
    'k_ezh2_sirt1_inhibition': 0.5,  # SIRT1 inhibition of EZH2
    'k_ezh2_ampk_inhibition': 0.70,  # AMPK inhibition of EZH2
    'k_dnmt_base': 0.28,           # Basal DNMT activity
    'k_sasp_methylation': 0.42,    # SASP-induced methylation drift

    # Erasers
    'k_tet_base': 0.02,            # Basal TET activity
    'k_tet_akg_boost': 0.30,       # Alpha-ketoglutarate boost to TET
    'k_gut_dnmt_inhibition': 0.45, # Gut metabolite DNMT inhibition

    # OSK reprogramming
    'k_osk_demethylation': 3.5,    # OSK demethylation strength
    'osk_pulse_duration': 0.000252,  # OSK pulse duration (normalized)
    'osk_pulse_interval': 0.00329,   # OSK pulse interval
    'k_osk_efficiency': 0.60,      # OSK reprogramming efficiency

    # Reference values
    'meth_young': 0.25,            # Youthful methylation level
}

# =============================================================================
# MITOCHONDRIA / HETEROPLASMY
# =============================================================================

HETEROPLASMY_PARAMS: Dict[str, float] = {
    # Initial conditions
    'H_0': 0.05,                   # Initial heteroplasmy
    'H_floor': 0.04,               # Minimum heteroplasmy

    # Dynamics (Kowald-Kirkwood model)
    'Neff': 100,                   # Effective population size
    'transcription_advantage': 0.20,  # Mutant mtDNA advantage
    'k_selection_base': 0.001,     # Baseline selection
    'k_selection_foxo3': 0.10,     # FOXO3 selection boost
    'k_selection_mitophagy': 0.60, # Mitophagy selection boost
    'time_scale': 12.0,            # Time scaling factor

    # Mitophagy intervention
    'mitophagy_pulse_interval': 0.00329,
    'mitophagy_pulse_duration': 0.000252,

}

# =============================================================================
# UPRmt (MITOCHONDRIAL UNFOLDED PROTEIN RESPONSE)
# =============================================================================

UPRMT_PARAMS: Dict[str, float] = {
    # Activation
    'k_uprmt_hetero': 0.40,        # Heteroplasmy activation
    'k_uprmt_ros': 0.30,           # ROS activation
    'k_uprmt_nad': 0.25,           # NAD+ modulation
    'Km_uprmt': 0.35,              # Half-max activation
    'Vmax_uprmt': 1.0,             # Maximum UPRmt

    # Effects
    'k_uprmt_proteostasis': 0.25,  # Proteostasis boost
    'k_uprmt_sen_reduction': 0.15, # Senescence reduction
    'k_uprmt_autophagy': 0.20,     # Autophagy boost
}

# =============================================================================
# SENESCENCE
# =============================================================================

SENESCENCE_PARAMS: Dict[str, float] = {
    # Entry
    'k_sen_from_damage': 0.40,     # Damage-induced senescence
    'k_sen_from_telomere': 0.12,   # Replicative senescence
    'k_sen_from_oncogene': 0.08,   # Oncogene-induced senescence
    'k_sen_from_sasp': 0.20,       # SASP-induced (paracrine) senescence

    # Clearance
    'k_sen_natural_clear': 0.02,   # Natural immune clearance
    'k_sen_autophagy_clear': 0.70, # Autophagy-mediated clearance
    'k_senolytic_clear': 0.55,     # Senolytic drug clearance

    # Senolytic intervention (biweekly pulsing — Baker 2016, PMID 26840489)
    'senolytic_pulse_interval': 0.0154,   # ~14 days normalized (biweekly)
    'senolytic_pulse_duration': 0.000072,
    'k_senolytic_kill_fraction': 0.30,    # Fraction killed per pulse (Zhu 2015, PMID 25754370)

    # Dynamics
    'sen_regrowth_rate': 0.10,     # Regrowth after clearance
    'sen_min_residual': 0.015,     # Minimum residual burden

    # Modifiers
    'k_ampk_sen_reduction': 0.45,  # AMPK protective effect
    'k_gut_sen_reduction': 0.40,   # Gut microbiome effect

}

# =============================================================================
# SASP (SENESCENCE-ASSOCIATED SECRETORY PHENOTYPE)
# =============================================================================

SASP_PARAMS: Dict[str, float] = {
    # Production
    'k_sasp_production': 0.45,     # SASP production rate
    'k_sasp_delay': 0.02,          # Delay in SASP onset

    # Clearance
    'k_sasp_decay': 0.60,          # SASP decay rate
    'k_sasp_liver_clear': 0.35,    # Liver-mediated clearance

    # Effects

    # Interventions
    'k_antiinflam_sasp_reduction': 0.40,  # Anti-inflammatory effect
    'k_senolytic_sasp_reduction': 0.50,   # Senolytic SASP reduction

    # Reference values
    'SASP_young': 0.005,           # Youthful SASP level (equilibrium with SenCells[0]=0.01)
    'SASP_old': 0.40,              # Aged SASP level
    # Effects (sender couplings — quantitatively negligible per perturbation
    # analysis, but structurally required by get_sasp_effects() in model.py)
    'k_sasp_cd38': 0.30,           # SASP induction of CD38
    'k_sasp_paracrine_sen': 0.12,  # Paracrine senescence induction
    'k_sasp_ros': 0.18,            # SASP-induced ROS
    'k_sasp_immune_suppress': 0.25,  # Immune suppression
    'k_sasp_meth_drift': 0.15,     # Methylation drift acceleration
    'k_sasp_toxicity': 0.20,       # Aggregate toxicity (infrastructure for future collapse)
}

# =============================================================================
# CANCER
# =============================================================================

CANCER_PARAMS: Dict[str, float] = {
    # Mutation accumulation
    'k_mut_from_damage': 0.015,    # Damage-to-mutation rate
    'k_mut_from_replication': 0.005,  # Replication errors
    'k_mut_from_ros': 0.010,       # ROS-induced mutations
    'k_mut_proliferation': 0.020,  # Proliferation-associated
    'k_osk_proliferation': 0.30,   # OSK-induced proliferation risk

    # Defense
    'k_immune_surveillance': 0.10, # Immune surveillance rate
    'k_immune_age_decline': 0.20,  # Age-related immune decline
    'k_p53_apoptosis': 0.15,       # p53-mediated apoptosis
    'k_senescence_tumor_suppression': 0.20,  # Senescence as tumor suppressor

    # Cancer formation
    'cancer_threshold': 50,        # Mutation threshold
    'cancer_hits': 2.5,            # Multi-hit exponent
    'k_cancer_rate': 0.001,        # Baseline cancer rate
}

# =============================================================================
# TSC2 / AMPK-mTORC1 CANONICAL PATHWAY
# =============================================================================

TSC2_PARAMS: Dict[str, float] = {
    # TSC2 activation by AMPK
    # AMPK phosphorylates TSC2 at Thr1271, Ser1387 (Inoki et al. 2003, PMID 14651849)
    # AMPK priming enables GSK3-beta sequential phosphorylation (Inoki et al. 2006, PMID 16959574)
    'TSC2_basal': 0.30,            # Basal TSC2 GAP activity without AMPK
    'TSC2_min': 0.10,              # Minimum residual TSC2 activity
    'k_ampk_tsc2_max': 0.70,      # Maximum AMPK-induced TSC2 activation
    'Km_ampk_tsc2': 0.35,         # Half-max AMPK for TSC2 activation
    'n_hill_ampk_tsc2': 1.5,      # Hill coefficient (AMPK/GSK3 cooperativity)

    # TSC2 inhibition of mTORC1 via Rheb inactivation
    # TSC2 GAP converts Rheb-GTP to Rheb-GDP (Huang & Manning 2008, PMID 18466115)
    # Both TSC2 and Raptor arms required for full effect (Kazyken et al. 2020, PMID 32912901)
    'k_tsc2_mtorc1_max': 0.65,    # Maximum TSC2-mediated mTORC1 inhibition
    'Km_tsc2_mtorc1': 0.50,       # Half-max TSC2 for mTORC1 inhibition
    'n_hill_tsc2_mtorc1': 1.2,    # Near-MM (single GTPase switch)
}

# =============================================================================
# BIOLOGICAL AGE INTEGRATION
# =============================================================================

BIOAGE_PARAMS: Dict[str, float] = {
    # Component weights (sum to 1.0)
    'w_meth': 0.30,                # Methylation weight
    'w_damage': 0.30,              # DNA damage weight
    'w_hetero': 0.25,              # Heteroplasmy weight (coupled with damage via PARP/NAD)
    'w_sen': 0.15,                 # Senescence weight
    'w_sasp': 0.00,                # SASP weight (currently 0, captured via senescence)

    # Death threshold
    'death_threshold': 1.0,        # BioAge at death (normalized)

    # Normalization constants — derived from control simulation at t=1.0
    # so that BioAge(t=1.0) = 1.0 by construction. ODE dynamics are
    # independent of these norms (no circularity).
    # Derivation: model.derive_normalization_constants()
    'meth_norm': 1.093528,         # Methylation at t=1.0 in control
    'damage_norm': 0.665621,       # DNA damage at t=1.0 in control
    'H_norm': 0.179975,            # Heteroplasmy at t=1.0 in control
    'sen_norm': 0.217642,          # SenCells at t=1.0 in control
    'SASP_norm': 0.038512,         # SASP at t=1.0 in control (w_sasp=0)
}

# =============================================================================
# CROSS-SPECIES TRANSLATION
# =============================================================================

# Pathway-specific dampening for mouse-to-human translation
CROSS_SPECIES_DAMPENING: Dict[str, float] = {
    'mtor_pathway': 0.35,          # mTOR effects dampen in humans
    'antioxidant': 0.25,           # Antioxidants less effective in humans
    'senolytic': 0.40,             # Senolytics may translate better
    'metabolic': 0.30,             # Metabolic interventions
}

# Intervention-type specific dampening
INTERVENTION_DAMPENING: Dict[str, float] = {
    'cr_like': 0.20,               # Caloric restriction mimetics
    'mtor_pathway': 0.35,          # mTOR inhibitors
    'ampk_pathway': 0.30,          # AMPK activators
    'nad_pathway': 0.40,           # NAD+ boosters
    'senolytic': 0.50,             # Senolytics
    'epigenetic': 0.35,            # Epigenetic modifiers
    'antioxidant': 0.25,           # Antioxidants
    'antiinflam': 0.30,            # Anti-inflammatories
}

# Global mouse-to-human scaling factor
HUMAN_SCALING_FACTOR: float = 0.32


# =============================================================================
# CONVENIENCE EXPORTS
# =============================================================================

# Collect all pathway parameters for easy iteration
ALL_PATHWAY_PARAMS = {
    'nad': NAD_PARAMS,
    'nmn_saturation': NMN_SATURATION_PARAMS,
    'cd38_saturation': CD38_SATURATION_PARAMS,
    'ampk': AMPK_PARAMS,
    'mtorc1': MTORC1_PARAMS,
    'sirtuin': SIRTUIN_PARAMS,
    'autophagy': AUTOPHAGY_PARAMS,
    'dna_repair': DNA_REPAIR_PARAMS,
    'methylation': METHYLATION_PARAMS,
    'heteroplasmy': HETEROPLASMY_PARAMS,
    'uprmt': UPRMT_PARAMS,
    'senescence': SENESCENCE_PARAMS,
    'sasp': SASP_PARAMS,
    'cancer': CANCER_PARAMS,
    'bioage': BIOAGE_PARAMS,
    'tsc2': TSC2_PARAMS,
}
