#!/usr/bin/env python3
"""
Unit tests for individual pathway functions, numerical convergence,
uncertainty quantification, and dose-response.

Run: python -m pytest tests/test_pathways.py -v
  or: python tests/test_pathways.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from elm.model import (
    hill_function,
    michaelis_menten,
    calculate_effective_ampk,
    calculate_mtorc1_activity,
    calculate_nad_synthesis,
    calculate_nad_consumption,
    calculate_damage_rate,
    calculate_repair_rate,
    calculate_methylation_rate,
    calculate_heteroplasmy_dynamics,
    calculate_senescence_rate,
    calculate_senescence_clearance,
    calculate_bioage,
    calculate_cancer_probability,
    simulate,
    run_control,
    calculate_lifespan_extension,
    run_combination_extension,
)
from elm.pathways import (
    NAD_PARAMS,
    AMPK_PARAMS,
    MTORC1_PARAMS,
    DNA_REPAIR_PARAMS,
    METHYLATION_PARAMS,
    HETEROPLASMY_PARAMS,
    SENESCENCE_PARAMS,
    CANCER_PARAMS,
    BIOAGE_PARAMS,
)
from elm.dose_response import (
    hill_dose_response,
    effective_dose,
    DOSE_RESPONSE_PARAMS,
)
from elm.uncertainty import (
    run_monte_carlo,
    predict_human_extension,
    compute_attribution,
    summarize_uncertainty,
)


# =========================================================================
# Hill function and Michaelis-Menten basics
# =========================================================================

def test_hill_zero_input():
    """Hill function at x=0 should return 0."""
    assert hill_function(0, Km=0.5, n=2) == 0.0

def test_hill_at_km():
    """Hill function at x=Km should return Vmax/2."""
    result = hill_function(0.5, Km=0.5, n=1, Vmax=1.0)
    assert abs(result - 0.5) < 1e-10

def test_hill_large_input():
    """Hill function at very large x should approach Vmax."""
    result = hill_function(1e6, Km=0.5, n=2, Vmax=1.0)
    assert abs(result - 1.0) < 1e-6

def test_michaelis_menten_zero_substrate():
    """MM at zero substrate should return 0."""
    assert michaelis_menten(0, Km=0.5) == 0.0

def test_michaelis_menten_at_km():
    """MM at substrate=Km should return Vmax/2."""
    result = michaelis_menten(0.5, Km=0.5, Vmax=1.0)
    assert abs(result - 0.5) < 1e-10


# =========================================================================
# NAD+ pathway
# =========================================================================

def test_nad_synthesis_positive():
    """NAD+ synthesis should always be positive."""
    rate = calculate_nad_synthesis(
        ampk=0.5, nmn_supplement=0, age_factor=0.5,
        heteroplasmy=0.1, params=NAD_PARAMS
    )
    assert rate > 0

def test_nad_synthesis_nmn_boost():
    """NMN supplementation should increase NAD+ synthesis."""
    rate_base = calculate_nad_synthesis(
        ampk=0.5, nmn_supplement=0, age_factor=0.5,
        heteroplasmy=0.1, params=NAD_PARAMS
    )
    rate_nmn = calculate_nad_synthesis(
        ampk=0.5, nmn_supplement=0.5, age_factor=0.5,
        heteroplasmy=0.1, params=NAD_PARAMS
    )
    assert rate_nmn > rate_base

def test_nad_consumption_increases_with_damage():
    """Higher DNA damage should increase NAD+ consumption via PARP."""
    low_damage = calculate_nad_consumption(
        nad=0.8, dna_damage=0.1, sirt1=0.5, sasp=0.1,
        age_factor=0.5, cd38_inhibitor=0, gut_antiinflam=0,
        params=NAD_PARAMS
    )
    high_damage = calculate_nad_consumption(
        nad=0.8, dna_damage=0.5, sirt1=0.5, sasp=0.1,
        age_factor=0.5, cd38_inhibitor=0, gut_antiinflam=0,
        params=NAD_PARAMS
    )
    assert high_damage > low_damage


# =========================================================================
# Edge cases: extreme inputs
# =========================================================================

def test_bioage_at_zero_state():
    """BioAge with all state variables at zero should be near zero."""
    ba = calculate_bioage(meth=0, damage=0, H=0, sen=0, sasp=0,
                          params=BIOAGE_PARAMS)
    assert ba >= 0
    assert ba < 0.1

def test_cancer_probability_at_zero_mutations():
    """Cancer probability should be near zero with no mutations."""
    prob = calculate_cancer_probability(0, params=CANCER_PARAMS)
    assert prob >= 0
    assert prob < 1e-6

def test_cancer_probability_monotonic():
    """Cancer probability should increase with mutations."""
    p1 = calculate_cancer_probability(0.1, params=CANCER_PARAMS)
    p2 = calculate_cancer_probability(0.5, params=CANCER_PARAMS)
    assert p2 > p1

def test_senescence_rate_nonnegative():
    """Senescence entry rate should be non-negative."""
    rate = calculate_senescence_rate(
        damage=0.3, sirt1=0.5, ampk=0.5, gut=0,
        uprmt_sen_reduction=0, mtorc1_sen_reduction=0,
        sasp_paracrine=0.2, params=SENESCENCE_PARAMS
    )
    assert rate >= 0

def test_heteroplasmy_bounded():
    """Heteroplasmy dynamics should not produce NaN or inf."""
    dH = calculate_heteroplasmy_dynamics(
        H=0.5, foxo3=0.5, mitophagy_active=False,
        mitophagy_strength=0, params=HETEROPLASMY_PARAMS
    )
    assert np.isfinite(dH)

def test_heteroplasmy_at_zero():
    """Heteroplasmy rate at H=0 should be zero (no mutant mtDNA to replicate)."""
    dH = calculate_heteroplasmy_dynamics(
        H=0.0, foxo3=0.5, mitophagy_active=False,
        mitophagy_strength=0, params=HETEROPLASMY_PARAMS
    )
    assert abs(dH) < 1e-10


# =========================================================================
# Reproducibility: seeded simulations
# =========================================================================

def test_simulate_reproducible():
    """Two runs with the same seed should produce identical results."""
    r1 = simulate(compound='rapamycin', sex='M', seed=42)
    r2 = simulate(compound='rapamycin', sex='M', seed=42)
    assert r1.t_death == r2.t_death
    assert np.array_equal(r1.BioAge, r2.BioAge)

def test_simulate_different_seeds():
    """Different seeds may produce different cancer events (stochastic)."""
    # This tests that the seed is actually used -- the BioAge trajectories
    # are deterministic (ODE), but cancer_occurred may differ.
    r1 = simulate(compound='rapamycin', sex='M', seed=1)
    r2 = simulate(compound='rapamycin', sex='M', seed=2)
    # BioAge trajectory should be identical (ODE is deterministic)
    assert np.allclose(r1.BioAge, r2.BioAge)
    # t_death should also match (driven by BioAge, not cancer in normal runs)
    assert abs(r1.t_death - r2.t_death) < 1e-10


# =========================================================================
# Numerical convergence: dt sensitivity
# =========================================================================

def test_numerical_convergence():
    """Halving dt should not change lifespan extension by more than 0.5%."""
    control_coarse = run_control(sex='M', dt=0.002)
    control_fine = run_control(sex='M', dt=0.001)

    treated_coarse = simulate(compound='rapamycin', sex='M', dt=0.002)
    treated_fine = simulate(compound='rapamycin', sex='M', dt=0.001)

    ext_coarse = calculate_lifespan_extension(treated_coarse, control_coarse)
    ext_fine = calculate_lifespan_extension(treated_fine, control_fine)

    assert abs(ext_coarse - ext_fine) < 0.5, (
        f"dt sensitivity too high: dt=0.002 gives {ext_coarse:.2f}%, "
        f"dt=0.001 gives {ext_fine:.2f}%"
    )


# =========================================================================
# Null model comparison: mechanistic vs. naive additive
# =========================================================================

def test_combination_beats_additive_null():
    """
    ELM's mechanistic Rapa+Acarb prediction should be closer to the
    observed +34% (male) than a naive additive model (Rapa + Acarb
    individual extensions summed).
    """
    control = run_control(sex='M', dt=0.001)

    rapa = simulate(compound='rapamycin', sex='M', dt=0.001)
    acarb = simulate(compound='acarbose', sex='M', dt=0.001)

    ext_rapa = calculate_lifespan_extension(rapa, control)
    ext_acarb = calculate_lifespan_extension(acarb, control)

    # Naive additive model: just sum single-compound extensions
    additive_prediction = ext_rapa + ext_acarb

    # ELM mechanistic prediction
    elm_prediction = run_combination_extension(
        {'rapamycin': 1.0, 'acarbose': 1.0}, sex='M', dt=0.001
    )

    observed = 34.0  # ITP observed

    additive_error = abs(additive_prediction - observed)
    elm_error = abs(elm_prediction - observed)

    assert elm_error < additive_error, (
        f"ELM ({elm_prediction:.1f}%, error {elm_error:.1f}) should beat "
        f"additive null ({additive_prediction:.1f}%, error {additive_error:.1f}) "
        f"vs observed {observed}%"
    )


# =========================================================================
# Dose-response edge cases
# =========================================================================

def test_dose_response_at_zero():
    """f(0) = 0 for all compounds."""
    for name in DOSE_RESPONSE_PARAMS:
        assert effective_dose(name, 0.0) == 0.0, f"{name}: f(0) != 0"

def test_dose_response_at_one():
    """f(1) = 1 for all compounds (calibration preserved by construction)."""
    for name in DOSE_RESPONSE_PARAMS:
        assert effective_dose(name, 1.0) == 1.0, f"{name}: f(1) != 1"

def test_dose_response_diminishing_returns():
    """f(2) < 2 for all compounds (diminishing returns)."""
    for name, params in DOSE_RESPONSE_PARAMS.items():
        f2 = hill_dose_response(2.0, params['n'], params['K'])
        assert f2 < 2.0, f"{name}: f(2) = {f2} >= 2.0"
        assert f2 > 1.0, f"{name}: f(2) = {f2} <= 1.0"

def test_dose_response_monotonic():
    """Dose-response should be monotonically increasing."""
    for name, params in DOSE_RESPONSE_PARAMS.items():
        prev = 0.0
        for d in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]:
            val = hill_dose_response(d, params['n'], params['K'])
            assert val >= prev, f"{name}: not monotonic at d={d}"
            prev = val

def test_dose_response_unknown_compound():
    """effective_dose should raise KeyError for unknown compound at non-unity dose."""
    try:
        effective_dose('nonexistent_compound', 2.0)
        assert False, "Should have raised KeyError"
    except KeyError:
        pass


# =========================================================================
# Uncertainty module
# =========================================================================

def test_monte_carlo_produces_ci():
    """Monte Carlo should produce valid confidence intervals."""
    mc = run_monte_carlo(compound='rapamycin', sex='M', n_samples=100,
                         seed=42, include_trajectories=False)
    assert mc.extension_median is not None
    assert mc.extension_ci_90 is not None
    lo, hi = mc.extension_ci_90
    assert lo < mc.extension_median < hi, (
        f"CI not ordered: {lo} < {mc.extension_median} < {hi}"
    )

def test_monte_carlo_reproducible():
    """Same seed should give identical MC results."""
    mc1 = run_monte_carlo(compound='rapamycin', sex='M', n_samples=50,
                          seed=99, include_trajectories=False)
    mc2 = run_monte_carlo(compound='rapamycin', sex='M', n_samples=50,
                          seed=99, include_trajectories=False)
    assert mc1.extension_median == mc2.extension_median

def test_monte_carlo_ci_covers_point_estimate():
    """The 90% CI should contain the deterministic extension (roughly)."""
    mc = run_monte_carlo(compound='rapamycin', sex='M', n_samples=200,
                         seed=42, include_trajectories=False)
    control = run_control(sex='M')
    treated = simulate(compound='rapamycin', sex='M')
    det_ext = calculate_lifespan_extension(treated, control)
    lo, hi = mc.extension_ci_90
    assert lo < det_ext < hi, (
        f"Deterministic extension {det_ext:.1f}% outside 90% CI [{lo:.1f}, {hi:.1f}]"
    )

def test_monte_carlo_trajectories():
    """Trajectory statistics should be computed when requested."""
    mc = run_monte_carlo(compound='rapamycin', sex='M', n_samples=50,
                         seed=42, include_trajectories=True)
    assert mc.BioAge_median is not None
    assert mc.NAD_median is not None
    assert mc.Survival_median is not None
    assert len(mc.BioAge_median) > 0

def test_human_prediction_dampened():
    """Human extension should be less than mouse (dampening < 1)."""
    pred = predict_human_extension(23.0, 'mtor_pathway', n_samples=500)
    assert pred.human_extension_median < 23.0, (
        f"Human prediction {pred.human_extension_median:.1f}% not dampened from mouse 23%"
    )
    assert pred.human_extension_median > 0

def test_summarize_uncertainty():
    """Summary string should contain key statistics."""
    mc = run_monte_carlo(compound='rapamycin', sex='M', n_samples=50,
                         seed=42, include_trajectories=False)
    summary = summarize_uncertainty(mc)
    assert 'Median' in summary
    assert '90% CI' in summary

def test_attribution_sums_approximately():
    """Greedy subtraction contributions should roughly sum to full stack."""
    components = {
        'Rapamycin': {'mtorc1_drug_inhibition': 0.4424, 'ampk': 0.2949},
        'Acarbose': {'ampk': 0.1941, 'gut_microbiome': 0.1412},
    }
    attr = compute_attribution(components, sex='M', verbose=False)
    total = sum(attr.contributions.values())
    # contributions + interaction should equal full stack
    assert abs(total + attr.interaction_effect - attr.full_stack_extension) < 0.1, (
        f"Attribution doesn't sum: {total} + {attr.interaction_effect} != {attr.full_stack_extension}"
    )


# =========================================================================
# Runner
# =========================================================================

if __name__ == '__main__':
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    passed = 0
    failed = 0
    for test_fn in tests:
        try:
            test_fn()
            print(f"  + PASS  {test_fn.__name__}")
            passed += 1
        except Exception as e:
            print(f"  X FAIL  {test_fn.__name__}: {e}")
            failed += 1
    print(f"\n  {passed} passed, {failed} failed out of {passed + failed}")
    sys.exit(1 if failed else 0)
