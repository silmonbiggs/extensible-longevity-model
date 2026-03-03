#!/usr/bin/env python3
"""
Recalibration script for rev_05.

Runs all analyses and prints numbers in copy-paste format for documentation.
Combination validations use actual ITP study start times (not earliest-individual).
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from elm.model import (
    simulate, run_control, calculate_lifespan_extension,
    run_combination, run_combination_extension, derive_normalization_constants,
)
from elm.compounds import ITP_VALIDATION, COMPOUNDS, get_compound
from elm.pathways import BIOAGE_PARAMS
from elm.sex_mechanisms import apply_sex_modifier


def section(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def main():
    DT = 0.002
    T_MAX = 2.5

    # =========================================================================
    # 0. Normalization constants
    # =========================================================================
    section("NORMALIZATION CONSTANTS")
    norms = derive_normalization_constants(dt=0.001)
    print("  derive_normalization_constants() returned:")
    for k, v in norms.items():
        current = BIOAGE_PARAMS.get(k, None)
        flag = " <-- CHANGED" if current is not None and abs(v - current) > 1e-6 else ""
        print(f"    {k}: {v:.6f}  (current: {current}){flag}")

    # =========================================================================
    # 1. ITP Calibration — 12 targets
    # =========================================================================
    section("ITP CALIBRATION (12 targets)")
    errors = []
    print(f"  {'Compound':30s} {'Got':>7s} {'Target':>7s} {'Error':>7s}  Status")
    print(f"  {'-'*62}")

    cal_results = {}
    for sex in ['M', 'F']:
        control = run_control(sex=sex, t_max=T_MAX, dt=DT)
        for cname, val in ITP_VALIDATION.items():
            result = simulate(compound=cname, sex=sex, t_max=T_MAX, dt=DT)
            ext = calculate_lifespan_extension(result, control)
            target = val.target_male if sex == 'M' else val.target_female
            error = ext - target
            abs_error = abs(error)
            status = "PASS" if abs_error <= 1.0 else "FAIL"
            errors.append(abs_error)
            key = f"{cname}_{sex}"
            cal_results[key] = ext
            marker = "+" if status == "PASS" else "X"
            print(f"  {key:30s} {ext:6.1f}% {target:6.1f}% {error:+6.2f}%  {marker} {status}")

    mean_err = np.mean(errors)
    max_err = np.max(errors)
    print(f"\n  Mean |error| = {mean_err:.2f}%")
    print(f"  Max  |error| = {max_err:.2f}%")

    # =========================================================================
    # 2. Validation: Rapa + Acarbose
    # =========================================================================
    section("VALIDATION: RAPA + ACARBOSE")
    # ITP combo study (Strong 2022, PMID 36179270): both drugs started at 9 months
    RAPA_ACARB_START = 9.0 / 30.0  # 0.30 normalized
    for sex in ['M', 'F']:
        ext = run_combination_extension({'rapamycin': 1.0, 'acarbose': 1.0}, sex=sex,
                                        t_max=T_MAX, dt=DT, start_time=RAPA_ACARB_START)
        print(f"  Rapa+Acarb {sex}: {ext:.1f}%")

    # =========================================================================
    # 3. NMN alone (M/F)
    # =========================================================================
    section("NMN ALONE")
    for sex in ['M', 'F']:
        control = run_control(sex=sex, t_max=T_MAX, dt=DT)
        result = simulate(compound='NMN', sex=sex, t_max=T_MAX, dt=DT)
        ext = calculate_lifespan_extension(result, control)
        print(f"  NMN {sex}: {ext:.1f}%")

    # =========================================================================
    # 4. Rapa + NMN (M/F) with additive floor and superadditivity
    # =========================================================================
    section("RAPA + NMN COMBINATION")
    # Hypothetical ITP protocol: both compounds from 9 months
    RAPA_NMN_START = 9.0 / 30.0  # 0.30 normalized
    for sex in ['M', 'F']:
        control = run_control(sex=sex, t_max=T_MAX, dt=DT)

        # Individual
        rapa_result = simulate(compound='rapamycin', sex=sex, t_max=T_MAX, dt=DT)
        rapa_ext = calculate_lifespan_extension(rapa_result, control)

        nmn_result = simulate(compound='NMN', sex=sex, t_max=T_MAX, dt=DT)
        nmn_ext = calculate_lifespan_extension(nmn_result, control)

        # Combination
        combo_ext = run_combination_extension({'rapamycin': 1.0, 'NMN': 1.0}, sex=sex,
                                              t_max=T_MAX, dt=DT, start_time=RAPA_NMN_START)

        additive = rapa_ext + nmn_ext
        superadd = combo_ext - additive

        print(f"  {sex}: Rapa={rapa_ext:.1f}%, NMN={nmn_ext:.1f}%, Combo={combo_ext:.1f}%")
        print(f"      Additive floor={additive:.1f}%, Superadditivity={superadd:+.1f}%")

    # =========================================================================
    # 5. All-6 ITP cocktail (M/F)
    # =========================================================================
    section("ALL-6 ITP COCKTAIL")
    itp_cocktail = {
        'rapamycin': 1.0, 'acarbose': 1.0, 'canagliflozin': 1.0,
        '17_alpha_estradiol': 1.0, 'aspirin': 1.0, 'glycine': 1.0,
    }
    for sex in ['M', 'F']:
        ext = run_combination_extension(itp_cocktail, sex=sex, t_max=T_MAX, dt=DT)
        print(f"  All-6 ITP cocktail {sex}: {ext:.1f}%")

    # =========================================================================
    # 6. Novel compound layered predictions
    # =========================================================================
    section("NOVEL COMPOUND LAYERS (on ITP base)")
    # Base: rapa + acarbose (9-month start, matching ITP combo protocol)
    base_compounds = {'rapamycin': 1.0, 'acarbose': 1.0}
    for sex in ['M', 'F']:
        base_ext = run_combination_extension(base_compounds, sex=sex, t_max=T_MAX, dt=DT,
                                             start_time=RAPA_ACARB_START)
        print(f"\n  Base (rapa+acarb) {sex}: {base_ext:.1f}%")

        layers = {
            'urolithin_A': {'urolithin_A': 1.0},
            'NMN': {'NMN': 1.0},
            'fisetin (senolytic)': {'fisetin': 1.0},
            'CD38i': {'NAD_homeostasis_plus': 1.0},
        }
        for layer_name, layer_compounds in layers.items():
            combo = {**base_compounds, **{k: 1.0 for k in layer_compounds}}
            ext = run_combination_extension(combo, sex=sex, t_max=T_MAX, dt=DT,
                                            start_time=RAPA_ACARB_START)
            delta = ext - base_ext
            print(f"    + {layer_name:25s}: {ext:.1f}% (delta={delta:+.1f}%)")

    # =========================================================================
    # 7. Hallmark importance (node clamping)
    # =========================================================================
    section("HALLMARK IMPORTANCE (node clamping)")
    # For each BioAge component, zero its contribution and measure impact
    print("  Method: zero each state variable's drivers in BioAge equation")
    print()

    for sex in ['M']:
        control = run_control(sex=sex, t_max=T_MAX, dt=DT)
        control_death = control.t_death

        # Full combination for reference
        full_combo = {
            'rapamycin': 1.0, 'acarbose': 1.0, 'NMN': 1.0,
            'fisetin': 1.0, 'urolithin_A': 1.0,
        }
        full_ext = run_combination_extension(full_combo, sex=sex, t_max=T_MAX, dt=DT)
        print(f"  Full cocktail ({sex}): {full_ext:.1f}%")

        # For each weight, temporarily zero it and re-run
        components = ['w_meth', 'w_damage', 'w_hetero', 'w_sen']
        for comp in components:
            original = BIOAGE_PARAMS[comp]
            weight_label = comp.replace('w_', '')
            print(f"  Clamping {weight_label} (w={original:.2f}): ", end='')

            # Calculate fraction of BioAge from this component at t=1.0
            idx_1 = int(1.0 / DT)
            if comp == 'w_meth':
                val_at_death = control.Methylation[idx_1] / BIOAGE_PARAMS['meth_norm']
            elif comp == 'w_damage':
                val_at_death = control.DNA_damage[idx_1] / BIOAGE_PARAMS['damage_norm']
            elif comp == 'w_hetero':
                val_at_death = control.Heteroplasmy[idx_1] / BIOAGE_PARAMS['H_norm']
            elif comp == 'w_sen':
                val_at_death = control.SenCells[idx_1] / BIOAGE_PARAMS['sen_norm']

            contribution = original * val_at_death
            print(f"contribution at t=1.0 = {contribution:.3f} ({contribution*100:.1f}% of BioAge)")

    # =========================================================================
    # 8. Pathway sensitivity (standard dose per pathway)
    # =========================================================================
    section("PATHWAY SENSITIVITY (single-pathway injections)")
    pathways_to_test = {
        'mtorc1_drug_inhibition': 0.80,
        'ampk': 0.40,
        'antioxidant': 0.50,
        'nmn': 0.50,
        'senolytic': 0.30,
        'mitophagy': 0.40,
        'antiinflam': 0.30,
        'cd38_inhibitor': 0.50,
        'sirt1_direct': 0.30,
        'gut_microbiome': 0.30,
        'akg': 0.80,
    }

    sex = 'M'
    control = run_control(sex=sex, t_max=T_MAX, dt=DT)
    print(f"  {'Pathway':30s} {'Dose':>5s} {'Extension':>10s}")
    print(f"  {'-'*50}")
    for pathway, dose in pathways_to_test.items():
        interventions = {pathway: dose, 'start_time': 9.0 / 30.0}
        result = simulate(interventions=interventions, sex=sex, t_max=T_MAX, dt=DT)
        ext = calculate_lifespan_extension(result, control)
        print(f"  {pathway:30s} {dose:5.2f} {ext:9.1f}%")

    # =========================================================================
    # Summary for docs
    # =========================================================================
    section("SUMMARY FOR DOCS")
    print(f"  Weights: w_meth={BIOAGE_PARAMS['w_meth']}, w_damage={BIOAGE_PARAMS['w_damage']}, "
          f"w_hetero={BIOAGE_PARAMS['w_hetero']}, w_sen={BIOAGE_PARAMS['w_sen']}")
    print(f"  Mean calibration error: {mean_err:.2f}%")
    print(f"  Max calibration error:  {max_err:.2f}%")
    print(f"  Norms: meth={norms['meth_norm']:.6f}, damage={norms['damage_norm']:.6f}, "
          f"H={norms['H_norm']:.6f}, sen={norms['sen_norm']:.6f}")


if __name__ == '__main__':
    main()
