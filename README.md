# Extensible Longevity Model (ELM)

[![Smoke Test](https://github.com/silmonbiggs/extensible-longevity-model/actions/workflows/smoke-test.yml/badge.svg)](https://github.com/silmonbiggs/extensible-longevity-model/actions/workflows/smoke-test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Release](https://img.shields.io/github/v/release/silmonbiggs/extensible-longevity-model)](https://github.com/silmonbiggs/extensible-longevity-model/releases/tag/v1.0.0)

A mechanistic ODE model of mouse aging calibrated to the NIA
Interventions Testing Program (ITP). Five state variables track
four kinds of biological damage — mitochondrial heteroplasmy,
DNA damage, epigenetic drift, and senescent cell accumulation —
coupled through NAD+ metabolism and SASP signaling. Six free
parameters (one scale factor per ITP compound) are fitted to male
lifespan data. Female predictions use four sex-specific mechanisms
from prior literature with no additional fitting.

## Key Results

- **Sex difference predictions:** 1.3 percentage point mean absolute
  error across all six female ITP outcomes, using zero female-fitted
  parameters. Four independently motivated sex mechanisms (AMPK
  saturation, CYP3A pharmacokinetics, testosterone dependency,
  COX-2/prostacyclin interference) are sufficient.
- **Combination validation:** Rapamycin + acarbose predicted +33.5%,
  ITP observed +34%. Not fitted — predicted from single-compound
  calibration alone. The model correctly identifies sub-additivity
  from shared pathway saturation.
- **Sensitivity:** One-at-a-time perturbation of 141 frozen parameters
  (±10%) shifts no prediction by more than 0.65 pp. Results are
  dominated by model structure, not individual parameter values.

## Installation

Requires Python 3.10 or later.

```bash
git clone https://github.com/silmonbiggs/extensible-longevity-model.git
cd extensible-longevity-model
pip install -e .             # installs elm package + dependencies (~10 s)
python tests/smoke_test.py   # verify all 12 predictions (< 1 min)
```

A single simulation (one compound, one sex) runs in about 1 second.
Regenerating all paper figures takes about 3 minutes.
The full sensitivity analysis takes about 19 minutes.

To view the interactive slide deck locally, open
`docs/ITPSexDiff.html` in any browser — no server needed.

## Quick Start

```python
from elm.model import run_control, simulate, calculate_lifespan_extension

control = run_control(sex='M')
treated = simulate(compound='rapamycin', sex='M')
ext = calculate_lifespan_extension(treated, control)
print(f"Rapamycin male extension: {ext:.1f}%")  # -> 23.0%
```

## How It Works

The model tracks five state variables across the mouse lifespan:

| Variable | What it tracks |
|----------|---------------|
| Heteroplasmy (*H*) | Fraction of mutant mitochondrial DNA |
| NAD+ level (*N*) | Cellular energy currency, normalized to young adult |
| DNA damage (*D*) | Double-strand break equivalents |
| Senescent cells (*S*) | p16-positive cell burden |
| Methylation drift (*M*) | Deviation from youthful CpG set points |

These feed into BioAge, a weighted composite. The mouse dies when
BioAge crosses a threshold. Longevity compounds work by perturbing
rate constants in the ODEs — rapamycin inhibits mTORC1, acarbose
activates AMPK, and so on.

Each compound has one free parameter: a scale factor found by
root-finding so that the predicted male lifespan extension matches
the ITP observation exactly. All other parameters (162 of 168) are
frozen priors set before calibration from published measurements,
biologically plausible ranges, or structural constants.

## Papers

- **Paper 1:** Four Biological Mechanisms Predict ITP Sex Differences
  from Male Data Alone
  ([source](docs/paper/ITPSexDiff/paper1_sex_differences.tex),
  [PDF](docs/paper/ITPSexDiff/paper1_sex_differences.pdf))
- **Supplementary Material:**
  Full parameter tables with provenance for all 168 parameters,
  identifiability analysis, and sensitivity results
  ([source](docs/paper/ITPSexDiff/supplementary.tex),
  [PDF](docs/paper/ITPSexDiff/supplementary.pdf))

## Reproducibility

**Determinism.** The simulation uses a fixed-step Euler loop
(`dt=0.002`) with no random seeding, so results are deterministic
on a given platform. Floating-point differences across OS/architecture
combinations can shift predictions by up to ~0.01 pp, well within
the smoke test tolerances.

**Dependencies.** `requirements.txt` uses range constraints (for
example, `numpy>=1.24,<3`), not exact pins — any compatible version
should work.

**Reproduce the paper from scratch:**

```bash
python tests/smoke_test.py              # verify 12 ITP targets (Table 1)
python scripts/generate_figures.py      # regenerate all figures -> docs/figures/
python scripts/oat_sensitivity.py       # OAT analysis -> oat_sensitivity.tsv
```

Expected outputs:
- `docs/figures/` — all paper and supplementary figures (PNGs)
- `oat_sensitivity.tsv` — sensitivity of all 141 frozen parameters

## Reviewer Guide

A six-step walkthrough for validating the model and results.

**1. Install and run the smoke test (< 1 min)**

```bash
pip install -e .
python tests/smoke_test.py
```

This checks all 12 ITP predictions (6 compounds × 2 sexes) against
published targets. Male predictions are calibrated by root-finding;
female predictions are genuine out-of-sample.

**2. Read the core model (~1200 lines across 3 files)**

- `elm/model.py` — ODE system, simulation loop, lifespan calculation
- `elm/compounds.py` — compound definitions, ITP targets, scale factors
- `elm/sex_mechanisms.py` — four sex-specific mechanisms

**3. Check parameter provenance**

Every parameter has a source annotation in `elm/pathways.py`. The
supplementary material (Table S1) lists all 168 parameters with
provenance categories: published measurement, biologically plausible
range, structural constant, or calibrated to control lifespan.
Zero parameters are calibrated to ITP intervention data.

**4. Regenerate all figures (~3 min)**

```bash
python scripts/generate_figures.py
```

Reproduces every figure in the paper from the model code.

**5. Browse the interactive deck**

[View the deck](https://silmonbiggs.github.io/extensible-longevity-model/ITPSexDiff.html)
— model architecture, calibration, sex difference mechanisms,
and validation results with interactive figures.

**6. Run sensitivity analysis (~19 min)**

```bash
python scripts/oat_sensitivity.py
```

One-at-a-time perturbation of all 141 frozen parameters (±10%).
Pre-computed results are in `oat_sensitivity.tsv` for quick inspection.

## Project Structure

```
elm/                  # Core model package
  model.py            #   ODE engine, simulation, lifespan calculation
  compounds.py        #   Compound definitions and ITP targets
  pathways.py         #   Pathway parameters (NAD+, AMPK, mTORC1, etc.)
  sex_mechanisms.py   #   Sex-specific modifiers
docs/
  paper/              #   LaTeX sources and PDFs
  figures/            #   All generated figures
  ITPSexDiff.html     #   Paper 1 slide deck (reviewer-ready)
  working_slides.html #   Master slide pool (all 61 slides)
scripts/
  generate_figures.py #   Regenerate all figures
  oat_sensitivity.py  #   One-at-a-time sensitivity analysis
  sweep_ampk_k.py     #   AMPK saturation parameter sweep
tests/
  smoke_test.py       #   Validates 12 ITP targets within ±1%
  test_pathways.py    #   Pathway unit tests
```

## Citation and Contact

**License:** MIT ([LICENSE](LICENSE))

**Citation:** See [CITATION.cff](CITATION.cff), or cite as:

> Biggs, S. J. (2026). Extensible Longevity Model (ELM) v1.0.0.
> https://github.com/silmonbiggs/extensible-longevity-model

**Contact:** james@longrun.bio
