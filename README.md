# Extensible Longevity Model (ELM)

[![Smoke Test](https://github.com/silmonbiggs/extensible-longevity-model/actions/workflows/smoke-test.yml/badge.svg)](https://github.com/silmonbiggs/extensible-longevity-model/actions/workflows/smoke-test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Release](https://img.shields.io/github/v/release/silmonbiggs/extensible-longevity-model)](https://github.com/silmonbiggs/extensible-longevity-model/releases/tag/v0.43.0)

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

### Observed vs predicted

Male predictions are calibrated by root-finding (error < 0.2 pp).
Female predictions are genuine out-of-sample.

| Compound | Male observed | Male predicted | Female observed | Female predicted |
|----------|:---:|:---:|:---:|:---:|
| Rapamycin | 23.0% | 23.2% | 26.0% | 26.8% |
| Acarbose | 22.0% | 22.2% | 5.0% | 7.8% |
| 17-alpha-estradiol | 19.0% | 19.2% | 0.0% | 0.0% |
| Canagliflozin | 14.0% | 14.0% | 9.0% | 8.8% |
| Aspirin | 8.0% | 8.2% | 0.0% | 0.0% |
| Glycine | 6.0% | 6.2% | 4.0% | 3.9% |

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

Five aging hallmarks feed into biological age:

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

## Degrees of Freedom

The model has 168 parameters. Six are free: one scale factor per
ITP compound, found by root-finding so that each male prediction
matches the ITP observation exactly. The remaining 162 are frozen
before calibration — set from published measurements, biologically
plausible ranges, or structural constants. None are fitted to ITP
intervention data.

The female predictions add four sex-specific mechanisms drawn from
independent experimental literature. These mechanisms introduce
zero additional free parameters. The 1.3 pp female error is not a
fit — it is a test.

The rapamycin + acarbose combination prediction (+33.5% vs +34%
observed) uses the same single-compound scale factors with no
re-fitting. Sub-additivity emerges from shared pathway saturation
in the ODE system.

The base model assigns each compound a single scale factor that
scales a fixed pathway profile. For five of six compounds this is
sufficient. Canagliflozin's residual pattern suggested that its
default pathway ratio (the balance between SGLT2-mediated and
AMPK-mediated effects) was inadequate. The ratio was then optimized,
adding a seventh fitted parameter. The hypothesis is data-suggested;
the parameter is fitted. This makes it stronger than ad hoc but
weaker than a pre-registered prediction.

In short: 6 numbers are fitted (to male data), with a seventh
(canagliflozin pathway split) motivated by the residual pattern.
Everything else — female predictions, combination predictions,
sensitivity rankings — is a consequence of the model structure.

## Limitations and Failure Modes

What would make this model wrong, and what does it not attempt:

- **BioAge threshold death rule.** The mouse dies when a weighted
  composite crosses a fixed threshold. Real mice die of specific
  causes (cancer, organ failure). The model cannot distinguish
  cause of death and would miss interventions that shift mortality
  from one cause to another without changing overall biological age.
- **Equal BioAge weights.** The four hallmark weights are set to
  0.25 each (equal prior). Weight sweeps show predictions are
  robust within the calibration-acceptable range, but the true
  weights are unknown. This is the largest structural assumption.
- **Six compounds, six parameters.** The system is exactly
  determined for males — there are no residual degrees of freedom
  to detect male model error. Male predictions pass by
  construction. All falsifiability comes from female and
  combination predictions.
- **ITP compounds only.** The model has been tested against six
  ITP compounds that act through AMPK, mTORC1, and related
  pathways. It says nothing about compounds that act through
  mechanisms not in the ODE system (for example, senolytics,
  gene therapy, or direct epigenetic reprogramming).
- **Mouse only.** All predictions are in mouse normalized lifespan
  units. No mouse-to-human translation is applied or claimed.
- **Frozen prior uncertainty.** The 162 frozen parameters are not
  known exactly. OAT sensitivity shows none shift predictions by
  more than 0.65 pp under ±10% perturbation, but correlated
  errors across multiple parameters could have larger effects.
- **No stochastic variation.** The ODE system is deterministic.
  Real aging has stochastic components (for example, mutation
  timing, immune stochasticity) that the model does not capture.
- **Acarbose female prediction is the worst.** Predicted 7.8%,
  observed 5.0% (2.8 pp error). This is within tolerance but is
  the largest single error, likely reflecting incomplete modeling
  of gut microbiome sex differences.

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
should work. CI is tested with Python 3.10 and 3.12. The release
was developed and tested against numpy 2.3, scipy 1.16, and
matplotlib 3.10.

**Reproduce the paper from scratch** (one command):

```bash
python scripts/reproduce_all.py        # runs all three steps below (~23 min)
```

Or run each step individually:

```bash
python tests/smoke_test.py              # verify 12 ITP targets (Table 1)
python scripts/generate_figures.py      # regenerate all figures -> docs/figures/
python scripts/oat_sensitivity.py       # OAT analysis -> oat_sensitivity.tsv
```

**Expected outputs:**
- `docs/figures/` — all paper and supplementary figures (PNGs)
- `oat_sensitivity.tsv` — sensitivity of all 141 frozen parameters

A successful smoke test prints this table:

```
=================================================================
  ELM Smoke Test
=================================================================
  Compound                           Got  Target   Error  Status
-----------------------------------------------------------------
  17_alpha_estradiol_F              0.0%    0.0%   0.00%  + PASS
  17_alpha_estradiol_M             19.2%   19.0%   0.20%  + PASS
  acarbose_F                        7.8%    5.0%   2.80%  + PASS
  acarbose_M                       22.2%   22.0%   0.18%  + PASS
  aspirin_F                         0.0%    0.0%   0.00%  + PASS
  aspirin_M                         8.2%    8.0%   0.19%  + PASS
  canagliflozin_F                   8.8%    9.0%   0.18%  + PASS
  canagliflozin_M                  14.0%   14.0%   0.01%  + PASS
  glycine_F                         3.9%    4.0%   0.06%  + PASS
  glycine_M                         6.2%    6.0%   0.20%  + PASS
  rapamycin_F                      26.8%   26.0%   0.84%  + PASS
  rapamycin_M                      23.2%   23.0%   0.19%  + PASS
-----------------------------------------------------------------
  All 12 targets within tolerance (M: ±1.0pp, F: ±5.0pp). PASS.
```

Values in the "Got" column may vary by up to ~0.01 pp across
platforms due to floating-point differences.

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
  reproduce_all.py    #   One-command full replication
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

> Biggs, S. J. (2026). Extensible Longevity Model (ELM) v0.43.0.
> https://github.com/silmonbiggs/extensible-longevity-model

**Contact:** james@longrun.bio
