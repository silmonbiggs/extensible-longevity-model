# Extensible Longevity Model (ELM)

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
  deck.html           #   Interactive slide deck
scripts/
  generate_figures.py #   Regenerate all figures
  oat_sensitivity.py  #   One-at-a-time sensitivity analysis
  sweep_ampk_k.py     #   AMPK saturation parameter sweep
tests/
  smoke_test.py       #   Validates 12 ITP targets within ±1%
  test_pathways.py    #   Pathway unit tests
```

## Reproducibility

Python 3.10+.

```bash
pip install -e .                    # install as editable package
python tests/smoke_test.py          # verify calibration
python scripts/generate_figures.py  # regenerate all figures
```

## Quick Start

```python
from elm.model import run_control, simulate, calculate_lifespan_extension

control = run_control(sex='M')
treated = simulate(compound='rapamycin', sex='M')
ext = calculate_lifespan_extension(treated, control)
print(f"Rapamycin male extension: {ext:.1f}%")  # -> 23.0%
```

## Interactive Deck

**[View the deck](https://silmonbiggs.github.io/extensible-longevity-model/deck.html)**
— model architecture, calibration, validation, and combination
predictions with interactive figures.

## License

MIT
