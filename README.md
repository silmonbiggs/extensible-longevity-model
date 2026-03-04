# Extensible Longevity Model (ELM)

A mechanistic ODE model of aging interventions, calibrated to the
NIA Interventions Testing Program (ITP) data. The model captures
8 interconnected aging subsystems (heteroplasmy, NAD+ metabolism,
DNA damage, AMPK/mTORC1 signaling, autophagy, senescence, SASP,
methylation) and predicts combination effects -- including
interactions that single-compound testing cannot reveal.

## LLM-Readable Content

For AI-assisted review, point your LLM at one of these plain-text
files directly:

- **[llms.txt](https://silmonbiggs.github.io/extensible-longevity-model/llms.txt)**
  -- brief summary (model, key results, experimental design).
- **[llms-full.txt](https://silmonbiggs.github.io/extensible-longevity-model/llms-full.txt)**
  -- complete slide-by-slide transcript with supplementary analyses.

## Interactive Slide Deck

**[View the deck](https://silmonbiggs.github.io/extensible-longevity-model/deck.html)**
-- model architecture, calibration, validation, and predictions
with interactive network diagrams, equations, and figures.

## Key Results

- **Calibration:** 12 ITP targets (6 compounds x 2 sexes). Mean
  absolute error 0.08 percentage points in predicted vs. observed
  median lifespan extension (max single-target error 0.20 pp).
  10 effective degrees of freedom from 12 data points, determined
  by Jacobian SVD rank (see Slide A1 in the deck).
- **Out-of-sample validation:** Rapamycin + Acarbose -- predicted
  +33%, ITP observed +34%. Not fit -- predicted from single-compound
  calibration only.
- **Testable prediction:** Rapamycin + NMN should be superadditive
  (+31.8% M / +36.0% F vs. additive floor of +26.4% / +29.6%).
- **Full-stack prediction:** 8 Taguchi compounds at 1x calibrated
  dose: +168% M / +140% F. Under full-stack conditions, all tracked
  variables (NAD+, heteroplasmy, DNA damage, senescent cells,
  methylation, SASP) stay within the ranges the model was calibrated
  on. Details: Slide A17.
- **Experimental design:** Taguchi L18 orthogonal array for 8
  longevity compounds in turquoise killifish (N. furzeri) -- $200K,
  ~8 months, 1,000 fish (20 groups x 50). Primary endpoint: median
  lifespan by Cox proportional hazards. Screens main effects,
  dose-response, and dominant interaction signals from one
  experiment. Compare: a single ITP mouse study costs ~$1.5M
  over 3+ years.

## What Is Validated vs. Predicted

| Status | Claim |
|--------|-------|
| Validated | 12 single-compound ITP lifespan extensions |
| Validated | Rapa+Acarb combination: +33% predicted, +34% observed |
| Predicted (calibrated) | All ITP pairwise/higher-order combinations |
| Predicted (novel) | NMN, CD38i, D+Q, Urolithin A -- literature-estimated, not calibrated to lifespan data |
| Predicted (full stack) | +168% M / +140% F from 8 compounds at 1x dose |

## Project Structure

```
elm/                  # Python model package (ODE engine, pathways, compounds)
docs/                 # Interactive slide deck + figures (GitHub Pages)
scripts/              # Figure generation, calibration, analyses
  clock_pca/          # Methylation weight constraint analysis (7 PMIDs)
tests/                # Smoke tests (12 ITP targets within +/-1%)
```

## Reproducibility

Requires Python 3.10+. Install dependencies:

```bash
pip install -r requirements.txt
```

Verify calibration (all 12 ITP targets within +/-1%):

```bash
python tests/smoke_test.py
```

Regenerate all figures used in the slide deck:

```bash
python scripts/generate_figures.py
```

Regenerate a single figure set:

```bash
python scripts/generate_figures.py stage_16
```

Regenerate the LLM-readable text files from the slide deck:

```bash
python scripts/generate_llms_txt.py
```

## Quick Start

```python
from elm.model import run_control, simulate, calculate_lifespan_extension

# Run a control simulation
control = run_control(sex='M')

# Simulate rapamycin treatment
treated = simulate(compound='rapamycin', sex='M')

# Calculate lifespan extension
ext = calculate_lifespan_extension(treated, control)
print(f"Rapamycin male extension: {ext:.1f}%")  # -> 23.0%
```

## License

MIT
