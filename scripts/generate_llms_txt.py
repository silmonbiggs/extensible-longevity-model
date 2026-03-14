#!/usr/bin/env python3
"""
Extract all slide content from the interactive HTML deck and generate
llms.txt (brief) and llms-full.txt (complete) for LLM accessibility.

Usage:
    python scripts/generate_llms_txt.py
"""

import re
import os
import textwrap
from pathlib import Path
from html import unescape

DOCS_DIR = Path(os.path.dirname(__file__)) / '..' / 'docs'
HTML_FILE = DOCS_DIR / 'ITPSexDiff.html'


def _safe_chr(m):
    """Convert a \\uXXXX escape to a character, replacing surrogates with ?."""
    cp = int(m.group(1), 16)
    if 0xD800 <= cp <= 0xDFFF:
        return ''  # skip surrogates (emoji halves)
    return chr(cp)


def ascii_normalize(text):
    """Replace Unicode punctuation and symbols with ASCII equivalents."""
    replacements = {
        '\u2014': ' -- ',   # em dash
        '\u2013': '-',      # en dash
        '\u2019': "'",      # right single quote
        '\u2018': "'",      # left single quote
        '\u201c': '"',      # left double quote
        '\u201d': '"',      # right double quote
        '\u2192': '->',     # rightwards arrow
        '\u2190': '<-',     # leftwards arrow
        '\u2265': '>=',     # greater than or equal
        '\u2264': '<=',     # less than or equal
        '\u2248': '~',      # approximately equal
        '\u00d7': 'x',      # multiplication sign
        '\u00b1': '+/-',    # plus-minus
        '\u0394': 'Delta',  # Greek capital delta
        '\u03b1': 'alpha',  # Greek alpha
        '\u03b2': 'beta',   # Greek beta
        '\u03ba': 'kappa',  # Greek kappa
        '\u03c3': 'sigma',  # Greek sigma
        '\u03b3': 'gamma',  # Greek gamma
        '\u2208': 'in',     # element of
        '\u2212': '-',      # minus sign
        '\u2260': '!=',     # not equal
        '\u221e': 'inf',    # infinity
        '\u00b2': '^2',     # superscript 2
        '\u00b3': '^3',     # superscript 3
        '\u2076': '^6',     # superscript 6
        '\u2019': "'",      # apostrophe
        '\u2026': '...',    # ellipsis
        '\u2714': '[done]', # check mark
        '\u2713': '[done]', # check mark variant
        '\ufe0f': '',       # variation selector
        '\u00b7': '*',      # middle dot
        '\u00b9': '^1',     # superscript 1
        '\u00bd': '1/2',    # one half
        '\u2077': '^7',     # superscript 7
        '\u2080': '_0',     # subscript 0
        '\u2082': '_2',     # subscript 2
        '\u2191': '^',      # up arrow
        '\u2193': 'v',      # down arrow
        '\u2205': '{}',     # empty set
        '\u25cf': '*',      # black circle
        '\u26a0': '[!]',    # warning
        '\u26a1': '[!]',    # lightning bolt
        '\u2705': '[done]', # green check
        '\u2753': '[?]',    # question mark
        '\U0001f527': '[wrench]',  # wrench emoji
        '\U0001f9ea': '[lab]',     # test tube emoji
    }
    for uni, ascii_eq in replacements.items():
        text = text.replace(uni, ascii_eq)
    return text


def strip_html(text):
    """Remove HTML tags, decode entities, clean up whitespace, ASCII-normalize."""
    # Decode surrogate pairs (emoji) first: \uD83E\uDDEA -> single codepoint
    def decode_surrogate_pair(m):
        hi = int(m.group(1), 16)
        lo = int(m.group(2), 16)
        cp = 0x10000 + (hi - 0xD800) * 0x400 + (lo - 0xDC00)
        return chr(cp)
    text = re.sub(r'\\u([dD][89abAB][0-9a-fA-F]{2})\\u([dD][c-fC-F][0-9a-fA-F]{2})',
                  decode_surrogate_pair, text)
    # Decode remaining unicode escapes like \u2014, \u2192, etc.
    text = re.sub(r'\\u([0-9a-fA-F]{4})', _safe_chr, text)
    # Remove HTML tags
    text = re.sub(r'<[^>]+>', ' ', text)
    # Decode HTML entities
    text = unescape(text)
    # Unescape JS string backslashes: \\frac -> \frac, \\text -> \text, etc.
    text = text.replace('\\\\', '\x00BACKSLASH\x00')  # protect literal \\
    text = text.replace('\\n', '\n')
    text = text.replace('\\t', ' ')
    text = text.replace("\\'", "'")
    text = text.replace('\x00BACKSLASH\x00', '\\')  # restore single backslash
    # ASCII-normalize Unicode punctuation
    text = ascii_normalize(text)
    # Collapse whitespace
    text = re.sub(r'\s+', ' ', text).strip()
    return text


def extract_field(stage_text, field):
    """Extract a JS string field from a stage object."""
    # Match field:'...' or field:"..." with string concatenation
    pattern = rf"{field}\s*:\s*'((?:[^'\\]|\\.)*)'"
    match = re.search(pattern, stage_text)
    if match:
        return match.group(1)

    # Try to match concatenated strings: field:'...' + '...' + '...'
    # First find where the field starts
    start_pattern = rf"{field}\s*:\s*'"
    start_match = re.search(start_pattern, stage_text)
    if not start_match:
        return None

    # Collect all concatenated string fragments
    pos = start_match.end()
    fragments = []
    remaining = stage_text[start_match.start():]

    # Extract all '...' fragments connected by +
    frag_pattern = r"'((?:[^'\\]|\\.)*)'"
    for m in re.finditer(frag_pattern, remaining):
        fragments.append(m.group(1))
        # Check if there's a + after this fragment (continuing concatenation)
        after = remaining[m.end():m.end() + 20].strip()
        if not after.startswith('+'):
            break

    return ''.join(fragments) if fragments else None


def extract_stages(html_content):
    """Extract stage objects from the STAGES array in the HTML."""
    # Find the STAGES array
    stages_match = re.search(r'const STAGES\s*=\s*\[', html_content)
    if not stages_match:
        return []

    # Find each stage object by looking for { // comment or { title:
    # Split on stage boundaries
    stages_section = html_content[stages_match.end():]
    # Find the closing ];
    bracket_depth = 1
    end_pos = 0
    for i, ch in enumerate(stages_section):
        if ch == '[':
            bracket_depth += 1
        elif ch == ']':
            bracket_depth -= 1
            if bracket_depth == 0:
                end_pos = i
                break
    stages_section = stages_section[:end_pos]

    # Split into individual stage blocks
    # Each stage starts with { and we need to find matching }
    stages = []
    i = 0
    while i < len(stages_section):
        if stages_section[i] == '{':
            depth = 1
            j = i + 1
            while j < len(stages_section) and depth > 0:
                if stages_section[j] == '{':
                    depth += 1
                elif stages_section[j] == '}':
                    depth -= 1
                j += 1
            stages.append(stages_section[i:j])
            i = j
        else:
            i += 1

    return stages


def parse_stage(stage_text):
    """Parse a stage object into a dict of fields."""
    result = {}

    # Title
    m = re.search(r"title\s*:\s*'([^']*)'", stage_text)
    result['title'] = m.group(1) if m else 'Untitled'

    # Sub identifier
    m = re.search(r"sub\s*:\s*'([^']*)'", stage_text)
    result['sub'] = m.group(1) if m else ''

    # mainGraph
    m = re.search(r"mainGraph\s*:\s*'([^']*)'", stage_text)
    result['mainGraph'] = m.group(1) if m else None

    # graph (single or array)
    m = re.search(r"graph\s*:\s*\[([^\]]*)\]", stage_text)
    if m:
        result['graph'] = re.findall(r"'([^']*)'", m.group(1))
    else:
        m = re.search(r"graph\s*:\s*'([^']*)'", stage_text)
        result['graph'] = [m.group(1)] if m else []

    # mainHtml - extract all concatenated strings
    result['mainHtml'] = extract_concat_strings(stage_text, 'mainHtml')

    # narrative - extract all concatenated strings
    result['narrative'] = extract_concat_strings(stage_text, 'narrative')

    # equation
    result['equation'] = extract_concat_strings(stage_text, 'equation')

    return result


def extract_concat_strings(text, field):
    """Extract a field that may be composed of concatenated JS strings."""
    # Find the field assignment
    pattern = rf"{field}\s*:\s*"
    match = re.search(pattern, text)
    if not match:
        return ''

    pos = match.end()
    remaining = text[pos:]

    # Check for null
    if remaining.strip().startswith('null'):
        return ''

    # Collect all string fragments
    fragments = []
    i = 0
    while i < len(remaining):
        # Skip whitespace
        while i < len(remaining) and remaining[i] in ' \t\n\r':
            i += 1
        if i >= len(remaining):
            break

        if remaining[i] == "'":
            # Find the end of this string (handling escapes)
            j = i + 1
            while j < len(remaining):
                if remaining[j] == '\\':
                    j += 2
                    continue
                if remaining[j] == "'":
                    break
                j += 1
            fragments.append(remaining[i+1:j])
            i = j + 1
            # Skip whitespace and check for +
            while i < len(remaining) and remaining[i] in ' \t\n\r':
                i += 1
            if i < len(remaining) and remaining[i] == '+':
                i += 1
                continue
            else:
                break
        else:
            break

    return ''.join(fragments)


def format_slide(idx, stage, verbose=True):
    """Format a single slide for text output."""
    lines = []
    label = stage['sub'] if stage['sub'] else f'Intro {idx}'
    lines.append(f"## SLIDE {label}: {stage['title']}")

    if stage.get('mainGraph'):
        lines.append(f"\n[Figure: {stage['mainGraph']}]")

    if stage.get('mainHtml'):
        cleaned = strip_html(stage['mainHtml'])
        if cleaned:
            lines.append('')
            lines.append(textwrap.fill(cleaned, width=80))

    if verbose and stage.get('narrative'):
        cleaned = strip_html(stage['narrative'])
        if cleaned:
            lines.append(f"\n--- Narrative ---")
            lines.append(textwrap.fill(cleaned, width=80))

    if verbose and stage.get('equation'):
        eq = stage['equation']
        # Keep LaTeX mostly intact but strip HTML wrappers
        eq = re.sub(r'<[^>]+>', ' ', eq)
        eq = re.sub(r'\\u([0-9a-fA-F]{4})', lambda m: chr(int(m.group(1), 16)), eq)
        # Unescape JS string backslashes for clean LaTeX: \\frac -> \frac
        eq = eq.replace('\\\\', '\x00BS\x00')
        eq = eq.replace('\\n', '\n')
        eq = eq.replace("\\'", "'")
        eq = eq.replace('\x00BS\x00', '\\')
        eq = ascii_normalize(eq)
        eq = re.sub(r'\s+', ' ', eq).strip()
        if eq:
            lines.append(f"\n--- Equation ---\n{eq}")

    if verbose and stage.get('graph'):
        figs = [g for g in stage['graph'] if g]
        if figs:
            lines.append(f"\n[Additional figures: {', '.join(figs)}]")

    lines.append('')
    return '\n'.join(lines)


def generate_brief(stages):
    """Generate llms.txt — brief summary."""
    lines = []
    lines.append("# Extensible Longevity Model (ELM)")
    lines.append("# https://silmonbiggs.github.io/extensible-longevity-model/")
    lines.append("#")
    lines.append("# This file provides an LLM-readable summary of the interactive slide deck.")
    lines.append("# For complete slide-by-slide content, fetch llms-full.txt from the same URL.")
    lines.append("")
    lines.append("## What This Is")
    lines.append("")
    lines.append("A mechanistic ODE model of aging interventions, calibrated to the NIA")
    lines.append("Interventions Testing Program (ITP) data. The model tracks 5 state")
    lines.append("variables governed by ODEs:")
    lines.append("")
    lines.append("  1. Heteroplasmy (H) -- mitochondrial DNA damage fraction")
    lines.append("  2. NAD+ (N) -- nicotinamide adenine dinucleotide level")
    lines.append("  3. DNA damage (D) -- genomic lesion accumulation")
    lines.append("  4. Senescent cell burden (S)")
    lines.append("  5. Methylation drift (M) -- epigenetic age")
    lines.append("")
    lines.append("Plus 10 signaling intermediates (AMPK, TSC2, mTORC1, SIRT1, FOXO3,")
    lines.append("Autophagy, SASP, ROS, EZH2, PGC-1alpha) computed at each time step")
    lines.append("without their own differential equations.")
    lines.append("")
    lines.append("These are NOT independent -- they form a coupled ODE system with")
    lines.append("cross-pathway feedback (e.g., NAD+ -> sirtuins -> FOXO3 -> heteroplasmy).")
    lines.append("The model predicts combination effects, including interactions that")
    lines.append("single-compound testing cannot reveal.")
    lines.append("")
    lines.append("## Scope of Evaluation")
    lines.append("")
    lines.append("The active model -- the part that makes predictions and should be")
    lines.append("evaluated -- consists of:")
    lines.append("")
    lines.append("  - 5 state variables + 10 signaling intermediates (listed above)")
    lines.append("  - 6 ITP compound calibrations (12 sex-specific targets)")
    lines.append("  - 4 sex-difference mechanisms (AMPK saturation, PK, testosterone, COX-2 interference)")
    lines.append("  - Dose-response Hill curves (3 fitted, 3 estimated)")
    lines.append("  - BioAge composite -> death threshold")
    lines.append("  - Combination prediction engine: compound pathway injections are")
    lines.append("    summed at the INPUT level, but the ODE dynamics produce nonlinear")
    lines.append("    interactions through Hill-function saturation, Michaelis-Menten")
    lines.append("    substrate dependence, and cross-pathway feedback. This is how the")
    lines.append("    model predicts super- and sub-additive combinations from single-")
    lines.append("    compound calibration alone.")
    lines.append("")
    lines.append("The codebase also contains two modules that are INACTIVE in all current")
    lines.append("predictions and should not factor into evaluation:")
    lines.append("")
    lines.append("  - Cancer module (elm/model.py: stochastic mutation accumulation).")
    lines.append("    Present for structural completeness. Cancer probability is computed")
    lines.append("    but does not influence BioAge or lifespan in any reported result.")
    lines.append("    It is scaffolding for future calibration when cancer-specific")
    lines.append("    endpoint data become available.")
    lines.append("")
    lines.append("  - Cross-species translation (elm/uncertainty.py: human scaling).")
    lines.append("    All current predictions are in mouse units (normalized lifespan).")
    lines.append("    No mouse-to-human translation is applied or claimed. The module")
    lines.append("    exists as a placeholder for future work; its parameters are")
    lines.append("    literature estimates with no anchor in the current calibration.")
    lines.append("")
    lines.append("Reviewers who evaluate these modules as if they contribute to the")
    lines.append("model's current claims will be evaluating code that the authors")
    lines.append("explicitly flag as not yet in service.")
    lines.append("")
    lines.append("## Why a Quantitative Model")
    lines.append("")
    lines.append("The biology of aging has accumulated decades of qualitative pathway")
    lines.append("diagrams: arrows from protein to protein, marked \"activates\" or")
    lines.append("\"inhibits,\" with no magnitudes attached. This represents the beginning")
    lines.append("of knowledge but not, in Lord Kelvin's sense, knowledge of a")
    lines.append("satisfactory kind. You cannot predict whether rapamycin + acarbose")
    lines.append("will be additive, sub-additive, or superadditive by staring at arrows.")
    lines.append("You need rate constants, saturation curves, and a system of equations")
    lines.append("that can be solved.")
    lines.append("")
    lines.append("This model is offered as a productive starting point, not a ground")
    lines.append("truth. The key properties of a starting point:")
    lines.append("")
    lines.append("  1. It must be quantitative -- otherwise it cannot make falsifiable")
    lines.append("     predictions about combinations.")
    lines.append("  2. It must be calibrated to existing data -- the ITP provides 12")
    lines.append("     sex-specific lifespan endpoints from controlled experiments.")
    lines.append("  3. It must be honest about what the data constrain -- one-at-a-time")
    lines.append("     perturbation of 141 frozen parameters (+/-10%) shifts no prediction")
    lines.append("     by more than 0.65 pp. Results are robust to parameter choices.")
    lines.append("  4. It must be extensible -- new compounds, tissue-specific endpoints,")
    lines.append("     and combination data can recalibrate and extend the model without")
    lines.append("     starting over.")
    lines.append("")
    lines.append("The 168 parameters are not all free. 6 compound scale factors are")
    lines.append("calibrated to male ITP data; the remaining 162 are frozen priors")
    lines.append("(115 literature, 18 set within plausible ranges, 6 calibrated to")
    lines.append("control biology, 23 structural). One-at-a-time perturbation of")
    lines.append("141 frozen parameters (+/-10%) shifts no prediction by more than")
    lines.append("0.65 pp (Slide A9). Could alternative parameter sets in distant")
    lines.append("basins fit the same 12 targets and make different combination")
    lines.append("predictions? In principle, yes. In practice, the Rapa + Acarbose")
    lines.append("out-of-sample prediction (+34% vs +34% observed) provides evidence")
    lines.append("that the current basin is not arbitrary.")
    lines.append("")
    lines.append("## Key Results")
    lines.append("")
    lines.append("- Calibration: 12 ITP targets, +/-0.08% mean error, 6 free parameters (one scale factor per compound)")
    lines.append("- Out-of-sample validation: Rapa+Acarb predicted +34%, observed +34% (not fit)")
    lines.append("- Testable prediction: Rapa+NMN superadditive (+31.8% M / +36.0% F)")
    lines.append("- Serial reliability framework explains why single interventions yield +10-26%")
    lines.append("")
    lines.append("## Individual Variation")
    lines.append("")
    lines.append("The model predicts median lifespan extension. Individual mice")
    lines.append("vary -- ITP reports ~12% CV in lifespan across genetically")
    lines.append("identical animals at multiple sites. A Monte Carlo resampling")
    lines.append("(1,000 log-normal draws with that CV) gives the 5th-95th")
    lines.append("percentile range of individual outcomes around the median:")
    lines.append("")
    lines.append("  - Rapamycin male: median +23%, individual range 0-49%")
    lines.append("  - Rapa+Acarb male: median +34%, individual range 8-61%")
    lines.append("  - Rapa+NMN male: median +32%, individual range 7-60%")
    lines.append("")
    lines.append("This is NOT uncertainty about whether the model is right -- it")
    lines.append("is the expected spread of individual lifespans if it is. Parameter")
    lines.append("sensitivity is a separate question, addressed by the tornado")
    lines.append("charts (Slides A5, A5b). Reproducible: run_monte_carlo(seed=42).")
    lines.append("")
    lines.append("## Compounds Modeled (ITP-validated)")
    lines.append("")
    lines.append("- Rapamycin (mTORC1 inhibitor): +23% M, +26% F")
    lines.append("- Acarbose (α-glucosidase inhibitor): +22% M, +5% F")
    lines.append("- 17α-Estradiol (non-feminizing estrogen): +19% M, 0% F")
    lines.append("- Canagliflozin (SGLT2 inhibitor): +14% M, +9% F")
    lines.append("- Aspirin (COX inhibitor): +8% M, 0% F")
    lines.append("- Glycine (amino acid supplement): +6% M, +4% F")
    lines.append("")
    lines.append("## Slide Deck Structure")
    lines.append("")
    lines.append(f"Total slides: {len(stages)}")
    lines.append("")

    # Group slides by the 10-section navigation structure
    SECTIONS = [
        ('Introduction',             0,  2),
        ('Architecture',             2, 11),
        ('Interventions',           11, 17),
        ('Calibration & Prediction',17, 19),
        ('Sex Differences',         19, 26),
        ('Limitations',             26, 27),
        ('Validation',              27, 29),
        ('Sensitivity',             29, 34),
        ('Specification',           34, len(stages)),
    ]

    for section_name, start, end in SECTIONS:
        entries = []
        for i in range(start, min(end, len(stages))):
            s = stages[i]
            sub = s['sub']
            title = s['title']
            entries.append(f"  [{sub or f'{i}'}] {title}")
        if entries:
            lines.append(f"### {section_name}")
            lines.extend(entries)
            lines.append("")

    lines.append("## Sensitivity Analysis")
    lines.append("")
    lines.append("The deck includes extensive sensitivity and robustness analysis:")
    lines.append("- Slide A3: Leave-one-out vulnerability analysis")
    lines.append("- Slide A6: BioAge weight (w_meth) sweep with clock PCA constraint")
    lines.append("- Slide A6b: heteroplasmy and senescence weight sweeps")
    lines.append("- Slide A7: Feedback loop analysis")
    lines.append("- Slide A8: Hallmark importance ranking")
    lines.append("- Slide A9: pathway intervention-point sensitivity ranking")
    lines.append("")
    lines.append("## Source Code and Reproducibility")
    lines.append("")
    lines.append("GitHub: https://github.com/silmonbiggs/extensible-longevity-model")
    lines.append("Model package: elm/ (Python, ODE engine)")
    lines.append("Dependencies: requirements.txt (numpy, matplotlib)")
    lines.append("Smoke test: python tests/smoke_test.py (all 12 ITP targets within +/-1%)")
    lines.append("Figure generation: scripts/generate_figures.py")
    lines.append("Slide deck source: docs/ITPSexDiff.html")
    lines.append("")
    lines.append("## For Complete Content")
    lines.append("")
    lines.append("Fetch: https://silmonbiggs.github.io/extensible-longevity-model/llms-full.txt")
    lines.append("This contains the full text of every slide including narratives and equations.")
    lines.append("")

    return '\n'.join(lines)


def generate_full(stages):
    """Generate llms-full.txt — complete slide-by-slide content."""
    lines = []
    lines.append("# Extensible Longevity Model (ELM) -- Complete Slide Content")
    lines.append("# https://silmonbiggs.github.io/extensible-longevity-model/")
    lines.append("#")
    lines.append("# This file contains the full text content of every slide in the")
    lines.append("# interactive presentation, extracted for LLM accessibility.")
    lines.append("# The interactive version includes Cytoscape network diagrams,")
    lines.append("# clickable nodes, and dynamically rendered figures.")
    lines.append("#")
    lines.append(f"# Total slides: {len(stages)}")
    lines.append(f"# Generated from: docs/ITPSexDiff.html")
    lines.append("")
    lines.append("## Scope of Evaluation")
    lines.append("")
    lines.append("The active model -- the part that makes predictions and should be")
    lines.append("evaluated -- consists of:")
    lines.append("")
    lines.append("  - 5 state variables (heteroplasmy, NAD+, DNA damage,")
    lines.append("    senescence, methylation) + 10 signaling intermediates")
    lines.append("  - 6 ITP compound calibrations (12 sex-specific targets)")
    lines.append("  - 4 sex-difference mechanisms (AMPK saturation, PK, testosterone, COX-2 interference)")
    lines.append("  - Dose-response Hill curves (3 fitted, 3 estimated)")
    lines.append("  - BioAge composite -> death threshold")
    lines.append("  - Combination prediction engine: compound pathway injections are")
    lines.append("    summed at the INPUT level, but the ODE dynamics produce nonlinear")
    lines.append("    interactions through Hill-function saturation, Michaelis-Menten")
    lines.append("    substrate dependence, and cross-pathway feedback. This is how the")
    lines.append("    model predicts super- and sub-additive combinations from single-")
    lines.append("    compound calibration alone.")
    lines.append("")
    lines.append("The codebase also contains two modules that are INACTIVE in all current")
    lines.append("predictions and should not factor into evaluation:")
    lines.append("")
    lines.append("  - Cancer module (elm/model.py: stochastic mutation accumulation).")
    lines.append("    Cancer_prob is computed but NOT fed to BioAge; death is determined")
    lines.append("    solely by BioAge >= death_threshold. The module is scaffolding for")
    lines.append("    future calibration when cancer-specific endpoint data become available.")
    lines.append("")
    lines.append("  - Cross-species translation (elm/uncertainty.py: human scaling).")
    lines.append("    All current predictions are in mouse units (normalized lifespan).")
    lines.append("    No mouse-to-human translation is applied or claimed.")
    lines.append("")
    lines.append("Parameter census: 168 total parameters. 6 free (compound scale")
    lines.append("factors calibrated to male ITP data). 162 frozen priors: 115")
    lines.append("literature, 18 set within plausible ranges, 6 calibrated to control")
    lines.append("biology, 23 structural defaults.")
    lines.append("")
    lines.append("## Why a Quantitative Model")
    lines.append("")
    lines.append("The biology of aging has accumulated decades of qualitative pathway")
    lines.append("diagrams: arrows from protein to protein, marked \"activates\" or")
    lines.append("\"inhibits,\" with no magnitudes attached. You cannot predict whether")
    lines.append("rapamycin + acarbose will be additive, sub-additive, or superadditive")
    lines.append("by staring at arrows. You need rate constants, saturation curves, and")
    lines.append("a system of equations that can be solved.")
    lines.append("")
    lines.append("This model is offered as a productive starting point, not a ground")
    lines.append("truth. It is designed to anneal toward agreement with biology through")
    lines.append("an explicit experiment -> calibrate -> extend cycle.")
    lines.append("")

    for i, stage in enumerate(stages):
        lines.append(format_slide(i, stage, verbose=True))

    # Append supplementary analyses
    clock_pca_path = Path(os.path.dirname(__file__)) / 'clock_pca' / 'clock_methylation_fraction.txt'
    if clock_pca_path.exists():
        lines.append("## SUPPLEMENTARY: Methylation Fraction of Biological Age")
        lines.append("(Referenced by Slide A6: BioAge Weight Sensitivity)")
        lines.append("")
        lines.append(ascii_normalize(clock_pca_path.read_text(encoding='utf-8')))

    return '\n'.join(lines)


def main():
    repo_root = DOCS_DIR / '..'

    html = HTML_FILE.read_text(encoding='utf-8')
    raw_stages = extract_stages(html)
    print(f"Extracted {len(raw_stages)} stages from HTML")

    stages = [parse_stage(s) for s in raw_stages]

    # Write brief (ASCII-normalized) to docs/ and repo root
    brief = ascii_normalize(generate_brief(stages))
    for dest in [DOCS_DIR / 'llms.txt', repo_root / 'llms.txt']:
        dest.write_text(brief, encoding='utf-8')
        print(f"Written: {dest.resolve()} ({len(brief)} bytes)")

    # Write full to docs/ and repo root
    full = generate_full(stages)
    for dest in [DOCS_DIR / 'llms-full.txt', repo_root / 'llms-full.txt']:
        dest.write_text(full, encoding='utf-8')
        print(f"Written: {dest.resolve()} ({len(full)} bytes)")


if __name__ == '__main__':
    main()
