#!/usr/bin/env python
"""Build the Paper 1 (ITP Sex Differences) slide deck from working_slides.html.

Reads docs/working_slides.html, keeps only stages with paper1:true,
injects a Paper 1 intro slide, and writes docs/ITPSexDiff.html.

Usage:
    python scripts/build_paper1_deck.py
"""

from pathlib import Path
import re

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "docs" / "working_slides.html"
DST = REPO / "docs" / "ITPSexDiff.html"

TITLE = "ELM: ITP Sex Differences (r1)"

INTRO_STAGE = (
    "{ // Paper 1 intro\n"
    "    title:'Paper 1: ITP Sex Differences', sub:'intro',\n"
    "    paper1:true,\n"
    "    nodes:[], edges:[],\n"
    "    mainHtml:'<div class=\"splash-content\">' +\n"
    "        '<h2>Four Biological Mechanisms Predict ITP Sex Differences from Male Data Alone</h2>' +\n"
    "        '<p class=\"subtitle\">We calibrate on male ITP data (6 compounds, 6 free parameters). Four sex-specific mechanisms from prior literature predict female outcomes with 1.3 pp mean error \\u2014 no female fitting.</p>' +\n"
    "        '<table class=\"splash-table\">' +\n"
    "        '<tr><th>Mechanism</th><th>What it does</th><th>Compounds affected</th></tr>' +\n"
    "        '<tr><td>AMPK saturation</td><td>Female AMPK already near ceiling</td><td>Acarbose, Canagliflozin, Glycine</td></tr>' +\n"
    "        '<tr><td>CYP3A pharmacokinetics</td><td>Faster drug clearance in females</td><td>Rapamycin</td></tr>' +\n"
    "        '<tr><td>Testosterone dependency</td><td>Male-only activation pathway</td><td>17\\u03B1-Estradiol</td></tr>' +\n"
    "        '<tr><td>COX-2 interference</td><td>Prostaglandin cross-talk in females</td><td>Aspirin</td></tr>' +\n"
    "        '</table>' +\n"
    "        '<p style=\"font-size:12px;color:#888;margin-top:12px\">Use \\u2192 arrow key or Next button to continue</p>' +\n"
    "        '</div>',\n"
    "    narrative:'',\n"
    "    equation:'',\n"
    "    graph:null,\n"
    "}"
)

# Paper 1 section nav: (label, sub-prefix to match)
# Index 0 is always Introduction (the injected intro stage).
SECTIONS = [
    ("Introduction",               "intro"),
    ("Architecture",               "0"),
    ("Interventions",              "8."),
    ("Calibration &amp; Prediction", "9"),
    ("Sex Differences",            "11"),
    ("Limitations",                "18"),
    ("Validation",                 "A3"),
    ("Sensitivity",                "A6"),
    ("Specification",              "A16"),
]


def extract_stages_span(html: str):
    """Return (start, end) character offsets of the STAGES array contents.

    start points to the first char after 'const STAGES = ['
    end   points to the ']' that closes the array.

    Uses a simple search for the standalone '];' line that ends the array,
    rather than bracket-depth tracking (which can be tripped by \\u escapes
    inside JS string literals).
    """
    marker = "const STAGES = ["
    start = html.index(marker) + len(marker)
    # The STAGES array ends with '];\n' on its own line.
    # Search for '\n];' after the array start.
    end = html.index("\n];", start) + 1  # +1 to point at the ']'
    return start, end


def split_stage_objects(array_body: str):
    """Split the text between [ and ] into individual stage object strings.

    Each top-level { ... } at brace-depth 0 is one stage.
    """
    stages = []
    depth = 0
    obj_start = None
    i = 0
    while i < len(array_body):
        ch = array_body[i]
        if ch == "{" and depth == 0:
            obj_start = i
            depth = 1
        elif ch == "{":
            depth += 1
        elif ch == "}" and depth == 1:
            depth = 0
            stages.append(array_body[obj_start : i + 1])
            obj_start = None
        elif ch == "}":
            depth -= 1
        # Skip string literals
        elif ch in ("'", '"', '`') and depth > 0:
            quote = ch
            i += 1
            while i < len(array_body) and array_body[i] != quote:
                if array_body[i] == "\\":
                    i += 1
                i += 1
        i += 1
    return stages


def has_paper1(stage_text: str) -> bool:
    """Check whether a stage object string contains paper1:true."""
    return "paper1:true" in stage_text or "paper1: true" in stage_text


def get_sub(stage_text: str) -> str:
    """Extract the sub:'...' value from a stage object string."""
    m = re.search(r"sub:\s*'([^']*)'", stage_text)
    return m.group(1) if m else ""


def build_section_nav(stages_text_list: list[str]) -> str:
    """Build the <select> dropdown HTML for Paper 1 sections."""
    lines = []
    lines.append(
        '        <select id="section-select" onchange="window.jumpToSection(this.value)">'
    )
    for label, sub_prefix in SECTIONS:
        if sub_prefix == "intro":
            idx = 0
        else:
            idx = None
            for j, stg in enumerate(stages_text_list):
                s = get_sub(stg)
                if sub_prefix.endswith("."):
                    # Match sub values that start with this prefix
                    if s.startswith(sub_prefix):
                        idx = j
                        break
                else:
                    if s == sub_prefix:
                        idx = j
                        break
            if idx is None:
                print(f"  WARNING: no stage found for section '{label}' (sub prefix '{sub_prefix}')")
                continue
        lines.append(f'            <option value="{idx}">{label}</option>')
    lines.append("        </select>")
    return "\n".join(lines)


def main():
    print(f"Reading {SRC}")
    html = SRC.read_text(encoding="utf-8")

    # --- 1. Extract and filter stages ---
    arr_start, arr_end = extract_stages_span(html)
    array_body = html[arr_start:arr_end]
    all_stages = split_stage_objects(array_body)
    print(f"  Total stages in source: {len(all_stages)}")

    paper1_stages = [s for s in all_stages if has_paper1(s)]
    print(f"  Stages with paper1:true: {len(paper1_stages)}")

    # Inject intro at position 0
    paper1_stages.insert(0, INTRO_STAGE)
    print(f"  Stages after injecting intro: {len(paper1_stages)}")

    # --- 2. Rebuild STAGES array ---
    new_array_body = "\n" + ",\n".join(paper1_stages) + ",\n"
    html = html[:arr_start] + new_array_body + html[arr_end:]

    # --- 3. Update <title> ---
    html = re.sub(r"<title>[^<]*</title>", f"<title>{TITLE}</title>", html)

    # --- 4. Update <h1> ---
    html = re.sub(r"(<h1>)[^<]*(</h1>)", rf"\g<1>{TITLE}\g<2>", html)

    # --- 5. Rebuild section nav ---
    new_nav = build_section_nav(paper1_stages)
    html = re.sub(
        r'<select id="section-select"[^>]*>.*?</select>',
        new_nav,
        html,
        flags=re.DOTALL,
    )

    # --- 6. Write output ---
    DST.write_text(html, encoding="utf-8")
    print(f"\nWrote {DST}")
    print(f"  Total slides: {len(paper1_stages)}")


if __name__ == "__main__":
    main()
