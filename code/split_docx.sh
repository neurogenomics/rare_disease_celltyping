#!/usr/bin/env bash
#
# code/split_docx.sh
#
# Split the rendered Word manuscript into the two files Genome Medicine wants
# as separate uploads:
#   - main_manuscript.docx   (title page through the Additional files subsection)
#   - additional_file_1.docx (the Supplementary Materials section only)
#
# Unlike split_pdf.sh (which splits the rendered PDF by page number), Word
# documents have no fixed pages, so we split at the "Supplementary Materials"
# heading instead. We do this by round-tripping through pandoc markdown:
#
#   index.docx --(pandoc)--> markdown (+ extracted media)
#              --(split at the "Supplementary Materials" heading)-->
#   main.md / supp.md --(pandoc, --reference-doc=index.docx)--> two .docx
#
# Reusing the rendered index.docx as the --reference-doc carries the same
# Word styles (fonts, heading styles, table styles) into both halves so they
# match the combined document.
#
# Why post-render split (rather than two Quarto documents)? Same reason as
# split_pdf.sh: the single index.qmd keeps Quarto's @fig-/@tbl- cross-references
# resolving and avoids duplicating the ~50 R setup chunks that the supplementary
# section depends on.
#
# Dependencies: pandoc (bundled with Quarto, or `brew install pandoc`).
#
# Usage (after `quarto render index.qmd --to docx`):
#   ./code/split_docx.sh

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
MAN_DIR="$PROJECT_ROOT/manuscript/_manuscript"

DOCX_IN="$MAN_DIR/index.docx"
DOCX_MAIN="$MAN_DIR/main_manuscript.docx"
DOCX_SUPP="$MAN_DIR/additional_file_1.docx"

# --- Pre-flight checks ------------------------------------------------------

if [[ ! -f "$DOCX_IN" ]]; then
  echo "ERROR: $DOCX_IN does not exist." >&2
  echo "Run \`quarto render index.qmd --to docx\` from manuscript/ first." >&2
  exit 1
fi

if ! command -v pandoc >/dev/null 2>&1; then
  echo "ERROR: required dependency 'pandoc' not found in PATH." >&2
  echo "Install on macOS via:  brew install pandoc   (or use Quarto's bundled pandoc)." >&2
  exit 1
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --- 1) docx -> markdown (extract embedded images) --------------------------
# Run pandoc from $WORK so the extracted-media paths are relative and resolve
# the same way on the way back out to docx.
(
  cd "$WORK"

  pandoc "$DOCX_IN" -f docx -t markdown --wrap=none --extract-media=media -o all.md

  # --- 2) locate the "Supplementary Materials" heading (any heading level) --
  SPLIT_LINE=$(grep -nE '^#{1,6}[[:space:]]+.*Supplementary Materials' all.md | head -1 | cut -d: -f1 || true)
  if [[ -z "${SPLIT_LINE:-}" ]]; then
    echo "ERROR: could not find a 'Supplementary Materials' heading in the extracted markdown." >&2
    echo "Inspect $WORK/all.md and update this script's heuristic." >&2
    exit 1
  fi

  TOTAL_LINES=$(wc -l < all.md)
  echo "Markdown has $TOTAL_LINES lines; 'Supplementary Materials' heading at line $SPLIT_LINE."

  # --- 3) split: main = before the heading; supp = heading onward -----------
  head -n "$((SPLIT_LINE - 1))" all.md > main.md
  tail -n "+$SPLIT_LINE"        all.md > supp.md

  # --- 4) markdown -> docx, reusing index.docx as the style reference -------
  pandoc main.md -f markdown -t docx --reference-doc="$DOCX_IN" -o "$DOCX_MAIN"
  pandoc supp.md -f markdown -t docx --reference-doc="$DOCX_IN" -o "$DOCX_SUPP"
)

echo
echo "Wrote $DOCX_MAIN"
echo "Wrote $DOCX_SUPP"
echo
echo "Word submission package (in manuscript/_manuscript/):"
echo "  main_manuscript.docx   - main article only"
echo "  additional_file_1.docx - supplementary materials"
