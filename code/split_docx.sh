#!/usr/bin/env bash
#
# code/split_docx.sh
#
# Split the rendered Word manuscript into the two files Genome Medicine wants
# as separate uploads:
#   - main_manuscript.docx   (title page through the Additional files subsection)
#   - additional_file_1.docx (the Supplementary Materials section only)
#
# This is a thin wrapper around code/split_docx.py, which performs a *surgical*
# XML split of the valid Quarto-generated index.docx: it trims the document
# body at the "Supplementary Materials" heading and keeps every other part of
# the package (styles, numbering, media, relationships) exactly as Quarto wrote
# them. It also balances bookmarks across the cut and renumbers supplementary
# floats to "S" form (Fig. S1, Table S1, ...) to match the PDF.
#
# An earlier version round-tripped through pandoc markdown
# (docx -> md -> docx); that re-created the figure tables and produced files
# Word flagged as "unreadable content". The surgical split avoids regenerating
# any OOXML, so the output is exactly as valid as Quarto's index.docx.
#
# Why split post-render (rather than two Quarto documents)? The single
# index.qmd keeps Quarto's @fig-/@tbl- cross-references resolving and avoids
# duplicating the ~50 R setup chunks the supplement depends on.
#
# Dependencies: python3 with lxml (`pip install lxml`).
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

if [[ ! -f "$DOCX_IN" ]]; then
  echo "ERROR: $DOCX_IN does not exist." >&2
  echo "Run \`quarto render index.qmd --to docx\` from manuscript/ first." >&2
  exit 1
fi

if ! python3 -c "import lxml" >/dev/null 2>&1; then
  echo "ERROR: python3 with the 'lxml' package is required (pip install lxml)." >&2
  exit 1
fi

python3 "$SCRIPT_DIR/split_docx.py" "$DOCX_IN" "$DOCX_MAIN" "$DOCX_SUPP"

echo
echo "Word submission package (in manuscript/_manuscript/):"
echo "  main_manuscript.docx   - main article only"
echo "  additional_file_1.docx - supplementary materials"
