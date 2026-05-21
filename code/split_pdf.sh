#!/usr/bin/env bash
#
# code/split_pdf.sh
#
# Split the rendered manuscript PDF into:
#   - main_manuscript.pdf   (title page through the Additional files subsection)
#   - additional_file_1.pdf (the Supplementary Materials section only)
#
# This script is intended to be run ONCE at submission time, after the
# manuscript has been rendered (manuscript/_manuscript/index.pdf).
# It is NOT part of the normal Quarto render — we deliberately keep all
# content in a single index.qmd so that Quarto's cross-references
# (`@fig-X`, `@tbl-X`) resolve naturally during rendering. The split happens
# only on the final PDF, leaving the source as a single document.
#
# Why post-render split (rather than two Quarto notebooks)?
#   * Quarto's `@fig-X` cross-references do NOT auto-resolve across
#     separate documents when targeting PDF output (only HTML).
#   * A multi-notebook setup would require duplicating ~50 R setup chunks
#     to recreate the analysis state in the supplementary notebook.
#   * The in-text "Additional file 1: Fig. S1" labels are plain text and
#     match the figure titles in the extracted Additional file 1 PDF, so
#     readers navigate by label rather than by PDF hyperlink.
#
# Dependencies (install once):
#   brew install poppler                  # macOS (provides pdftotext, pdfseparate, pdfunite)
#   apt install poppler-utils             # Debian / Ubuntu
#   choco install xpdf-utils              # Windows
# Optionally also install qpdf:
#   brew install qpdf                     # qpdf is faster for single-step splits
#
# Usage:
#   ./code/split_pdf.sh
#
# After running, the two output PDFs live next to the source:
#   manuscript/_manuscript/main_manuscript.pdf
#   manuscript/_manuscript/additional_file_1.pdf
#
# Together with manuscript/_manuscript/additional_file_2.xlsx (generated
# automatically by the Quarto render at render time), these form the full
# Genome Medicine submission package.

set -euo pipefail

# Resolve project root (one level above this script's directory)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

PDF_IN="$PROJECT_ROOT/manuscript/_manuscript/index.pdf"
PDF_MAIN="$PROJECT_ROOT/manuscript/_manuscript/main_manuscript.pdf"
PDF_SUPP="$PROJECT_ROOT/manuscript/_manuscript/additional_file_1.pdf"

# --- Pre-flight checks ------------------------------------------------------

if [[ ! -f "$PDF_IN" ]]; then
  echo "ERROR: $PDF_IN does not exist." >&2
  echo "Run \`quarto render index.qmd --to pdf\` from manuscript/ first." >&2
  exit 1
fi

for cmd in pdftotext; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: required dependency '$cmd' not found in PATH." >&2
    echo "Install on macOS via:    brew install poppler" >&2
    echo "Install on Debian via:   apt install poppler-utils" >&2
    exit 1
  fi
done
# Pick a split backend: prefer qpdf (single-call), fall back to pdfseparate+pdfunite (poppler).
if command -v qpdf >/dev/null 2>&1; then
  BACKEND=qpdf
elif command -v pdfseparate >/dev/null 2>&1 && command -v pdfunite >/dev/null 2>&1; then
  BACKEND=poppler
else
  echo "ERROR: need either qpdf, or pdfseparate+pdfunite, on PATH." >&2
  echo "Install with: brew install poppler   (or:  brew install qpdf)" >&2
  exit 1
fi

# --- Find the page where 'Supplementary Materials' starts -------------------
#
# We extract the PDF to plain text with -layout (to preserve heading position)
# and walk page by page until we hit the "Supplementary Materials" heading.
# pdftotext page separator is the form-feed character (\f).

# Run the pipeline with pipefail temporarily off, because awk's early-exit
# closes the pipe to pdftotext (SIGPIPE) and pipefail would otherwise abort
# the script even when we got the line we wanted.
set +o pipefail
SUPP_START_PAGE=$(pdftotext -layout "$PDF_IN" - | awk '
  BEGIN { page = 1 }
  /^\x0c/ { page++ }
  /Supplementary Materials/ {
    # First line where the heading appears on its own; skip any in-body
    # references that just mention "supplementary materials" in prose.
    # Heading lines tend to have leading whitespace then the bold text.
    if (match($0, /^[[:space:]]*[0-9]*[[:space:]]+Supplementary Materials[[:space:]]*$/)) {
      print page
      exit
    }
  }
')
set -o pipefail

if [[ -z "${SUPP_START_PAGE:-}" ]]; then
  echo "ERROR: could not find 'Supplementary Materials' heading in the PDF." >&2
  echo "Inspect the PDF and update this script's heuristic." >&2
  exit 1
fi

# Total page count (pdfinfo from poppler if available; otherwise qpdf)
if command -v pdfinfo >/dev/null 2>&1; then
  TOTAL_PAGES=$(pdfinfo "$PDF_IN" | awk '/^Pages:/ {print $2}')
else
  TOTAL_PAGES=$(qpdf --show-npages "$PDF_IN")
fi
MAIN_LAST_PAGE=$((SUPP_START_PAGE - 1))

if (( MAIN_LAST_PAGE < 1 || SUPP_START_PAGE > TOTAL_PAGES )); then
  echo "ERROR: nonsensical page split (supp starts at $SUPP_START_PAGE / $TOTAL_PAGES total)." >&2
  exit 1
fi

echo "PDF has $TOTAL_PAGES pages total. (split backend: $BACKEND)"
echo "Main manuscript    : pages 1 through $MAIN_LAST_PAGE"
echo "Additional file 1  : pages $SUPP_START_PAGE through $TOTAL_PAGES"

# --- Perform the split ------------------------------------------------------

if [[ "$BACKEND" == "qpdf" ]]; then
  qpdf --empty --pages "$PDF_IN" 1-$MAIN_LAST_PAGE -- "$PDF_MAIN"
  qpdf --empty --pages "$PDF_IN" $SUPP_START_PAGE-$TOTAL_PAGES -- "$PDF_SUPP"
else
  # poppler backend: pdfseparate writes one PDF per page; pdfunite then merges.
  # NOTE: pdfunite is strictly positional — files are concatenated in the
  # order given on the command line. We must build the file list in
  # numerical page order; alphabetic order (e.g. macOS BSD `ls -v` which
  # ignores `-v`, or default `ls`) puts page 10 before page 2 and produces
  # a scrambled output.
  TMPDIR_MAIN=$(mktemp -d)
  TMPDIR_SUPP=$(mktemp -d)
  trap 'rm -rf "$TMPDIR_MAIN" "$TMPDIR_SUPP"' EXIT
  pdfseparate -f 1 -l $MAIN_LAST_PAGE "$PDF_IN" "$TMPDIR_MAIN/p-%d.pdf"
  pdfseparate -f $SUPP_START_PAGE -l $TOTAL_PAGES "$PDF_IN" "$TMPDIR_SUPP/p-%d.pdf"

  main_files=()
  for i in $(seq 1 $MAIN_LAST_PAGE); do
    main_files+=("$TMPDIR_MAIN/p-$i.pdf")
  done
  pdfunite "${main_files[@]}" "$PDF_MAIN"

  supp_files=()
  for i in $(seq $SUPP_START_PAGE $TOTAL_PAGES); do
    supp_files+=("$TMPDIR_SUPP/p-$i.pdf")
  done
  pdfunite "${supp_files[@]}" "$PDF_SUPP"
fi

echo
echo "Wrote $PDF_MAIN"
echo "Wrote $PDF_SUPP"
echo
echo "Submission package (in manuscript/_manuscript/):"
echo "  main_manuscript.pdf    - main article only"
echo "  additional_file_1.pdf  - supplementary materials (figures + tables in PDF form)"
echo "  additional_file_2.xlsx - supplementary tables in machine-readable XLSX (generated by Quarto render)"
