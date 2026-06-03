#!/usr/bin/env bash
#
# code/export_figures.sh
#
# Assemble a per-figure PDF for every manuscript figure into a single folder,
# named to match the manuscript labels so the editorial office can identify
# each at a glance:
#
#   manuscript/_manuscript/figures/
#     Figure_1.pdf  ... Figure_9.pdf      (main figures)
#     Figure_S1.pdf ... Figure_S16.pdf    (supplementary figures)
#     figures_manifest.csv                (filename -> label -> short title)
#
# Most figures are produced by code chunks and are already vector PDFs in the
# Quarto render output (manuscript/_manuscript/_tex/index_files/figure-pdf/).
# A handful are static images (img/*.png); those are wrapped into single-page
# PDFs with ImageMagick.
#
# Run AFTER `quarto render index.qmd --to pdf` (so the figure-pdf cache is
# current). Dependencies: ImageMagick (`magick`).
#
# Usage:  ./code/export_figures.sh

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
MAN="$ROOT/manuscript"
PDFDIR="$MAN/_manuscript/_tex/index_files/figure-pdf"
IMG="$MAN/img"
# Output folder; override by passing a path as the first argument.
OUT="${1:-$HOME/Downloads/manuscript_figures}"

command -v magick >/dev/null 2>&1 || { echo "ERROR: ImageMagick (magick) not found." >&2; exit 1; }
[[ -d "$PDFDIR" ]] || { echo "ERROR: $PDFDIR missing — run 'quarto render index.qmd --to pdf' first." >&2; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"
MAN_CSV="$OUT/figures_manifest.csv"
echo "file,label,title" > "$MAN_CSV"

# emit: copy a render-cache PDF (P:stem) or wrap a static image (I:path) into
# Figure_<label>.pdf, and append a manifest row. Order of the emit calls below
# matches the figures' order of appearance, i.e. the manuscript numbering.
emit () {  # $1=label  $2=kind:src  $3=title
  local label="$1" spec="$2" title="$3"
  local kind="${spec%%:*}" src="${spec#*:}"
  local dest="$OUT/Figure_${label}.pdf"
  if [[ "$kind" == "P" ]]; then
    cp "$PDFDIR/${src}.pdf" "$dest"
  else
    magick "$IMG/${src}" "$dest"
  fi
  printf '%s,%s,%s\n' "Figure_${label}.pdf" "Figure ${label}" "$title" >> "$MAN_CSV"
  echo "  Figure_${label}.pdf  <-  ${src}"
}

echo "Main figures:"
emit 1  "I:study_design.png"            "Multi-modal data fusion (study design)"
emit 2  "P:fig-summary-1"               "Cell types underlying thousands of phenotypes"
emit 3  "P:fig-ontology-lvl-1"          "Specific phenotypes: fewer genes and cell types"
emit 4  "P:fig-rni-1"                   "Recurrent bacterial infection subtypes"
emit 5  "P:fig-network-rni-1"           "Causal network of recurrent Neisserial infections"
emit 6  "P:fig-celltype-severity-dot-1" "Cell types and phenotype severity"
emit 7  "P:fig-congenital-1"            "Congenital phenotypes and foetal cell types"
emit 8  "P:fig-therapy-validate-1"      "Prioritised targets recapitulate gene therapies"
emit 9  "P:fig-therapy-examples-1"      "Evidence-based gene therapy targets"

echo "Supplementary figures:"
emit S1  "P:fig-evidence-histograms-1"  "GenCC evidence score distributions"
emit S2  "I:fig-diagram.png"            "Multi-scale disease investigation strategy"
emit S3  "P:fig-ctd-correlation-1"      "Inter- and intra-dataset validation"
emit S4  "P:fig-monarch-recall-1"       "Monarch Knowledge Graph recall"
emit S5  "P:fig-celltype-severity-bar-1" "Cell types ordered by mean severity"
emit S6  "P:fig-therapy-filter-1"       "Prioritised target filtering steps"
emit S7  "P:fig-therapy-validate-all-1" "Validation of prioritised targets (all therapies)"
emit S8  "P:fig-animal-models-1"        "Experimental model translatability"
emit S9  "P:fig-therapy-examples2-1"    "Causal multi-scale networks (examples)"
emit S10 "I:fig-therapy-examples-supp/respiratory_failure.png"     "Respiratory failure network"
emit S11 "I:fig-therapy-examples-supp/dementia.png"                "Dementia network"
emit S12 "I:fig-therapy-examples-supp/lethal_skeletal_dysplasia.png" "Lethal skeletal dysplasia network"
emit S13 "I:fig-therapy-examples-supp/small_vessel_disease.png"    "Small vessel disease network"
emit S14 "I:fig-therapy-examples-supp/parkinson.png"               "Parkinson's disease networks"
emit S15 "I:fig-therapy-examples-supp/alzheimer.png"               "Alzheimer's disease networks"
emit S16 "I:fig-therapy-examples-supp/als.png"                     "Amyotrophic lateral sclerosis network"

echo
echo "Wrote $(ls "$OUT"/Figure_*.pdf | wc -l | tr -d ' ') figure PDFs to: $OUT"
echo "Manifest: $MAN_CSV"
