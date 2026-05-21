#!/usr/bin/env bash
#
# code/zip_tex.sh
#
# Bundle the rendered LaTeX submission package into
# manuscript/_manuscript/manuscript_tex.zip.
#
# Run as a Quarto post-render hook so the zip always matches the latest
# rendered PDF. The hook is configured in manuscript/_quarto.yml under
# project.post-render. To rebuild manually:
#
#   ./code/zip_tex.sh
#
# Bundle layout (preserves _tex/ as the top-level directory so a recipient
# can drop the contents straight into a journal submission portal):
#   _tex/index.tex
#   _tex/references.bib
#   _tex/agujournal2019.cls
#   _tex/trackchanges.sty
#   _tex/index_files/figure-pdf/*.pdf
#   _tex/img/*

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
TEX_DIR="$PROJECT_ROOT/manuscript/_manuscript/_tex"
ZIP_OUT="$PROJECT_ROOT/manuscript/_manuscript/manuscript_tex.zip"

if [[ ! -d "$TEX_DIR" ]]; then
  echo "[zip_tex] skipping: $TEX_DIR not present (run quarto render first)" >&2
  exit 0
fi

# Build the zip from inside _manuscript/ so the archive paths start with _tex/.
cd "$PROJECT_ROOT/manuscript/_manuscript"

# Remove any stale zip so we rebuild from scratch
rm -f "$ZIP_OUT"

zip -r -q "$ZIP_OUT" _tex
echo "[zip_tex] wrote $ZIP_OUT ($(du -h "$ZIP_OUT" | awk '{print $1}'), $(unzip -l "$ZIP_OUT" | tail -1 | awk '{print $2}') files)"
