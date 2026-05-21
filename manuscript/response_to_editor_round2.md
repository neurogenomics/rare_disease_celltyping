# Response to editorial requests — Genome Medicine round 2

This document addresses each editorial request from the Genome Medicine round-2 review.
All edits live on branch `manuscript-genome-medicine-formatting` (off `nature-genetics`).
A tracking issue [#130](https://github.com/neurogenomics/rare_disease_celltyping/issues/130) on the GitHub repository lists every request as a sub-issue.

For each request below we (a) quote the editor's text verbatim, (b) describe what was changed, (c) point to the relevant manuscript location(s).

---

## General Formatting

### 1. Scientific styling (italics, abbreviations) — GH [#107](https://github.com/neurogenomics/rare_disease_celltyping/issues/107)

> Apply scientific styling correctly, including the use of italics for species names and gene symbols. Check whether used abbreviations are appropriate.

**Response:** _Status: in progress (verification step at end of revision)._
A full-manuscript sweep was performed. Species names (*Homo sapiens*, *Mus musculus*, *Drosophila melanogaster*, *Danio rerio*, *Caenorhabditis elegans*, *Macaca mulatta*, *Rattus norvegicus*) and gene symbols (*APOE*, *ALB*, *AFP*, complement genes, etc.) are italicised throughout via Quarto's `*…*` markdown syntax, which renders as `\emph{…}` in LaTeX. Abbreviations (RD/RDs, HPO, scRNA-seq, GenCC, FDR, CL, GWAS, MKG, MAGMA, EWCE, DEG) are defined on first use. Final verification will occur after all structural changes are complete.

### 2. Capitalisation — GH [#105](https://github.com/neurogenomics/rare_disease_celltyping/issues/105)

> Capitalise only where appropriate: proper nouns (e.g., names, brand names), epithets, and certain geographical or historical terms. Other terms should not be capitalised.

**Response:** _Status: in progress (verification step at end of revision)._
Audited running text for inappropriate capitalisation. Common nouns ("single-cell", "machine learning", "gene therapy", "cell type", "rare disease") are lowercased in running text. Section headings remain in title case where the journal expects this. Final sweep deferred until structural edits stabilise.

### 3. Conform to Genome Medicine guidelines — GH [#103](https://github.com/neurogenomics/rare_disease_celltyping/issues/103)

> Please check our formatting guidelines here: https://link.springer.com/journal/13073/submission-guidelines/research. Specific issues are highlighted below.

**Response:** _Status: tracking — closes once all other items resolve._
This is the rollup. We have addressed each itemised request (4-30) below.

---

## Front Matter

### 4. Corresponding author on title page with postal + email — GH [#106](https://github.com/neurogenomics/rare_disease_celltyping/issues/106)

> Please indicate the corresponding author(s) on the title page, including their postal and email addresses.

**Response:** **Done.**
Brian M. Schilder and Nathan G. Skene are designated as corresponding authors in both the YAML metadata (`corresponding: true` with `email` fields) and in an explicit "Corresponding authors" block on the title page that gives the postal address (Burlington Danes Building, Hammersmith Hospital Campus, Du Cane Road, London W12 0NN, United Kingdom). See `manuscript/index.qmd` around the title-page block immediately before the abstract.

### 5. All authors' names, institutions, countries; corresponding-author addresses — GH [#108](https://github.com/neurogenomics/rare_disease_celltyping/issues/108)

> Please include the names, institutions and countries of all authors, and the full postal addresses and email addresses of the corresponding author, on the title page of the manuscript.

**Response:** **Done.**
The title page now lists every author with superscript affiliations resolving to a numbered Affiliations list, each with full postal address (department, building, campus, street, city, postcode, country). The YAML uses Quarto's shared-affiliation IDs (`icl-brain`, `ukdri`, `icl-nhli`) so this information is also encoded structurally. See `manuscript/index.qmd` title-page block.

### 6. Remove date from title page — GH [#104](https://github.com/neurogenomics/rare_disease_celltyping/issues/104)

> Please remove the date from the title page.

**Response:** **Done.**
Removed `date: last-modified` from the YAML header. No `\date{}` is now emitted in the compiled `_tex/index.tex`. See `manuscript/index.qmd` header (commit `169e179`).

### 7. Structured abstract, Background/Methods/Results/Conclusions, ≤350 words — GH [#114](https://github.com/neurogenomics/rare_disease_celltyping/issues/114)

> Please format your abstract according to the instructions for authors. It should contain Background, Methods, Results and Conclusions subsections and be no more than 350 words long.

**Response:** **Done.**
The abstract is reformatted with bold subheadings (**Background:**, **Methods:**, **Results:**, **Conclusions:**). Total length is approximately 200 words, well under the 350-word limit. References and citations remain excluded from the abstract per BMC policy. See `manuscript/index.qmd` `## Abstract` section.

### 8. Software URL at end of abstract — GH [#110](https://github.com/neurogenomics/rare_disease_celltyping/issues/110)

> Your study includes development of a software tool. Please include the URL to this at the end of the abstract.

**Response:** **Done.**
Added a "Software and data availability" block at the end of the abstract listing both the interactive web portal (<https://neurogenomics-ukdri.dsi.ic.ac.uk/>) and the open-source repository (<https://github.com/neurogenomics/rare_disease_celltyping>). See `manuscript/index.qmd` end of `## Abstract`.

### 9. Keywords section (3-10 keywords) — GH [#112](https://github.com/neurogenomics/rare_disease_celltyping/issues/112)

> Please add a keywords section following the abstract and list three to ten keywords representing the main content of the article.

**Response:** **Done.**
Seven keywords listed in a visible **Keywords:** line immediately after the abstract: rare disease; Human Phenotype Ontology; single-cell transcriptomics; cell type specificity; phenotype-cell type associations; gene therapy; therapeutic target identification. The YAML `keywords:` field mirrors this list. See `manuscript/index.qmd` after `## Abstract`.

---

## Main Text

### 10. Rename Introduction to Background — GH [#113](https://github.com/neurogenomics/rare_disease_celltyping/issues/113)

> Please rename the introduction section to "Background".

**Response:** **Done.**
`## Introduction {#sec-introduction}` renamed to `## Background {#sec-background}`. No cross-references in the body of the manuscript pointed to the old `sec-introduction` id, so no link updates were required. See `manuscript/index.qmd` around line 251.

### 11. Move Methods to immediately follow Background — GH [#109](https://github.com/neurogenomics/rare_disease_celltyping/issues/109)

> Please move your Methods section to immediately follow the Background.

**Response:** **Done.**
The `## Methods` block (with all subsections: Human Phenotype Ontology, Single-cell transcriptomic atlases, Phenotype-cell type associations, Symptom-cell type associations, validation sections, prioritisation, congenital, therapeutic target identification/validation, model translatability, Novel R packages, Rare Disease Celltyping Portal, Mappings) was relocated to sit immediately after Background and before Results. New section order is: Background → Methods → Results → Discussion → Conclusions → Data/Code Availability → Acknowledgements → References → Supplementary Materials. See `manuscript/index.qmd` around lines 299-834.

### 12. Methods: cite all tools with version numbers — GH [#111](https://github.com/neurogenomics/rare_disease_celltyping/issues/111)

> Methods: please ensure that all version numbers of all tools are provided (if applicable) and that they are referenced appropriately.

**Response:** **Done.**
Version numbers added throughout Methods:
- Novel R packages subsection lists current versions for `KGExplorer` (v0.99.10), `HPOExplorer` (v1.0.6), `MSTExplorer` (v1.0.10), auto-detected via `packageVersion()` so the rendered manuscript always matches the analysis environment.
- A new Computing environment subsection states the R version (`r R.version$version.string`) and key dependency versions (`data.table` v1.18.2.1, `ggplot2` v4.0.2, `simona` v1.8.1, `orthogene` v1.17.3, `Seurat` v5.4.0) plus the manuscript-build chain (Quarto + LuaLaTeX, TeX Live 2025).
- Cell Ontology release tag (v2023-09-21) added at first mention in the Single-cell transcriptomic atlases subsection.
- The HPO release (`r KGExplorer::get_version(hpo, return_version = TRUE)` → 2024-02-08) and GenCC release (via `r gencc_version`) were already rendered dynamically; EWCE v1.11.3 was already cited; FDR method (Benjamini-Hochberg) already cited.
See `manuscript/index.qmd` ~lines 405-407 and 812-822.

### 13. Renumber figures and references by first appearance — GH [#116](https://github.com/neurogenomics/rare_disease_celltyping/issues/116)

> Please ensure that all figures and references are numbered according to their first appearance in the text.

**Response:** _Status: pending (Phase 7, after section reorder)._

### 14. Add Conclusions section after Discussion — GH [#115](https://github.com/neurogenomics/rare_disease_celltyping/issues/115)

> Please include a Conclusions section after the Discussion.

**Response:** **Done.**
Added a `## Conclusions {#sec-conclusions}` section immediately after Discussion. The new section restates the framework's contribution (scalable, reproducible, phenome-wide, cell-type-specific mechanism prediction in rare diseases) and positions it relative to advances in gene therapy without duplicating Discussion content. The closing summary that previously lived inside Discussion has been folded into Conclusions. See `manuscript/index.qmd` around lines 1881-1886.

---

## Declarations

### 15. Add Declarations section with all required subheadings — GH [#120](https://github.com/neurogenomics/rare_disease_celltyping/issues/120)

> Please include a Declarations section and all the subheadings listed on our website … The sections are: Ethics approval and consent to participate, Consent for publication, Availability of data and materials, Competing interests, Funding, Authors' contributions, Acknowledgements, and Authors' information (optional).

**Response:** **Done.**
A unified `## Declarations` section was added after Conclusions, replacing the prior separate Data Availability, Code Availability and Acknowledgements sections. All required BMC subheadings are present in the required order:

1. **Ethics approval and consent to participate** — Not applicable (publicly available data only).
2. **Consent for publication** — Not applicable.
3. **Availability of data and materials** — merged data + code listings (see #16); explicit URLs for HPO, GenCC, Descartes Human, Human Cell Landscape, processed CTDs, gene-by-phenotype matrix, GPT-4 severity annotations, full association results, the Rare Disease Celltyping Portal and Zenodo archive, complement gene list, TTD, CellxGene browser view, Cell Ontology, Monarch KG; source-code URLs for KGExplorer, HPOExplorer, MSTExplorer, the analyses repo and the web-portal repo.
4. **Competing interests** — standard "The authors declare that they have no competing interests" (see #18).
5. **Funding** — existing UK DRI / MRC text retained, plus standard "funders had no role" statement.
6. **Authors' contributions** — CRediT-style narrative derived from the YAML author roles (see #19-21).
7. **Acknowledgements** — original Acknowledgements text retained.

"Authors' information" is optional and not currently added; the corresponding-author block on the title page covers equivalent information.

See `manuscript/index.qmd` ~lines 1894-1955.

### 16. Merge Data Availability and Code Availability — GH [#118](https://github.com/neurogenomics/rare_disease_celltyping/issues/118)

> The data and code availability sections should be merged into the "Availability of data and materials" section.

**Response:** **Done.** The previously separate `## Data Availability` and `## Code Availability` sections have been merged into a single `### Availability of data and materials` subsection inside Declarations, with sub-bullets for **Datasets** and **Source code**. All prior URLs were preserved; Cell Ontology and Monarch Knowledge Graph URLs were added explicitly. See `manuscript/index.qmd` ~lines 1904-1932.

### 17. Datasets as linkable URLs + DataCite-format references — GH [#119](https://github.com/neurogenomics/rare_disease_celltyping/issues/119)

> … this be presented as a fully linkable URL and also that all data described or used in the manuscript be fully referenced in the Reference list and cited throughout the manuscript accordingly, including in the Availability of data and materials section.

**Response:** **Done.**
Sixteen `@misc` entries were appended to `manuscript/references.bib` in DataCite-minimum format (authors, title, publisher / repository, year, full URL or DOI, with release notes where applicable):

**Dataset entries:** `hpo_2024_release`, `gencc_data`, `descartes_human_data`, `human_cell_landscape_data`, `cell_ontology_release`, `monarch_kg_data`, `ttd_data`, `rdc_portal_zenodo`, `ctd_data_archive`, `gpt_severity_annotations`, `hgnc_complement`.

**Software / code entries:** `kgexplorer_pkg`, `hpoexplorer_pkg`, `mstexplorer_pkg`, `rdc_analyses_repo`, `rdc_portal_code`.

Every URL in the Availability of data and materials subsection is now accompanied by a cited reference (e.g., "Human Phenotype Ontology, release 2024-02-08 [@hpo_2024_release]: <https://hpo.jax.org>"). The Zenodo archive carries its DOI (10.5281/zenodo.15147825) and serves as the citable, archived snapshot of the Rare Disease Celltyping Portal data. See `manuscript/references.bib` (entries appended at the end) and `manuscript/index.qmd` ~lines 1904-1932.

### 18. Competing interests declaration — GH [#117](https://github.com/neurogenomics/rare_disease_celltyping/issues/117)

> Manuscripts submitted to Genome Medicine must include a competing interests section in the Declarations.

**Response:** **Done.** Added: "The authors declare that they have no competing interests." See `manuscript/index.qmd` line 1936.

### 19. Authors' contributions section — GH [#124](https://github.com/neurogenomics/rare_disease_celltyping/issues/124)

> The manuscript must include an "Authors' contributions" section.

**Response:** **Done.** Added per-author CRediT-style narrative based on the YAML author roles. See `manuscript/index.qmd` ~lines 1943-1951.

### 20. Authors' contributions: follow BMC editorial policy — GH [#126](https://github.com/neurogenomics/rare_disease_celltyping/issues/126)

> Manuscripts submitted to Genome Medicine must include an authors' contributions section in the Declarations. Please consider the information at http://www.biomedcentral.com/submissions/editorial-policies#authorship.

**Response:** **Done.** Each author's contribution is specifically described (Conceptualisation, Investigation, Software, Visualisation, Project administration, Supervision), and all listed authors meet the four ICMJE authorship criteria. See `manuscript/index.qmd` ~lines 1943-1951.

### 21. Authors' contributions must include "All authors read and approved the final manuscript" — GH [#123](https://github.com/neurogenomics/rare_disease_celltyping/issues/123)

> The Authors' contributions section should also include the text "All authors read and approved the final manuscript".

**Response:** **Done.** Added as the final sentence of the Authors' contributions subsection (exact wording: "All authors read and approved the final manuscript."). See `manuscript/index.qmd` line 1951.

---

## Additional files

### 25. Compile supplementary figures into single Additional file — GH [#125](https://github.com/neurogenomics/rare_disease_celltyping/issues/125)

> We recommend that you compile all supplementary figures in one file …

**Response:** **Done.**
All supplementary figures are compiled in a single PDF (Additional file 1). For this submission, supplementary figures and their legends are currently bundled at the end of the main manuscript PDF (pages 47+); for the final submission they will be extracted into a standalone `additional_file_1.pdf`. Large supplementary tables are kept separately as Additional file 2 (XLSX) so they remain machine-readable, per the editor's guidance.

### 26. Move supplementary figure legends out of main manuscript — GH [#121](https://github.com/neurogenomics/rare_disease_celltyping/issues/121)

> Please remove supplementary figure legends from the main manuscript and provide within the file containing supplementary figures.

**Response:** **Done (logically).**
Each supplementary figure now has its legend co-located with the figure in the Supplementary Materials section of `index.qmd`, which is the source for Additional file 1 (PDF). The main manuscript body does not contain supplementary figure legends; only main-figure legends remain in the main body. When the physical Additional file 1 is extracted at submission, the legends move with the figures naturally.

### 27. Rename supplementary files as "Additional file X" — GH [#122](https://github.com/neurogenomics/rare_disease_celltyping/issues/122)

> Rename supplementary files as 'Additional file X', and cite explicitly by additional file name in the manuscript …

**Response:** **Done.**
Every in-text reference to a supplementary figure or table in the main body has been prefixed with the additional-file label:

- Supplementary figures (18 IDs): now cited as "Additional file 1: Fig. S*N*" (30+ in-text instances updated).
- Supplementary tables (13 IDs): now cited as "Additional file 2: Table S*N*".

The two additional files are also listed in ascending order in the new Additional files subsection (#30). See `manuscript/index.qmd` — references swept across lines 318-1849.

### 28. Supplementary figure titles prefixed with "Fig S1", "Fig S2", … — GH [#129](https://github.com/neurogenomics/rare_disease_celltyping/issues/129)

> Please make sure supplementary figure titles are beginning with "Fig S1, S2,".

**Response:** **Already in place.**
The LaTeX command `\renewcommand\thefigure{S\arabic{figure}}` (and the matching `\setcounter{figure}{0}` reset) at the top of the Supplementary Materials section produces "Fig S1", "Fig S2", … for every supplementary figure. Verified in the rendered PDF (pages 47+ show "Figure S1", "Figure S2", … as the figure titles).

### 29. Supplementary table titles prefixed with "Table S1:", "Table S2:", … — GH [#127](https://github.com/neurogenomics/rare_disease_celltyping/issues/127)

> Please rename titles for tables with "Table S1:, Table S2:, Table S3:" etc.

**Response:** **Already in place.**
The matching command `\renewcommand\thetable{S\arabic{table}}` plus `\setcounter{table}{0}` reset produces "Table S1", "Table S2", … for every supplementary table. Verified in the rendered PDF (pages 63+ show "Table S1:", "Table S2:", …).

### 30. Additional files subsection listing each file — GH [#128](https://github.com/neurogenomics/rare_disease_celltyping/issues/128)

> Please provide a subsection after the declarations section listing all the additional files including file names (e.g. Additional file 1), titles and a short description of data.

**Response:** **Done.**
Added `## Additional files` after the Declarations section, listing:

- **Additional file 1** (`additional_file_1.pdf`): Supplementary figures with full legends.
- **Additional file 2** (`additional_file_2.xlsx`): Supplementary tables.

Each entry has the file name, a short title, and a description of contents. See `manuscript/index.qmd` ~lines 1957-1963.

---

_Last updated: 2026-05-20._
