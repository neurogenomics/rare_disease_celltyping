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

**Response:** _Status: pending (Phase 3)._

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

**Response:** _Status: pending (Phase 4)._

### 16. Merge Data Availability and Code Availability — GH [#118](https://github.com/neurogenomics/rare_disease_celltyping/issues/118)

> The data and code availability sections should be merged into the "Availability of data and materials" section.

**Response:** _Status: pending (Phase 4)._

### 17. Datasets as linkable URLs + DataCite-format references — GH [#119](https://github.com/neurogenomics/rare_disease_celltyping/issues/119)

> … this be presented as a fully linkable URL and also that all data described or used in the manuscript be fully referenced in the Reference list and cited throughout the manuscript accordingly, including in the Availability of data and materials section.

**Response:** _Status: pending (Phase 4)._

### 18. Competing interests declaration — GH [#117](https://github.com/neurogenomics/rare_disease_celltyping/issues/117)

> Manuscripts submitted to Genome Medicine must include a competing interests section in the Declarations.

**Response:** _Status: pending (Phase 4)._

### 19. Authors' contributions section — GH [#124](https://github.com/neurogenomics/rare_disease_celltyping/issues/124)

> The manuscript must include an "Authors' contributions" section.

**Response:** _Status: pending (Phase 4)._

### 20. Authors' contributions: follow BMC editorial policy — GH [#126](https://github.com/neurogenomics/rare_disease_celltyping/issues/126)

> Manuscripts submitted to Genome Medicine must include an authors' contributions section in the Declarations. Please consider the information at http://www.biomedcentral.com/submissions/editorial-policies#authorship.

**Response:** _Status: pending (Phase 4)._

### 21. Authors' contributions must include "All authors read and approved the final manuscript" — GH [#123](https://github.com/neurogenomics/rare_disease_celltyping/issues/123)

> The Authors' contributions section should also include the text "All authors read and approved the final manuscript".

**Response:** _Status: pending (Phase 4)._

---

## Additional files

### 25. Compile supplementary figures into single Additional file — GH [#125](https://github.com/neurogenomics/rare_disease_celltyping/issues/125)

> We recommend that you compile all supplementary figures in one file …

**Response:** _Status: pending (Phase 5)._

### 26. Move supplementary figure legends out of main manuscript — GH [#121](https://github.com/neurogenomics/rare_disease_celltyping/issues/121)

> Please remove supplementary figure legends from the main manuscript and provide within the file containing supplementary figures.

**Response:** _Status: pending (Phase 5)._

### 27. Rename supplementary files as "Additional file X" — GH [#122](https://github.com/neurogenomics/rare_disease_celltyping/issues/122)

> Rename supplementary files as 'Additional file X', and cite explicitly by additional file name in the manuscript …

**Response:** _Status: pending (Phase 5)._

### 28. Supplementary figure titles prefixed with "Fig S1", "Fig S2", … — GH [#129](https://github.com/neurogenomics/rare_disease_celltyping/issues/129)

> Please make sure supplementary figure titles are beginning with "Fig S1, S2,".

**Response:** _Status: pending (Phase 5)._

### 29. Supplementary table titles prefixed with "Table S1:", "Table S2:", … — GH [#127](https://github.com/neurogenomics/rare_disease_celltyping/issues/127)

> Please rename titles for tables with "Table S1:, Table S2:, Table S3:" etc.

**Response:** _Status: pending (Phase 5)._

### 30. Additional files subsection listing each file — GH [#128](https://github.com/neurogenomics/rare_disease_celltyping/issues/128)

> Please provide a subsection after the declarations section listing all the additional files including file names (e.g. Additional file 1), titles and a short description of data.

**Response:** _Status: pending (Phase 5)._

---

_Last updated: 2026-05-20._
