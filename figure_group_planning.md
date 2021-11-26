# Plot groupings
Individual plots need to be grouped into larger figures. It can be done with
[this](https://patchwork.data-imaginist.com/) package. This is where to plan
what to include in each figure. I will be referring to plots from the masters
writeup.

## 1 - Introduction and examples
These are plots that were used in the introduction section of the masters report.
A lot of it may not be relevant for the paper.
Group:
* Figure 1 ADHD ancestor terms - From introduction section, used to help describe
the DAG structure.
* Figure 2 Rare disease EWCE homepage - Screen shot of the EWCE homepage (not
needed in the paper?)
* Figure 3 Interactive web ap - Screen shot of the interactive web app (Pos not
needed as can just include a link to the app?)
* Figure 4 Print friendly version of cell select - Example of network plot produced
by the non interactive version of cell select. This is the modified version of the
plotting funciton included in ontologyIndex (I changed it to allow mapping colour
to a variable/ heatmap).

## 2 - Overview results
These are figures that were made to get a sense of the overall number of significant
results and that they seemed to be sensible.
Group:
* Figure 5 Number of significant enrichments per cell - Shows how many significant
phenotypes were found for each cell type
* Figure 6 Main branches of HPO - shows the main branches (child terms of
Phenotypic abnormality) and three of them are highlighted as they are to be used as
examples in the following plots
* Figure 7 Number of signif enrichemnts per cell in highlighted branches

## 3 - HPO patterns and relationships
These plots show general observations about the HPO such as relationship between
ontology level and number of genes.
Group:
* Figure 8 Relationship between significance trehshold and proportion of enrichments
found in the expeced HPO branch - e.g neuronal cells are typically associated with
abnormalities of the nervous system, and these enrichments are typically the most
significant. (Maybe this should be in the overview results plot above as it
references  the 3 example branches again - cardio, neuro, immune).
* Figure 9 Ontology level relationships - Replace this with the new 3 plot version
that uses violin plots (see figures_demo.rmd).

## 4 - Infectoin related plots
These plots show infection related results at different levels of resolution (ont lvl).
Group:
* Figure 10 - Significant enrichments per cell for descendants of recurrent infections
* Figure 11 - Cell-phenotype assocaitions for child terms of recurrent bacterial infection.
* Figure 12 - Enrichments for child terms of recurrent gram-negative bacterial infection.


## 5 - Social iteraction related plots
Plots related to imparied social interaction results at differnt levels of resolution.
Not sure if these are wanted in the actual paper as its showing the amacrine cells
poor eye contact relatioship which may not really be significant? Would need to do
conditional enrichemnt tests maybe?
Group:
* Figure 13 - Cell types enriched for imparied social interaction
* Figure NA - Possibly add a ontology network plot here showing impaired social
interactoin with its child terms branching off it?
* Figure 14 - Poor eye contact plot - shows that the result for amacrine cells of
the retina was specific to poor eye contact.

