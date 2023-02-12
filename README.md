## Identification of cell types underlying thousands of rare diseases and subtraits

Kitty B. Murphy, Bobby Gordon-Smith, Jai Chapman, Momoko Otani, Brian M. Schilder, Nathan G. Skene

### Introduction
Rare diseases (RDs) are uncommon as individual diagnoses, but as a group contribute to an enormous disease burden globally. However, partly due the low prevalence and high diversity of individual RDs, this category of diseases is understudied and under-resourced. The advent of large, standardised genetics databases has enabled high-throughput, comprehensive approaches that uncover new insights into the multi-scale aetiology of thousands of diseases. Here, using the Human Phenotype Ontology (9,677 annotated phenotypes) and multiple single-cell transcriptomic atlases (77 human cell types and 38 mouse cell types), we conducted >688,000 enrichment tests (x100,000 bootstrap iterations each) to identify >13,888 genetically supported cell type-phenotype associations. 

This repository contains the data and code needed to replicate the analyses in our preprint [insert link to preprint], as well as links to the R packages required (see below). 

- [HPOExplorer](https://github.com/neurogenomics/HPOExplorer): contains useful functions for working with the Human Phenotype Ontology (HPO). It allows you to create interactive phenotype network plots, as well as many other useful functions.

- [MultiEWCE](https://github.com/neurogenomics/MultiEWCE): allows you to run [EWCE](https://github.com/NathanSkene/EWCE) on multiple gene lists in parallel. 
