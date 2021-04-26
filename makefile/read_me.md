# Make File rare disease EWCE
## Table of Contents
1. [Makefile.R](#1)
2. [source ](#2)
3. [Writeup.Rmd](#3)
4. [hpo_to_genelists-wholebody_2019hpo_format.r](#4)
5. [Issues ](#5)


## Makefile.R<a name = 1></a>
This is the main code that calls the other scripts, produces the CTD dataset, generates results and outputs the report. 

## source <a name = 2></a>
This directory contains the CTD gen script and the Results gen script.

## Writeup.Rmd<a name = 3></a>
This is the template for the writeup that is used to output the final report.

## hpo_to_genelists-wholebody_2019hpo_format.r<a name = 4></a>
This is a slightly modified version of `process_hpo_to_genelists-wholebody.r`. It now downlaods the genes_to_phenotype.csv data if its not already present from [here](https://ewce.page.link/gene2phenotype "stored on fig share"). 
### missing phenotypes
The new dataset from hpo linked above seems to be missing some phenotypes, and also has a few extra

|New total  |Original total  |Missing  |Extra  |
|:--        |:--             |:--      |:--    |
|1791       |2854            |1314     |251    |

## Issues <a name = 5></a>
* I dont have the full results so i cannot get the markdown file to knit as its expecting certain cardiac cells. I almost got it to work by replacing them with some bladder cells, but I think I need to test it with the full dataset. 
* Also, This is still using the tabula muris data from figshare (FACS directory and annotations_facs.csv)
