# KRN-GWAS

This is a research repo for our project entitled "**Empirical Comparisons of Different Statistical Modelsto Identify and Validate Kernel RowNumber-Associated Variants using NestedAssociation Mapping and Related Populations in Maize**". 

## Introduction
In this study, we compared different statistical models for doing GWAS and cross-validated the KRN-associated variants identified in the initial GWAS using three unrelated populations. The KRN-associated variants identified in this study have the potential to enhance our understanding of the developmental steps involved in ear development.

## Architecture about this Repo
This project contains ~400 commits. A `largedata` directory was intentionally ignored by adding to `gitignore` because of the large size of the files within the folder. To guide the visitors having a better idea about the repo, here we briefly introduce the functions or sepecific purposes of the directory system. The layout of directories is based on the idea from [ProjectTemplate](http://projecttemplate.net/architecture.html). 

1. **cache**: Here we store intermediate data sets that are generated during a preprocessing step.
2. **data**: Here we store our raw data of small size. Data of large size, i.e. > 100M, store in a `largedata` folder that has been ignored using `gitignore`.
3. **doc**: Documentation codes (i.e. Rmd files) for generating the figures.
4. **graphs**: Graphs produced during the study.
5. **lib**: Some functions for our work.
6. **profilling**: Analysis scripts for the project. It contains some sub-directories.
7. **table**: Table produced during the study.

## Figures
Rmd code to generate some of the Figures in the paper.


## License
This repo is free and open source for research usage, licensed under [GPLv2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
