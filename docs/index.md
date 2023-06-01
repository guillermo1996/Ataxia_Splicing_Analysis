# Ataxia Splicing Noise Analysis

This repository contains the code and reports related to the splicing analysis on Ataxia RNAseq samples.

The code presented here is based on the [splicing-accuracy-manuscript](https://github.com/SoniaRuiz/splicing-accuracy-manuscript) repository from Sonia García Ruiz, and was developed in order to provide additional information regarding the effects ataxia on splicing accuracy.

## Analyses

### Median MSR studies

Analysis focused on the differences in mis-splicing ratio medians between cases and controls for different levels of study. Provided [report](https://guillermo1996.github.io/Ataxia_Splicing_Analysis/RMarkdown/Ataxia_Pseudobulk_MedianMSR.html) and [code](https://github.com/guillermo1996/Ataxia_Splicing_Analysis/blob/main/RMarkdown/Ataxia_Pseudobulk_MedianMSR.Rmd).

### Distances, deltaMES and MSR distribution studies

Analysis focused on the differences in distances, deltaMES and MSR distributions between cases and controls for different levels of study. Provided [report](https://guillermo1996.github.io/Ataxia_Splicing_Analysis/RMarkdown/Ataxia_Pseudobulk_Distances.html) and [code](https://github.com/guillermo1996/Ataxia_Splicing_Analysis/blob/main/RMarkdown/Ataxia_Pseudobulk_Distances.Rmd).

### Junction category distribution

A series of reports studying the junction category distributions. These reports are based on the research by Regina H. Reynolds ([GitHub](https://github.com/RHReynolds/LBD-seq-bulk-analyses)):

-   Leafcutter junction categories: provided [report](https://guillermo1996.github.io/Ataxia_Splicing_Analysis/RMarkdown/Leafcutter_Junction_categories.html) and [code](https://github.com/guillermo1996/Ataxia_Splicing_Analysis/blob/main/RMarkdown/Leafcutter_Junction_categories.Rmd).

-   STAR junction categories: provided [report](https://guillermo1996.github.io/Ataxia_Splicing_Analysis/RMarkdown/STAR_Junction_categories.html) and [code](https://github.com/guillermo1996/Ataxia_Splicing_Analysis/blob/main/RMarkdown/STAR_Junction_categories.Rmd).

-   Targeted gene junction categories: provided [report](https://guillermo1996.github.io/Ataxia_Splicing_Analysis/RMarkdown/Gene_Junction_categories.html) and [code](https://github.com/guillermo1996/Ataxia_Splicing_Analysis/blob/main/RMarkdown/Gene_Junction_categories.Rmd).

## Repository structure

Within this repository, the following structure was followed:

-   [**R**](R/): helper functions employed across all the relevant scripts and RMarkdowns.

-   [**RMarkdown**](RMarkdown/): directory containing all the reports for the different studies. Provided are both the `.Rmd` and their `.html` output, which can be viewed via [github pages](https://guillermo1996.github.io/Ataxia_Splicing_Analysis/).

-   [**Scripts**](Scripts/): scripts employed for the different analyses:

    -   [**download_and_extraction.R**](Scripts/download_and_extraction.R): script to download the BAM files from AWS S3. It requires login credentials which must be manually updated every 12 hours. Other helper scripts related to the credential generation can also be found.

    -   [**splicing_noise_analysis.R**](Scripts/splicing_noise_analysis.R): main script to run the splicing noise analysis.

    -   [**splicing_analysis_levels**](Scritps/splicing_analysis_levels): directory containing the scripts to finalize the splicing noise analysis for each of the levels and tissues studied.

## Issues

The code in this repository was executed using R version 4.2.1 in Ubuntu 20.04.4 LTS under [Docker bioconductor](https://registry.hub.docker.com/r/bioconductor/bioconductor_docker) image `RELEASE_3_15` with additional installations (e.g. regtools v1.0.0, samtools v1.17, bedtools v2.30.0 and others).

The software in this repository is provided "as is", without warranty of any kind under an MIT license. For any bugs or mistakes found in the code, please raise an issue on the GitHub repository or contact [Guillermo Rocamora Pérez](mailto:guillermorocamora@gmail.com).
