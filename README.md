# Assessment-of-log-transformation-and-z-scores-for-scRNA-seq-data-analysis

## Team
* Roshan Sharma @ NYGC
* Heather Geiger @ NYGC
* Ravneet Kaur @ Emory
* Vincent Liu @ MSKCC

## Problem
Numerous methods have been proposed for the normalization of single-cell RNA-seq (scRNA-seq) data. Yet, these methods have not been thoroughly benchmarked to assess robustness and influence on downstreasm analysis. A critical step in the analysis pipeline is to account for unwanted biological and technical effects that mask the signal of interest.

## Proposed Solution and Methods
Here we benchmark the performance of eight commonly used scRNA-seq normalization methods, listed in the **Methods Tested** section, on four 10x datasets. Reasoning that each normalization method has its own objectives that may not be shared by other methods, we believe there is no single metric that can fully represent the quality of any given method. Therefore, we explore three propertieis of the resulting normalization: robustness, correlation to library size (# transcripts per cell), and influence on gene expression.

To assess the robustness of a normalization method, we first downsample the counts from data to leave all cells with a pre‐specified percentage of counts or fewer. For each normalization method, we then assess correlation of the normalized, downsampled data to normalized full data. Downsampling can deliver a more realistic representation of what cellular expression profiles would look like at similar count depths. We will look at the correlation of highly variable genes and differentially expressed genes on various different data sets and compare them to library size. Correlation of genes before and after mean and standard deviation. This gives us the set of differentially expressed genes.



## Methods Tested
* Median
* Median + Log
* Median + Log + Z-score
* Median + Log + Linear Regression
* GLMPCA
* scVI
* scRAN
* scTransform

## Datasets Used
[Dropbox](https://www.dropbox.com/sh/4uk9bdtk5t8ud7z/AAAZBJmMAw6IPt_qcaLVcCtYa?dl=0)
* 10x PBMC 3k

## Issues
We had no coffee.

## Lessons Learned
Always a pain to set up environment for various methods.
