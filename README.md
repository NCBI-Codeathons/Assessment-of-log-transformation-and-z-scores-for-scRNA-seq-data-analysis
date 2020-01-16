# Assessment-of-log-transformation-and-z-scores-for-scRNA-seq-data-analysis

## Team
* Roshan Sharma
* Heather Geiger
* Ravneet Kaur
* Vincent Liu

## Problem
Numerous methods have been proposed for the normalization of single-cell RNA-seq (scRNA-seq) data. Yet, these methods have not been thoroughly benchmarked and assessed, partly owing to the lack of a metric that appropriately measures the performancee of any normalization method.

## Proposed Solution and Methods
Here we benchmark the performance of eight commonly used scRNA-seq normalization methods, listed in the **Methods Tested** section, on various datasets. We use blah blah as a measure of normalization quality, based on the rationale that blah blah.

We will normalize the PBMC_10X_10K data using the list, listed in the **Methods Tested** section. Downsampling the counts from the data to leave all cells with a pre‐specified number of counts or fewer. Downsampling can deliver a more realistic representation of what cellular expression profiles would look like at similar count depths. 

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
