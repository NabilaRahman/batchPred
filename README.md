# Objective Prediction of Confounders
batchPred is an R package that helps identify the most relevant set of covariates to correct for when you have many variables that could potentially confound your data, e.g known variables such as age, sex, read quality, or unknown variables such as principle components. The package requires an expression dataset, at least 1 dataframe of covariates, and a reference gene-gene association network with a confidence score. In our examples, we use 'full network' gold standards from the [GIANT database](http://giant.princeton.edu/download).

## Installation
You can install the most recent development version from github
using devtools:

``` r
# install.packages("devtools")
devtools::install_github("NabilaRahman/batchPred")
```

## Quick Start
Step 1: Prepare your reference network using **refNet** function.  
  
Step 2: Take none overlapping subset(s) from this reference network using **refNetSubsets** function.
  
Step 3: Compute the most relevant set of covariates to correct for these subset(s) using **batchPred** function.  
It is recommended to run batchPred using at least 3 distinct reference subsets and add covariates to your model if they show at least 0.1% improvement in AUC score in each run. 
  
Step 4: Visualise selected covariates using the **plotBceF2** function.

### References 
Somekh, J., Shen-Orr, S.S. & Kohane, I.S. Batch correction evaluation framework using a-priori gene-gene associations: applied to the GTEx dataset. BMC Bioinformatics 20, 268 (2019). https://doi.org/10.1186/s12859-019-2855-9
