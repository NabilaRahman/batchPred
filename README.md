# Objective Prediction of Confounders
batchPred tool helps identify the optimal set of covariates to correct for when you have many that could potentially confound your data. e.g. known variables: age, sex, alignment qc; as well as unknown variables: principle components. The package requires an expression dataset, at least 1 dataframe of covariates, and a reference gene-gene association network with a confidence score. (More Details coming soon)

## Installation
You can install the most recent development version from github
using devtools:

``` r
# install.packages("devtools")
devtools::install_github("NabilaRahman/batchPred")
```

### References 
Somekh, J., Shen-Orr, S.S. & Kohane, I.S. Batch correction evaluation framework using a-priori gene-gene associations: applied to the GTEx dataset. BMC Bioinformatics 20, 268 (2019). https://doi.org/10.1186/s12859-019-2855-9
