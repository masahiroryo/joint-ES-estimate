# Effect Size Analysis for Joint Factors

This repository offers the R script used in Rillig et al. (2019) in Science.

It takes the following outcomes:
1) The effect size of each factor with 95% confidence intervals (bootstrap resampling)
2) The effect size of each joint factor level with 95%CI, compared with the effect size expected from single factor level in three different null expectations
3) The nonlinear regression fit to the relationship between the number of factors and the measured response for each response variable
4) The variability explained using a random forest algorithm (R2) quantified in three different modeling assumptions


Note that the dataset stored in this repository is essentially identical to the dataset in figshare (https://figshare.com/s/30b479473025f366013d); however, for the R script, some factors are named differently and the information of all single factor level treatments are duplicated in the dataset with the variable "remark" as "1". This data duplication does not affect to the result at all but needed for simplifying the script for visualization.
