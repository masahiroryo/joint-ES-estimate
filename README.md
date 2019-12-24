# Effect Size Analysis for Joint Factors  

## This repository offers the R script used in Rillig et al. (2019) in Science.  
![image](https://user-images.githubusercontent.com/30728887/71410433-81c52000-2645-11ea-98a1-59d9d3e5754d.png)  


## The R script brings the following outcomes:  
1) The effect size of each factor with 95% confidence intervals (bootstrap resampling)
2) The effect size of each joint factor level with 95%CI, compared with the effect size expected from single factor level in three different null expectations
3) The nonlinear regression fit to the relationship between the number of factors and the measured response for each response variable
4) The variability explained using a random forest algorithm (R2) quantified in three different modeling assumptions

![Rplot_example](https://user-images.githubusercontent.com/30728887/71410345-301c9580-2645-11ea-9d04-97416a2ac019.jpeg)  

## Note  
The dataset stored in this repository is essentially identical to the dataset in figshare (https://figshare.com/s/30b479473025f366013d); however, for the R script, some factors are named differently and the information of all single factor level treatments are duplicated in the dataset with the variable "remark" as "1". This data duplication does not affect to the result at all but needed for simplifying the script for visualization.
