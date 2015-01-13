	
#Software for Hierarchical Modeling of Genome-wide Association Studies

## Purpose
This page provides example code for the manuscript 'Enriching the analysis of genome-wide association studies with hierarchical modeling' by Chen and Witte (American Journal of Human Genetics 2007;81:397-404).  Specifically, the R script generates hierarchical modeling estimates for regression based association tests through the use of SNP level covariate data. 

## Use
The first stage estimate file and covariate file are provided as example input files.  Note that header rows must not be included and that both files must share the same number of rows. Please refer to the comments in hm.R for more information.

1. Load the R script in R by typing: source("hm.r")
 
1. Invoke the hm function in our example by typing: hm("first_stage.txt","covariates.txt",16,0.1,0.01)
where,

*. the first parameter is the name of the file containing first stage regression coefficients and standard errors,
*. the second parameter is the name of the file containing the SNP level covariates considered as the prior level information,
*. the third parameter is the column number in the covariates.txt file that is to be used in the weighting function,
*. the fourth column (0.1 here) is the baseline residual standard error tau,
*. and the final column (0.01 here) is the minimum standard deviation rho used in the weighting function. 
 
1. Two output files are generated. The first lists the estimates and standard error for the second stage coefficients pi, which measure the effect each covariate has on first stage coefficients beta hat. The second lists the posterior hierarchical modeling estimates and standard errors followed by the second stage hierarchical modeling estimates and standard errors.
