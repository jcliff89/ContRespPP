
# ContRespPP

ContRespPP is an implementation of the Bayesian approach to using 
predictive probability in an ANOVA construct with a continuous normal response, 
when threshold values must be obtained for the question of interest to be 
evaluated as successful ([Sieck and Christensen (2021)](https://doi.org/10.1002/qre.2802).
In this package, the Bayesian Mission Mean (BMM) is used to evaluate a question 
of interest (that is, a mean that randomly selects combination of factor levels 
based on their probability of occurring instead of averaging over the factor 
levels, as in the grand mean). Under this construct, in contrast to a Gibbs 
sampler (or Metropolis-within-Gibbs sampler), a two-stage sampling method is 
required. The nested sampler determines the conditional posterior distribution 
of the model parameters, given Y, and the outside sampler determines the marginal 
posterior distribution of Y (also commonly called the predictive distribution for Y). 
This approach provides a sample from the joint posterior distribution of Y and 
the model parameters, while also accounting for the threshold value that must be 
obtained in order for the question of interest to be evaluated as successful.

## Installation

You can install the released version of ContRespPP from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ContRespPP")
```

or from [Github](https://github.com/jcliff89/ContRespPP) with:

``` r
devtools::install_github("jcliff89/ContRespPP")
```

## How to Use

For details, please refer to the package vignette `vignette("gibbs-sampler")`.
