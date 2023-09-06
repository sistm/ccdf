
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `citcdf` <a><img src='man/figures/logo.svg' align="right" height="139" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/citcdf)](https://CRAN.R-project.org/package=citcdf)
[![R-CMD-check](https://github.com/sistm/citcdf/workflows/R-CMD-check/badge.svg)](https://github.com/sistm/citcdf/actions)

## Overview

`citcdf` is a package to perform conditional independence testing using
empirical conditional cumulative distribution function estimations.

The main function of the package is `cit()`. It uses an asymptotic test
(for large sample size , or a permutation test for small sample size
with the argument `method`) to perform conditional independence testing.

The approach implemented in this package is detailed in the following
article:

> Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2020).
> Distribution-free complex hypothesis testing for single-cell RNA-seq
> differential expression analysis, *BioRxiv*
> [doi:10.1101/2021.05.21.445165](https://doi.org/10.1101/2021.05.21.445165)

## Installation

**`citcdf` is available from
[GitHub](https://github.com/sistm/citcdf):**

``` r
#install.packages("devtools")
remotes::install_github("sistm/citcdf")
```

## Example

Here is a basic example which shows how to use `citcdf` with simple
generated data.

``` r
## Data Generation
X <- data.frame("X1" = as.factor(rbinom(n=100, size = 1, prob = 0.5)))
Y <- data.frame("Y1" = t(replicate(10, ((X$X1==1)*rnorm(n = 50,0,1)) + ((X$X1==0)*rnorm(n = 50,0.5,1)))))
```

``` r
# Hypothesis testing
res_asymp <- cit(exprmat=Y, variable2test=X, test="asymptotic") # asymptotic test
res_perm <- cit(exprmat=Y, variable2test=X, test="permutation",
                         adaptive=TRUE) # adaptive permutation test
```

– Marine Gauthier, Denis Agniel, Kalidou Ba, Rodolphe Thiébaut & Boris
Hejblum

*hex illustration by Jérôme Dubois.*
