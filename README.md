
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ccdf`

[![CRAN
status](https://www.r-pkg.org/badges/version/ccdf)](https://CRAN.R-project.org/package=ccdf)
[![R-CMD-check](https://github.com/sistm/ccdf/workflows/R-CMD-check/badge.svg)](https://github.com/sistm/ccdf/actions)

## Overview

`ccdf` is a package for performing single-cell RNA-seq differential
expression analysis and more generally **complex hypothesis testing**.

The main function of the package is `ccdf_testing()`. It allows to use
either an asymptotic test for large sample size or a permutation test
for small sample size with the argument `method`.

The methods implemented in this package are detailed in the following
article:

> Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2020).
> Distribution-free complex hypothesis testing for single-cell RNA-seq
> differential expression analysis, *BioRxiv*
> [doi:10.1101/2021.05.21.445165](https://doi.org/10.1101/2021.05.21.445165)

## Installation

**`ccdf` is available from CRAN.**

``` r
install.packages("ccdf")
```

**To install `ccdf`, you can download the development version on
[GitHub](https://github.com/Mgauth/ccdf).**

``` r
#install.packages("devtools")
devtools::install_github("sistm/ccdf")
```

## Example

Here is a basic example which shows how to use `ccdf` with simple
generated data.

``` r
## Data Generation
X <- data.frame("X1" = as.factor(rbinom(n=100, size = 1, prob = 0.5)))
Y <- data.frame("gene1" = t(replicate(10, ((X$X1==1)*rnorm(n = 50,0,1)) + ((X$X1==0)*rnorm(n = 50,0.5,1)))))
```

``` r
# Hypothesis testing
res_asymp <- ccdf_testing(exprmat=Y, variable2test=X, test="asymptotic") # asymptotic test
res_perm <- ccdf_testing(exprmat=Y, variable2test=X, test="permutations",
                         adaptive=TRUE) # adaptive permutation test
```

– Marine Gauthier, Denis Agniel, Rodolphe Thiébaut & Boris Hejblum
