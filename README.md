
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dexterMML

DexterMML is currently under development. Interfaces to funtions are
likely to change.

## Overview

DexterMML is an extension of the R package
[dexter](https://CRAN.R-project.org/package=dexter). Dexter fits IRT
models using Conditional Maximum Likelihood (CML) which is our preferred
method. DexterMML, as the name implies, offers Marginal Maximum
Likelihood (MML) estimation of IRT models.

The main difference is that in MML one needs to assume a population
distribution. Misspecification of that distribution can result in biased
estimation of item parameters. Since CML does not need *a priori*
specification of a population distribution, it is naturally robust
against misspecification. MML also supposes a random sample of the
population which is not needed for CML.

Since assumptions buy information, MML can be used in some circumstances
where CML cannot and it can be used to fit much more elaborate models
(which we do not but other packages do). Specifically MML can be used -
with caution\! - when your dataset includes adaptive or multi stage
tests.

DexterMML distinguishes itself from other MML R-packages by:

  - including far fewer models and options
  - being considerably faster
  - support for the dexter data(base) structure

If the user desires more flexibility in model choice and estimation
options, they are better of using R packages
[TAM](https://CRAN.R-project.org/package=TAM) or
[mirt](https://CRAN.R-project.org/package=mirt)

## Installation

DexterMML is not currently on CRAN and can only be installed from github

``` r
# install.packages("devtools")
devtools::install_github("jessekps/dexterMML")
```

## Usage

``` r
# to do: examples
```
