
<!-- README.md is generated from README.Rmd. Please edit that file -->

# passt

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/johannes-titz/passt.svg?branch=master)](https://travis-ci.org/johannes-titz/passt)
<!-- badges: end -->

The *passt* package is an R implementation of the PASS-T model, an
artificial neural network designed to explain how humans make judgments
of frequency and duration.

## Installation

<!-- 
You can install the released version of passt from [CRAN](https://CRAN.R-project.org) with:


```r
#install.packages("passt")
```
-->

You can install the development version from
[GitHub](https://github.com/) with devtools (and vignettes build, this
takes a couple of seconds):

``` r
devtools::install_github("johannes-titz/passt", build_vignettes = TRUE)
```

## Example

Run a simulation with 10 orthogonal stimuli presented for 1, 2, …, 10
times and durations of 10 s, 9 s, …, 1 s.

``` r
library(passt)
sim1 <- run_sim(patterns = diag(10), frequency = 1:10, duration = 10:1,
                lrate_onset = 0.05, lrate_drop_time = 2, lrate_drop_perc = 0)
head(sim1$output)
#>           [,1]      [,2]      [,3]      [,4]      [,5]     [,6]     [,7]
#> [1,] 0.8120493 0.8627214 0.9062597 0.9497354 0.9896714 1.025086 1.062488
#> [2,] 0.8141704 0.8568780 0.9065295 0.9491241 0.9914349 1.031216 1.062683
#> [3,] 0.8127721 0.8577301 0.9103283 0.9497560 0.9878989 1.027558 1.061171
#> [4,] 0.8110864 0.8590330 0.9068912 0.9490866 0.9862814 1.030735 1.066927
#> [5,] 0.8077268 0.8613726 0.9039062 0.9476095 0.9876164 1.028995 1.064632
#> [6,] 0.8135064 0.8561094 0.9083850 0.9451558 0.9929584 1.029055 1.065343
#>          [,8]     [,9]    [,10]
#> [1,] 1.095977 1.129194 1.166817
#> [2,] 1.092363 1.133150 1.162450
#> [3,] 1.098592 1.131645 1.162548
#> [4,] 1.097970 1.131607 1.160383
#> [5,] 1.100242 1.131482 1.166417
#> [6,] 1.095027 1.133914 1.160546
```

The output is the *memory strength* for each stimulus. One row
corresponds with one simulation, while one column is one stimulus
pattern. You can see that stimuli presented more often (going from the
first column to the last column) produce higher activations.

For more details please check out the vignette:

``` r
vignette("passt")
```

<!-- 
# Citation

Titz, J. (2019). passt: An R implementation of the PASS-T model. R package version 0.1.0. https://CRAN.R-project.org/package=passt

A BibTex entry for LaTeX users is:

@Manual{titz2019, title = {passt: an R implementation of the PASS-T model}, author = {Titz, Johannes}, year = {2019}, note = {R package version 0.1.0}, url = {https://CRAN.R-project.org/package=passt}, }
-->
