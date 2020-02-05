
<!-- README.md is generated from README.Rmd. Please edit that file -->

# passt

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/johannes-titz/passt.svg?branch=master)](https://travis-ci.org/johannes-titz/passt)
[![Codecov test
coverage](https://codecov.io/gh/johannes-titz/passt/branch/master/graph/badge.svg)](https://codecov.io/gh/johannes-titz/passt?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/passt)](https://CRAN.R-project.org/package=passt)
<!-- badges: end -->

The *passt* package is an R implementation of the Probability ASSociator
Time (PASS-T) model, an artificial neural network designed to explain
how humans make judgments of frequency and duration (Titz & Sedlmeier,
2019). The package was developed with two purposes in mind: (1) to
provide a simple way to reproduce simulation results in judgments of
frequency and duration as described in Titz & Sedlmeier (2019), and (2)
to explore the PASS-T model by allowing users to run their own
individual simulations. The package is targeted at cognitive
psychologists who study judgments of frequency and duration, memory,
heuristics and artificial neural networks.

The general idea of the original PASS model (Sedlmeier, 1999, 2002) is
that information about the frequency of events is naturally incorporated
in artificial neural networks. If an object is presented repeatedly, the
weights of the network change systematically. This way, the network is
able to deduce the frequency of occurrence of events based on the final
weights. Put simply, the artificial neural network produces a higher
activation the more often a stimulus is presented.

As an extension of the PASS model, PASS-T is also able to process the
presentation time of the stimuli. One can define how often and how long
different objects are presented and test whether the network is able to
make valid judgments of frequency and duration in retrospect. The
principal aim of PASS-T is to explain empirical results that are found
in studies with humans (Titz et al., 2019). Specifically, two empirical
results are quite robust: (1) the correlation between exposure
*frequency* and judgment is usually larger than between exposure
*duration* and judgment (sometimes termed *frequency primacy*); and (2)
the amount of attention that participants pay to the stimuli moderates
effect sizes. These findings have been of some interest in cognitive
psychology in the last decade (Betsch, Glauer, Renkewitz, Winkler, &
Sedlmeier, 2010; Titz et al., 2019; Titz & Sedlmeier, 2019; Winkler,
2009; Winkler, Glauer, Betsch, & Sedlmeier, 2015), although its roots go
back half a century (Cooper & Pantle, 1967; Hintzman, 1970; Hintzman,
Summers, & Block, 1975) and touch upon different research areas such as
memory, heuristics, and perception of magnitudes. PASS-T can be seen as
a summary of some of this research, cast into a formal model. The R
package *passt* is the first concrete implementation of PASS-T to make
the model available to a wider audience.

## Installation

You can install the released version of passt from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("passt")
```

You can install the development version from
[GitHub](https://github.com/) with devtools (and vignettes build, this
takes a couple of seconds):

``` r
devtools::install_github("johannes-titz/passt", build_vignettes = TRUE)
```

## Using *passt*

There are only two functions of *passt* that you will need to use. The
first function is *run\_sim*, which runs several simulations with
specific parameters and returns the final output activation for each
input pattern. The second is *run\_exp*, which aggregates the data and
gives effect sizes for each simulation run. Let us first look at
*run\_sim*.

Since PASS-T extends the PASS model, PASS-T should be able to produce
all the results that PASS can produce. The simplest PASS simulation is a
type of counter: an artificial neural network sensitive only to
frequency information.

### A simple counter

``` r
library(passt)
set.seed(20191015)
```

``` r
sim1 <- run_sim(patterns = diag(10), frequency = 1:10, duration = 10:1,
                lrate_onset = 0.05, lrate_drop_time = 2, lrate_drop_perc = 0)
```

Some explanation is necessary: *patterns* is a matrix of input patterns
that the network will process. Here we will use orthogonal stimuli to
avoid any interference on the input level. The function *diag* will
create an identity matrix so that every input pattern has only one
active unit. This active unit will be unique for each pattern. You can
also use other stimulus patterns as long as the activation for each unit
is either 0 or 1.

Next, we specify how often the ten patterns should be shown and the
duration of each presentation. The first pattern will be shown once for
10 s, the second pattern will be shown twice for 9 s each (a total
duration of 18 s), and so on. You can also use other arrangements for a
simple counter.

Finally, we set the learning rate parameters: the initial learning rate
(*lrate\_onset*) is 0.05, but after 2 s (*lrate\_drop\_time*) it will
drop to 0 (*lrate\_drop\_perc*). This parameter selection is crucial for
a counter because only the first two seconds of a stimulus presentation
are processed by the network.

The function *run\_sim* returns a list with three elements, of which
*output* is the most interesting. Note that, by default, 100 simulations
are run (you can change this with the parameter *n\_runs*). We will only
look at the first couple of simulations:

``` r
head(sim1$output)
#>           [,1]      [,2]      [,3]      [,4]      [,5]     [,6]     [,7]
#> [1,] 0.8134030 0.8623727 0.9090286 0.9526607 0.9893114 1.027876 1.065307
#> [2,] 0.8161057 0.8618859 0.9028135 0.9479106 0.9905695 1.024340 1.062587
#> [3,] 0.8117953 0.8602871 0.9049667 0.9490695 0.9895403 1.027940 1.061659
#> [4,] 0.8148240 0.8633805 0.9009140 0.9490142 0.9884823 1.026991 1.061192
#> [5,] 0.8114915 0.8611520 0.9011725 0.9449056 0.9886680 1.027973 1.067783
#> [6,] 0.8177491 0.8598771 0.9061970 0.9493045 0.9831946 1.031035 1.059464
#>          [,8]     [,9]    [,10]
#> [1,] 1.093939 1.125975 1.160127
#> [2,] 1.098260 1.132112 1.163416
#> [3,] 1.097431 1.134466 1.162845
#> [4,] 1.100671 1.128499 1.166032
#> [5,] 1.103455 1.132242 1.161158
#> [6,] 1.102026 1.132638 1.158515
```

The *output* gives the sum of the activity of all output units for the
specific patterns for each simulation run. Each row corresponds with one
simulation and each column with one input pattern. Sometimes the output
activation is referred to as the *memory strength* in an abstract sense.
One can easily see that increasing the presentation frequency
(i.e. going from the first to the tenth column) results in a higher
*memory strength* (activation).

### Simulating the moderating role of attention

The PASS-T model was used to explain how attention moderates effect
sizes in judgments of frequency and duration (Titz & Sedlmeier, 2019). A
short example for a setup could look like this:

``` r
duration <- c(4, 2, 1, 8, 4, 2, 12, 6, 3)
frequency <- c(2, 4, 8, 2, 4, 8, 2, 4, 8)
lrate_drop_perc <- seq(0, 1, 0.04)
sim4 <- lapply(lrate_drop_perc, function(x) 
  run_exp(frequency, duration, 0.05, 2, x, diag(9), 30, 0.1))
```

The *lapply* function goes through each *lrate\_drop\_perc* value and
runs a simulation. Note that we reduced the number of simulation runs to
30 to keep the computation time within reasonable limits. Since we now
have a list of simulations, we need to convert it into a data frame and
then add the information about the drop of the learning parameter:

``` r
sim4 <- plyr::ldply(sim4, "data.frame")
sim4 <- cbind(sim4, lrate_drop_perc)
sim4
#>            f_dv      td_dv        d_dv lrate_drop_perc
#> 1   0.701028160 0.01351331 -0.47651489            0.00
#> 2   0.663278975 0.07333753 -0.39285858            0.04
#> 3   0.650119216 0.13604230 -0.37299350            0.08
#> 4   0.633004783 0.28672642 -0.30706420            0.12
#> 5   0.556395109 0.31918010 -0.19419062            0.16
#> 6   0.551342130 0.37495150 -0.18118624            0.20
#> 7   0.481881025 0.42253064 -0.11725408            0.24
#> 8   0.416859327 0.48672492  0.02817775            0.28
#> 9   0.442182027 0.54100892  0.02780947            0.32
#> 10  0.352404794 0.55631562  0.07249250            0.36
#> 11  0.334024242 0.60051914  0.13500556            0.40
#> 12  0.295262006 0.62575890  0.19079468            0.44
#> 13  0.307331334 0.70923977  0.21864830            0.48
#> 14  0.271276504 0.66693212  0.20536546            0.52
#> 15  0.191859228 0.70464540  0.26520307            0.56
#> 16  0.188205704 0.73051860  0.28517694            0.60
#> 17  0.172752955 0.74948246  0.34016942            0.64
#> 18  0.100509302 0.76208425  0.38346817            0.68
#> 19  0.129906979 0.77737432  0.34669211            0.72
#> 20  0.151104550 0.76956575  0.33782937            0.76
#> 21  0.073552426 0.81581281  0.42232720            0.80
#> 22  0.030263447 0.79239836  0.46008706            0.84
#> 23  0.023810247 0.81660532  0.44032598            0.88
#> 24  0.005663264 0.78959218  0.43559317            0.92
#> 25  0.064438336 0.83468908  0.44817420            0.96
#> 26 -0.045779705 0.82140950  0.51397210            1.00
```

If the learning rate does not drop much (i.e. *lrate\_drop\_perc*
approaches 1), the influence of total duration (*td\_dv*) on the output
activation is strong, while the influence of frequency (*f\_dv*) is
weak. But if the learning rate drops substantially
(i.e. *lrate\_drop\_perc* approaches 0), the opposite is true.

For more details please check out the vignette in R:

``` r
vignette("passt")
```

Or on CRAN:
<https://cran.r-project.org/web/packages/passt/vignettes/passt.html>

## Issues and Support

If you find any bugs, please use the issue tracker at:

<https://github.com/johannes-titz/passt/issues>

If you need assistance in how to use the package, drop me an e-mail at
johannes at titz.science or johannes.titz at gmail.com

## Contributing

Contributions of any kind are very welcome\! I will sincerely consider
every suggestion on how to improve the code, the documentation and the
presented examples. Even minor things, such as suggestions for better
wording or improving grammar in any part of the package are considered
as valuable contributions.

If you want to make a pull request, please check that you can still
build the package without any errors, warnings or notes. Overall, simply
stick to the R packages book: <https://r-pkgs.org/> and follow the code
style described here: <http://r-pkgs.had.co.nz/r.html#style>

<!-- 
# Citation

Titz, J. (2019). passt: An R implementation of the PASS-T model. R package version 0.1.0. https://CRAN.R-project.org/package=passt

A BibTex entry for LaTeX users is:

@Manual{titz2019, title = {passt: an R implementation of the PASS-T model}, author = {Titz, Johannes}, year = {2019}, note = {R package version 0.1.0}, url = {https://CRAN.R-project.org/package=passt}, }
-->

## References

<div id="refs" class="references">

<div id="ref-Betsch2010">

Betsch, T., Glauer, M., Renkewitz, F., Winkler, I., & Sedlmeier, P.
(2010). Encoding, storage and judgment of experienced frequency and
duration. *Judgment and Decision Making*, *5*, 347–364.

</div>

<div id="ref-Cooper1967">

Cooper, E. H., & Pantle, A. J. (1967). The total-time hypothesis in
verbal learning. *Psychological Bulletin*, *68*, 221–234.
<https://doi.org/10.1037/h0025052>

</div>

<div id="ref-Hintzman1970a">

Hintzman, D. L. (1970). Effects of repetition and exposure duration on
memory. *Journal of Experimental Psychology*, *83*, 435–444.
<https://doi.org/10.1037/h0028865>

</div>

<div id="ref-Hintzman1975">

Hintzman, D. L., Summers, J. J., & Block, R. A. (1975). What causes the
spacing effect? Some effects of repetition, duration, and spacing on
memory for pictures. *Memory & Cognition*, *3*, 287–294. Retrieved from
<https://link.springer.com/article/10.3758/BF03212913>

</div>

<div id="ref-Sedlmeier1999">

Sedlmeier, P. (1999). *Improving statistical reasoning: Theoretical
models and practical implications*. Mawah, NJ: Erlbaum.

</div>

<div id="ref-Sedlmeier2002a">

Sedlmeier, P. (2002). Associative learning and frequency judgments: The
PASS model. In P. Sedlmeier & T. Betsch (Eds.), *Frequency processing
and cognition* (pp. 137–152). Oxford, England: Oxford University Press.

</div>

<div id="ref-Titz2018">

Titz, J., Burkhardt, M., & Sedlmeier, P. (2019). *Asymmetries in
judgments of frequency and duration: A meta-analysis*. Manuscript
submitted for publication.

</div>

<div id="ref-Titz2019">

Titz, J., & Sedlmeier, P. (2019). *Simulating asymmetries in judgments
of frequency and duration: The PASS-T model*. Manuscript submitted for
publication.

</div>

<div id="ref-Winkler2009a">

Winkler, I. (2009). *The processing of frequency and duration* (Doctoral
dissertation, Chemnitz University of Technology). Retrieved from
<https://nbn-resolving.org/urn:nbn:de:bsz:ch1-200900917>

</div>

<div id="ref-Winkler2015">

Winkler, I., Glauer, M., Betsch, T., & Sedlmeier, P. (2015). The impact
of attention on judgments of frequency and duration. *PLoS ONE*, *10*,
e0126974. <https://doi.org/10.1371/journal.pone.0126974>

</div>

</div>
