---
title: 'passt: An R implementation of the Probability Associator Time (PASS-T) model'
tags:
  - R
  - judgments of frequency
  - judgments of duration
  - PASS-T
  - artificial neural network
authors:
  - name: Johannes Titz
    orcid: 0000-0002-1102-5719
    affiliation: 1
affiliations:
 - name: Department of Psychology, TU Chemnitz, Germany
   index: 1
date: 12 November 2019
bibliography: library.bib
---

# Summary
The *passt* package is an R implementation of the Probability ASSociator Time (PASS-T) model, an artificial neural network designed to explain how humans make judgments of frequency and duration [@Titz2019]. The package was developed with two purposes in mind: (1) to provide a simple way to reproduce simulation results in judgments of frequency and duration as described in @Titz2019, and (2) to explore the PASS-T model by allowing users to run their own individual simulations.

The general idea of the original PASS model [@Sedlmeier1999;@Sedlmeier2002a] is that information about the frequency of events is naturally incorporated in artificial neural networks. If an object is presented repeatedly, the weights of the network change systematically. This way, the network is able to deduce the frequency of occurrence of events based on the final weights. Put simply, the artificial neural network produces a higher activation the more often a stimulus is presented.

As an extension of the PASS model, PASS-T is also able to process the presentation time of the stimuli. One can define how often and how long different objects are presented and test whether the network is able to make valid judgments of frequency and duration in retrospect. The principal aim of PASS-T is to explain empirical results that are found in studies with humans [@Titz2018]. Specifically, two empirical results are quite robust: (1) the correlation between exposure *frequency* and judgment is usually larger than between exposure *duration* and judgment (sometimes termed *frequency primacy*); and (2) the amount of attention that participants pay to the stimuli moderates effect sizes. These findings have been of some interest in cognitive psychology in the last decade [@Betsch2010;@Titz2019;@Titz2018;@Winkler2015;@Winkler2009a], although its roots go back half a century [@Hintzman1970a;@Hintzman1975;@Cooper1967] and touch upon different research areas such as memory, heuristics, and perception of magnitudes. PASS-T can be seen as a summary of some of this research, cast into a formal model. The R package *passt* is the first concrete implementation of PASS-T to make the model available to a wider audience.

# References
