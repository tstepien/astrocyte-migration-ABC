# astrocyte-migration-ABC

<a href="https://github.com/tstepien/astrocyte-migration-ABC/"><img src="https://img.shields.io/badge/github-tstepien%2Fastrocyte--migration--ABC-blue" /></a> <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" /></a>

The MATLAB code contained in the astrocyte-migration-ABC project was developed to numerically simulate astrocyte migration and differentiation during embryonic development of the retina and its corresponding approximate Bayesian computation (ABC) analysis for parameter estimation and model selection. It is described in:
>[Tracy L. Stepien](https://github.com/tstepien/), An approximate Bayesian computation approach for embryonic astrocyte migration model selection, Submitted to ---

## Programs
+ Files beginning with ``program_``: run these programs to solve the astrocyte model equations
+ Files beginning with ``uq_abc``: run these programs to run approximate Bayesian computation analysis

The number in the file corresponds with the number of nonzero parameters in the model.

### Additional MATLAB Packages
The following MATLAB packages by others were used:
+ [latin_random.m](https://people.sc.fsu.edu/~jburkardt/m_src/latin_random/latin_random.html): Latin hypercube sampling
+ [plot_uq/ws_distance.m](https://github.com/nklb/wasserstein-distance): Computes 1- and 2- Wasserstein distance between two discrete probability measures 

## Description of Folders
+ [ABC_results](ABC_results): Contains .mat files with the parameter sets and corresponding errors (comparing the model to data)
+ [LHpts](LHpts): Contains .mat files with Latin hypercube points (scaled between 0 and 1) for different numbers of model parameters
+ [plot_simulations](plot_simulations): Code to create Figure 1 in the paper
+ [plot_uq](plot_uq): Code to create Figures 2-6 and Table 3 in the paper

## Licensing
Copyright 2017-2024 [Tracy Stepien](https://github.com/tstepien/). This is free software made available under the MIT License. For details see the [LICENSE](LICENSE) file.
