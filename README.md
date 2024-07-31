# astrocyte-migration-ABC

<a href="https://github.com/tstepien/astrocyte-migration-ABC/"><img src="https://img.shields.io/badge/github-tstepien%2Fastrocyte--migration--ABC-blue" /></a> <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" /></a>

The MATLAB code contained in the astrocyte-migration-ABC project was developed to numerically simulate astrocyte migration and differentiation during embryonic development of the retina and its corresponding approximate Bayesian computation (ABC) analysis for parameter estimation and model reduction. It is described in:
>[Tracy L. Stepien](https://github.com/tstepien/), An approximate Bayesian computation approach for embryonic astrocyte migration model reduction, Submitted to the _Bulletin of Mathematical Biology_.

## Programs
+ Files beginning with ``program_``: Run these programs to solve the astrocyte model equations
+ Files beginning with ``uq_``: Run these programs to run approximate Bayesian computation
+ Files beginning with ``plot_uq/plot_abc``: Run these programs to analyze the results of the approximate Bayesian computation
+ [plot_uq/analyze_modelselect.m](plot_uq/analyze_modelselect.m): Run this program to compute the Bayes factors to compare models

The number in the file corresponds with the number of nonzero parameters in the model.

The other files in the main directory are subfunctions necessary for the simulations.

## Description of Folders
+ [parameter_analysis](parameter_analysis): Code to run gradient descent (fminsearch in MATLAB) on the 5 smallest parameter sets found using approximate Bayesian computation
+ [plot_simulations](plot_simulations): Code to create Figure 2 in the paper
+ [plot_uq](plot_uq): Code to create Figures 3-7 and Table 4 in the paper

### Folders with MATLAB data hosted elsewhere
The following two folders are available for download on OSF at [https://osf.io/rm5qj/](https://osf.io/rm5qj/)

+ *ABC_results*: Contains .mat files with the parameter sets for running approximate Bayesian computation and the corresponding errors (comparing the model to data)
+ *LHpts*: Contains .mat files with Latin hypercube points (scaled between 0 and 1) for different numbers of model parameters

## Additional MATLAB Packages
The following MATLAB packages were used:
+ [latin_random](https://people.sc.fsu.edu/~jburkardt/m_src/latin_random/latin_random.html): Latin hypercube sampling
+ [wasserstein-distance](https://github.com/nklb/wasserstein-distance): Computes 1- and 2- Wasserstein distance (Earth mover's distance) between two discrete probability measures 

## Licensing
Copyright 2017-2024 [Tracy Stepien](https://github.com/tstepien/). This is free software made available under the MIT License. For details see the [LICENSE](LICENSE) file.
