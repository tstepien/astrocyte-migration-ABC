%% INVERSION: ASTROCYTE MIGRATION

%% 1 - INITIALIZE UQLAB
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

filename = 'parameter_analysis/inversion18.mat';

%%
% Load the measured population size stored in |Data|:
load('uq_Example_BayesianPreyPred.mat')

%% 2 - FORWARD MODEL
% Specify the forward model as a UQLab MODEL object:
ModelOpts.mFile = 'uq_eqns_and_error18';

myForwardModel = uq_createModel(ModelOpts);

%% 3 - PRIOR DISTRIBUTION OF THE MODEL PARAMETERS
% Specify prior distributions (given available information about the model 
% parameters before experimental observations) as a UQLab INPUT object:

PriorOpts.Marginals(1).Name = '$\mu$';  % adhesion constant
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [0.01 10];  % (mN h/mm^3)
PriorOpts.Marginals(1).Bounds = [0.01 10];  % (mN h/mm^3)

PriorOpts.Marginals(2).Name = '$\alpha_{10}$';  % base proliferation rate APC
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(2).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(3).Name = '$\alpha_{11}$';  % proliferation rate APC wrt PDGFA
PriorOpts.Marginals(3).Type = 'Uniform';
PriorOpts.Marginals(3).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(3).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(4).Name = '$\alpha_{12}$';  % proliferation rate APC wrt choroid oxygen
PriorOpts.Marginals(4).Type = 'Uniform';
PriorOpts.Marginals(4).Parameters = [0 1];  % (/hr)
PriorOpts.Marginals(4).Bounds = [0 1];  % (/hr)

PriorOpts.Marginals(5).Name = '$\alpha_{13}$';  % proliferation rate APC wrt hyaloid oxygen
PriorOpts.Marginals(5).Type = 'Uniform';
PriorOpts.Marginals(5).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(5).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(6).Name = '$\alpha_{20}$';  % base proliferation rate IPA
PriorOpts.Marginals(6).Type = 'Uniform';
PriorOpts.Marginals(6).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(6).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(7).Name = '$\alpha_{21}$';  % proliferation rate IPA wrt PDGFA
PriorOpts.Marginals(7).Type = 'Uniform';
PriorOpts.Marginals(7).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(7).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(8).Name = '$\alpha_{22}$';  % proliferation rate IPA wrt choroid oxygen
PriorOpts.Marginals(8).Type = 'Uniform';
PriorOpts.Marginals(8).Parameters = [0 1];  % (/hr)
PriorOpts.Marginals(8).Bounds = [0 1];  % (/hr)

PriorOpts.Marginals(9).Name = '$\alpha_{23}$';  % proliferation rate IPA wrt hyaloid oxygen
PriorOpts.Marginals(9).Type = 'Uniform';
PriorOpts.Marginals(9).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(9).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(10).Name = '$\beta_0$';  % base differentiation rate
PriorOpts.Marginals(10).Type = 'Uniform';
PriorOpts.Marginals(10).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(10).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(11).Name = '$\beta_1$';  % mass action rate
PriorOpts.Marginals(11).Type = 'Uniform';
PriorOpts.Marginals(11).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(11).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(12).Name = '$\beta_2$';  % differentiation rate wrt LIF
PriorOpts.Marginals(12).Type = 'Uniform';
PriorOpts.Marginals(12).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(12).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(13).Name = '$\beta_3$';  % differentiation rate wrt choroid oxygen
PriorOpts.Marginals(13).Type = 'Uniform';
PriorOpts.Marginals(13).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(13).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(14).Name = '$\beta_4$';  % differentiation rate wrt hyaloid oxygen
PriorOpts.Marginals(14).Type = 'Uniform';
PriorOpts.Marginals(14).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(14).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(15).Name = '$\eta_1$';  % apoptosis rate APC
PriorOpts.Marginals(15).Type = 'Uniform';
PriorOpts.Marginals(15).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(15).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(16).Name = '$\eta_2$';  % apoptosis rate IPA
PriorOpts.Marginals(16).Type = 'Uniform';
PriorOpts.Marginals(16).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(16).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(17).Name = '$P_\mathrm{hy}$';  % partial pressure of oxygen due to hyaloid artery
PriorOpts.Marginals(17).Type = 'Uniform';
PriorOpts.Marginals(17).Parameters = [0 20];  % (dimensionless)
PriorOpts.Marginals(17).Bounds = [0 20];  % (dimensionless)

PriorOpts.Marginals(18).Name = '$r_\mathrm{hy}$';  % radius at half-maximum of Hill function for hyaloid
PriorOpts.Marginals(18).Type = 'Uniform';
PriorOpts.Marginals(18).Parameters = [0.001 2];  % (mm)
PriorOpts.Marginals(18).Bounds = [0.001 2];  % (mm)

myPriorDist = uq_createInput(PriorOpts);

%% 4 - MEASUREMENT DATA
%
% Because the lynx and hare populations have different discrepancy options, 
% the measurement data is stored in two different data structures:  
% myData(1).y = Data.hare.'/1000; %in 1000
% myData(1).Name = 'Hare data';
% myData(1).MOMap = 1:21; % Output ID
% 
% myData(2).y = Data.lynx.'/1000; %in 1000
% myData(2).Name = 'Lynx data';
% myData(2).MOMap = 22:42; % Output ID

% ????????????????????

%% 5 - DISCREPANCY MODEL
% To infer the discrepancy variance, lognormal priors are put on the
% discrepancy parameters:
%
% * $\sigma^2_{\mathrm{APC}} \sim \mathcal{LN}(\lambda_{\sigma^2_\mathrm{APC}} 
%       = -1, \zeta_{\sigma^2_\mathrm{APC}} = 1)$
% * $\sigma^2_{\mathrm{IPA}} \sim \mathcal{LN}(\lambda_{\sigma^2_\mathrm{IPA}} 
%       = -1, \zeta_{\sigma^2_\mathrm{IPA}} = 1)$
%
% Specify these distributions in UQLab separately as two INPUT objects:
SigmaOpts.Marginals(1).Name = 'Sigma2APC';
SigmaOpts.Marginals(1).Type = 'Lognormal';
SigmaOpts.Marginals(1).Parameters = [-1 1];

SigmaDist1 = uq_createInput(SigmaOpts);

SigmaOpts.Marginals(1).Name = 'Sigma2IPA';
SigmaOpts.Marginals(1).Type = 'Lognormal';
SigmaOpts.Marginals(1).Parameters = [-1 1];

SigmaDist2 = uq_createInput(SigmaOpts);

% Assign these distributions to the discrepancy model options:
DiscrepancyOpts(1).Type = 'Gaussian';
DiscrepancyOpts(1).Prior = SigmaDist1;
DiscrepancyOpts(2).Type = 'Gaussian';
DiscrepancyOpts(2).Prior = SigmaDist2;
% Save variables to file
save(filename)

%% 6 - BAYESIAN ANALYSIS
%
%% 6.1 MCMC solver options
% To sample from the posterior distribution, the affine invariant ensemble
% algorithm is chosen, with $400$ iterations and $100$ parallel chains:
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.Steps = 400;
Solver.MCMC.NChains = 100;

% Enable progress visualization during iteration for specific parameters.
% Update the plots every $40$ iterations:
Solver.MCMC.Visualize.Parameters = [1 2];
Solver.MCMC.Visualize.Interval = 40;

%% 6.2 Posterior sample generation
% The options of the Bayesian analysis are gathered within a single
% structure with fields: 
BayesOpts.Type = 'Inversion';
BayesOpts.Name = 'Bayesian model';
BayesOpts.Prior = myPriorDist;
%BayesOpts.Data = myData;
BayesOpts.Discrepancy = DiscrepancyOpts;
BayesOpts.Solver = Solver;

% Perform and store in UQLab the Bayesian inversion analysis:
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

% Print out a report of the results:
uq_print(myBayesianAnalysis)

%% 6.3 Posterior sample post-processing
% Diagnose the quality of the results,
% create a trace plot of the first parameter:
uq_display(myBayesianAnalysis, 'trace', 1)

% From the plots, one can see that several chains have not converged yet.
% From the trace plot, the non-converged chains are all characterized by a
% final value $x_1^{(T)}>0.8$:
badChainsIndex = squeeze(myBayesianAnalysis.Results.Sample(end,1,:) > 0.8);

% These chains can be removed from the sample through post-processing. 
% Additionally, draw a sample of size $10^3$ from the prior and posterior
% predictive distributions: 
uq_postProcessInversionMCMC(myBayesianAnalysis,...
                        'badChains', badChainsIndex,...
                        'prior', 1000,...
                        'priorPredictive', 1000,...
                        'posteriorPredictive', 1000);

% *Note*: sampling prior predictive samples requires new
% model evaluations

% Display the post processed results:
uq_display(myBayesianAnalysis)

% Save variables to file
save(filename)