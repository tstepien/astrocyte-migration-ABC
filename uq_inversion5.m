%% INVERSION: ASTROCYTE MIGRATION

%% 1 - INITIALIZE UQLAB
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

num_param = 5;
num_steps = 5000;
num_chains = 100;

savefiles = 'yes';

% if strcmp(savefiles,'yes')==1
%     doublecheck = input('Are you sure you would like to save the output files? (it may overwrite): ');
%     if strcmp(doublecheck,'y')==1
         diary(strcat('parameter_analysis/diary_inversion',num2str(num_param),...
             '_',num2str(num_chains),'chains_',num2str(num_steps),'steps.txt'));
         filename = strcat('parameter_analysis/inversion',num2str(num_param),...
             '_',num2str(num_chains),'chains_',num2str(num_steps),'steps.mat');
%     else
%         return;
%     end
% end

%% 2 - FORWARD MODEL
% Specify the forward model as a UQLab MODEL object:
ModelOpts.mFile = strcat('uq_eqns_and_error',num2str(num_param));

myForwardModel = uq_createModel(ModelOpts);

%% 3 - PRIOR DISTRIBUTION OF THE MODEL PARAMETERS
% Specify prior distributions (given available information about the model 
% parameters before experimental observations) as a UQLab INPUT object:

PriorOpts.Marginals(1).Name = '$\mu$';  % adhesion constant
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [0.0001 0.1];  % (mN h/mm^3)
PriorOpts.Marginals(1).Bounds = [0.0001 0.1];  % (mN h/mm^3)

PriorOpts.Marginals(2).Name = '$\alpha_{11}$';  % proliferation rate APC wrt PDGFA
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(2).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(3).Name = '$\alpha_{21}$';  % proliferation rate IPA wrt PDGFA
PriorOpts.Marginals(3).Type = 'Uniform';
PriorOpts.Marginals(3).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(3).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(4).Name = '$\beta_1$';  % mass action rate
PriorOpts.Marginals(4).Type = 'Uniform';
PriorOpts.Marginals(4).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(4).Bounds = [0 2];  % (/hr)

PriorOpts.Marginals(5).Name = '$\beta_4$';  % differentiation rate wrt hyaloid oxygen
PriorOpts.Marginals(5).Type = 'Uniform';
PriorOpts.Marginals(5).Parameters = [0 2];  % (/hr)
PriorOpts.Marginals(5).Bounds = [0 2];  % (/hr)

myPriorDist = uq_createInput(PriorOpts);

%% 4 - MEASUREMENT DATA

myData.y = 0;
myData.Name = 'Minimal error is zero';

%% 5 - DISCREPANCY MODEL
% To infer the discrepancy variance, lognormal priors are put on the
% discrepancy parameters:
%
% * $\sigma^2_{\mathrm{error}} \sim \mathcal{LN}(\lambda_{\sigma^2_\mathrm{error}} 
%       = -1, \zeta_{\sigma^2_\mathrm{error}} = 1)$
%
SigmaOpts.Marginals(1).Name = 'Sigma2Error';
SigmaOpts.Marginals(1).Type = 'Lognormal';
SigmaOpts.Marginals(1).Parameters = [-1 1];

SigmaDist = uq_createInput(SigmaOpts);

% Assign these distributions to the discrepancy model options:
DiscrepancyOpts(1).Type = 'Gaussian';
DiscrepancyOpts(1).Prior = SigmaDist;

% Save variables to file
if strcmp(savefiles,'yes')==1
    save(filename)
end

%% 6 - BAYESIAN ANALYSIS
%
%% 6.1 MCMC solver options
% To sample from the posterior distribution, the affine invariant ensemble
% algorithm is chosen, with $400$ iterations and $100$ parallel chains:
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.Steps = num_steps;
Solver.MCMC.NChains = num_chains;

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
BayesOpts.Data = myData;
BayesOpts.Discrepancy = DiscrepancyOpts;
BayesOpts.Solver = Solver;

% Perform and store in UQLab the Bayesian inversion analysis:
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

% Print out a report of the results:
uq_print(myBayesianAnalysis)

%% 6.3 Posterior sample post-processing
% Diagnose the quality of the results,
% create a trace plot of the first parameter:
if ~isdeployed
    for i=1:num_param
        uq_display(myBayesianAnalysis, 'trace', i)
    end
end

% From the plots, one can see that several chains have not converged yet.
% From the trace plot, the non-converged chains are all characterized by a
% final value $x_1^{(T)}>0.8$:
%badChainsIndex = squeeze(myBayesianAnalysis.Results.Sample(end,1,:) > 0.8);

% These chains can be removed from the sample through post-processing. 
% Additionally, draw a sample of size $10^3$ from the prior and posterior
% predictive distributions: 
% uq_postProcessInversionMCMC(myBayesianAnalysis,...
%                         'badChains', badChainsIndex,...
%                         'prior', 1000,...
%                         'priorPredictive', 1000,...
%                         'posteriorPredictive', 1000);

% *Note*: sampling prior predictive samples requires new
% model evaluations

% Display the post processed results:
if ~isdeployed
    uq_display(myBayesianAnalysis)
end

% Save variables to file
if strcmp(savefiles,'yes')==1
    save(filename)
    diary off
    beep
end