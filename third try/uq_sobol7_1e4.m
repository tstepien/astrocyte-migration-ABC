%% SENSITIVITY: ASTROCYTE MIGRATION
%% 1 - INITIALIZE UQLAB
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clc;
clearvars
rng(100,'twister')
uqlab

num_param = 7;
sampleN = 1e4;

filename = strcat('parameter_analysis/sobol',num2str(num_param),'_1e4.mat');
load(strcat('plot_uq/distributions',num2str(num_param),'.mat'));

%% 2 - COMPUTATIONAL MODEL
% Create a MODEL object from the function file:
ModelOpts.mFile = strcat('uq_eqns_and_error',num2str(num_param));

myModel = uq_createModel(ModelOpts);

%% 3 - PROBABILISTIC INPUT MODEL
% The probabilistic input model consists of 18 independent random variables.
% Specify the marginals as follows:
InputOpts.Marginals(1).Name = '$\mu$';  % adhesion constant
InputOpts.Marginals(1).Type = 'Gamma';
InputOpts.Marginals(1).Parameters = [bestfitdist_param{1}.a bestfitdist_param{1}.b];  % (mN h/mm^3)
InputOpts.Marginals(1).Bounds = [0.0001 0.1];  % (mN h/mm^3)

InputOpts.Marginals(2).Name = '$\alpha_{11}$';  % proliferation rate APC wrt PDGFA
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(2).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(3).Name = '$\alpha_{12}$';  % proliferation rate APC wrt choroid oxygen
InputOpts.Marginals(3).Type = 'Weibull';
InputOpts.Marginals(3).Parameters = [bestfitdist_param{3}.A bestfitdist_param{3}.B];  % (/hr)
InputOpts.Marginals(3).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(4).Name = '$\alpha_{21}$';  % proliferation rate IPA wrt PDGFA
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(4).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(5).Name = '$\beta_1$';  % differentiation rate wrt LIF
InputOpts.Marginals(5).Type = 'Uniform';
InputOpts.Marginals(5).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(5).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(6).Name = '$\beta_2$';  % differentiation rate wrt choroid oxygen
InputOpts.Marginals(6).Type = 'Gamma';
InputOpts.Marginals(6).Parameters = [bestfitdist_param{6}.a bestfitdist_param{6}.b];  % (/hr)
InputOpts.Marginals(6).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(7).Name = '$\beta_4$';  % mass action rate
InputOpts.Marginals(7).Type = 'Uniform';
InputOpts.Marginals(7).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(7).Bounds = [0 2];  % (/hr)

% Create an INPUT object based on the specified marginals:
myInput = uq_createInput(InputOpts);

%% 4 - SENSITIVITY ANALYSIS
% Sensitivity analysis is performed with the following methods:
%
% * Sobol' sensitivity indices

%% 4.6 Sobol' indices
% Select the sensitivity tool and the Sobol' method:
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';

% Specify the maximum order of the Sobol' indices calculation:
SobolOpts.Sobol.Order = 1;

% Specify the sample size for the indices estimation of each variable
SobolOpts.Sobol.SampleSize = sampleN;

% Note that the total cost of computation is $(M+2)\times N$,
% where $M$ is the input dimension (number of parameters) and $N$ is the 
% sample size.
% Therefore, if M=8 and N=sampleN=10^4, the total cost for the current 
% setup is $(8+2)\times 10^4 = 10^5$ evaluations of the full computational 
% model.

% Add boostrap-based confidence intervals (0.025 and 0.975 quantiles)
SobolSensOpts.Bootstrap.Replications = 100;
SobolSensOpts.Alpha = 0.05;

% Run the sensitivity analysis:
SobolAnalysis = uq_createAnalysis(SobolOpts);

% Print the results of the analysis:
uq_print(SobolAnalysis)

% Create a graphical representation of the results:
if ~isdeployed
    uq_display(SobolAnalysis)
end

% Save variables to file
save(filename)