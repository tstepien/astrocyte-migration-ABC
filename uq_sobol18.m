%% SENSITIVITY: ASTROCYTE MIGRATION
%% 1 - INITIALIZE UQLAB
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clc;
clearvars
rng(100,'twister')
uqlab

sampleN = 1e5;
filename = '../parameter_analysis/sobol18.mat';
load('distributions18.mat');

%% 2 - COMPUTATIONAL MODEL
% Create a MODEL object from the function file:
ModelOpts.mFile = 'uq_eqns_and_error18';

myModel = uq_createModel(ModelOpts);

%% 3 - PROBABILISTIC INPUT MODEL
% The probabilistic input model consists of 18 independent random variables.
% Specify the marginals as follows:
InputOpts.Marginals(1).Name = '$\mu$';  % adhesion constant
InputOpts.Marginals(1).Type = 'Gaussian';
InputOpts.Marginals(1).Parameters = [bestfitdist_param{1}.mu bestfitdist_param{1}.sigma];  % (mN h/mm^3)
InputOpts.Marginals(1).Bounds = [0.01 10];  % (mN h/mm^3)

InputOpts.Marginals(2).Name = '$\alpha_{10}$';  % base proliferation rate APC
InputOpts.Marginals(2).Type = 'Gaussian';
InputOpts.Marginals(2).Parameters = [bestfitdist_param{2}.mu bestfitdist_param{2}.sigma];  % (/hr)
InputOpts.Marginals(2).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(3).Name = '$\alpha_{11}$';  % proliferation rate APC wrt PDGFA
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(3).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(4).Name = '$\alpha_{12}$';  % proliferation rate APC wrt choroid oxygen
InputOpts.Marginals(4).Type = 'Gaussian';
InputOpts.Marginals(4).Parameters = [bestfitdist_param{4}.mu bestfitdist_param{4}.sigma];  % (/hr)
InputOpts.Marginals(4).Bounds = [0 1];  % (/hr)

InputOpts.Marginals(5).Name = '$\alpha_{13}$';  % proliferation rate APC wrt hyaloid oxygen
InputOpts.Marginals(5).Type = 'Weibull';
InputOpts.Marginals(5).Parameters = [bestfitdist_param{5}.A bestfitdist_param{5}.B];  % (/hr)
InputOpts.Marginals(5).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(6).Name = '$\alpha_{20}$';  % base proliferation rate IPA
InputOpts.Marginals(6).Type = 'Uniform';
InputOpts.Marginals(6).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(6).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(7).Name = '$\alpha_{21}$';  % proliferation rate IPA wrt PDGFA
InputOpts.Marginals(7).Type = 'Uniform';
InputOpts.Marginals(7).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(7).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(8).Name = '$\alpha_{22}$';  % proliferation rate IPA wrt choroid oxygen
InputOpts.Marginals(8).Type = 'Uniform';
InputOpts.Marginals(8).Parameters = [0 1];  % (/hr)
InputOpts.Marginals(8).Bounds = [0 1];  % (/hr)

InputOpts.Marginals(9).Name = '$\alpha_{23}$';  % proliferation rate IPA wrt hyaloid oxygen
InputOpts.Marginals(9).Type = 'Uniform';
InputOpts.Marginals(9).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(9).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(10).Name = '$\beta_0$';  % base differentiation rate
InputOpts.Marginals(10).Type = 'Weibull';
InputOpts.Marginals(10).Parameters = [bestfitdist_param{10}.A bestfitdist_param{10}.B];  % (/hr)
InputOpts.Marginals(10).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(11).Name = '$\beta_1$';  % mass action rate
InputOpts.Marginals(11).Type = 'Uniform';
InputOpts.Marginals(11).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(11).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(12).Name = '$\beta_2$';  % differentiation rate wrt LIF
InputOpts.Marginals(12).Type = 'Weibull';
InputOpts.Marginals(12).Parameters = [bestfitdist_param{12}.A bestfitdist_param{12}.B];  % (/hr)
InputOpts.Marginals(12).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(13).Name = '$\beta_3$';  % differentiation rate wrt choroid oxygen
InputOpts.Marginals(13).Type = 'Weibull';
InputOpts.Marginals(13).Parameters = [bestfitdist_param{13}.A bestfitdist_param{13}.B];  % (/hr)
InputOpts.Marginals(13).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(14).Name = '$\beta_4$';  % differentiation rate wrt hyaloid oxygen
InputOpts.Marginals(14).Type = 'Uniform';
InputOpts.Marginals(14).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(14).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(15).Name = '$\eta_1$';  % apoptosis rate APC
InputOpts.Marginals(15).Type = 'Gaussian';
InputOpts.Marginals(15).Parameters = [bestfitdist_param{15}.mu bestfitdist_param{15}.sigma];  % (/hr)
InputOpts.Marginals(15).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(16).Name = '$\eta_2$';  % apoptosis rate IPA
InputOpts.Marginals(16).Type = 'Uniform';
InputOpts.Marginals(16).Parameters = [0 2];  % (/hr)
InputOpts.Marginals(16).Bounds = [0 2];  % (/hr)

InputOpts.Marginals(17).Name = '$P_\mathrm{hy}$';  % partial pressure of oxygen due to hyaloid artery
InputOpts.Marginals(17).Type = 'Weibull';
InputOpts.Marginals(17).Parameters = [bestfitdist_param{17}.A bestfitdist_param{17}.B];  % (dimensionless)
InputOpts.Marginals(17).Bounds = [0 20];  % (dimensionless)

InputOpts.Marginals(18).Name = '$r_\mathrm{hy}$';  % radius at half-maximum of Hill function for hyaloid
InputOpts.Marginals(18).Type = 'Uniform';
InputOpts.Marginals(18).Parameters = [0.001 2];  % (mm)
InputOpts.Marginals(18).Bounds = [0.001 2];  % (mm)

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
uq_display(SobolAnalysis)

% Save variables to file
save(filename)