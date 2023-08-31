clear variables global;
clc;
close all;
addpath plot_simulations

%%% Note: the subfunction oxygen_michmen uses a parfor to speed up
%%% calculations

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 15; %%% adhesion constant (Pa h/mm)
alpha10 = 0.08; %%% (/hr) basal proliferation rate APC
alpha11 = 0.09; %%% (/hr) proliferation rate APC wrt PDGFA
alpha12 = 0.08; %%% (/hr) proliferation rate APC wrt choroid oxygen
alpha13 = 0; %%% (/hr) proliferation rate APC wrt hyaloid oxygen
alpha20 = 0; %%% (/hr) basal proliferation rate IPA
alpha21 = 0.01; %%% (/hr) proliferation rate IPA wrt PDGFA
alpha22 = 0.005; %%% (/hr) proliferation rate IPA wrt choroid oxygen
alpha23 = 0; %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
beta0 = 0.07; %%% (/hr) basal differentiation rate
beta1 = 0.02; %%% (/hr) differentiation rate wrt LIF
beta2 = 0.03; %%% (/hr) differentiation rate wrt choroid oxygen
beta3 = 0; %%% (/hr) differentiation rate wrt hyaloid oxygen
beta4 = 0; %%% (/hr) mass action rate
eta1 = 0; %%% (/hr) apoptosis rate APC
eta2 = 0; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyaloid artery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_hy = 0; %%% partial pressure of oxygen due to hyaloid artery
r_hy = 1; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%% loading a parameter set file %%%%%%%%%%%%%%%%%%%%%%%
% load('uq/parameter_analysis/problem_parameter_set12726.mat');
% mu = param_ind(1);
% alpha10 = param_ind(2);
% alpha11 = param_ind(3);
% alpha12 = param_ind(4);
% alpha13 = param_ind(5);
% alpha20 = param_ind(6);
% alpha21 = param_ind(7);
% alpha22 = param_ind(8);
% alpha23 = param_ind(9);
% beta0 = param_ind(10);
% beta1 = param_ind(11);
% beta2 = param_ind(12);
% beta3 = param_ind(13);
% beta4 = param_ind(14);
% eta1 = param_ind(15);
% eta2 = param_ind(16);
% P_hy = param_ind(17);
% r_hy = param_ind(18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad,choroidPO2] = eqnsolver(mu,alpha10,...
    alpha11,alpha12,alpha13,alpha20,alpha21,alpha22,alpha23,beta0,beta1,...
    beta2,beta3,beta4,eta1,eta2,P_hy,r_hy,m);
toc


%% error
%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_tot,err_time,err_rad,err_dens] = errorfunction(t,r,mvgbdy,c1,c2);


%% area under curve - trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%% testing: conservation of mass %%%%%%%%%%%%%%%%%%%%%%
areaundercurve = zeros(length(t),1);
dr = m.dr;
k = c1+c2;
for i = 1:length(t)
    areaundercurve(i) = dr*sum( k(1:mvgbdy(1)/dr+(i-1)) ...
        + k(2:mvgbdy(1)/dr+1+(i-1)) ) / 2;
end

%% plots
plot_the_plots
% plot_the_plots_APCIPA
% plot_movingbdy