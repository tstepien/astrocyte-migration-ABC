clear variables global;
clc;
close all;
addpath plot_simulations

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%% loading a parameter set file %%%%%%%%%%%%%%%%%%%%%%%
load('uq/param_inversion13.mat');

mu = modeParam(1)/10^3; %%% adhesion constant
alpha10 = modeParam(2); %%% (/hr) basal proliferation rate APC
alpha11 = modeParam(3); %%% (/hr) proliferation rate APC wrt PDGFA
alpha12 = 0; %%% (/hr) proliferation rate APC wrt choroid oxygen
alpha13 = 0; %%% (/hr) proliferation rate APC wrt hylaoid oxygen
alpha20 = modeParam(4); %%% (/hr) basal proliferation rate IPA
alpha21 = modeParam(5); %%% (/hr) proliferation rate IPA wrt PDGFA
alpha22 = modeParam(6); %%% (/hr) proliferation rate IPA wrt choroid oxygen
alpha23 = modeParam(7); %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
beta0 = 0; %%% (/hr) basal differentiation rate
beta1 = modeParam(8); %%% (/hr) differentiation rate wrt LIF
beta2 = 0; %%% (/hr) differentiation rate wrt choroid oxygen
beta3 = 0; %%% (/hr) differentiation rate wrt hyaloid oxygen
beta4 = modeParam(9); %%% (/hr) mass action rate
eta1 = modeParam(10); %%% (/hr) apoptosis rate APC
eta2 = modeParam(11); %%% (/hr) apoptosis rate IPA
P_hy = modeParam(12); %%% partial pressure of oxygen due to hyaloid artery
r_hy = modeParam(13); %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(mu,alpha10,alpha11,...
    alpha12,alpha13,alpha20,alpha21,alpha22,alpha23,beta0,beta1,beta2,...
    beta3,beta4,eta1,eta2,P_hy,r_hy,m);
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
plot_the_plots_APCIPA
plot_movingbdy