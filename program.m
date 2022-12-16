clear variables global;
clc;
close all;
addpath plot_simulations

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 0.001; %%% adhesion constant
alpha10 = 0; %%% (/hr) basal proliferation rate APC
alpha11 = 0.4491; %%% (/hr) proliferation rate APC wrt oxygen
alpha12 = 0.9004; %%% (/hr) proliferation rate APC wrt PDGFA
alpha20 = 0; %%% (/hr) basal proliferation rate IPA
alpha21 = 0.3199; %%% (/hr) proliferation rate IPA wrt oxygen
alpha22 = 0.5114; %%% (/hr) proliferation rate IPA wrt PDGFA
beta0 = 0; %%% (/hr) basal differentiation rate
beta1 = 0.0094; %%% (/hr) mass action rate
beta2 = 0.0462; %%% (/hr) differentiation rate wrt oxygen
beta3 = 0.5908; %%% (/hr) differentiation rate wrt LIF
eta1 = 0.0091; %%% (/hr) apoptosis rate APC
eta2 = 0.4670; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyaloid artery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_hy = 4.6343; %%% partial pressure of oxygen due to hyaloid artery
r_hy = 0.6963; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(mu,alpha10,alpha11,...
    alpha12,alpha20,alpha21,alpha22,beta0,beta1,beta2,beta3,eta1,eta2,...
    P_hy,r_hy,m);
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