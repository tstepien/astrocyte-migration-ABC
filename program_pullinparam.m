clear variables global;
clc;
close all;
addpath plot_simulations

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%% loading a parameter set file %%%%%%%%%%%%%%%%%%%%%%%
% load('parameter_analysis/minparam_18.mat');
load('parameter_analysis/pointestimate_inversion18.mat');
ind = 1;

mu = param_val(ind,1); %%% adhesion constant
alpha10 = param_val(ind,2); %%% (/hr) basal proliferation rate APC
alpha11 = param_val(ind,3); %%% (/hr) proliferation rate APC wrt PDGFA
alpha12 = param_val(ind,4); %%% (/hr) proliferation rate APC wrt choroid oxygen
alpha13 = param_val(ind,5); %%% (/hr) proliferation rate APC wrt hylaoid oxygen
alpha20 = param_val(ind,6); %%% (/hr) basal proliferation rate IPA
alpha21 = param_val(ind,7); %%% (/hr) proliferation rate IPA wrt PDGFA
alpha22 = param_val(ind,8); %%% (/hr) proliferation rate IPA wrt choroid oxygen
alpha23 = param_val(ind,9); %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
beta0 = param_val(ind,10); %%% (/hr) basal differentiation rate
beta1 = param_val(ind,11); %%% (/hr) differentiation rate wrt LIF
beta2 = param_val(ind,12); %%% (/hr) differentiation rate wrt choroid oxygen
beta3 = param_val(ind,13); %%% (/hr) differentiation rate wrt hyaloid oxygen
beta4 = param_val(ind,14); %%% (/hr) mass action rate
eta1 = param_val(ind,15); %%% (/hr) apoptosis rate APC
eta2 = param_val(ind,16); %%% (/hr) apoptosis rate IPA
P_hy = param_val(ind,17); %%% partial pressure of oxygen due to hyaloid artery
r_hy = param_val(ind,18); %%% radius at half-maximum of Hill function for hyaloid

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