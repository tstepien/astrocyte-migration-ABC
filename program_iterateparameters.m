clear variables global;
clc;
close all;
addpath plot_simulations

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
bnd.mu = [0,1]; %%% adhesion constant
bnd.alpha10 = [0,0.1]; %%% (/hr) basal proliferation rate APC
bnd.alpha11 = [0,0.1]; %%% (/hr) proliferation rate APC wrt oxygen
bnd.alpha12 = [0,0.1]; %%% (/hr) proliferation rate APC wrt PDGFA
bnd.alpha20 = [0,0.1]; %%% (/hr) basal proliferation rate IPA
bnd.alpha21 = [0,0.1]; %%% (/hr) proliferation rate IPA wrt oxygen
bnd.alpha22 = [0,0.1]; %%% (/hr) proliferation rate IPA wrt PDGFA
bnd.beta0 = [0,0.1]; %%% (/hr) basal differentiation rate
bnd.beta1 = [0,0.1]; %%% (/hr) mass action rate
bnd.beta2 = [0,0.1]; %%% (/hr) differentiation rate wrt oxygen
bnd.beta3 = [0,0.1]; %%% (/hr) differentiation rate wrt LIF
bnd.eta1 = [0,0.1]; %%% (/hr) apoptosis rate APC
bnd.eta2 = [0,0.1]; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyaloid artery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnd.P_hy = [0,5]; %%% partial pressure of oxygen due to hyaloid artery
bnd.r_hy = [0,1]; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%% initialize number of parameters %%%%%%%%%%%%%%%%%%%%%
numparamval = 1000;

param = zeros(numparamval,15);
err = zeros(numparamval,4);

for i=1:numparamval
    mu = bnd.mu(1) + (bnd.mu(2)-bnd.mu(1))*rand(1);
    alpha10 = bnd.alpha10(1) + (bnd.alpha10(2)-bnd.alpha10(1))*rand(1);
    alpha11 = bnd.alpha11(1) + (bnd.alpha11(2)-bnd.alpha11(1))*rand(1);
    alpha12 = bnd.alpha12(1) + (bnd.alpha12(2)-bnd.alpha12(1))*rand(1);
    alpha20 = bnd.alpha20(1) + (bnd.alpha20(2)-bnd.alpha20(1))*rand(1);
    alpha21 = bnd.alpha21(1) + (bnd.alpha21(2)-bnd.alpha21(1))*rand(1);
    alpha22 = bnd.alpha22(1) + (bnd.alpha22(2)-bnd.alpha22(1))*rand(1);
    beta0 = bnd.beta0(1) + (bnd.beta0(2)-bnd.beta0(1))*rand(1);
    beta1 = bnd.beta1(1) + (bnd.beta1(2)-bnd.beta1(1))*rand(1);
    beta2 = bnd.beta2(1) + (bnd.beta2(2)-bnd.beta2(1))*rand(1);
    beta3 = bnd.beta3(1) + (bnd.beta3(2)-bnd.beta3(1))*rand(1);
    eta1 = bnd.eta1(1) + (bnd.eta1(2)-bnd.eta1(1))*rand(1);
    eta2 = bnd.eta2(1) + (bnd.eta2(2)-bnd.eta2(1))*rand(1);
    P_hy = bnd.P_hy(1) + (bnd.P_hy(2)-bnd.P_hy(1))*rand(1);
    r_hy = bnd.r_hy(1) + (bnd.r_hy(2)-bnd.r_hy(1))*rand(1);

    param(i,:) = [mu,alpha10,alpha11,alpha12,alpha20,alpha21,alpha22,...
        beta0,beta1,beta2,beta3,eta1,eta2,P_hy,r_hy];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(mu,alpha10,...
        alpha11,alpha12,alpha20,alpha21,alpha22,beta0,beta1,beta2,beta3,...
        eta1,eta2,P_hy,r_hy,m);
    toc

    % print variable values
    parameter_table = table(mu,alpha10,alpha11,alpha12,alpha20,...
        alpha21,alpha22,beta0,beta1,beta2,beta3,eta1,eta2,P_hy,r_hy)


%% error
%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [err_tot,err_time,err_rad,err_dens] = errorfunction(t,r,mvgbdy,c1,c2);

    err(i,:) = [err_tot,err_time,err_rad,err_dens];

    error_table = table(err_tot,err_time,err_rad,err_dens)

%% plots
%     if t(end)>1 && err_tot<1000
%         plot_the_plots_APCIPA
%     end

end