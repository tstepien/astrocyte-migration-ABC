clear variables global;
clc;
close all;
addpath plot_simulations

%%% time unit: hr
%%% space unit: mm

oxyfunc = 'oxygen_zeroorder';

%%%%%%%%%%%%%%%%%%%%%% loading a parameter set file %%%%%%%%%%%%%%%%%%%%%%%
num_param = 9;
ninetype = 'bio';
load(strcat('parameter_analysis/bestfitparam',num2str(num_param),ninetype,'.mat'));

if num_param==18
    mu = param_new(1); %%% adhesion constant
    alpha10 = param_new(2); %%% (/hr) basal proliferation rate APC
    alpha11 = param_new(3); %%% (/hr) proliferation rate APC wrt PDGFA
    alpha12 = param_new(4); %%% (/hr) proliferation rate APC wrt choroid oxygen
    alpha13 = param_new(5); %%% (/hr) proliferation rate APC wrt hylaoid oxygen
    alpha20 = param_new(6); %%% (/hr) basal proliferation rate IPA
    alpha21 = param_new(7); %%% (/hr) proliferation rate IPA wrt PDGFA
    alpha22 = param_new(8); %%% (/hr) proliferation rate IPA wrt choroid oxygen
    alpha23 = param_new(9); %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
    beta0 = param_new(10); %%% (/hr) basal differentiation rate
    beta1 = param_new(11); %%% (/hr) differentiation rate wrt LIF
    beta2 = param_new(12); %%% (/hr) differentiation rate wrt choroid oxygen
    beta3 = param_new(13); %%% (/hr) differentiation rate wrt hyaloid oxygen
    beta4 = param_new(14); %%% (/hr) mass action rate
    eta1 = param_new(15); %%% (/hr) apoptosis rate APC
    eta2 = param_new(16); %%% (/hr) apoptosis rate IPA
    P_hy = param_new(17); %%% partial pressure of oxygen due to hyaloid artery
    r_hy = param_new(18); %%% radius at half-maximum of Hill function for hyaloid
elseif num_param==13
    mu = param_new(1); %%% adhesion constant
    alpha10 = param_new(2); %%% (/hr) basal proliferation rate APC
    alpha11 = param_new(3); %%% (/hr) proliferation rate APC wrt PDGFA
    alpha12 = param_new(4); %%% (/hr) proliferation rate APC wrt choroid oxygen
    alpha13 = 0; %%% (/hr) proliferation rate APC wrt hylaoid oxygen
    alpha20 = param_new(5); %%% (/hr) basal proliferation rate IPA
    alpha21 = param_new(6); %%% (/hr) proliferation rate IPA wrt PDGFA
    alpha22 = param_new(7); %%% (/hr) proliferation rate IPA wrt choroid oxygen
    alpha23 = 0; %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
    beta0 = param_new(8); %%% (/hr) basal differentiation rate
    beta1 = param_new(9); %%% (/hr) differentiation rate wrt LIF
    beta2 = param_new(10); %%% (/hr) differentiation rate wrt choroid oxygen
    beta3 = 0; %%% (/hr) differentiation rate wrt hyaloid oxygen
    beta4 = param_new(11); %%% (/hr) mass action rate
    eta1 = param_new(12); %%% (/hr) apoptosis rate APC
    eta2 = param_new(13); %%% (/hr) apoptosis rate IPA
    P_hy = 0; %%% partial pressure of oxygen due to hyaloid artery
    r_hy = 1; %%% radius at half-maximum of Hill function for hyaloid
elseif num_param==11
    mu = param_new(1); %%% adhesion constant
    alpha10 = param_new(2); %%% (/hr) basal proliferation rate APC
    alpha11 = param_new(3); %%% (/hr) proliferation rate APC wrt PDGFA
    alpha12 = param_new(4); %%% (/hr) proliferation rate APC wrt choroid oxygen
    alpha13 = 0; %%% (/hr) proliferation rate APC wrt hylaoid oxygen
    alpha20 = param_new(5); %%% (/hr) basal proliferation rate IPA
    alpha21 = param_new(6); %%% (/hr) proliferation rate IPA wrt PDGFA
    alpha22 = param_new(7); %%% (/hr) proliferation rate IPA wrt choroid oxygen
    alpha23 = 0; %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
    beta0 = 0; %%% (/hr) basal differentiation rate
    beta1 = param_new(8); %%% (/hr) differentiation rate wrt LIF
    beta2 = param_new(9); %%% (/hr) differentiation rate wrt choroid oxygen
    beta3 = 0; %%% (/hr) differentiation rate wrt hyaloid oxygen
    beta4 = param_new(10); %%% (/hr) mass action rate
    eta1 = 0; %%% (/hr) apoptosis rate APC
    eta2 = param_new(11); %%% (/hr) apoptosis rate IPA
    P_hy = 0; %%% partial pressure of oxygen due to hyaloid artery
    r_hy = 1; %%% radius at half-maximum of Hill function for hyaloid
elseif num_param==9 && strcmp(ninetype,'bio')==1
    mu = param_new(1); %%% adhesion constant
    alpha10 = param_new(2); %%% (/hr) basal proliferation rate APC
    alpha11 = param_new(3); %%% (/hr) proliferation rate APC wrt PDGFA
    alpha12 = param_new(4); %%% (/hr) proliferation rate APC wrt choroid oxygen
    alpha13 = 0; %%% (/hr) proliferation rate APC wrt hylaoid oxygen
    alpha20 = param_new(5); %%% (/hr) basal proliferation rate IPA
    alpha21 = param_new(6); %%% (/hr) proliferation rate IPA wrt PDGFA
    alpha22 = 0; %%% (/hr) proliferation rate IPA wrt choroid oxygen
    alpha23 = 0; %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
    beta0 = 0; %%% (/hr) basal differentiation rate
    beta1 = param_new(7); %%% (/hr) differentiation rate wrt LIF
    beta2 = 0; %%% (/hr) differentiation rate wrt choroid oxygen
    beta3 = 0; %%% (/hr) differentiation rate wrt hyaloid oxygen
    beta4 = param_new(8); %%% (/hr) mass action rate
    eta1 = 0; %%% (/hr) apoptosis rate APC
    eta2 = param_new(9); %%% (/hr) apoptosis rate IPA
    P_hy = 0; %%% partial pressure of oxygen due to hyaloid artery
    r_hy = 1; %%% radius at half-maximum of Hill function for hyaloid
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad,choroidPO2] = eqnsolver(mu,...
    alpha10,alpha11,alpha12,alpha13,alpha20,alpha21,alpha22,alpha23,...
    beta0,beta1,beta2,beta3,beta4,eta1,eta2,P_hy,r_hy,m,oxyfunc);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_tot,err_time,err_rad,err_dens,err_flag] = errorfunction(t,r,mvgbdy,c1,c2);
disp(['total error: ',num2str(err_tot)])

%% plots
plot_the_plots_4panels