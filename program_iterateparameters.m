clear variables global;
clc;
close all;
addpath plot_simulations

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 0.01; %%% adhesion constant
alpha10 = [0,0.01]; %%% (/hr) basal proliferation rate APC
alpha11 = [0,0.01]; %%% (/hr) proliferation rate APC wrt oxygen
alpha12 = [0,0.01]; %%% (/hr) proliferation rate APC wrt PDGFA
alpha20 = [0,0.01]; %%% (/hr) basal proliferation rate IPA
alpha21 = [0,0.01]; %%% (/hr) proliferation rate IPA wrt oxygen
alpha22 = [0,0.01]; %%% (/hr) proliferation rate IPA wrt PDGFA
beta0 = [0,0.001]; %%% (/hr) basal differentiation rate
beta1 = [0,0.001]; %%% (/hr) mass action rate
beta2 = [0,0.001]; %%% (/hr) differentiation rate wrt oxygen
beta3 = [0,0.001]; %%% (/hr) differentiation rate wrt LIF
eta1 = [0,0.001]; %%% (/hr) apoptosis rate APC
eta2 = [0,0.001]; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyaloid artery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_hy = 4.6343; %%% partial pressure of oxygen due to hyaloid artery
r_hy = 0.6963; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

numparamval = 10;

for i=1:numparamval
    a10 = alpha10(1) + (alpha10(2)-alpha10(1))*rand(1);
    a11 = alpha11(1) + (alpha11(2)-alpha11(1))*rand(1);
    a12 = alpha12(1) + (alpha12(2)-alpha12(1))*rand(1);
    a20 = alpha20(1) + (alpha20(2)-alpha20(1))*rand(1);
    a21 = alpha21(1) + (alpha21(2)-alpha21(1))*rand(1);
    a22 = alpha22(1) + (alpha22(2)-alpha22(1))*rand(1);
    b0 = beta0(1) + (beta0(2)-beta0(1))*rand(1);
    b1 = beta1(1) + (beta1(2)-beta1(1))*rand(1);
    b2 = beta2(1) + (beta2(2)-beta2(1))*rand(1);
    b3 = beta3(1) + (beta3(2)-beta3(1))*rand(1);
    e1 = eta1(1) + (eta1(2)-eta1(1))*rand(1);
    e2 = eta2(1) + (eta2(2)-eta2(1))*rand(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(mu,a10,a11,a12,a20,...
    a21,a22,b0,b1,b2,b3,e1,e2,P_hy,r_hy,m);
toc


%% error
%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_tot,err_time,err_rad,err_dens] = errorfunction(t,r,mvgbdy,c1,c2)

%% plots
if t(end)>1
plot_the_plots_APCIPA
end

end