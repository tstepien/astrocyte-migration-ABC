clear variables global;
clc;
close all;
addpath plot_simulations

%%% time unit: hr
%%% space unit: mm

oxyfunc = 'oxygen_zeroorder';

%%%%%%%%%%%%%%%%%%%%%% loading a parameter set file %%%%%%%%%%%%%%%%%%%%%%%
load('ABC_results/abc18_5e5.mat');
smallestval = 1;
errorsmall = mink(err_tot,smallestval);
ind = find(err_tot == errorsmall(end));

mu = mu(ind); %%% adhesion constant
alpha10 = alpha10(ind); %%% (/hr) basal proliferation rate APC
if length(alpha11)>1
    alpha11 = alpha11(ind); %%% (/hr) proliferation rate APC wrt PDGFA
end
alpha12 = alpha12(ind); %%% (/hr) proliferation rate APC wrt choroid oxygen
if length(alpha13)>1
    alpha13 = alpha13(ind); %%% (/hr) proliferation rate APC wrt hylaoid oxygen
end
alpha20 = alpha20(ind); %%% (/hr) basal proliferation rate IPA
alpha21 = alpha21(ind); %%% (/hr) proliferation rate IPA wrt PDGFA
if length(alpha22)>1
    alpha22 = alpha22(ind); %%% (/hr) proliferation rate IPA wrt choroid oxygen
end
if length(alpha23)>1
    alpha23 = alpha23(ind); %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
end
if length(beta0)>1
    beta0 = beta0(ind); %%% (/hr) basal differentiation rate
end
beta1 = beta1(ind); %%% (/hr) differentiation rate wrt LIF
if length(beta2)>1
    beta2 = beta2(ind); %%% (/hr) differentiation rate wrt choroid oxygen
end
if length(beta3)>1
    beta3 = beta3(ind); %%% (/hr) differentiation rate wrt hyaloid oxygen
end
if length(beta4)>1
    beta4 = beta4(ind); %%% (/hr) mass action rate
end
if length(eta1)>1
    eta1 = eta1(ind); %%% (/hr) apoptosis rate APC
end
eta2 = eta2(ind); %%% (/hr) apoptosis rate IPA
if length(P_hy)>1
    P_hy = P_hy(ind); %%% partial pressure of oxygen due to hyaloid artery
end
if length(r_hy)>1
    r_hy = r_hy(ind); %%% radius at half-maximum of Hill function for hyaloid
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


%% error
%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_tot,err_time,err_rad,err_dens,err_flag] = errorfunction(t,r,mvgbdy,c1,c2);
disp(['total error: ',num2str(err_tot)])

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
plotletter = '';
plot_the_plots_4panels