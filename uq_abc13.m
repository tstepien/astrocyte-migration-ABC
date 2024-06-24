clear variables global;
clc;

multiplier = 5;
power = 5;
N = (multiplier)*10^(power);
num_param = 13;

rng(100,'twister');

oxyfunc = 'oxygen_zeroorder';

diary(strcat(pwd,'/ABC_results/diary_abc',num2str(num_param),'.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%% parameters to investigate %%%%%%%%%%%%%%%%%%%%%%%%
bound = [5 100; %mu - adhesion constant
    0 5; %alpha10 - (/hr) base proliferation rate APC
    0 5; %alpha11 - (/hr) proliferation rate APC wrt PDGFA
    0 5; %alpha12 - (/hr) proliferation rate APC wrt choroid oxygen
    % 0 5; %alpha13 - (/hr) proliferation rate APC wrt hyaloid oxygen
    0 5; %alpha20 - (/hr) base proliferation rate IPA
    0 5; %alpha21 - (/hr) proliferation rate IPA wrt PDGFA
    0 5; %alpha22 - (/hr) proliferation rate IPA wrt choroid oxygen
    % 0 5; %alpha23 - (/hr) proliferation rate IPA wrt hyaloid oxygen
    0 0.5; %beta0 - (/hr) base differentiation rate
    0 0.5; %beta1 - (/hr) differentiation rate wrt LIF
    0 0.5; %beta2 - (/hr) differentiation rate wrt choroid oxygen
    % 0 0.5; %beta3 - (/hr) differentiation rate wrt hyaloid oxygen
    0 0.5; %beta4 - (/hr) mass action rate
    0 5; %eta1 - (/hr) apoptosis rate APC
    0 10]; %eta2 - (/hr) apoptosis rate IPA
    % 0 20; %P_hy - partial pressure of oxygen due to hyaloid artery
    % 0.1 4]; %r_hy - radius at half-maximum of Hill function for hyaloid

numpar = length(bound);

mu      = (bound(1,2) - bound(1,1))*rand(N,1) + bound(1,1);
alpha10 = (bound(2,2) - bound(2,1))*rand(N,1) + bound(2,1);
alpha11 = (bound(3,2) - bound(3,1))*rand(N,1) + bound(3,1);
alpha12 = (bound(4,2) - bound(4,1))*rand(N,1) + bound(4,1);
alpha13 = 0;
alpha20 = (bound(5,2) - bound(5,1))*rand(N,1) + bound(5,1);
alpha21 = (bound(6,2) - bound(6,1))*rand(N,1) + bound(6,1);
alpha22 = (bound(7,2) - bound(7,1))*rand(N,1) + bound(7,1);
alpha23 = 0;
beta0   = (bound(8,2) - bound(8,1))*rand(N,1) + bound(8,1);
beta1   = (bound(9,2) - bound(9,1))*rand(N,1) + bound(9,1);
beta2   = (bound(10,2) - bound(10,1))*rand(N,1) + bound(10,1);
beta4   = (bound(11,2) - bound(11,1))*rand(N,1) + bound(11,1);
beta3   = 0;
eta1    = (bound(12,2) - bound(12,1))*rand(N,1) + bound(12,1);
eta2    = (bound(13,2) - bound(13,1))*rand(N,1) + bound(13,1);
P_hy    = 0;
r_hy    = 1;

%%%%%%%%%%%%%%%%%%%%%%%% preallocate error vectors %%%%%%%%%%%%%%%%%%%%%%%%
err_tot = zeros(N,1);
err_time = zeros(N,1);
err_rad = zeros(N,1);
err_dens = zeros(N,1);
err_flag = cell(N,1);

save(strcat(pwd,'/ABC_results/abc',num2str(num_param),'_setup.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor i=1:N
    disp(['iteration i: ',num2str(i)])
    %%% solve equation
    [t,r,c1,c2,~,~,mvgbdy,~,~] = eqnsolver(mu(i),alpha10(i),alpha11(i),...
        alpha12(i),alpha13,alpha20(i),alpha21(i),alpha22(i),alpha23,...
        beta0(i),beta1(i),beta2(i),beta3,beta4(i),eta1(i),eta2(i),...
        P_hy,r_hy,m,oxyfunc);
    
    %%% error calculation
    [err_tot(i),err_time(i),err_rad(i),err_dens(i),err_flag{i}] ...
        = errorfunction(t,r,mvgbdy,c1,c2);

    %%% save data for each iteration
    par_save13(sprintf(strcat(pwd,'/ABC_results/abc13/output%d.mat'), i), ...
        mu(i),alpha10(i),alpha11(i),alpha12(i),alpha20(i),alpha21(i),...
        alpha22(i),beta0(i),beta1(i),beta2(i),beta4(i),eta1(i),eta2(i),...
        err_tot(i),err_time(i),err_rad(i),err_dens(i),err_flag{i});
end

save(strcat(pwd,'/ABC_results/abc',num2str(num_param),'_allresults.mat'));
diary off