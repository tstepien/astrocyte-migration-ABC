clear variables global;
clc;

% savefiles = 'yes';

multiplier = 5;
power = 5;
N = (multiplier)*10^(power);
num_param = 10;

iternum = 3;

%%%%%%%%%%%%%%%%%%%%%% create Latin Hypercube points %%%%%%%%%%%%%%%%%%%%%%
% %%% MATLAB function (painstakenly slow for large N)
% % LHpts = lhsdesign(N,13);
% %%% LATIN_RANDOM function
% % https://people.sc.fsu.edu/~jburkardt/m_src/latin_random/latin_random.html
% LHpts = latin_random(num_param, N)';
% save(strcat(pwd,'/LHpts/LHpts',num2str(num_param),'_',num2str(multiplier),...
%     'e',num2str(power),'.mat'),'LHpts')
% stop


%%%%%%%%%%%%%%%%%%%%%%% load Latin Hypercube points %%%%%%%%%%%%%%%%%%%%%%%
load(strcat(pwd,'/LHpts/LHpts',num2str(num_param),'_',num2str(multiplier),...
    'e',num2str(power),'.mat'))

% if strcmp(savefiles,'yes')==1
%     doublecheck = input('Are you sure you would like to save the output files? (it may overwrite): ');
%     if strcmp(doublecheck,'y')==1
        diary(strcat(pwd,'/parameter_analysis/diary_abc',num2str(num_param),'_',...
            num2str(multiplier),'e',num2str(power),'_',num2str(iternum),'.txt'));
%     else
%         return;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%% parameters to investigate %%%%%%%%%%%%%%%%%%%%%%%%
% set equal to 0: alpha13, alpha20, alpha23, beta3, P_hy, r_hy
% note: set r_hy=1 to avoid singularity at r=0, but essentially r_hy=0
% because P_hy=0

bound = [0.0001 0.1; %mu - adhesion constant
    0 2; %alpha10 - (/hr) base proliferation rate APC
    0 2; %alpha11 - (/hr) proliferation rate APC wrt PDGFA
    0 2; %alpha12 - (/hr) proliferation rate APC wrt choroid oxygen
    0 2; %alpha21 - (/hr) proliferation rate IPA wrt PDGFA
    0 2; %alpha22 - (/hr) proliferation rate IPA wrt choroid oxygen
    0 2; %beta1 - (/hr) differentiation rate wrt LIF
    0 2; %beta2 -  (/hr) differentiation rate wrt choroid oxygen
    0 2; %beta4 -  (/hr) mass action rate
    0 2]; %eta2 - (/hr) apoptosis rate IPA
numpar = length(bound);

mu      = (bound(1,2) - bound(1,1))*LHpts(:,1) + bound(1,1);
alpha10 = (bound(2,2) - bound(2,1))*LHpts(:,2) + bound(2,1);
alpha11 = (bound(3,2) - bound(3,1))*LHpts(:,3) + bound(3,1);
alpha12 = (bound(4,2) - bound(4,1))*LHpts(:,4) + bound(4,1);
alpha21 = (bound(5,2) - bound(5,1))*LHpts(:,5) + bound(5,1);
alpha22 = (bound(6,2) - bound(6,1))*LHpts(:,6) + bound(6,1);
beta1   = (bound(7,2) - bound(7,1))*LHpts(:,7) + bound(7,1);
beta2   = (bound(8,2) - bound(8,1))*LHpts(:,8) + bound(8,1);
beta4   = (bound(9,2) - bound(9,1))*LHpts(:,9) + bound(9,1);
eta2    = (bound(10,2) - bound(10,1))*LHpts(:,10) + bound(10,1);

%%%%%%%%%%%%%%%%%%%%%%%% preallocate error vectors %%%%%%%%%%%%%%%%%%%%%%%%
err_tot = zeros(N,1);
err_time = zeros(N,1);
err_rad = zeros(N,1);
err_dens = zeros(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1+2e5:3e5
    disp(['iteration i: ',num2str(i)])
    %%% solve equation
    [t,r,c1,c2,~,~,mvgbdy,~,~] = eqnsolver(mu(i),alpha10(i),alpha11(i),...
        alpha12(i),0,0,alpha21(i),alpha22(i),0,...
        0,beta1(i),beta2(i),0,beta4(i),0,eta2(i),...
        0,1,m);
    
    %%% error calculation
    [err_tot(i),err_time(i),err_rad(i),err_dens(i)]  = errorfunction(t,r,mvgbdy,c1,c2);
end

% if strcmp(savefiles,'yes')==1
    save(strcat(pwd,'/parameter_analysis/abc',num2str(num_param),'_',...
        num2str(multiplier),'e',num2str(power),'_',num2str(iternum),'.mat'));
    diary off
% end
