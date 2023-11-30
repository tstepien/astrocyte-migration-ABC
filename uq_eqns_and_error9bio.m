function Y = uq_eqns_and_error9bio(X)
% Y = uq_eqns_and_error9bio(X)
%
% This function returns the computed values of the moving boundary location
% astrocyte migration for a single population of cells and the final
% simulated time (days)
%
% inputs:
%   X = [mu, alpha10, alpha11, alpha12, alpha20, alpha21, ...
%           beta1, beta4, eta2]
%
% output:
%   Y = total error

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters to examine %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APCs and IPAs
mu = X(:,1);
alpha10 = X(:,2);
alpha11 = X(:,3);
alpha12 = X(:,4);
alpha13 = 0;
alpha20 = X(:,5);
alpha21 = X(:,6);
alpha22 = 0;
alpha23 = 0;
beta0 = 0;
beta1 = X(:,7);
beta2 = 0;
beta3 = 0;
beta4 = X(:,8);
eta1 = 0;
eta2 = X(:,9);

%%% hyaloid artery
P_hy = 0;
r_hy = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mesh parameters
m.dr = 0.01;
m.rmax = 5;
m.tmax = 7*24;

%%%%%%%%%%%%%%%%%%% solve equation and calculate error %%%%%%%%%%%%%%%%%%%%
%%% initialize
N = size(X,1);
Y = zeros(N,1);

parfor i=1:N
    [t,r,c1,c2,~,~,mvgbdy,~,~,~] = eqnsolver(mu(i),alpha10(i),alpha11(i),...
        alpha12(i),alpha13,alpha20(i),alpha21(i),alpha22,alpha23,...
        beta0,beta1(i),beta2,beta3,beta4(i),eta1,eta2(i),...
        P_hy,r_hy,m,'oxygen_zeroorder');
    
    [Y(i),~,~,~,~] = errorfunction(t,r,mvgbdy,c1,c2);
end