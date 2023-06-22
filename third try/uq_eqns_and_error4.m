function Y = uq_eqns_and_error4(X)
% Y = uq_eqns_and_error4(X)
%
% This function returns the computed values of the moving boundary location
% astrocyte migration for a single population of cells and the final
% simulated time (days)
%
% inputs:
%   X = [alpha11, alpha21, beta1, beta4]
%
% note: removed alpha13, alpha20, alpha23, beta3, Phy, rhy
%   beta0, eta1
%   eta2
%   alpha10, alpha22
%   beta2
%   alpha12 = 0.125
%   mu = 0.033
%
% output:
%   Y = total error

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters to examine %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APCs and IPAs
mu = 0.033;
alpha10 = 0;
alpha11 = X(:,1);
alpha12 = 0.125;
alpha13 = 0;
alpha20 = 0;
alpha21 = X(:,2);
alpha22 = 0;
alpha23 = 0;
beta0 = 0;
beta1 = X(:,3);
beta2 = 0;
beta3 = 0;
beta4 = X(:,4);
eta1 = 0;
eta2 = 0;

%%% hyaloid artery
P_hy = 0;
r_hy = 1; % set r_hy=1 to avoid singularity at r=0, but essentially r_hy=0
          % because P_hy=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mesh parameters
m.dr = 0.01;
m.rmax = 5;
m.tmax = 7*24;

%%%%%%%%%%%%%%%%%%% solve equation and calculate error %%%%%%%%%%%%%%%%%%%%
%%% initialize
N = size(X,1);
Y = zeros(N,1);

for i=1:N
    [t,r,c1,c2,~,~,mvgbdy,~,~] = eqnsolver(mu,alpha10,alpha11(i),...
        alpha12,alpha13,alpha20,alpha21(i),alpha22,alpha23,...
        beta0,beta1(i),beta2,beta3,beta4(i),eta1,eta2,...
        P_hy,r_hy,m);
    
    [Y(i),~,~,~] = errorfunction(t,r,mvgbdy,c1,c2);
end