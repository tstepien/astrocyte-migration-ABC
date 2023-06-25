function Y = uq_eqns_and_error13(X)
% Y = uq_eqns_and_error13(X)
%
% This function returns the computed values of the moving boundary location
% astrocyte migration for a single population of cells and the final
% simulated time (days)
%
% inputs:
%   X = [mu, alpha10, alpha11, alpha12, alpha20 , alpha21, ...
%           alpha22, beta0, beta1, beta2, beta4, ...
%           eta1, eta2]
%
% note: removed alpha13, alpha23, beta3, Phy, rhy
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
alpha22 = X(:,7);
alpha23 = 0;
beta0 = X(:,8);
beta1 = X(:,9);
beta2 = X(:,10);
beta3 = 0;
beta4 = X(:,11);
eta1 = X(:,12);
eta2 = X(:,13);

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
    [t,r,c1,c2,~,~,mvgbdy,~,~] = eqnsolver(mu(i),alpha10(i),alpha11(i),...
        alpha12(i),alpha13,alpha20(i),alpha21(i),alpha22(i),alpha23,...
        beta0(i),beta1(i),beta2(i),beta3,beta4(i),eta1(i),eta2(i),...
        P_hy,r_hy,m);
    
    [Y(i),~,~,~] = errorfunction(t,r,mvgbdy,c1,c2);
end