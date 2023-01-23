function Y = uq_eqns_and_error18(X)
% Y = uq_eqns_and_error18(X)
%
% This function returns the computed values of the moving boundary location
% astrocyte migration for a single population of cells and the final
% simulated time (days)
%
% inputs:
%   X = [mu, alpha10, alpha11, alpha12, alpha13, alpha20, alpha21, ...
%           alpha22, alpha23, beta0, beta1, beta2, beta3, beta4, ...
%           eta1, eta2, P_hy, r_hy]
%
% output:
%   Y = total error

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters to examine %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APCs and IPAs
mu = X(:,1);
alpha10 = X(:,2);
alpha11 = X(:,3);
alpha12 = X(:,4);
alpha13 = X(:,5);
alpha20 = X(:,6);
alpha21 = X(:,7);
alpha22 = X(:,8);
alpha23 = X(:,9);
beta0 = X(:,10);
beta1 = X(:,11);
beta2 = X(:,12);
beta3 = X(:,13);
beta4 = X(:,14);
eta1 = X(:,15);
eta2 = X(:,16);

%%% hyaloid artery
P_hy = X(:,17);
r_hy = X(:,18);

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
        alpha12(i),alpha13(i),alpha20(i),alpha21(i),alpha22(i),alpha23(i),...
        beta0(i),beta1(i),beta2(i),beta3(i),beta4(i),eta1(i),eta2(i),...
        P_hy(i),r_hy(i),m);
    
    [Y(i),~,~,~] = errorfunction(t,r,mvgbdy,c1,c2);
end