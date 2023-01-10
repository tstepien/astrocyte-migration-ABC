function [models,logP] = myemceerun(mu_input,numruns)

global time centerdensincrease;

%%% initialize indicators for prior
time.simulation_end_early = 0;
centerdensincrease = 0;

%%% range used for previous parameter estimation
mrange = [0.1 , 1.5 ; ... % F/k
    500 , 5000 ; ... % k/b
    0.2 , 1 ; ... % alpha
    1000 , 3000]; % rho0

%%% make a set of starting points for the entire ensemble of walkers
% from a truncated normal distribution (to guarantee positive parameters)
if isnumeric(mu_input)==1
    mu = mu_input;
    
    %%% some sort of "variance"
    mvar = max(abs(mrange-mu),[],2)/5;
    
    %%% number of walkers (at least 2 times number of parameters)
    numwalkers = 2*length(mu);

    minit = [rmvnrnd(mu(1),mvar(1),numwalkers,-1, 0)' ; ... % F/k
        rmvnrnd(mu(2),mvar(2),numwalkers,-1, 0)' ; ... % k/b
        rmvnrnd(mu(3),mvar(3),numwalkers,-1, 0)' ; ... % alpha
        rmvnrnd(mu(4),mvar(4),numwalkers,-1, 0)']; % rho0
elseif ischar(mu_input)==1
    load(mu_input,'models')
    minit = models(:,:,end);
    numwalkers = size(minit,2);
end

%%% prior - all 4 parameters need to be positive
%%%       - k/b needs to be larger than 100
%%%       - alpha needs to be less than 10
%%%       - cell boundary cannot go beyond the computational domain
%%%       - density in center needs to decrease
logprior = @(m) (m(1)>0) && (m(2)>0) && (m(2)>100) ...
    && (m(3)>0) && (m(3)<10) && (m(4)>0) ...
    && (time.simulation_end_early==0) && (centerdensincrease==0);

%%% likelihood
% this is based on the error function
loglike = @likelihoodfunc;

%%% prior and likelihood
logPfuns = {logprior loglike};

%%% apply the MCMC hammer
[models,logP] = gwmcmc(minit,logPfuns,numruns*numwalkers,'ThinChain',1);