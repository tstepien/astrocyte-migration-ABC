function [models,logP] = myemceerun_par(mu_input,numruns,filestring,growthfunction,...
    time,number_nodes,exp_px,exp_scale)

%%% initialize indicators for prior
time_simulation_end_early = 0;
centerdensincrease = 0;

%%% range for uniform distribution
mrange = [0 , 1.5 ; ... % F/k
    500 , 5000 ; ... % k/b
    0 , 1 ; ... % alpha
    1000 , 2000]; % rho0

%%% make a set of starting points for the entire ensemble of walkers
% from a uniform distribution
if isnumeric(mu_input)==1
    mu = mu_input;
    
    %%% number of walkers (at least 2 times number of parameters)
    numwalkers = 2*length(mu);
    
    %%% uniform distribution with values in interval (a,b): a + (b-a).*rand(N,1)
    minit = [mrange(1,1) + (mrange(1,2) - mrange(1,1))*rand(1,numwalkers) ; ... % F/k
        mrange(2,1) + (mrange(2,2) - mrange(2,1))*rand(1,numwalkers) ; ... % k/b
        mrange(3,1) + (mrange(3,2) - mrange(3,1))*rand(1,numwalkers) ; ... % alpha
        mrange(4,1) + (mrange(4,2) - mrange(4,1))*rand(1,numwalkers)]; % rho0
elseif ischar(mu_input)==1
    load(mu_input,'models')
    minit = models(:,:,end);
    numwalkers = size(minit,2);
end

%%% prior - all 4 parameters need to be positive
%%%       - k/b needs to be larger than 100 --- remove?
%%%       - alpha needs to be less than 10 --- remove?
%%%       - cell boundary cannot go beyond the computational domain
%%%       - density in center needs to decrease
logprior = @(m) (m(1)>0) && (m(2)>0) ...%&& (m(2)>100) ...
    && (m(3)>0) ...%&& (m(3)<10) 
    && (m(4)>0) ...
    && (time_simulation_end_early==0) && (centerdensincrease==0);

%%% likelihood
% this is based on the error function
loglike = @(m) likelihoodfunc(m,filestring,growthfunction,time,...
    number_nodes,exp_px,exp_scale);

%%% prior and likelihood
logPfuns = {logprior loglike};

%%% apply the MCMC hammer
[models,logP] = gwmcmc(minit,logPfuns,numruns*numwalkers,'ThinChain',1,...
    'Parallel',true);