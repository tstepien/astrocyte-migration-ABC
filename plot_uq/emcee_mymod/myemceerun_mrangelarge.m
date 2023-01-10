function [models,logP] = myemceerun_mrangelarge(mu_input,numruns,filestring,growthfunction,...
    time,number_nodes,exp_px,exp_scale)

%%% initialize indicators for prior
time_simulation_end_early = 0;
centerdensincrease = 0;

%%% large range
mrange = [0 , 2 ; ... % F/k
    0 , 5000 ; ... % k/b
    0 , 1 ; ... % alpha
    1 , 3000]; % rho0

%%% make a set of starting points for the entire ensemble of walkers
% from a truncated normal distribution (to guarantee positive parameters)
if isnumeric(mu_input)==1
    mu = mu_input;
    
    %%% some sort of "variance"
    mvar = zeros(4,1);
    for kk = 1:4
        if kk==1 || kk==3
            mvar(kk) = max(abs(mrange(kk,:)-mu(kk)),[],2)/5;
        elseif kk==2 || kk==4
            mvar(kk) = max(abs(mrange(kk,:)-mu(kk)),[],2)*5;
        end
    end
    
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

keyboard

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