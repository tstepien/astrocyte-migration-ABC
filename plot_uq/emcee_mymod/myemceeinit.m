function [paramval,minquant] = myemceeinit(numwalkers,threshold,filestring,...
    nameofsavedfile,growthfunction,time,number_nodes,exp_px,exp_scale)

%%% range for uniform distribution
mrange = [0 , 1.5 ; ... % F/k
    500 , 5000 ; ... % k/b
    0 , 1 ; ... % alpha
    1000 , 2000]; % rho0

paramval = zeros(4,numwalkers);
minquant = zeros(1,numwalkers);

count_numparametersets = 0;

%%% make a set of starting points for the entire ensemble of walkers
% from a uniform distribution
for ii=1:numwalkers
    disp(['Walker #',num2str(ii)])
    min_quantity = NaN;
    
    while min_quantity>threshold || isnan(min_quantity)==1
        count_numparametersets = count_numparametersets + 1;
        
        %%% uniform distribution with values in interval (a,b):
        %%% a + (b-a).*rand(N,1)
        minit = [mrange(1,1) + (mrange(1,2) - mrange(1,1))*rand(1,1) ; ... % F/k
            mrange(2,1) + (mrange(2,2) - mrange(2,1))*rand(1,1) ; ... % k/b
            mrange(3,1) + (mrange(3,2) - mrange(3,1))*rand(1,1) ; ... % alpha
            mrange(4,1) + (mrange(4,2) - mrange(4,1))*rand(1,1)]; % rho0
        
        min_quantity = likelihoodfunc(minit,filestring,growthfunction,...
            time,number_nodes,exp_px,exp_scale);
    end
    
    paramval(:,ii) = minit;
    minquant(ii) = min_quantity;
    
    save(strcat('emceeinit_',nameofsavedfile,'.mat'),'paramval','minquant',...
        'count_numparametersets')
end