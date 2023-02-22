clear variables
clc

num_param = 13;
num_steps = 5e3;
num_chains = 500;

filename = strcat('../parameter_analysis/inversion',num2str(num_param),...
            '_',num2str(num_chains),'chains_',num2str(num_steps),'steps.mat');
load(filename);

PostSample = cell(num_param,1);
ksd = cell(num_param,1);
bounds = cell(num_param,1);
dom = cell(num_param,1);
modeParam = zeros(num_param,1);

for i=1:num_param
    PostSample{i} = squeeze(myBayesianAnalysis.Results.PostProc.PostSample(end,i,:));
    ksd{i} = ksdensity(PostSample{i});
    bounds{i} = PriorOpts.Marginals(i).Bounds;

    ind = (max(ksd{i}) == ksd{i});
    dom{i} = linspace(bounds{i}(1),bounds{i}(2),length(ksd{i}));
    modeParam(i) = dom{i}(ind);

    clear ind
end