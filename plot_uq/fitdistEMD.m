function [bestfitdist,bestfitdist_param] = fitdistEMD(num_param,...
    num_hold,param_sort_hold,bound)
% [bestfitdist,bestfitdist_param] = fitdistEMD(num_param,num_hold,param_sort_hold,bound)
%
% Fit data to different probability distributions, then create synthetic
% data based on the fitted distributions, and then calculate difference 
% between the synthetic and data probability distributions using the 
% Weisserstein metric / Earth mover's distance. The smallest difference
% indicates the best fitting probability distribution.
%
% inputs:
%   num_param       = number of parameters in the model
%   num_hold        = number of accepted parameter sets from ABC that are
%                     being analyzed
%   param_sort_hold = accepted parameter sets that are sorted according to
%                     increasing error values
%   bound           = lower and upper bounds for the parameter values
%
% outputs:
%   bestfitdist       = the best fitting probability distribution for each
%                       parameter
%   bestfitdist_param = the parameters for the best fitting probability
%                       distributions (e.g. mu and sigma for a normal
%                       distribution)

%% fit the data to different probability distributions

dist_type = {'Normal';'Lognormal';'Gamma';'Exponential';'Weibull';...
    'Logistic';'Uniform'};
%%% didn't use these distributions:
%%% 'beta';'birnbaumsaunders';'burr';'negative binomial';'extreme value';'kernel';
%%% 'generalized extreme value';'generalized pareto';'inversegaussian';
%%% 'nakagami';'loglogistic';'poisson';'rayleigh';'rician';'tlocationscale';
num_dist = length(dist_type);

dist_param = cell(num_dist,num_param);

% all distributions but uniform
for i=1:num_dist-1
    for j=1:num_param
        dist_param{i,j} = fitdist(param_sort_hold(:,j),dist_type{i});
    end
end

% uniform
for j=1:num_param
    dist_param{num_dist,j}.Lower = bound(j,1);
    dist_param{num_dist,j}.Upper = bound(j,2);
end

%% create synthetic data based on the fitted distributions
rng(100,'twister')

dist_synth = cell(num_dist,num_param);

for i=1:num_dist
    for j=1:num_param
        if strcmp(dist_type{i},'Normal')==1 || strcmp(dist_type{i},'Lognormal')==1 ...
                || strcmp(dist_type{i},'Logistic')==1
            dist_synth{i,j} = random(dist_type{i},dist_param{i,j}.mu,...
                dist_param{i,j}.sigma,num_hold,1);

         elseif strcmp(dist_type{i},'Gamma')==1
            dist_synth{i,j} = random(dist_type{i},dist_param{i,j}.a,...
                dist_param{i,j}.b,num_hold,1);

        elseif strcmp(dist_type{i},'Exponential')==1
            dist_synth{i,j} = random(dist_type{i},dist_param{i,j}.mu,...
                num_hold,1);

        elseif strcmp(dist_type{i},'Weibull')==1
            dist_synth{i,j} = random(dist_type{i},dist_param{i,j}.A,...
                dist_param{i,j}.B,num_hold,1);

        elseif strcmp(dist_type{i},'Uniform')==1
            dist_synth{i,j} = dist_param{i,j}.Lower ...
                + (dist_param{i,j}.Upper - dist_param{i,j}.Lower).*rand(num_hold,1);
        end
    end
end


%% calculate difference between synthetic and data probability distributions 
%%% using Weisserstein metric / Earth mover's distance

wsd1 = zeros(num_dist,num_param);
wsd2 = zeros(num_dist,num_param);

for i=1:num_dist
    for j=1:num_param
        wsd1(i,j) = ws_distance(dist_synth{i,j},param_sort_hold(:,j),1);
        wsd2(i,j) = ws_distance(dist_synth{i,j},param_sort_hold(:,j),2);
    end
end

bestfitdist = cell(1,num_param);
bestfitdist_param = cell(1,num_param);

for j=1:num_param
    ind = min(wsd1(:,j))==wsd1(:,j);
    bestfitdist{j} = dist_type{ind};
    bestfitdist_param(j) = dist_param(ind,j);
end
disp(bestfitdist);