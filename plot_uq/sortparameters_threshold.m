function [num_hold,param_sort_hold] = sortparameters_threshold(param_original,...
    err_original,threshold)
% [num_hold,param_sort_hold] = sortparameters_threshold(param_original,...
%   err_original,threshold)
%
% Sort and hold onto parameter sets smaller than 'threshold'
%
% inputs:
%   param_original = ABC parameter sets
%   err_original   = corresponding error of the ABC parameter sets
%   threshold      = error threshold upper bound
%
% outputs:
%   num_hold        = number of accepted parameter sets from ABC to analyze
%   param_sort_hold = accepted parameter sets that are sorted according to
%                     increasing error values

[~,ind_sort] = sort(err_original(:,4));

err_sort = err_original(ind_sort,:);
param_sort = param_original(ind_sort,:);

ind_threshold = (err_sort(:,4) < threshold);

num_hold = sum(ind_threshold);

param_sort_hold = param_sort(ind_threshold,:);