function [num_hold,param_sort_hold] = sortparameters_percent(num_param,N,...
    param_original,err_original,err_names,percentholdon)
% [num_hold,param_sort_hold] = sortparameters_percent(num_param,N,...
%   param_original,err_original,err_names,percentholdon)
%
% Sort and hold onto 'percentholdon' smallest parameter sets
%
% inputs:
%   num_param      = number of parameters in the model
%   N              = number of ABC parameter sets
%   param_original = ABC parameter sets
%   err_original   = corresponding error of the ABC parameter sets
%   err_names      = names of the types of errors calculated (for plot)
%   percentholdon  = percent of accepted parameter sets to hold on to
%
% outputs:
%   num_hold        = number of accepted parameter sets from ABC to analyze
%   param_sort_hold = accepted parameter sets that are sorted according to
%                     increasing error values

%% remove errors that were set to 10^4
maxthreshold = 10^4;

ind = (1:N)';
ind_maxthreshold = ind(err_original(:,4) < maxthreshold);
num_maxthreshold = length(ind_maxthreshold);

err_maxthreshold = err_original(ind_maxthreshold,:);

[~,ind_sort_maxthreshold] = sort(err_maxthreshold(:,4));
err_maxthreshold_sort = err_maxthreshold(ind_sort_maxthreshold,:);

fig1 = figure;
tiledlayout(2,2)
for i=1:4
    nexttile
    scatter(1:num_maxthreshold,err_maxthreshold_sort(:,i))
    xlim([0,num_maxthreshold])
    xlabel(err_names{i})
end
sgtitle(strcat(['Errors <10^4 (',num2str(num_maxthreshold),' parameter sets)']))

modes_error = zeros(1,4);
for i=1:4
    modes_error(i) = mode(err_maxthreshold(:,i));
end

%% look at errors that are smaller than the mode errors for density, radius, and time

% ind_maxmode = ind( err_original(:,1) <= modes_error(1) ...
%     & err_original(:,2) <= modes_error(2) ...
%     & err_original(:,3) <= modes_error(3) );
% num_maxmode = length(ind_maxmode);
% 
% err_maxmode = err_original(ind_maxmode,:);
% 
% [~,ind_sort_maxmode] = sort(err_maxmode(:,4));
% err_maxmode_sort = err_maxmode(ind_sort_maxmode,:);
% 
% fig2 = figure;
% tiledlayout(2,2)
% for i=1:4
%     nexttile
%     scatter(1:num_maxmode,err_maxmode_sort(:,i))
%     xlim([0,num_maxmode])
%     xlabel(err_names{i})
% end
% sgtitle(strcat(['Errors < modes for density/radius/time (',num2str(num_maxmode),' parameter sets)']))

%% sort and hold onto 'percentholdon' smallest parameter sets

num_parametersets = num_maxthreshold;
ind_parametersets = ind_maxthreshold;
ind_sort = ind_sort_maxthreshold;

num_hold = ceil(percentholdon * num_parametersets);

param_sort = zeros(num_parametersets,num_param);
param_sort_hold = zeros(num_hold,num_param);
for i = 1:num_param
    param_sort(:,i) = param_original(ind_parametersets,i);
    param_sort_hold(:,i) = param_sort(ind_sort(1:num_hold),i);
end