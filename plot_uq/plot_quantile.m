function [errorlevels,percentaccepted,levelp1] = plot_quantile(N,err_original)
% [num_hold,param_sort_hold,levelp1] = sortparameters(N,err_original)
%
% Sort and determine the percent of accepted parameters based on various
% error levels
%
% inputs:
%   N              = number of ABC parameter sets
%   err_original   = corresponding error of the ABC parameter sets
%
% outputs:
%   errorlevels     = levels of error thresholds
%   percentaccepted = percent of parameter sets accepted at different error
%                     thresholds
%   levelp1         = error threshold that results in approximately 10% of
%                     the parameter sets chosen

%% remove errors that were set to 10^4
maxthreshold = 10^4;

ind = (1:N)';
ind_parametersets = ind(err_original(:,4) < maxthreshold);
num_parametersets = length(ind_parametersets);

err_maxthreshold = err_original(ind_parametersets,:);

[~,ind_sort] = sort(err_maxthreshold(:,4));
err_sort = err_maxthreshold(ind_sort,:);

%% determine percent accepted sets based on different error levels

errorlevels = (floor(min(err_sort(:,4))):ceil(max(err_sort(:,4))))';
EL = length(errorlevels);

percentaccepted = zeros(EL,1);

for i=1:EL
    percentaccepted(i) = sum( err_sort(:,4) < errorlevels(i)) ...
        / num_parametersets;
end

ind1 = find(percentaccepted<0.1,1,'last');
ind2 = find(percentaccepted>=0.1,1,'first');
levelp1 = interp1([percentaccepted(ind1),percentaccepted(ind2)],...
    [errorlevels(ind1),errorlevels(ind2)],0.1);

fig = figure;
set(gca,'FontSize',18)
hold on
plot(percentaccepted,errorlevels,'-o','Color','k','LineWidth',1.5)
plot(percentaccepted,levelp1*ones(size(percentaccepted)),'--','Color',...
    [0.5 0.5 0.5],'LineWidth',1.5)
ylimval = get(gca,'YLim');
line([0.1,0.1],ylimval,'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
hold off
box on
annotation(fig,'textbox',...
    [0.65 0.23 0.2 0.05],...
    'String',['$\mathcal{E}$=',num2str(levelp1)],'EdgeColor','none',...
    'Interpreter','latex','FontSize',24);
xlabel('Percent of parameter sets accepted','Interpreter','latex','FontSize',24)
ylabel('Error $\mathcal{E}$','Interpreter','latex','FontSize',24)
set(gca,'XLim',[0,0.5])