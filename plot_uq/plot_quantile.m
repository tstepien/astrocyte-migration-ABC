function [errorlevels,percentaccepted,levelp1] = plot_quantile(N,...
    err_original,num_subplot)
% [num_hold,param_sort_hold,levelp1] = sortparameters(N,...
%   err_original,num_subplot)
%
% Sort and determine the percent of accepted parameters based on various
% error levels
%
% inputs:
%   N            = number of ABC parameter sets
%   err_original = corresponding error of the ABC parameter sets
%   num_subplot  = number of subplot (corresponds to letter label)
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

percentcutoff = 0.01;
percentplotupperbound = 0.2;

ind1 = find(percentaccepted<percentcutoff,1,'last');
ind2 = find(percentaccepted>=percentcutoff,1,'first');
levelp1 = interp1([percentaccepted(ind1),percentaccepted(ind2)],...
    [errorlevels(ind1),errorlevels(ind2)],percentcutoff);

ind3 = find(percentaccepted>=percentplotupperbound,1,'first');

fig = figure;
set(gca,'FontSize',18)
hold on
plot(percentaccepted(1:ind3),errorlevels(1:ind3),'-o','Color','k',...
    'LineWidth',1.5)
plot(percentaccepted(1:ind3),levelp1*ones(size(percentaccepted(1:ind3))),...
    '--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
ylimval = get(gca,'YLim');
line([percentcutoff,percentcutoff],ylimval,'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
hold off
box on

if num_subplot==1
    title('A','FontSize',24);
    yannoteloc = levelp1/ylimval(2) + 0.1;
elseif num_subplot==2
    title('B','FontSize',24);
    yannoteloc = levelp1/ylimval(2) + 0.1;
elseif num_subplot==3
    title('C','FontSize',24);
    yannoteloc = levelp1/ylimval(2) - 0.01;
elseif num_subplot==4
    title('D','FontSize',24);
    yannoteloc = levelp1/ylimval(2) - 0.01;
elseif num_subplot==5
    title('E','FontSize',24);
    yannoteloc = levelp1/ylimval(2) - 0.01;
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';

annotation(fig,'textbox',...
    [0.65 yannoteloc 0.2 0.05],...
    'String',['$\mathcal{E}$=',num2str(levelp1)],'EdgeColor','none',...
    'Interpreter','latex','FontSize',24);
xlabel('Percent of parameter sets accepted','Interpreter','latex','FontSize',24)
ylabel('Error $\mathcal{E}$','Interpreter','latex','FontSize',24)
set(gca,'XLim',[0,percentplotupperbound])
