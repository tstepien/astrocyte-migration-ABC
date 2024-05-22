function plot_histograms(pos_tiled,pos_tiled_ylabel,num_param,percentholdon,...
    num_hold,param_sort_hold,bestfitdist,bound,fit_dist_plot,titles_on,...
    param_names,param_names_words)
% plot_histograms(pos_tiled,pos_tiled_ylabel,num_param,percentholdon,...
%    num_hold,param_sort_hold,bestfitdist,bound,fit_dist_plot,titles_on,...
%    param_names,param_names_words)
%
% Plot histograms of the probability density distributions of each
% parameter (the posterior distribution from the ABC rejection method) with
% the best fitting pdf curve from minimizing the Earth mover's distance
%
% inputs:
%   pos_tiled         = indices where the plots will be tiled
%   pos_tiled_ylabel  = indices where the ylabel 'Probability' will be printed
%   num_param         = number of parameters in the model
%   percentholdon     = percent of accepted parameter sets to hold on to
%   num_hold          = number of accepted parameter sets from ABC that are
%                       being analyzed
%   param_sort_hold   = accepted parameter sets that are sorted according to
%                       increasing error values
%   bestfitdist       = the best fitting probability distribution for each
%                       parameter
%   bound             = lower and upper bounds for the parameter values
%   fit_dist_plot     = include the fitted pdf curve (1=yes, 0=no)
%   titles_on         = include overall title with percenthold on listed
%                       (1=yes, 0=no)
%   param_names       = cell of strings of parameter variables
%   param_names_words = cell of strings of parameter meanings in words

curvecolor = [0.8500 0.3250 0.0980];
bincolor = [0.3 0.3 0.3];

fig3 = figure;
tiledlayout(4,5,'TileSpacing','compact','Padding','compact')
for i=1:num_param
    nexttile(pos_tiled(i))
    
    if strcmp(fit_dist_plot,'no')==1
        histogram(param_sort_hold(:,i),'Normalization','probability',...
            'BinMethod','sturges','FaceColor','none','EdgeColor',bincolor);
    elseif strcmp(fit_dist_plot,'yes')==1
        if strcmp(bestfitdist{i},'Uniform')==0
            h = histfit(param_sort_hold(:,i),[],bestfitdist{i});
            h(1).FaceColor = 'none';
            h(1).EdgeColor = bincolor;
            h(2).Color = curvecolor;
            box on
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',round(yt/num_hold,2));
        else
            hold on
            histogram(param_sort_hold(:,i),'Normalization','probability',...
                'BinMethod','sqrt','FaceColor','none','EdgeColor',bincolor);
            hh = get(gca,'YLim');
            plot(linspace(bound(i,1),bound(i,2),100),hh(2)/2*ones(1,100),...
                'Color',curvecolor,'LineWidth',2.5)
            box on
            hold off
        end
    end
    
    xlabel(param_names{i},'Interpreter','latex')
    if strcmp(titles_on,'yes')==1
        title(param_names_words{i},'FontWeight','normal')
    end
    if ismember(i,pos_tiled_ylabel)
        ylabel('Probability','Interpreter','latex')
    end
    xlim(bound(i,:))
    
    set(gca,'FontSize',14)
end

if strcmp(titles_on,'yes')==1
    sgtitle(strcat(['Smallest ',num2str(percentholdon*100),'% Error (',num2str(num_hold),' parameter sets)']))
end

set(fig3,'Units','inches','Position',[0,0,15,8],'PaperPositionMode','auto')