function plot_cornerplot(num_param,param_names,param_sort_hold)
% plot_cornerplot(num_param,param_names,param_sort_hold)
%
% Corner plot showing the 2D projected histograms of error values for pairs
% of parameter values
%
% inputs:
%   num_param       = number of parameters in the model
%   param_names     = cell of strings of parameter variables
%   param_sort_hold = accepted parameter sets that are sorted according to
%                     increasing error values

num_keepscatter = 3; %number of smallest errors to keep for scatter plot
param_min = param_sort_hold(1:num_keepscatter,:);

param_mean = zeros(1,num_param);
param_mode = zeros(1,num_param);
for i=1:num_param
    param_mean(i) = mean(param_sort_hold(:,i));
    param_mode(i) = mode(param_sort_hold(:,i));
end

fig4 = figure;

tiledlayout(num_param,num_param,'TileSpacing','compact','Padding','compact')
pos_tiled = 1;
for j=1:num_param-1
    pos_tiled = [pos_tiled num_param*j+1:num_param*j+1+j];
end

for i=1:length(pos_tiled)
    nexttile(pos_tiled(i))
    [cc,rr] = ind2sub([num_param,num_param],pos_tiled(i));
    if rr==cc % on the diagonal
        [f,xi] = ksdensity(param_sort_hold(:,cc));
        plot(xi,f);
    else % below the diagonal
        [pdfx,xi] = ksdensity(param_sort_hold(:,cc));
        [pdfy,yi] = ksdensity(param_sort_hold(:,rr));
        [xxi,yyi] = meshgrid(xi,yi);
        [pdfxx,pdfyy] = meshgrid(pdfx,pdfy);
        pdfxy = pdfxx.*pdfyy;
        contourf(xxi,yyi,pdfxy,'LineColor','none');
    end
    if cc==1
        ylabel(param_names{rr},'Interpreter','latex');
    end
    if rr==num_param
        xlabel(param_names{cc},'Interpreter','latex');
    end
    clear cc rr
end
set(fig4,'Units','inches','Position',[-1,-1,11,9],'PaperPositionMode','auto')