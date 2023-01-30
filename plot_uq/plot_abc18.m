clear variables global;
clc;

addpath emcee_mymod

multiplier = 1;
power = 6;
N = (multiplier)*10^(power);
num_param = 18;

percentholdon = 0.1;
what_set = 'maxthreshold'; %'maxthreshold' or 'maxmode'
fit_dist_plot = 'no'; % using percentholdon = 0.01 for distribution fits
titles_on = 'yes';

load(strcat('../parameter_analysis/abc',num2str(num_param),'_',...
    num2str(multiplier),'e',num2str(power),'.mat'))

err_original = [err_dens err_rad err_time err_tot];
err_names = {'Density Error','Radius Error','Time Error','Total Error'};

% errDensity.original = err_dens;
% errRadius.original = err_rad;
% errTime.original = err_time;
% errTotal.original = err_tot;

param_original = [mu, alpha10, alpha11, alpha12, alpha13, alpha20, alpha21, ...
    alpha22, alpha23, beta0, beta1, beta2, beta3, beta4, eta1 , eta2 , P_hy , r_hy];

clear err_dens err_rad err_time err_tot mu alpha11 alpha12 alpha13 ...
    alpha21 alpha22 alpha23 beta1 beta2 beta3 beta4 eta1 eta2 Te P_hy r_hy;

param_names = {'$\mu$','$\alpha_{10}$','$\alpha_{11}$','$\alpha_{12}$',...
    '$\alpha_{13}$','$\alpha_{20}$','$\alpha_{21}$','$\alpha_{22}$',...
    '$\alpha_{23}$','$\beta_0$','$\beta_1$','$\beta_2$','$\beta_3$',...
    '$\beta_4$','$\eta_1$','$\eta_2$','$P_\mathrm{hy}$','$r_\mathrm{hy}$'};
num_param = length(param_names);
param_names_words = {'Adhesion constant','APC base prolif rate',...
    'APC prolif wrt PDGFA','APC prolif wrt choroid O_2',...
    'APC prolif wrt hyaloid O_2','IPA base prolif rate',...
    'IPA prolif wrt PDGFA','APC prolif wrt choroid O_2',...
    'IPA prolif wrt hyaloid O_2','Base diff rate',...
    'Mass action rate','Diff rate wrt LIF',...
    'Diff rate wrt choroid O_2','Diff rate wrt hyaloid O_2',...
    'APC apoptosis rate','IPA apoptosis rate','Hyaloid max',...
    'Hyaloid half-max value'};

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

ind_maxmode = ind( err_original(:,1) <= modes_error(1) ...
    & err_original(:,2) <= modes_error(2) ...
    & err_original(:,3) <= modes_error(3) );
num_maxmode = length(ind_maxmode);

err_maxmode = err_original(ind_maxmode,:);

[~,ind_sort_maxmode] = sort(err_maxmode(:,4));
err_maxmode_sort = err_maxmode(ind_sort_maxmode,:);

fig2 = figure;
tiledlayout(2,2)
for i=1:4
    nexttile
    scatter(1:num_maxmode,err_maxmode_sort(:,i))
    xlim([0,num_maxmode])
    xlabel(err_names{i})
end
sgtitle(strcat(['Errors < modes for density/radius/time (',num2str(num_maxmode),' parameter sets)']))

%% histograms of parameters

if strcmp(what_set,'maxmode')==1
    num_parametersets = num_maxmode;
    ind_parametersets = ind_maxmode;
    ind_sort = ind_sort_maxmode;
elseif strcmp(what_set,'maxthreshold')==1
    num_parametersets = num_maxthreshold;
    ind_parametersets = ind_maxthreshold;
    ind_sort = ind_sort_maxthreshold;
end

num_hold = ceil(percentholdon * num_parametersets);

param_sort = zeros(num_parametersets,num_param);
param_sort_hold = zeros(num_hold,num_param);
for i = 1:num_param
    param_sort(:,i) = param_original(ind_parametersets,i);
    param_sort_hold(:,i) = param_sort(ind_sort(1:num_hold),i);
end

fig3 = figure;
tiledlayout(4,5,'TileSpacing','compact','Padding','compact')
tiledpos = [1:5,7:19];
for i=1:num_param
    nexttile(tiledpos(i))
    
    if strcmp(fit_dist_plot,'no')==1
        histogram(param_sort_hold(:,i),'Normalization','probability',...
            'BinMethod','sturges','FaceColor','none','LineWidth',1.5);
    elseif strcmp(fit_dist_plot,'yes')==1
        distributionfit = {'exponential','normal','uniform','exponential','uniform',...
            'uniform','exponential','uniform','exponential','normal',...
            'normal','normal','normal'};
        if strcmp(distributionfit{i},'uniform')==0
            h = histfit(param_sort_hold(:,i),[],distributionfit{i});
            h(1).FaceColor = 'none';
            h(2).Color = 'k';
            box on
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',round(yt/num_hold,2));
        else
            hold on
            histogram(param_sort_hold(:,i),'Normalization','probability',...
                'BinMethod','sqrt','FaceColor','none');
            plot(linspace(bound(i,1),bound(i,2),100),0.02*ones(1,100),'k',...
                'LineWidth',2.5)
            box on
            hold off
        end
    end
    
    xlabel(param_names{i},'Interpreter','latex')
    if strcmp(titles_on,'yes')==1
        title(param_names_words{i},'FontWeight','normal')
    end
    if i==1 || i==6 || i==10 || i==15
        ylabel('Percentage','Interpreter','latex')
    end
    xlim([0,bound(i,2)])
    
    set(gca,'FontSize',14)
end

if strcmp(titles_on,'yes')==1
    sgtitle(strcat(['Smallest ',num2str(percentholdon*100),'% Error (',num2str(num_hold),' parameter sets)']))
end

set(fig3,'Units','inches','Position',[2,2,15,8],'PaperPositionMode','auto')

%% determine type of distribution

disttype = {'Normal';'Lognormal';'Gamma';'Exponential';'Weibull';...
    'Logistic';'Uniform'};
%%% didn't use these distributions:
%%% 'beta';'birnbaumsaunders';'burr';'negative binomial';'extreme value';'kernel';
%%% 'generalized extreme value';'generalized pareto';'inversegaussian';
%%% 'nakagami';'loglogistic';'poisson';'rayleigh';'rician';'tlocationscale';
num_dist = length(disttype);

dist_param = cell(num_dist-1,num_param);
GoF_dist = zeros(num_dist,num_param);
pval = zeros(num_dist,num_param);

dist_normal = cell(1,num_param);
dist_lognormal = cell(1,num_param);
dist_gamma = cell(1,num_param);
dist_exponential = cell(1,num_param);
dist_weibull = cell(1,num_param);
dist_logistic = cell(1,num_param);
dist_uniform = cell(1,num_param);

hist_normal = cell(1,num_param);

for i=1:num_dist-1
    for j=1:num_param
        dist_param{i,j} = fitdist(param_sort_hold(:,j),disttype{i});
%         [GoF_dist(i,j),pval(i,j)] = chi2gof(param_sort_hold(:,j),'CDF',dist_param{i,j}); %param_dist{i,j}
    end
end


%% calculate difference between sample and standard distributions 
%%% using Weisserstein metric / Earth mover's distance

dist_create = cell(num_dist,num_param);
distN = 1e5;

for i=1:num_dist
    for j=1:num_param
        if strcmp(disttype{i},'Normal')==1 || strcmp(disttype{i},'Lognormal')==1 ...
                || strcmp(disttype{i},'Logistic')==1
            dist_create{i,j} = pdf(disttype{i},linspace(bound(j,1),bound(j,2),distN),...
                dist_param{i,j}.mu,dist_param{i,j}.sigma);

        elseif strcmp(disttype{i},'Gamma')==1
            dist_create{i,j} = pdf(disttype{i},linspace(bound(j,1),bound(j,2),distN),...
                dist_param{i,j}.a,dist_param{i,j}.b);

        elseif strcmp(disttype{i},'Exponential')==1
            dist_create{i,j} = pdf(disttype{i},linspace(bound(j,1),bound(j,2),distN),...
                dist_param{i,j}.mu);

        elseif strcmp(disttype{i},'Weibull')==1
            dist_create{i,j} = pdf(disttype{i},linspace(bound(j,1),bound(j,2),distN),...
                dist_param{i,j}.A,dist_param{i,j}.B);

        elseif strcmp(disttype{i},'Uniform')==1
            dist_create{i,j} = pdf(disttype{i},linspace(bound(j,1),bound(j,2),distN),...
                bound(j,1),bound(j,2));
        end
    end
end

wsd1 = zeros(num_dist,num_param);
wsd2 = zeros(num_dist,num_param);

for i=1:num_dist
    for j=1:num_param
        wsd1(i,j) = ws_distance(dist_create{i,j},param_sort_hold(:,j),1);
        wsd2(i,j) = ws_distance(dist_create{i,j},param_sort_hold(:,j),2);
    end
end


% h=histogram(param_sort_hold(:,1),100);
% h2=histogram(dist_create{1,1},100);
% h=histogram(param_sort_hold(:,1),100);
% counts1=h.BinCounts;
% binLoc1=h.BinEdges;
% h2=histogram(dist_create{1,1},100);
% counts2=h2.BinCounts;
% binLoc2=h2.BinEdges;
% [f,fval]=emd(binLoc1',binLoc2',counts1'/sum(counts1),counts2'/sum(counts2),@gdf)

%% corner plot

% num_keepscatter = 3; %number to keep for scatter plot
% param_min = param_sort_hold(1:num_keepscatter,:);
% 
% param_mean = zeros(1,num_param);
% param_mode = zeros(1,num_param);
% for i=1:num_param
%     param_mean(i) = mean(param_sort_hold(:,i));
%     param_mode(i) = mode(param_sort_hold(:,i));
% end
% 
% fig4 = figure;
% ecornerplot(param_sort_hold,param_min,param_mean,bound,'names',param_names,'ks',true);
% set(fig4,'Units','inches','Position',[2,2,10,8],'PaperPositionMode','auto')

%% correlation

corrmatrix = zeros(num_param,num_param);
for i=1:num_param
    for j=1:num_param
        if j<i
            corrmatrix(i,j) = corr(param_sort_hold(:,i),param_sort_hold(:,j));
        end
    end
end

max(abs(corrmatrix(:)))