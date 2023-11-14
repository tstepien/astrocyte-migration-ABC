clear variables global;
clc;

addpath emcee_mymod

multiplier = 5;
power = 5;
N = (multiplier)*10^(power);
num_param = 11;

percentholdon = 0.025;
what_set = 'maxthreshold'; %'maxthreshold' or 'maxmode'
fit_dist_plot = 'yes'; % using percentholdon for distribution fits
titles_on = 'yes';

load(strcat('../parameter_analysis/abc',num2str(num_param),'_',...
    num2str(multiplier),'e',num2str(power),'.mat'))

err_original = [err_dens err_rad err_time err_tot];
err_names = {'Density Error','Radius Error','Time Error','Total Error'};

param_original = [mu, alpha10, alpha11, alpha12, alpha20, alpha21, ...
    alpha22, beta1, beta2, beta4, eta2];

clear err_dens err_rad err_time err_tot mu alpha10 alpha11 alpha12 alpha20 ...
    alpha21 alpha22 beta1 beta2 beta4 eta2;

param_names = {'$\mu$','$\alpha_{10}$','$\alpha_{11}$','$\alpha_{12}$',...
    '$\alpha_{20}$','$\alpha_{21}$','$\alpha_{22}$',...
    '$\beta_1$','$\beta_2$','$\beta_4$','$\eta_2$'};
num_param = length(param_names);
param_names_words = {'Adhesion constant','APC base prolif rate',...
    'APC prolif wrt PDGFA','APC prolif wrt choroid O_2',...
    'IPA base prolif rate','IPA prolif wrt PDGFA',...
    'IPA prolif wrt choroid O_2','Diff rate wrt LIF',...
    'Diff rate wrt choroid O_2','Mass action rate','IPA apoptosis rate'};

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

%% fit the data to different probabiltiy distributions

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

% export distribution information
save(strcat('distributions',num2str(num_param),'.mat'),'bestfitdist','bestfitdist_param')

%% histograms of parameters

fig3 = figure;
tiledlayout(4,5,'TileSpacing','compact','Padding','compact')
pos_tiled = [1:4,7:9,12:13,15,17];
for i=1:num_param
    nexttile(pos_tiled(i))
    
    if strcmp(fit_dist_plot,'no')==1
        histogram(param_sort_hold(:,i),'Normalization','probability',...
            'BinMethod','sturges','FaceColor','none','LineWidth',1.5);
    elseif strcmp(fit_dist_plot,'yes')==1
        if strcmp(bestfitdist{i},'Uniform')==0
            h = histfit(param_sort_hold(:,i),[],bestfitdist{i});
            h(1).FaceColor = 'none';
            h(2).Color = 'k';
            box on
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',round(yt/num_hold,2));
        else
            hold on
            histogram(param_sort_hold(:,i),'Normalization','probability',...
                'BinMethod','sqrt','FaceColor','none');
            hh = get(gca,'YLim');
            plot(linspace(bound(i,1),bound(i,2),100),hh(2)/2*ones(1,100),'k',...
                'LineWidth',2.5)
            box on
            hold off
        end
    end
    
    xlabel(param_names{i},'Interpreter','latex')
    if strcmp(titles_on,'yes')==1
        title(param_names_words{i},'FontWeight','normal')
    end
    if i==1 || i==5 || i==8 || i==11
        ylabel('Percentage','Interpreter','latex')
    end
    xlim(bound(i,:))
    
    set(gca,'FontSize',14)
end

if strcmp(titles_on,'yes')==1
    sgtitle(strcat(['Smallest ',num2str(percentholdon*100),'% Error (',num2str(num_hold),' parameter sets)']))
end

set(fig3,'Units','inches','Position',[0,0,15,8],'PaperPositionMode','auto')

%% corner plot

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

%% correlation

corrmatrix = zeros(num_param,num_param);
for i=1:num_param
    for j=1:num_param
        if j<i
            corrmatrix(i,j) = corr(param_sort_hold(:,i),param_sort_hold(:,j));
        end
    end
end

maxcorrmatrix = max(abs(corrmatrix(:)))