clear variables global;
clc;

multiplier = 5;
power = 5;
N = (multiplier)*10^(power);
num_param = 9;

percentholdon = 0.025;
what_set = 'maxthreshold'; %'maxthreshold' or 'maxmode'
fit_dist_plot = 'yes'; % using percentholdon for distribution fits
titles_on = 'yes';

load(strcat('../ABC_results/abc',num2str(num_param),'uni_',...
    num2str(multiplier),'e',num2str(power),'.mat'))

err_original = [err_dens err_rad err_time err_tot];
err_names = {'Density Error','Radius Error','Time Error','Total Error'};

param_original = [mu, alpha10, alpha12, alpha20, alpha21, ...
    alpha22, beta1, beta2, eta2];

clear err_dens err_rad err_time err_tot mu alpha10 alpha11 alpha12 alpha20 ...
    alpha21 alpha22 beta1 beta2 beta4 eta2;

param_names = {'$\mu$','$\alpha_{10}$','$\alpha_{12}$',...
    '$\alpha_{20}$','$\alpha_{21}$','$\alpha_{22}$',...
    '$\beta_1$','$\beta_2$','$\eta_2$'};
num_param = length(param_names);
param_names_words = {'Adhesion constant','APC base prolif rate',...
    'APC prolif wrt choroid O_2',...
    'IPA base prolif rate','IPA prolif wrt PDGFA',...
    'IPA prolif wrt choroid O_2','Diff rate wrt LIF',...
    'Diff rate wrt choroid O_2','IPA apoptosis rate'};

%% sort and hold onto 'percentholdon' smallest parameter sets

[num_hold,param_sort_hold] = sortparameters(num_param,N,param_original,...
    err_original,err_names,percentholdon,what_set);

%% fit the data to probability distributions, calculate Earth mover's
% distance, and report best fitting distribution
[bestfitdist,bestfitdist_param] = fitdistEMD(num_param,num_hold,...
    param_sort_hold,bound);

% export distribution information
save(strcat('distributions',num2str(num_param),'uni.mat'),'bestfitdist',...
    'bestfitdist_param')

%% histograms of parameters

pos_tiled = [1:2,4,7:9,12:13,17];
pos_tiled_ylabel = [1,4,7,9];

plot_histograms(pos_tiled,pos_tiled_ylabel,num_param,percentholdon,num_hold,...
    param_sort_hold,bestfitdist,bound,fit_dist_plot,titles_on,param_names,...
    param_names_words)

%% corner plot

plot_cornerplot(num_param,param_names,param_sort_hold)

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