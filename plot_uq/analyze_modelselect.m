clear variables global;
clc;

load('../ABC_results/modelselect_allresults.mat');

threshold = 16.25;

err_original = [err_dens err_rad err_time err_tot];
err_names = {'Density Error','Radius Error','Time Error','Total Error'};

param_original = [mu, alpha10, alpha11, alpha12, alpha13, alpha20, alpha21, ...
    alpha22, alpha23, beta0, beta1, beta2, beta3, beta4, eta1 , eta2 , P_hy , r_hy];

clear err_dens err_rad err_time err_tot mu alpha10 alpha11 alpha12 alpha13 alpha20 ...
    alpha21 alpha22 alpha23 beta0 beta1 beta2 beta3 beta4 eta1 eta2 P_hy r_hy;

param_names = {'$\mu$','$\alpha_{10}$','$\alpha_{11}$','$\alpha_{12}$',...
    '$\alpha_{13}$','$\alpha_{20}$','$\alpha_{21}$','$\alpha_{22}$',...
    '$\alpha_{23}$','$\beta_0$','$\beta_1$','$\beta_2$','$\beta_3$',...
    '$\beta_4$','$\eta_1$','$\eta_2$','$P_\mathrm{hy}$','$r_\mathrm{hy}$'};
num_param = length(param_names);
param_names_words = {'Adhesion constant','APC base prolif rate',...
    'APC prolif wrt PDGFA','APC prolif wrt choroid O_2',...
    'APC prolif wrt hyaloid O_2','IPA base prolif rate',...
    'IPA prolif wrt PDGFA','IPA prolif wrt choroid O_2',...
    'IPA prolif wrt hyaloid O_2','Base diff rate',...
    'Diff rate wrt LIF','Diff rate wrt choroid O_2',...
    'Diff rate wrt hyaloid O_2','Mass action rate',...
    'APC apoptosis rate','IPA apoptosis rate','Hyaloid max',...
    'Hyaloid half-max value'};

ind_hold = (err_original(:,4) < threshold);
err_new = err_original(ind_hold,:);
param_new = param_original(ind_hold,:);
modelnumber_new = modelnumber(ind_hold);

N = length(modelnumber_new);

percent1 = sum(modelnumber_new==1)/N;
percent2 = sum(modelnumber_new==2)/N;
percent3 = sum(modelnumber_new==3)/N;
percent4 = sum(modelnumber_new==4)/N;
percent5 = sum(modelnumber_new==5)/N;

[percent1 percent2 percent3 percent4 percent5]