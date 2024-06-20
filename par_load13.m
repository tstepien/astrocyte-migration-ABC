function [params,errors] = par_load13(fname)

load(fname,'mu','alpha10','alpha11','alpha12','alpha20','alpha21',...
    'alpha22','beta0','beta1','beta2','beta4','eta1','eta2','err_tot',...
    'err_time','err_rad','err_dens','err_flag')
params = [mu,alpha10,alpha11,alpha12,alpha20,alpha21,alpha22,beta0,beta1,...
    beta2,beta4,eta1,eta2];
errors = {err_tot,err_time,err_rad,err_dens,err_flag};

end