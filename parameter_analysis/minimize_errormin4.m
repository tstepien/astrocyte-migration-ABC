num_param = 13;

load(strcat('abc',num2str(num_param),'_5e5.mat'));

totalnumber = 5;
numsmall = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_min,err_ind] = mink(err_tot,totalnumber);
ind = err_ind(numsmall);

mu = mu(ind); %%% adhesion constant
alpha10 = alpha10(ind); %%% (/hr) basal proliferation rate APC
alpha11 = alpha11(ind); %%% (/hr) proliferation rate APC wrt PDGFA
alpha12 = alpha12(ind); %%% (/hr) proliferation rate APC wrt choroid oxygen
if length(alpha13)>1
    alpha13 = alpha13(ind); %%% (/hr) proliferation rate APC wrt hylaoid oxygen
end
alpha20 = alpha20(ind); %%% (/hr) basal proliferation rate IPA
alpha21 = alpha21(ind); %%% (/hr) proliferation rate IPA wrt PDGFA
if length(alpha22)>1
    alpha22 = alpha22(ind); %%% (/hr) proliferation rate IPA wrt choroid oxygen
end
if length(alpha23)>1
    alpha23 = alpha23(ind); %%% (/hr) proliferation rate IPA wrt hyaloid oxygen
end
beta0 = beta0(ind); %%% (/hr) basal differentiation rate
beta1 = beta1(ind); %%% (/hr) differentiation rate wrt LIF
if length(beta2)>1
    beta2 = beta2(ind); %%% (/hr) differentiation rate wrt choroid oxygen
end
if length(beta3)>1
    beta3 = beta3(ind); %%% (/hr) differentiation rate wrt hyaloid oxygen
end
if length(beta4)>1
    beta4 = beta4(ind); %%% (/hr) mass action rate
end
if length(eta1)>1
    eta1 = eta1(ind); %%% (/hr) apoptosis rate APC
end
eta2 = eta2(ind); %%% (/hr) apoptosis rate IPA
if length(P_hy)>1
    P_hy = P_hy(ind); %%% partial pressure of oxygen due to hyaloid artery
end
if length(r_hy)>1
    r_hy = r_hy(ind); %%% radius at half-maximum of Hill function for hyaloid
end

if num_param==18
    param_init = [mu,alpha10,alpha11,alpha12,alpha13,alpha20,alpha21,alpha22,...
        alpha23,beta0,beta1,beta2,beta3,beta4,eta1,eta2,P_hy,r_hy];
elseif num_param==13
    param_init = [mu,alpha10,alpha11,alpha12,alpha20,alpha21,alpha22,...
        beta0,beta1,beta2,beta4,eta1,eta2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ..
options = optimset('Display','iter','MaxIter',50);
if num_param==18
    [param_new,err_output,exitflag,fmsoutput] = fminsearch(@errorfunc18,param_init);
elseif num_param==13
    [param_new,err_output,exitflag,fmsoutput] = fminsearch(@errorfunc13,param_init);
end

save(strcat('fminsearchresults',num2str(num_param),'_',num2str(numsmall)),...
    'param_new','param_init','err_output','exitflag','fmsoutput');