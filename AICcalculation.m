clear variables
clc

pn = [18;13;11;10];
P = length(pn);

% number of data points: 7 APC radii, 7 IPA radii, 1 time point = 15
numdata = 15;

param_val = zeros(P,18);
error_val = zeros(P,4);
AIC_val = zeros(P,1);
AICc_val = zeros(P,1);

for j=1:P
    load(strcat('parameter_analysis/abc',num2str(pn(j)),'_5e5.mat'));

    ind = find(err_tot==min(err_tot));

    if pn(j)==10
        P_alpha22 = 0;
    else
        P_alpha22 = alpha22(ind);
    end
    
    if pn(j)==10 || pn(j)==11
        P_beta2 = 0;
        P_eta1 = 0;
    else
        P_beta2 = beta2(ind);
        P_eta1 = eta1(ind);
    end

    if pn(j)==10 || pn(j)==11 || pn(j)==13
        P_alpha13 = 0;
        P_alpha23 = 0;
        P_beta3 = 0;
        P_Phy = 0;
        P_rhy = 1;
    else
        P_alpha13 = alpha13(ind);
        P_alpha23 = alpha23(ind);
        P_beta3 = beta3(ind);
        P_Phy = P_hy(ind);
        P_rhy = r_hy(ind);
    end

    P_mu = mu(ind);
    P_alpha10 = alpha10(ind);
    P_alpha11 = alpha11(ind);
    P_alpha12 = alpha12(ind);
    P_alpha20 = alpha20(ind);
    P_alpha21 = alpha21(ind);
    P_beta0 = beta0(ind);
    P_beta1 = beta1(ind);
    P_beta4 = beta4(ind);
    P_eta2 = eta2(ind);

    param_val(j,:) = [P_mu,P_alpha10,P_alpha11,P_alpha12,P_alpha13,...
        P_alpha20,P_alpha21,P_alpha22,P_alpha23,P_beta0,P_beta1,P_beta2,...
        P_beta3,P_beta4,P_eta1,P_eta2,P_Phy,P_rhy];
    error_val(j,:) = [err_dens(ind),err_rad(ind),err_time(ind),err_tot(ind)];

    AIC_val(j) = 2*pn(j) - 2*log(error_val(j,4));
    AICc_val(j) = AIC_val(j) + 2*pn(j)*(pn(j)+1)/(numdata-pn(j)-1);
    
    clear ind
end

min_error = error_val(:,4);
table(pn,min_error,AIC_val,AICc_val)

% param_val = [mu(ind),alpha10(ind),alpha11(ind),alpha12(ind),alpha13(ind),...
%         alpha20(ind),alpha21(ind),alpha22(ind),alpha23(ind),beta0(ind),beta1(ind),beta2(ind),...
%         beta3(ind),beta4(ind),eta1(ind),eta2(ind),P_hy(ind),r_hy(ind)];