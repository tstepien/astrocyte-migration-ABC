multiplier = 5;
power = 5;
N = (multiplier)*10^(power);

mu = zeros(N,1);
alpha10 = zeros(N,1);
alpha11 = zeros(N,1);
alpha12 = zeros(N,1);
% alpha13 = zeros(N,1);
alpha20 = zeros(N,1);
alpha21 = zeros(N,1);
alpha22 = zeros(N,1);
% alpha23 = zeros(N,1);
beta0 = zeros(N,1);
beta1 = zeros(N,1);
beta2 = zeros(N,1);
% beta3 = zeros(N,1);
beta4 = zeros(N,1);
eta1 = zeros(N,1);
eta2 = zeros(N,1);
% P_hy = zeros(N,1);
% r_hy = zeros(N,1);
err_tot = zeros(N,1);
err_time = zeros(N,1);
err_rad = zeros(N,1);
err_dens = zeros(N,1);
err_flag = cell(N,1);

parfor i = 1:N
    filename = sprintf(strcat(pwd,'/ABC_results/july2024/abc13/output%d.mat'),i);
    if isfile(filename)==1
    [params,errors] = par_load13(filename);
        mu(i) = params(1);
        alpha10(i) = params(2);
        alpha11(i) = params(3);
        alpha12(i) = params(4);
        % alpha13(i) = params(5);
        alpha20(i) = params(5);
        alpha21(i) = params(6);
        alpha22(i) = params(7);
        % alpha23(i) = params(9);
        beta0(i) = params(8);
        beta1(i) = params(9);
        beta2(i) = params(10);
        % beta3(i) = params(13);
        beta4(i) = params(11);
        eta1(i) = params(12);
        eta2(i) = params(13);
        % P_hy(i) = params(17);
        % r_hy(i) = params(18);
        err_tot(i) = errors{1};
        err_time(i) = errors{2};
        err_rad(i) = errors{3};
        err_dens(i) = errors{4};
        err_flag{i} = errors{5};
    end
end