function err_tot = errorfunc13(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = param(1);
alpha10 = param(2);
alpha11 = param(3);
alpha12 = param(4);
alpha13 = 0;
alpha20 = param(5);
alpha21 = param(6);
alpha22 = param(7);
alpha23 = 0;
beta0 = param(8);
beta1 = param(9);
beta2 = param(10);
beta3 = 0;
beta4 = param(11);
eta1 = param(12);
eta2 = param(13);
P_hy = 0;
r_hy = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

oxyfunc = 'oxygen_zeroorder';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,r,c1,c2,~,~,mvgbdy,~,~,~] = eqnsolver(mu,...
    alpha10,alpha11,alpha12,alpha13,alpha20,alpha21,alpha22,alpha23,...
    beta0,beta1,beta2,beta3,beta4,eta1,eta2,P_hy,r_hy,m,oxyfunc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_tot,~,~,~,~] = errorfunction(t,r,mvgbdy,c1,c2);