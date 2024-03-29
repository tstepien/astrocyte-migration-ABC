function [t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad,choroidPO2] = eqnsolver(mu,...
    alpha10,alpha11,alpha12,alpha13,alpha20,alpha21,alpha22,alpha23,beta0,...
    beta1,beta2,beta3,beta4,eta1,eta2,P_hy,r_hy,m,oxyfunc)
% [c1,c2] = [t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad,choroidPO2] = eqnsolver(mu,...
%     alpha10,alpha11,alpha12,alpha13,alpha20,alpha21,alpha22,alpha23,beta0,...
%     beta1,beta2,beta3,beta4,eta1,eta2,P_hy,r_hy,m,oxyfunc)
%
% this is the solver
%
% Inputs:
%   [...]   = parameter values
%   m       = structure of mesh values
%   oxyfunc = oxygen function - either 'oxygen_zeroorder' for zero-order 
% %           kinetics or 'oxygen_michmen' for Michaelis-Menten
%
% Outputs:
%   t       = time
%   r       = spatial mesh
%   c1      = density of APCs
%   c2      = density of IPAs
%   q1      = concentration of PDGFA
%   q2      = concentration of LIF
%   mvgbdy  = location of moving boundary
%   vel_cir = circumferential spreading (v/r)
%   vel_rad = radial spreading (partial v/partial r)
%   choroidPO2  = partial pressure of choroid oxygen

global whatstep tcurr;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%% rename inputted parameters from mesh structure %%%%%%%%%%%%%%
dr = m.dr;
rmax = m.rmax; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = m.tmax; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%% input all fixed parameters that are %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% known/derived from literature %%%%%%%%%%%%%%%%%%%%%%
parameters_fixed


%%%%%%%%%%%%%%%%%%%%% parameter scalings/calculations %%%%%%%%%%%%%%%%%%%%%
xi1 = xibar_PDGFA / phi;
xi2 = xibar_LIF / phi;

Tprimeatce = Tderivative(ce,kappa); % T'(ce)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0:dr:rmax;
R = length(r);
% rhalf = (r(1:R-1)+r(2:R))/2;

tol = 10^(-6); % tolerance for predictor-corrector scheme

if abs(s0/dr - round(s0/dr))>2*eps
    error('error specifying s(0): moving boundary must be located on grid node');
end

hy = hyaloid(r,P_hy,r_hy);


%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time
tcurr = 0;

%%% astrocytes
c1_init = zeros(1,R);
ce_prime = 1.021*ce;
for i=1:R
    if r(i)<=s0 %fitted parabola, c1_init=0 for r(i)>s0
        c1_init(i) = ce + (ce_prime - ce)*(1 - r(i)^2/s0^2);
    end
end
c2_init = zeros(1,R);

c1_old = c1_init;
c2_old = c2_init;
k_old = c1_old + c2_old;

%%% growth factors
q1_old = zeros(1,R);
q2_old = zeros(1,R);


%%%%%%%%%%%%%%%%%%%%%%% initialize final variables %%%%%%%%%%%%%%%%%%%%%%%%
mvgbdy = s0;
s = s0;
c1 = c1_init;
c2 = c2_init;
k = k_old;
q1 = q1_old;
q2 = q2_old;
t = tcurr;

%%% subscript i is for space, j is for time (write in the order (j,i))
j_init = s0/dr+1;
j = j_init;

%%% velocity
[vel_cir,vel_rad] = velocity(j,c1,c2,r,kappa,mu);

mvgbdy_vel = [];

%%% choroid oxygen PO2
[thickness_ret,~,~,~] = thick_rad(tcurr,r);
choroidPO2 = feval(oxyfunc,r,thickness_ret,P0,Pm,Dalpha,M0);

%%% concentration on the moving boundary (mb)
c1mb = c1(j);
c2mb = c2(j);

%%% parameters for predictor/corrector steps
aa = 1;
bb = 1/2;

while tcurr < tmax && j<R-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%% predictor step %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    whatstep = 'predictor';
    
    %%% should have s0 start further than 2 nodes in
    %%% lines commented so code runs faster, but left in case they're
    %%% desired
%     if s0==0
%         dt_p = dr;
%     elseif s0==dr
%         dt_p = 1/aa * mu/Tprimeatce * dr^2 / ( k_old(j) - k_old(j-1) );
%     else
        dt_p = 1/aa * mu/Tprimeatce * 2*dr^2 / ...
            ( 3*k_old(j) - 4*k_old(j-1) + k_old(j-2) );
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%% solve eqn's with dt_p %%%%%%%%%%%%%%%%%%%%%%%%
    dt_c = 0;
    numiter = 0;
    while abs(dt_p-dt_c)>=tol
        %%% cell layer thickness and radius
        [thickness_ret,thickness_RGC,radius_endo,~] = thick_rad(tcurr+dt_p,r);
    
        %%% oxygen
        PO2 = feval(oxyfunc,r,thickness_ret,P0,Pm,Dalpha,M0);
        
        %%% growth factors
        [q1_hat,~] = growthfactors_implicit(q1_old,q2_old,dt_p,tcurr,...
            r,dr,R,thickness_RGC,radius_endo,maxRGCthick,thickness_ret,D1,D2,...
            xi1,xi2,gamma1,gamma2);
        
        %%% cell sum
        k_hat = cellpops_sum(j,c1_old,c2_old,q1_hat,PO2,dt_p,r,Pm,kappa,...
            mu,alpha10,alpha11,alpha12,alpha13,alpha20,alpha21,alpha22,...
            alpha23,eta1,eta2,ce,cmax,hy);
        
        %%%%%%%%%%%%%%%%%%%%%%%%% corrector step %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% should have s0 start further than 2 nodes in
        %%% lines commented so code runs faster, but left in case they're
        %%% desired
%         if s0==0
%             dt_c = mu/Tprimeatce * dr / ( ...
%                 bb*( k_hat(j+1) - k_hat(j) )/dr ...
%                 + (1-bb)*( 1/aa*mu/Tprimeatce ) );
%         elseif s0==dr
%             dt_c = mu/Tprimeatce * dr^2 / ( ...
%                 bb*( 3*k_hat(j+1) - 4*k_hat(j) + k_hat(j-1) )/2 ...
%                 + (1-bb)*( k_old(j) - k_old(j-1) ) );
%         else
            dt_c = mu/Tprimeatce * 2*dr^2 / ( ...
                bb*( 3*k_hat(j+1) - 4*k_hat(j) + k_hat(j-1) ) ...
                + (1-bb)*( 3*k_old(j) - 4*k_old(j-1) + k_old(j-2) ) );
%         end
        
%         [tcurr/24 dt_p, dt_c] % max(q1_hat) max(q2_hat)]
%         keyboard
        
        if abs(dt_p-dt_c)<tol
            break;
        else
            dt_p = abs(real(dt_c));
            dt_c = 0;
            numiter = numiter + 1;
            if numiter >50
                disp('***stopping simulation since predictor-corrector stuck in a loop***')
                return;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% solve next time step %%%%%%%%%%%%%%%%%%%%%%%%%
    whatstep = 'corrector';
    
    %%% cell layer thickness and radius
    [thickness_ret,thickness_RGC,radius_endo,~] = thick_rad(tcurr+dt_c,r);
    
    %%% oxygen
    PO2 = feval(oxyfunc,r,thickness_ret,P0,Pm,Dalpha,M0);
    
    %%% growth factors
    [q1_new,q2_new] = growthfactors_implicit(q1_old,q2_old,dt_c,tcurr,...
        r,dr,R,thickness_RGC,radius_endo,maxRGCthick,thickness_ret,D1,D2,...
        xi1,xi2,gamma1,gamma2);
    
    %%% cell sum
    k_new = cellpops_sum(j,c1_old,c2_old,q1_new,PO2,dt_c,r,Pm,kappa,mu,...
        alpha10,alpha11,alpha12,alpha13,alpha20,alpha21,alpha22,alpha23,...
        eta1,eta2,ce,cmax,hy);
    
    %%% cells separate
    [c1_new,c2_new] = cellpops_separate(j,c1_old,c2_old,k_new,q1_new,...
        q2_new,PO2,dt_c,r,Pm,kappa,mu,alpha10,alpha11,alpha12,alpha13,...
        alpha20,alpha21,alpha22,alpha23,beta0,beta1,beta2,beta3,beta4,...
        eta1,eta2,ce,cmax,hy);

    %%%%%%%%%%%%%%%%%%%%%% reset for next time step %%%%%%%%%%%%%%%%%%%%%%%
    j = j+1;
    s = s + dr;
    tcurr = tcurr + dt_c;
%     if tcurr>67
%         keyboard
%     end
    c1_old = c1_new;
    c2_old = c2_new;
    k_old = k_new;
    q1_old = q1_new;
    q2_old = q2_new;
    
    %%% save variables
    mvgbdy = [mvgbdy ; s];
    c1 = [c1 ; c1_new];
    c2 = [c2 ; c2_new];
    k = [k ; k_new];
    q1 = [q1 ; q1_new];
    q2 = [q2 ; q2_new];
    t = [t ; tcurr];
    c1mb = [c1mb ; c1_new(j)];
    c2mb = [c2mb ; c2_new(j)];
    choroidPO2 = [choroidPO2 ; PO2];

    if sum(c1_new<0 & abs(c1_new)>10*eps )>0 ...
            || sum(c2_new<0 & abs(c2_new)>10*eps)>0
        disp('***stopping simulation since one of the cell densities went negative***')
        % note: but the negative value is larger than 10*machine epsilon
        return;
    end
    
    %%% velocity calculation
    [vel_cir_new,vel_rad_new] = velocity(j,c1_new,c2_new,r,kappa,mu);
    vel_cir = [vel_cir ; vel_cir_new];
    vel_rad = [vel_rad ; vel_rad_new];
    
    mvgbdy_vel = [mvgbdy_vel ; (mvgbdy(end)-mvgbdy(end-1))/(t(end)-t(end-1))];
    
    if mvgbdy_vel(end)<10^(-4) || isreal(mvgbdy_vel(end))==0
        disp('***stopping simulation since moving boundary velocity is <10^(-4)***')
        return;
    elseif mvgbdy_vel(end)<10^(-3)
        disp('**fyi moving boundary is moving slow**')
    end
    
    if tcurr/24>=12
        disp('***stopping simulation since ending time greater than 12 days***')
        return;
    end
end

disp(['location of moving boundary at last time step: ',num2str(mvgbdy(end))])
% disp('PO2 at last time step')
% PO2(end)
disp(['ending time in hours: ',num2str(t(end)/24)])
disp(['velocity of moving boundary at last time step: ',num2str(mvgbdy_vel(end))])