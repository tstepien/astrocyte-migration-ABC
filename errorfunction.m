function [err_tot,err_time,err_rad,err_dens,err_flag] = errorfunction(t,r,mvgbdy,c1,c2)
% [err_tot,err_time,err_rad,err_dens,err_flag] = errorfunction(t,mvgbdy,c1,c2)
%
% This is the error function for comparing the experimental data with
% simulations
%
% inputs:
%   t      = time vector
%   r      = spatial grid vector
%   mvgbdy = vector with location of moving boundary over time
%   c1     = density of APCs
%   c2     = density of IPAs
%
% outputs:
%   err_tot  = total error = (err_rad + err_dens + err_time)
%   err_time = error from simulation time end
%   err_rad  = error from astrocyte radius
%   err_dens = error from astrocyte density
%   err_flag = error flag

err_flag = [];

%%% penalties
if t(end)/24>8 || t(end)/24<6 || ~isreal(t(end)) ...
        || isequal(sum(c1(end,:)<1),length(r)) ...
        || isequal(sum(c2(end,:)<1),length(r)) ...
        || sum(c1(:)<0 & abs(c1(:))>10*eps)>0 ...
        || sum(c2(:)<0 & abs(c2(:))>10*eps)>0
    err_tot = 10^4; %NaN;
    err_time = 10^4; %NaN;
    err_rad = 10^4; %NaN;
    err_dens = 10^4; %NaN;
    if t(end)/24>8
        err_flag = [err_flag , 1];
    end
    if t(end)/24<6
        err_flag = [err_flag , 2];
    end
    if ~isreal(t(end))
        err_flag = [err_flag , 3];
    end
    if isequal(sum(c1(end,:)<1),length(r))
        err_flag = [err_flag , 4];
    end
    if isequal(sum(c2(end,:)<1),length(r))
        err_flag = [err_flag , 5];
    end
    if sum(c1(:)<0 & abs(c1(:))>10*eps)>0
        err_flag = [err_flag , 6];
    end
    if sum(c2(:)<0 & abs(c2(:))>10*eps)>0
        err_flag = [err_flag , 7];
    end
    return
end

%%% fixed parameters
parameters_fixed

%%% APC and IPA radius (mm) for E15-E16, E18-E22/P0
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67];
rad_IPA = [0; 0.04; 0.08; 0.33; 1; 1.5; 2];

%%% times: day 0, 1, 2, 3, 4, 5, 6, 7
%%% but note that we don't have data for day 2
dayswithdata = [1:2 4:8];
numdays = length(dayswithdata);

%%% pre-allocate arrays
ind = zeros(numdays,1);
comp_APClarger = zeros(numdays,length(r));
comp_IPAlarger = zeros(numdays,length(r));
nodes_APC = zeros(numdays,length(r));
nodes_IPA = zeros(numdays,length(r));
nodes_retina = zeros(numdays,length(r));
numnodes_APC = zeros(numdays,1);
numnodes_IPA = zeros(numdays,1);
numnodes_retina = zeros(numdays,1);
err_APC = zeros(numdays,1);
err_IPA = zeros(numdays,1);

%%% set |entries|<eps to 0
c1(abs(c1)<eps) = 0;
c2(abs(c2)<eps) = 0;

%%% calculate values on each day of data
for i=1:numdays
    %%% radius
    jj = dayswithdata(i);
    ind(i) = find(abs((t/24-(jj-1)))==min(abs(t/24-(jj-1))),1,'first');
    
    %%% density - data
    if i==1
        nodes_APC(i,:) = (r<=rad_APC(i)); %% annulus
    else
        nodes_APC(i,:) = (r<=rad_APC(i) & r>rad_IPA(i)); %% annulus
        nodes_IPA(i,:) = (r<=rad_IPA(i)); %% disc
    end
    nodes_retina(i,:) = (r<=rad_APC(i)); %% retina
    
    numnodes_APC(i) = sum(nodes_APC(i,:));
    numnodes_IPA(i) = sum(nodes_IPA(i,:));
    numnodes_retina(i) = sum(nodes_retina(i,:));
    
    %%% density - simulations
    comp_APClarger(i,:) = (c1(ind(i),:) > c2(ind(i),:));
    comp_IPAlarger(i,:) = (c1(ind(i),:) < c2(ind(i),:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time error
err_time = abs(7 - t(end)/24);

%%% radius error
%%% note that rad_APC(1)=mvgbdy(1) by initial condition
err_rad = sum( abs(rad_APC(2:end) - mvgbdy(ind(2:end))) );

%%% density error
%%% start at i=2 since APC and IPA density will match by initial condition
for i=2:numdays
    err_APC(i) = 1 - jaccard(comp_APClarger(i,:),nodes_APC(i,:));
    err_IPA(i) = 1 - jaccard(comp_IPAlarger(i,:),nodes_IPA(i,:));
end
err_dens = sum(err_APC) + sum(err_IPA);

%%% total error
err_tot = err_time + 2*err_rad + err_dens;

%%% figure
% figure
% tiledlayout(2,2)
% 
% nexttile
% plot(r,comp_APClarger)
% title('where c1>c2 (APC>IPA)')
% 
% nexttile
% plot(r,comp_IPAlarger)
% title('where c1<c2 (APC<IPA)')
% 
% nexttile
% plot(r,nodes_APC)
% title('nodes where APCs are (annulus)')
% 
% nexttile
% plot(r,nodes_IPA)
% title('nodes where IPAs are (disc)')