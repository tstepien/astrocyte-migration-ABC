function [Lvec, Pvec] = oxygen_jtb(r,thickness_ret,P0,Pm,Dalpha,M0)
% [Lvec, Pvec] = oxygen_setup(r,thickness_ret,P0,Pm,Dalpha,M0)
% 
% 1D oxygen diffusion with Michaelis-Menten uptake
% Calculate PO2 on inner surface of retina - TWS, June 2021
% dy1/dx = y2
% dy2/dx = M0 / Dalpha * y1/(Pm + y1)

    k = M0 / Dalpha;
    Lvec = [];
    Pvec = [];
    % Lvec = [0];
    % Pvec = [0]; %put PO2 = 0 when retina has zero thickness    %L is thickness of retina
    for i = 1:length(thickness_ret)
        if thickness_ret(i)>0
        L = thickness_ret(i);
        xmesh = linspace(0,L,21);
        solinit = bvpinit(xmesh, @(x) guess(x,L,P0));
        sol = bvp4c(@(x,y) bvpfcn(x,y,k,Pm), @(ya,yb) bcfcn(ya,yb,P0), solinit);
        Lvec = [Lvec; L];
        Pvec = [Pvec; sol.y(1,length(sol.x))]; 
        else
            Lvec = [Lvec; L];
            Pvec = [Pvec; 0]; 
        end
    end

Lvec = Lvec';
Pvec = Pvec';
end

function dydx = bvpfcn(~,y,k,Pm)
dydx = [y(2) (k * y(1) / (Pm + y(1)))];
end

function res = bcfcn(ya,yb,P0)
res = [(ya(1)-P0) yb(2)];
end

function g = guess(x,L,P0)
g = [(P0 * (1- x / L)^2) ( -2 * P0 * (1 - x / L) / L) ];
end