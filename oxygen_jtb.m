function Pvec = oxygen_jtb(thickness_ret,P0,Pm,Dalpha,M0)
% Pvec = oxygen_jtb(thickness_ret,P0,Pm,Dalpha,M0)
% 
% 1D oxygen diffusion with Michaelis-Menten uptake
% Calculate PO2 on inner surface of retina - TWS, June 2021
%   updated to preallocate arrays and add Jacobian - TLS, August 2023
% dy1/dx = y2
% dy2/dx = M0 / Dalpha * y1/(Pm + y1)

k = M0 / Dalpha;

N = length(thickness_ret);
Pvec = zeros(1,N); % put PO2 = 0 when retina has zero thickness

for i = 1:N
    if thickness_ret(i)>0
        xmesh = linspace(0,thickness_ret(i),21);
        solinit = bvpinit(xmesh, @(x) guess(x,thickness_ret(i),P0));
        opts = bvpset('FJacobian',@(x,y) jacbvp(x,y,k,Pm),'BCJacobian',@jacbc);
        sol = bvp4c(@(x,y) bvpfcn(x,y,k,Pm), @(ya,yb) bcfcn(ya,yb,P0), solinit, opts);
        Pvec(i) = sol.y(1,length(sol.x));
    end
end

end

%%%%%%%%%%%%%%%
function dydx = bvpfcn(~,y,k,Pm)
dydx = [y(2) 
    (k * y(1) / (Pm + y(1)))];
end

%%%%%%%%%%%%%%%
function res = bcfcn(ya,yb,P0)
res = [(ya(1)-P0) 
    yb(2)];
end

%%%%%%%%%%%%%%%
function g = guess(x,L,P0)
g = [(P0 * (1- x / L)^2) 
    ( -2 * P0 * (1 - x / L) / L) ];
end

%%%%%%%%%%%%%%%
function dfdy = jacbvp(~,y,k,Pm)
dfdy = [0 1 ;
       (k*Pm/(Pm+y(1))^2) 0];
end

%%%%%%%%%%%%%%%
function [dBCdya,dBCdyb] = jacbc(~,~)
dBCdya = [1 0
          0 0];
dBCdyb = [0 0
          0 1];
end