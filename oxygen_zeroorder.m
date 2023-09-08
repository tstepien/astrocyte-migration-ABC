function PO2 = oxygen_zeroorder(r,thickness_ret,P0,Pm,Dalpha,M0)
% PO2 = oxygen_zeroorder(r,thickness_ret,P0,Pm,Dalpha,M0)
% choroid PO2
%
% Equation from Dollery et al. (1969) PMID: 5395349
% Uniform oxygen demand throughout the retina
% a.k.a. Zero-order kinetics from Popel (1989) PMID: 2673661
%
% inputs:
%   r             = spatial mesh
%   thickness_ret = total retinal thickness (mm)
%   {others}      = parameters
%                   (note that Pm is extraneous)
%
% outputs:
%   PO2 = partial pressure of oxygen (O2)

%%% partial pressure of O2
ind = (0 < thickness_ret) & (thickness_ret <= sqrt(2*P0*Dalpha/M0));
if length(thickness_ret)==1
    PO2 = ( P0 - M0/(2*Dalpha)*thickness_ret.^2 ) .*ind .* ones(size(r));
else
    PO2 = ( P0 - M0/(2*Dalpha)*thickness_ret.^2 ) .*ind;
end