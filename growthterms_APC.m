function g1 = growthterms_APC(c1,c2,q1,q2,PO2,Pm,alpha10,alpha11,alpha12,...
    alpha13,beta0,beta1,beta2,beta3,beta4,eta1,cmax,hyaloid)
% g1 = growthterms_APC(c1,c2,q1,q2,PO2,Pm,alpha10,alpha11,alpha12,...
%     alpha13,beta0,beta1,beta2,beta3,beta4,eta1,cmax,hyaloid)
%
% growth terms just APCs (c1): g1

choroid = PO2./(Pm+PO2);
hyaloid = hyaloid(1:length(c1));

g1 = ( alpha10 + alpha11*q1 + alpha12*choroid + alpha13*hyaloid ).*c1.*(1 - (c1+c2)/cmax) ...
    - ( beta0 + beta1*q2 + beta2*choroid + beta3*hyaloid + beta4/cmax*c2 ).*c1 ...
    - eta1*c1;