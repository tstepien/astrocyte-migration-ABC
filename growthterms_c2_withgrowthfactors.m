function g = growthterms_c2_withgrowthfactors(c1,c2,q1,q2,PO2,Pm,alpha20,...
    alpha21,alpha22,beta0,beta1,beta2,beta3,eta2,cmax,hyaloid)
% g = growthterms_c2_withgrowthfactors(c1,c2,q1,q2,PO2,Pm,alpha20,...
%     alpha21,alpha22,beta0,beta1,beta2,beta3,eta2,cmax,hyaloid)
%
% growth terms just IPAs (c2): g2

hyaloid = hyaloid(1:length(c1));
choroid = PO2./(Pm+PO2);

g = ( alpha20 + alpha21*(hyaloid + choroid) + alpha22*q1 ).*c2.*(1 - (c1+c2)/cmax) ...
    + ( beta0 + beta1*c2 + beta2*(hyaloid + choroid) + beta3*q2 ).*c1 ...
    - eta2.*c2;