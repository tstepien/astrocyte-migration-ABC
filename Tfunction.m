function T = Tfunction(c,kappa,cmin,rbar)
% T = Tfunction(c,kappa,cmin,rbar)
%
% Tension function evaluated at c: T(c)

T = kappa * ( 1/sqrt(pi*c) - rbar);
