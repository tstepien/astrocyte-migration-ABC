function T = Tfunction(c,kappa,rbar)
% T = Tfunction(c,kappa,rbar)
%
% Tension function evaluated at c: T(c)

T = kappa * ( 1/sqrt(pi*c) - rbar);