function Tp = Tderivative(c,kappa,cmin,rbar)
% Tp = Tderivative(c,kappa,cmin,rbar)
%
% Derivative of tension function evaluated at c: T'(c)

Tp = -kappa/(2 * sqrt(pi) * c^(3/2));
