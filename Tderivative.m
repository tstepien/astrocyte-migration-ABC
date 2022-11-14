function Tp = Tderivative(c,kappa)
% Tp = Tderivative(c,kappa)
%
% Derivative of tension function evaluated at c: T'(c)

Tp = -kappa./(2 * sqrt(pi) * c.^(3/2));
