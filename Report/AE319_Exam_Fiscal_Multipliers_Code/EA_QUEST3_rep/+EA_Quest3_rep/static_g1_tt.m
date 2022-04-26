function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 39);

T = EA_Quest3_rep.static_resid_tt(T, y, x, params);

T(28) = getPowerDeriv(T(25),params(100),1);
T(29) = exp(y(43))*(-((-params(49))/((1+y(10)-params(47))*(1+y(10)-params(47)))));
T(30) = getPowerDeriv(T(1),(-params(99)),1);
T(31) = getPowerDeriv(T(1),1-params(99),1);
T(32) = getPowerDeriv(T(15)/(1+y(32)-params(47)),params(117),1);
T(33) = (-(params(62)*exp(y(94))*(exp(y(38))-exp(y(38))*params(50))*getPowerDeriv(exp(y(38))-exp(y(38))*params(50),params(58),1)));
T(34) = getPowerDeriv(T(4),1-params(99),1);
T(35) = T(33)*T(34);
T(36) = getPowerDeriv(T(4),(-params(99)),1);
T(37) = (exp(y(38))-exp(y(38))*params(50))*getPowerDeriv(exp(y(38))-exp(y(38))*params(50),params(58)-1,1);
T(38) = getPowerDeriv(exp(params(77)*(y(52)-y(53))+(y(52)-y(53))*(1-params(77))),params(101),1);
T(39) = getPowerDeriv(T(12),1-params(117),1);

end
