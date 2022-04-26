function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
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

assert(length(T) >= 27);

T(1) = exp(y(43))*(1-params(49)/(1+y(10)-params(47)));
T(2) = exp(y(87))*T(1)^(-params(99));
T(3) = params(62)*exp(y(94))*(exp(y(38))-exp(y(38))*params(50))^params(58);
T(4) = 1-T(3);
T(5) = T(4)^(1-params(99));
T(6) = exp(y(42))^(-params(99));
T(7) = exp(y(87))*T(1)^(1-params(99));
T(8) = T(4)^(-params(99));
T(9) = (exp(y(38))-exp(y(38))*params(50))^(params(58)-1);
T(10) = exp(y(42))^(1-params(99));
T(11) = params(62)*params(58)*T(8)*T(10);
T(12) = (y(81)*(1-params(102))+y(82)*params(102))/(exp(y(56))*(1-params(102))+exp(y(57))*params(102));
T(13) = (1+params(113))*T(12)^(1-params(117));
T(14) = (1-y(78)-params(103))/(1+params(113))*(params(106)-1)/params(106);
T(15) = T(14)/exp(y(60));
T(16) = (T(15)/(1+y(32)-params(47)))^params(117);
T(17) = (params(106)-1)/params(106)/exp(y(60));
T(18) = params(43)*params(27)/params(106)/exp(y(60));
T(19) = params(27)*params(6)*params(43)/params(106)/exp(y(60));
T(20) = params(94)+(1-params(94))*exp(y(53))^(1-params(101));
T(21) = exp(y(5))^params(3);
T(22) = exp(y(41))+exp(y(50))+exp(y(46))+exp(y(45));
T(23) = (1-params(94))*exp(params(77)*(y(52)-y(53))+(y(52)-y(53))*(1-params(77)))^params(101);
T(24) = T(23)*exp(y(53)-y(52));
T(25) = exp(params(78)*(y(5)*params(94)*params(3)-y(54))+(y(5)*params(94)*params(3)-y(54))*(1-params(78)));
T(26) = exp(y(54))*(1-params(94))*T(25)^params(100);
T(27) = exp(y(61))^params(3);

end
