function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 30);

T(1) = exp(y(98))*(1-params(49)/(1+y(65)-params(47)));
T(2) = exp(y(142))*T(1)^(-params(99));
T(3) = params(62)*exp(y(149))*(exp(y(93))-params(50)*exp(y(11)))^params(58);
T(4) = 1-T(3);
T(5) = T(4)^(1-params(99));
T(6) = exp(y(97))^(-params(99));
T(7) = exp(y(142))*T(1)^(1-params(99));
T(8) = T(4)^(-params(99));
T(9) = (exp(y(93))-params(50)*exp(y(11)))^(params(58)-1);
T(10) = exp(y(97))^(1-params(99));
T(11) = params(62)*params(58)*T(8)*T(10);
T(12) = (y(136)*(1-params(102))+y(137)*params(102))/(exp(y(111))*(1-params(102))+exp(y(112))*params(102));
T(13) = (1+params(113))*T(12)^(1-params(117));
T(14) = (1-y(133)-params(103))/(1+params(113))*(params(106)-1)/params(106);
T(15) = T(14)/exp(y(31));
T(16) = (T(15)/(1+y(87)-params(47)))^params(117);
T(17) = (params(106)-1)/params(106)/exp(y(115));
T(18) = params(43)*params(27)/params(106)/exp(y(115));
T(19) = params(27)*params(6)*params(43)/params(106)/exp(y(115));
T(20) = 1/exp(y(115))*params(23);
T(21) = params(23)*1/exp(y(168))*(1+y(166)-params(47))/(1+y(127));
T(22) = params(94)+(1-params(94))*exp(y(108))^(1-params(101));
T(23) = exp(y(60))^params(3);
T(24) = exp(y(96))+exp(y(105))+exp(y(101))+exp(y(100));
T(25) = exp(params(77)*(y(25)-y(26))+(1-params(77))*(y(107)-y(108)));
T(26) = (1-params(94))*T(25)^params(101);
T(27) = T(26)*exp(y(108)-y(107));
T(28) = exp(params(78)*(y(3)*params(94)*params(3)-y(27))+(1-params(78))*(y(60)*params(94)*params(3)-y(109)));
T(29) = exp(y(109))*(1-params(94))*T(28)^params(100);
T(30) = exp(y(116))^params(3);

end
