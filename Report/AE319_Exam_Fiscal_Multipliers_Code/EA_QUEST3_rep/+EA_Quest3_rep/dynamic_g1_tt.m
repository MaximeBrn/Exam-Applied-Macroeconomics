function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 41);

T = EA_Quest3_rep.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(31) = getPowerDeriv(T(28),params(100),1);
T(32) = exp(y(98))*(-((-params(49))/((1+y(65)-params(47))*(1+y(65)-params(47)))));
T(33) = getPowerDeriv(T(1),(-params(99)),1);
T(34) = getPowerDeriv(T(1),1-params(99),1);
T(35) = getPowerDeriv(T(15)/(1+y(87)-params(47)),params(117),1);
T(36) = getPowerDeriv(exp(y(93))-params(50)*exp(y(11)),params(58),1);
T(37) = getPowerDeriv(T(4),1-params(99),1);
T(38) = getPowerDeriv(T(4),(-params(99)),1);
T(39) = getPowerDeriv(exp(y(93))-params(50)*exp(y(11)),params(58)-1,1);
T(40) = getPowerDeriv(T(25),params(101),1);
T(41) = getPowerDeriv(T(12),1-params(117),1);

end
