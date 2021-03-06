function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = EA_Quest3_rep.static_g1_tt(T, y, x, params);
end
g1 = zeros(107, 107);
g1(1,36)=1-params(53);
g1(1,58)=(-((1-params(53))*params(111)));
g1(1,64)=(-((1-params(53))*params(107)));
g1(1,95)=(-1);
g1(1,101)=(-((1-params(53))*(1-params(107))));
g1(2,10)=(-(T(5)*exp(y(87))*T(29)*T(30)));
g1(2,38)=(-(T(2)*T(35)));
g1(2,43)=(-(T(5)*exp(y(87))*T(1)*T(30)));
g1(2,56)=exp(y(56));
g1(2,87)=(-(T(2)*T(5)));
g1(2,94)=(-(T(2)*T(34)*(-T(3))));
g1(3,38)=(-(T(6)*T(35)));
g1(3,42)=(-(T(5)*exp(y(42))*getPowerDeriv(exp(y(42)),(-params(99)),1)));
g1(3,57)=exp(y(57));
g1(3,94)=(-(T(6)*T(34)*(-T(3))));
g1(4,10)=(-(exp(y(94))*T(9)*params(62)*params(58)*T(8)*exp(y(87))*T(29)*T(31)));
g1(4,38)=(-(exp(y(94))*(T(9)*params(62)*params(58)*T(7)*T(33)*T(36)+params(62)*params(58)*T(7)*T(8)*T(37))));
g1(4,43)=(-(exp(y(94))*T(9)*params(62)*params(58)*T(8)*exp(y(87))*T(1)*T(31)));
g1(4,81)=1;
g1(4,87)=(-(exp(y(94))*params(62)*params(58)*T(7)*T(8)*T(9)));
g1(4,94)=(-(exp(y(94))*params(62)*params(58)*T(7)*T(8)*T(9)+exp(y(94))*T(9)*params(62)*params(58)*T(7)*T(36)*(-T(3))));
g1(5,38)=(-(exp(y(94))*(T(11)*T(37)+T(9)*params(62)*params(58)*T(10)*T(33)*T(36))));
g1(5,42)=(-(exp(y(94))*T(9)*params(62)*params(58)*T(8)*exp(y(42))*getPowerDeriv(exp(y(42)),1-params(99),1)));
g1(5,82)=1;
g1(5,94)=(-(exp(y(94))*T(9)*T(11)+exp(y(94))*T(9)*params(62)*params(58)*T(10)*T(36)*(-T(3))));
g1(6,29)=(-1);
g1(6,36)=(-1);
g1(6,64)=1;
g1(7,29)=(-1);
g1(7,32)=(-params(99));
g1(7,63)=(-params(99));
g1(7,64)=params(99);
g1(8,42)=exp(y(42))*(1+params(113));
g1(8,73)=1;
g1(8,76)=(-y(85));
g1(8,78)=y(85);
g1(8,85)=(-(1-y(78)-params(103)+y(76)));
g1(9,38)=(-exp(y(38)-y(60)));
g1(9,60)=exp(y(38)-y(60));
g1(9,85)=1;
g1(10,41)=exp(y(41));
g1(10,42)=(-(exp(y(42))*params(44)*params(102)));
g1(10,43)=(-(exp(y(43))*(1-params(44)*params(102))));
g1(11,32)=T(13)*(-T(15))/((1+y(32)-params(47))*(1+y(32)-params(47)))*T(32);
g1(11,56)=T(16)*(1+params(113))*(-((y(81)*(1-params(102))+y(82)*params(102))*exp(y(56))*(1-params(102))))/((exp(y(56))*(1-params(102))+exp(y(57))*params(102))*(exp(y(56))*(1-params(102))+exp(y(57))*params(102)))*T(39);
g1(11,57)=T(16)*(1+params(113))*T(39)*(-((y(81)*(1-params(102))+y(82)*params(102))*exp(y(57))*params(102)))/((exp(y(56))*(1-params(102))+exp(y(57))*params(102))*(exp(y(56))*(1-params(102))+exp(y(57))*params(102)));
g1(11,60)=T(13)*T(32)*(-(T(14)*exp(y(60))))/(exp(y(60))*exp(y(60)))/(1+y(32)-params(47))-((1-y(78)-params(103))*(-(exp(y(60))*(params(106)-1)/params(106)))/(exp(y(60))*exp(y(60)))+(y(83)-params(28)-params(47)-(1-params(98))*(y(63)-params(28)))*(-(exp(y(60))*params(43)*params(27)/params(106)))/(exp(y(60))*exp(y(60)))-(y(83)-params(28)-params(47)-(1-params(98))*(y(63)-params(28)))*(-(exp(y(60))*params(27)*params(6)*params(43)/params(106)))/(exp(y(60))*exp(y(60))));
g1(11,63)=(-(T(18)*(-(1-params(98)))-T(19)*(-(1-params(98)))));
g1(11,78)=T(17)+T(13)*T(32)*(params(106)-1)*(-1)/(1+params(113))/params(106)/exp(y(60))/(1+y(32)-params(47));
g1(11,81)=T(16)*(1+params(113))*T(39)*(1-params(102))/(exp(y(56))*(1-params(102))+exp(y(57))*params(102));
g1(11,82)=T(16)*(1+params(113))*T(39)*params(102)/(exp(y(56))*(1-params(102))+exp(y(57))*params(102));
g1(11,83)=(-(T(18)-T(19)));
g1(12,6)=(-((1+y(51))*params(4)/exp(y(38))));
g1(12,38)=(-((1+y(51))*(-(exp(y(38))*y(6)*params(4)))/(exp(y(38))*exp(y(38)))));
g1(12,51)=(-(y(6)*params(4)/exp(y(38))));
g1(12,60)=(-(exp(y(60))*(1+y(100))))/(exp(y(60))*exp(y(60)));
g1(12,100)=1/exp(y(60));
g1(13,16)=params(21)-params(21)/(1+y(36));
g1(13,36)=(-((y(16)-params(47)-params(29))*(-params(21))/((1+y(36))*(1+y(36)))));
g1(13,48)=params(22)*exp(y(48));
g1(13,71)=(-1);
g1(14,6)=exp(y(59))*(1-params(105))*(1-params(4));
g1(14,59)=y(6)*(1-params(105))*(1-params(4))*exp(y(59));
g1(14,65)=(-y(71));
g1(14,71)=(-(1-(params(29)+1-y(72)-params(10)-params(93)-y(98)-y(65))));
g1(14,72)=(-y(71));
g1(14,79)=(-((1-params(105))*(params(1)+params(2)*2*(y(79)-params(116)))));
g1(14,98)=(-y(71));
g1(15,6)=(1-params(4))*exp(y(59));
g1(15,59)=exp(y(59))*y(6)*(1-params(4));
g1(15,79)=(-(params(1)+(y(79)-params(116))*2*params(2)+y(79)*2*params(2)));
g1(16,6)=1;
g1(16,63)=params(43)*params(24)*(params(6)-1);
g1(16,88)=1;
g1(17,36)=(-(getPowerDeriv(1+y(36),(-params(118)),1)));
g1(17,62)=1;
g1(18,21)=1;
g1(18,48)=(-exp(y(48)));
g1(19,22)=1;
g1(19,49)=(-exp(y(49)));
g1(20,21)=1;
g1(20,48)=(-1);
g1(20,50)=1;
g1(20,59)=1;
g1(21,21)=(-(1-params(4)));
g1(21,22)=(-(1-params(5)));
g1(21,23)=(-(params(4)*(1+params(60))));
g1(21,26)=(-params(4));
g1(21,30)=(-(1-params(4)));
g1(21,32)=1;
g1(22,36)=(-1);
g1(22,63)=1;
g1(22,72)=1;
g1(23,38)=(-params(4));
g1(23,39)=params(4);
g1(23,58)=1;
g1(23,79)=(-((1-params(4))*1/y(79)));
g1(23,80)=(-((1-params(4))*(-(1/y(80)))));
g1(24,38)=(-(1-params(70)));
g1(24,39)=1-params(70);
g1(25,79)=(-(1-params(81)));
g1(25,80)=1-params(81);
g1(26,52)=exp(y(52));
g1(26,53)=(-((1-params(94))*exp(y(53))*getPowerDeriv(exp(y(53)),1-params(101),1)*getPowerDeriv(T(20),1/(1-params(101)),1)));
g1(27,5)=(-((1+y(89)+params(43)*params(25)*(params(28)+params(6)*(params(96)*y(66)+y(66)*(1-params(96))-params(28))-y(66)))*exp(y(5))*getPowerDeriv(exp(y(5)),params(3),1)));
g1(27,53)=exp(y(53));
g1(27,66)=(-(T(21)*params(43)*params(25)*(params(6)-1)));
g1(27,89)=(-T(21));
g1(28,54)=exp(y(54));
g1(28,69)=(-(params(43)*params(26)*(params(6)-1)));
g1(28,90)=(-1);
g1(29,41)=(-exp(y(41)));
g1(29,45)=(-exp(y(45)));
g1(29,46)=(-exp(y(46)));
g1(29,50)=(-exp(y(50)));
g1(29,74)=(-1);
g1(30,44)=(-exp(y(44)));
g1(30,47)=exp(y(47));
g1(30,74)=1;
g1(30,91)=(-1);
g1(31,41)=(-(exp(y(41))*T(24)));
g1(31,45)=(-(exp(y(45))*T(24)));
g1(31,46)=(-(exp(y(46))*T(24)));
g1(31,47)=exp(y(47));
g1(31,50)=(-(exp(y(50))*T(24)));
g1(31,52)=(-(T(22)*(exp(y(53)-y(52))*(1-params(94))*exp(params(77)*(y(52)-y(53))+(y(52)-y(53))*(1-params(77)))*T(38)+T(23)*(-exp(y(53)-y(52))))));
g1(31,53)=(-(T(22)*(T(24)+exp(y(53)-y(52))*(1-params(94))*T(38)*exp(params(77)*(y(52)-y(53))+(y(52)-y(53))*(1-params(77)))*((-params(77))-(1-params(77))))));
g1(32,5)=(-(T(27)*exp(y(54))*(1-params(94))*T(25)*(params(78)*params(94)*params(3)+params(94)*params(3)*(1-params(78)))*T(28)));
g1(32,44)=exp(y(44));
g1(32,54)=(-(T(27)*(T(26)+exp(y(54))*(1-params(94))*T(28)*T(25)*((-params(78))-(1-params(78))))));
g1(32,61)=(-(T(26)*exp(y(61))*getPowerDeriv(exp(y(61)),params(3),1)));
g1(33,2)=params(92);
g1(33,11)=(-1);
g1(33,36)=1;
g1(33,37)=(-1);
g1(33,97)=(-1);
g1(34,2)=1-(1+y(36)-y(63)-y(32)-params(30));
g1(34,32)=y(2);
g1(34,36)=(-y(2));
g1(34,63)=y(2);
g1(34,74)=(-1);
g1(35,32)=exp(y(40));
g1(35,38)=(-(y(76)*exp(y(38)-y(60))));
g1(35,40)=exp(y(40))-exp(y(40))*(1+y(72)-y(32)-params(30));
g1(35,41)=params(113)*exp(y(41));
g1(35,45)=(-exp(y(45)));
g1(35,46)=(-exp(y(46)));
g1(35,60)=(-(y(76)*(-exp(y(38)-y(60)))));
g1(35,72)=(-exp(y(40)));
g1(35,73)=1;
g1(35,76)=(-exp(y(38)-y(60)));
g1(35,78)=y(85);
g1(35,85)=(-(params(105)-(y(78)+params(103))));
g1(36,14)=1-params(32);
g1(36,45)=(-(params(34)*params(33)));
g1(36,92)=(-(params(34)*params(40)));
g1(37,17)=1-params(51);
g1(37,46)=(-(params(37)*params(52)));
g1(37,93)=(-(params(37)*params(41)));
g1(38,16)=(-1);
g1(38,17)=1;
g1(39,38)=(-(params(38)*params(108)*(-exp(y(38)))));
g1(39,76)=1;
g1(39,99)=(-1);
g1(40,40)=(-(exp(y(40))*params(7)));
g1(41,58)=(-(params(114)*params(115)*params(39)));
g1(41,78)=1;
g1(42,38)=(-(y(76)*exp(y(38)-y(60))));
g1(42,60)=(-(y(76)*(-exp(y(38)-y(60)))));
g1(42,76)=(-exp(y(38)-y(60)));
g1(42,77)=1;
g1(43,38)=(-(y(76)*exp(y(38)-y(60))));
g1(43,60)=(-(y(76)*(-exp(y(38)-y(60)))));
g1(43,73)=1;
g1(43,75)=1;
g1(43,76)=(-exp(y(38)-y(60)));
g1(44,78)=y(85);
g1(44,85)=(-(1-y(78)-params(103)));
g1(44,86)=1;
g1(45,35)=(-params(84));
g1(45,37)=1-params(82);
g1(45,68)=(-params(83));
g1(46,35)=(-params(87));
g1(46,37)=(-params(85));
g1(46,68)=1-params(86);
g1(47,35)=1-params(90);
g1(47,37)=(-params(88));
g1(47,61)=(-params(91));
g1(47,68)=(-params(89));
g1(48,32)=1;
g1(48,35)=(-1);
g1(49,26)=1;
g1(50,51)=1-params(72);
g1(51,65)=1;
g1(51,96)=(-1);
g1(52,87)=1-params(63);
g1(53,88)=1-params(64);
g1(54,89)=1-params(65);
g1(55,90)=1-params(66);
g1(56,91)=1-params(67);
g1(57,92)=1-params(68);
g1(58,93)=1-params(69);
g1(59,94)=1-params(71);
g1(60,95)=1;
g1(61,96)=1-(params(76)+params(75)+params(73)+params(74));
g1(62,97)=1-params(79);
g1(63,98)=1-params(80);
g1(64,100)=1;
g1(65,99)=1-params(110);
g1(66,101)=1;
g1(67,7)=1;
g1(67,32)=(-1);
g1(67,63)=(-1);
g1(67,64)=1;
g1(68,9)=1;
g1(68,32)=(-1);
g1(68,63)=(-1);
g1(68,64)=1;
g1(69,10)=1;
g1(69,32)=(-1);
g1(69,63)=(-1);
g1(69,64)=1;
g1(70,11)=1;
g1(70,63)=(-1);
g1(70,68)=1;
g1(71,12)=1;
g1(71,32)=(-1);
g1(71,63)=(-1);
g1(71,69)=1;
g1(72,14)=1;
g1(72,32)=(-1);
g1(72,63)=(-1);
g1(72,64)=1;
g1(73,16)=1;
g1(73,21)=(-1);
g1(74,17)=1;
g1(74,22)=(-1);
g1(75,19)=1;
g1(75,32)=(-1);
g1(75,63)=(-1);
g1(75,66)=1;
g1(76,21)=(-1);
g1(76,32)=1;
g1(76,65)=1;
g1(77,23)=1;
g1(78,25)=1;
g1(78,32)=(-1);
g1(78,63)=(-1);
g1(79,26)=(-params(4));
g1(79,27)=1;
g1(79,30)=(-(1-params(4)));
g1(80,23)=(-1);
g1(80,28)=1;
g1(80,84)=(-1);
g1(81,30)=1;
g1(82,31)=1;
g1(83,32)=1;
g1(83,34)=(-1);
g1(84,1)=1;
g1(84,40)=(-exp(y(40)));
g1(85,4)=1;
g1(86,3)=1;
g1(86,42)=(-exp(y(42)));
g1(87,24)=1;
g1(87,45)=(-exp(y(45)));
g1(88,55)=1;
g1(88,77)=(-(1/y(77)));
g1(89,63)=(-1);
g1(89,64)=1;
g1(90,63)=(-1);
g1(90,66)=1;
g1(91,63)=(-1);
g1(91,69)=1;
g1(92,32)=1;
g1(92,63)=1;
g1(92,83)=(-1);
g1(93,63)=(-1);
g1(93,83)=1;
g1(93,84)=(-1);
g1(94,32)=(-1);
g1(94,33)=1;
g1(95,7)=(-1);
g1(95,8)=1;
g1(96,16)=(-1);
g1(96,18)=1;
g1(97,14)=(-1);
g1(97,15)=1;
g1(98,12)=(-1);
g1(98,13)=1;
g1(99,19)=(-1);
g1(99,20)=1;
g1(100,66)=(-1);
g1(100,67)=1;
g1(101,69)=(-1);
g1(101,70)=1;
g1(102,41)=(-1);
g1(102,52)=1;
g1(102,102)=1;
g1(103,45)=(-1);
g1(103,52)=1;
g1(103,103)=1;
g1(104,38)=(-1);
g1(104,60)=1;
g1(104,104)=1;
g1(105,96)=(-1);
g1(105,105)=1;
g1(106,96)=(-1);
g1(106,106)=1;
g1(107,96)=(-1);
g1(107,107)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
