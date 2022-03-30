% Run and plot IRF of replication file 


% EA_QUEST3 model

clc;
close all;


% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run replication dynare file


dynare EA_Quest3_rep_rand.mod
stoch_simul (drop=200, order=1, irf = 41, ar=0, periods=0, nograph, noprint, conditional_variance_decomposition =[1 4],irf_shocks=(E_EPS_G,E_EPS_IG)) E_INOM E_PHIC E_LYGAP E_PHI E_GY E_GC E_GI E_GCNLC E_GCLC E_TBYN E_GG E_GIG E_GTR E_GL E_GWRY E_DBGYN E_R E_GE E_GTFP;
set_param_value('G1E',-0.08);
% Will have to use https://forum.dynare.org/t/loop-over-parameters/3335/2

oo_.irfs;

%%% GOVERNMENT CONSUMPTION %%%%

% Variables in Figure 1.a
yirf =  cumsum(E_GY_E_EPS_G); % Output growth
cirf =  cumsum(E_GC_E_EPS_G); % Private consumption growth
iirf =  cumsum(E_GI_E_EPS_G); % Private investment growth
cricirf =  cumsum(E_GCNLC_E_EPS_G); % Private consumption growth (ricardian HH)
cliqcirf =  cumsum(E_GCLC_E_EPS_G); % Private consumption growth (non-Ricardian HH)
tbynirf =  E_TBYN_E_EPS_G; % Nominal trade balance to GDP share
ggirf =  cumsum(E_GG_E_EPS_G);  % Gov consumption growth
gigirf =  cumsum(E_GIG_E_EPS_G); % Gov investment growth
gtrirf =  cumsum(E_GTR_E_EPS_G); % Transfers growth

% Variables in Figure 1.b
glirf =  cumsum(E_GL_E_EPS_G); % Employment rate growth
rwirf =  cumsum(E_GWRY_E_EPS_G+E_GY_E_EPS_G); % (real wages/real gdp) growth 
dbgynirf =  E_DBGYN_E_EPS_G; % Deficit
ygapirf =  E_LYGAP_E_EPS_G; % Log ouput gap
phiirf = E_PHI_E_EPS_G; % Inflation
infirf = E_PHIC_E_EPS_G; % Inflation of C deflator
inomirf = E_INOM_E_EPS_G; % Nominal interest rate
rirf =  E_R_E_EPS_G; % Real interest rate
geirf =  cumsum(E_GE_E_EPS_G); % Nominal exchange rate

fiscalmultiplier=yirf./ggirf/GSN; % Formula of fiscal multiplier

spending_fm_short=fiscalmultiplier(1);
spending_fm_long=fiscalmultiplier(4);

 % Go back to original path
cd('..');


% Plot replicated IRF
t = 1:1:length(infirf);
zeroline = ones(length(t),1)*0;

figure % replicate 1.a
subplot(3,3,1); % Output 
plot(t,yirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -5*10^(-4) 10*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$Y_t$$','interpreter','latex','FontSize',10);

subplot(3,3,2); % Consumption
plot(t,cirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -6*10^(-4) 0*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$C_t$$','interpreter','latex','FontSize',10);

subplot(3,3,3); % Investment
plot(t,iirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -6*10^(-4) 0*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$I_t$$','interpreter','latex','FontSize',10);

subplot(3,3,4); % Consumption (Ricardian)
plot(t,cricirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -10*10^(-4) 5*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$C_t^i$$','interpreter','latex','FontSize',10);

subplot(3,3,5); % Consumption (Non-Ricardian)
plot(t,cliqcirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -2*10^(-4) 4*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$C_t^k$$','interpreter','latex','FontSize',10);

subplot(3,3,6); % Trade balance
plot(t,tbynirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -2*10^(-4) 2*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$TBYN_t$$','interpreter','latex','FontSize',10);

subplot(3,3,7); % Government consumption
plot(t,ggirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -2*10^(-3) 6*10^(-3)]);
xlabel('quarters','FontSize',8);
title('$$G_t$$','interpreter','latex','FontSize',10);

subplot(3,3,8); % Government investment
plot(t,gigirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -2*10^(-4) 6*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$GI_t$$','interpreter','latex','FontSize',10);

subplot(3,3,9); % Government transfers
plot(t,gtrirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -6*10^(-4) 2*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$TR_t$$','interpreter','latex','FontSize',10);



figure % Replicat 1.b
subplot(3,3,1); % Employment rate growth
plot(t,glirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 0*10^(-4) 4*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$L_t$$','interpreter','latex','FontSize',10);

subplot(3,3,2); % (real wages/real gdp) growth 
plot(t,rwirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -4*10^(-4) 4*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$WR_t$$','interpreter','latex','FontSize',10);

subplot(3,3,3); % Deficit
plot(t,dbgynirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -2*10^(-3) 2*10^(-3)]);
xlabel('quarters','FontSize',8);
title('$$deficit_t$$','interpreter','latex','FontSize',10);

subplot(3,3,4); % Log ouput gap
plot(t,ygapirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -5*10^(-4) 10*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$log(YGAP_t)$$','interpreter','latex','FontSize',10);

subplot(3,3,5);% Inflation
plot(t,phiirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -2*10^(-4) 4*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$\pi_t$$','interpreter','latex','FontSize',10);

subplot(3,3,6);% Inflation of C deflator
plot(t,infirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 0*10^(-4) 1.5*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$\pi_t^c$$','interpreter','latex','FontSize',10);

subplot(3,3,7); % Nominal interest rate
plot(t,inomirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 0*10^(-4) 1.5*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$inomt$$','interpreter','latex','FontSize',10);

subplot(3,3,8); % Real interest rate
plot(t,rirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -2*10^(-4) 4*10^(-4)]);
xlabel('quarters','FontSize',8);
title('$$R_t$$','interpreter','latex','FontSize',10);

subplot(3,3,9); % Nominal exchange rate
plot(t,geirf,'LineWidth',2);hold on;
plot(t,zeroline,'LineWidth',1);
axis([0 42 -1*10^(-3) 2*10^(-3)]);
xlabel('quarters','FontSize',8);
title('$$e_t$$','interpreter','latex','FontSize',10);

figure
plot(t,fiscalmultiplier,'LineWidth',2);
xlabel('quarters','FontSize',8);
title('Government spending multipliers','FontSize',10);

GAMWE
G1E