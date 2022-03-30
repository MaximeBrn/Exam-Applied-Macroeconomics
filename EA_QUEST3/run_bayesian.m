% Run and plot IRF of replication file 


% EA_QUEST3 model

clear parameters;
clc;
close all;


% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run replication dynare file


dynare EA_Quest3_rep.mod

y_baseline =  cumsum(E_GY_E_EPS_G); % Output growth
phi_baseline = E_PHI_E_EPS_G; % Inflation
gg_baseline =  cumsum(E_GG_E_EPS_G);  % Gov consumption growth
fm_baseline=y_baseline./gg_baseline/GSN; % Formula of fiscal multiplier


% Plot baseline IRF
t = 1:1:length(y_baseline);

figure
plot(t,fm_baseline,'LineWidth',2);
xlabel('quarters','FontSize',8);
title('Government spending multipliers','FontSize',10);

% Make one parameter (GSLAG)

H=13; % Number of simulation
GSLAG_vec=[-0.6 -0.55 -0.50 -0.45 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0];
n=length(y_baseline);
res=nan(n,H);
for i=1:H
    set_param_value('GSLAG',GSLAG_vec(i));
    stoch_simul(M_,options_,oo_,{'E_GY'});
    res(:,i)=cumsum(E_GY_E_EPS_G)./cumsum(E_GG_E_EPS_G)/GSN;
end

fm_mean=mean(res,2);
fm_CI_lower=quantile(res,0.05,2); 
fm_CI_upper=quantile(res,0.95,2);
t = 1:1:n;
figure
plot(t,fm_mean,'LineWidth',2); hold on
plot(t,fm_CI_lower,'LineWidth',2); hold on
plot(t,fm_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Government spending multipliers','FontSize',10);

% Randomize one parameter (GAMIE)

H=20; % Number of simulation
GAMIE_vec=gamrnd(13.6871125370781,5.5553426476197,1,H);
n=length(y_baseline);
res=nan(n,H);
for i=1:H
    set_param_value('GAMIE',GAMIE_vec(i));
    stoch_simul(M_,options_,oo_,{'E_GY'});
    res(:,i)=cumsum(E_GY_E_EPS_G)./cumsum(E_GG_E_EPS_G)/GSN;
end

fm_mean=mean(res,2);
fm_CI_lower=quantile(res,0.05,2); 
fm_CI_upper=quantile(res,0.95,2);
t = 1:1:n;
figure
plot(t,fm_mean,'LineWidth',2); hold on
plot(t,fm_CI_lower,'LineWidth',2); hold on
plot(t,fm_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Government spending multipliers','FontSize',10);

figure
plot(t,fm_baseline,'LineWidth',2); hold on
plot(t,fm_mean,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Compare baseline and mean','FontSize',10);

% Randomize three parameters (GAMIE,GAMI2E,RHOGE)

H=100; % Number of simulations
GAMIE_vec=gamrnd(13.6871125370781,5.5553426476197,1,H);
GAMI2E_vec=gamrnd(4.67927403523814,0.239695301355207,1,H);
RHOGE_vec=betarnd(5.9456293913,13.9860816087,1,H);
n=length(y_baseline);
res=nan(n,H);
for i=1:H
    set_param_value('GAMIE',GAMIE_vec(i));
    set_param_value('GAMI2E',GAMI2E_vec(i));
    set_param_value('RHOGE',RHOGE_vec(i));
    stoch_simul(M_,options_,oo_,{'E_GY'});
    res(:,i)=cumsum(E_GY_E_EPS_G)./cumsum(E_GG_E_EPS_G)/GSN;
end

fm_mean=mean(res,2);
fm_CI_lower=quantile(res,0.05,2); 
fm_CI_upper=quantile(res,0.95,2);
t = 1:1:n;
figure
plot(t,fm_mean,'LineWidth',2); hold on
plot(t,fm_CI_lower,'LineWidth',2); hold on
plot(t,fm_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Government spending multipliers','FontSize',10);

figure
plot(t,fm_baseline,'LineWidth',2); hold on
plot(t,fm_mean,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Compare baseline and mean','FontSize',10);
