% Run and plot IRF of replication file 


% EA_QUEST3 model

s = 30; % Set seed
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
res_y=nan(n,H);
res_gg=nan(n,H);
res_phi=nan(n,H);
res_fm=nan(n,H);
for i=1:H
    set_param_value('GSLAG',GSLAG_vec(i));
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG'});
    res_y(:,i)=cumsum(E_GY_E_EPS_G);
    res_gg(:,i)=cumsum(E_GG_E_EPS_G);
    res_phi(:,i)=cumsum(E_PHI_E_EPS_G);
    res_fm(:,i)=res_y(:,i)./res_gg(:,i)/GSN;
end

y_mean=mean(res_y,2);
y_CI_lower=quantile(res_y,0.05,2); 
y_CI_upper=quantile(res_y,0.95,2);
t = 1:1:n;
figure
plot(t,y_mean,'LineWidth',2); hold on
plot(t,y_CI_lower,'LineWidth',2); hold on
plot(t,y_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Y response','FontSize',10);


fm_mean=mean(res_fm,2);
fm_CI_lower=quantile(res_fm,0.05,2); 
fm_CI_upper=quantile(res_fm,0.95,2);
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
res_y=nan(n,H);
res_gg=nan(n,H);
res_phi=nan(n,H);
res_fm=nan(n,H);
for i=1:H
    set_param_value('GAMIE',GAMIE_vec(i));
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG'});
    res_y(:,i)=cumsum(E_GY_E_EPS_G);
    res_gg(:,i)=cumsum(E_GG_E_EPS_G);
    res_phi(:,i)=cumsum(E_PHI_E_EPS_G);
    res_fm(:,i)=res_y(:,i)./res_gg(:,i)/GSN;
end

y_mean=mean(res_y,2);
y_CI_lower=quantile(res_y,0.05,2); 
y_CI_upper=quantile(res_y,0.95,2);
t = 1:1:n;
figure
plot(t,y_mean,'LineWidth',2); hold on
plot(t,y_CI_lower,'LineWidth',2); hold on
plot(t,y_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Y response','FontSize',10);


fm_mean=mean(res_fm,2);
fm_CI_lower=quantile(res_fm,0.05,2); 
fm_CI_upper=quantile(res_fm,0.95,2);
t = 1:1:n;
figure
plot(t,fm_mean,'LineWidth',2); hold on
plot(t,fm_CI_lower,'LineWidth',2); hold on
plot(t,fm_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Government spending multipliers','FontSize',10);



% figure
% plot(t,fm_baseline,'LineWidth',2); hold on
% plot(t,fm_mean,'LineWidth',2); hold on
% xlabel('quarters','FontSize',8);
% title('Compare baseline and mean','FontSize',10);
% 

% Randomize three parameters (GAMIE,GAMI2E,RHOGE)

H=20; % Number of simulations
rng(s);
GAMIE_vec=gamrnd(13.6871125370781,5.5553426476197,1,H);
GAMI2E_vec=gamrnd(4.67927403523814,0.239695301355207,1,H);
RHOGE_vec=betarnd(5.9456293913,13.9860816087,1,H);
n=length(y_baseline);
res_y=nan(n,H);
res_gg=nan(n,H);
res_phi=nan(n,H);
res_fm=nan(n,H);
for i=1:H
    set_param_value('GAMIE',GAMIE_vec(i));
    set_param_value('GAMI2E',GAMI2E_vec(i));
    set_param_value('RHOGE',RHOGE_vec(i));
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG'});
    res_y(:,i)=cumsum(E_GY_E_EPS_G);
    res_gg(:,i)=cumsum(E_GG_E_EPS_G);
    res_phi(:,i)=cumsum(E_PHI_E_EPS_G);
    res_fm(:,i)=res_y(:,i)./res_gg(:,i)/GSN;
end

y_mean=mean(res_y,2);
y_CI_lower=quantile(res_y,0.05,2); 
y_CI_upper=quantile(res_y,0.95,2);
t = 1:1:n;
figure
plot(t,y_mean,'LineWidth',2); hold on
plot(t,y_CI_lower,'LineWidth',2); hold on
plot(t,y_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Y response','FontSize',10);

gg_mean=mean(res_gg,2);
gg_CI_lower=quantile(res_gg,0.05,2); 
gg_CI_upper=quantile(res_gg,0.95,2);
t = 1:1:n;
figure
plot(t,gg_mean,'LineWidth',2); hold on
plot(t,gg_CI_lower,'LineWidth',2); hold on
plot(t,gg_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('GG response','FontSize',10);

fm_mean=mean(res_fm,2);
fm_CI_lower=quantile(res_fm,0.05,2); 
fm_CI_upper=quantile(res_fm,0.95,2);
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

figure
horizon=29;
t_small=t(1:horizon);
plot(t_small,fm_baseline(1:horizon),'LineWidth',2); hold on
plot(t_small,fm_mean(1:horizon),'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Compare baseline and mean','FontSize',10);

% Randomize (almost) all parameters 

H=10; % Number of simulations
rng(s);
A2E_vec=betarnd(11.9122825378418,251.052011895752,1,H);
GAMIE_vec=   gamrnd(13.6871125370781,5.5553426476197,1,H);
GAMI2E_vec   =   gamrnd(4.67927403523814,0.239695301355207,1,H);
GAMLE_vec   =   gamrnd(22.5285830778624,2.58375326130466,1,H);
GAMPE_vec   =   gamrnd(34.7633401512157,1.76742222504333,1,H);
GAMPME_vec   =   gamrnd(3.40697211803127,0.492578143248719,1,H);
GAMPXE_vec   =   gamrnd(2.40760377777554,10.8528655093496,1,H);
GAMWE_vec   =   gamrnd(2.44564013067148,0.528246156823283,1,H);
HABE_vec   =   betarnd(81.080357597323,62.8322401970025,1,H);
HABLE_vec   =   betarnd(19.8492375075171,4.68931794744286,1,H);
IG1E_vec   =   betarnd(1.77116614802584,10.0602710465355,1,H);
IGSLAG_vec   =   betarnd(13.365,16.5009217877095,1,H);
ILAGE_vec   =   betarnd(333.88212714256,36.7274045952133,1,H);
KAPPAE_vec   =   gamrnd(18.7634597430808,0.102454452767374,1,H);
RHOCE_vec   =   betarnd(81.3292104751508,7.61349564378053,1,H);
RHOETA_vec   =   betarnd(1.68669468122152,13.7169097134955,1,H);
RHOETAM_vec   =   betarnd(149.482871932629,6.92904805547293,1,H);
RHOETAX_vec   =   betarnd(27.0550083713561,6.30916522755388,1,H);
RHOGE_vec   =   betarnd(5.9456293913,13.9860816087,1,H);
RHOIG_vec   =   betarnd(15.9853198285919,2.75479720375499,1,H);
RHOLE_vec   =   betarnd(93.0308739765042,2.38540702503857,1,H);
RHOL0_vec   =   betarnd(163.236723064735,11.6472742191037,1,H);
RHOPCPM_vec   =   betarnd(6.79701820447584,3.42098871746619,1,H);
RHOPWPX_vec   =   betarnd(7.55065652002992,27.4222778015538,1,H);
RHORPE_vec   =   betarnd(143.276939711565,2.30011750400604,1,H);
RHORPK_vec   =   betarnd(130.420074851259,12.1466882130819,1,H);
RHOUCAP0_vec   =   betarnd(124.150243398439,6.30078465498012,1,H);
RPREME_vec   =   betarnd(7.13850986121256,349.786983199416,1,H);
RPREMK_vec   =   betarnd(86.5944164201183,3447.87155991124,1,H);
SE_vec   =   betarnd(270.226910453978,44.4294827155352,1,H);
SFPE_vec   =   betarnd(29.5031576539166,4.35403497164756,1,H);
SFPME_vec   =   betarnd(8.76171924619187,3.14117335833451,1,H);
SFPXE_vec   =   betarnd(29.9690857153074,2.67697715539782,1,H);
SFWE_vec   =   betarnd(4.75838075897478,1.39257678882095,1,H);
SIGC_vec   =   gamrnd(25.3852365844692,0.16136150578585,1,H);
SIGEXE_vec   =   gamrnd(62.795719140625,0.0403817335752031,1,H);
SIGIME_vec   =   gamrnd(30.126530741068,0.0389158648925282,1,H);
SLC_vec   =   betarnd(13.6959979217823,25.3573180798781,1,H);
RHOTR_vec   =   betarnd(54.6695232771421,8.63469543191545,1,H);
TYE1_vec   =   betarnd(7.60691935607018,10.1912073544356,1,H);
TYE2_vec   =   betarnd(7.28637478137341,85.7707744060264,1,H);
WRLAG_vec   =   betarnd(2.72512582214576,7.54673931975307,1,H);
n=length(y_baseline);
res_y=nan(n,H);
res_gg=nan(n,H);
res_phi=nan(n,H);
res_fm=nan(n,H);
for i=1:H
    set_param_value('A2E',A2E_vec(i));
    set_param_value('GAMIE',GAMIE_vec(i));
    set_param_value('GAMI2E',GAMI2E_vec (i));
    set_param_value('GAMLE',GAMLE_vec (i));
    set_param_value('GAMPE',GAMPE_vec (i));
    set_param_value('GAMPME',GAMPME_vec (i));
    set_param_value('GAMPXE',GAMPXE_vec (i));
    set_param_value('GAMWE',GAMWE_vec (i));
    set_param_value('HABE',HABE_vec (i));
    set_param_value('HABLE',HABLE_vec (i));
    set_param_value('IG1E',IG1E_vec (i));
    set_param_value('IGSLAG',IGSLAG_vec (i));
    set_param_value('ILAGE',ILAGE_vec (i));
    set_param_value('KAPPAE',KAPPAE_vec (i));
    set_param_value('RHOCE',RHOCE_vec (i));
    set_param_value('RHOETA',RHOETA_vec (i));
    set_param_value('RHOETAM',RHOETAM_vec (i));
    set_param_value('RHOETAX',RHOETAX_vec (i));
    set_param_value('RHOGE',RHOGE_vec (i));
    set_param_value('RHOIG',RHOIG_vec (i));
    set_param_value('RHOLE',RHOLE_vec (i));
    set_param_value('RHOL0',RHOL0_vec (i));
    set_param_value('RHOPCPM',RHOPCPM_vec (i));
    set_param_value('RHOPWPX',RHOPWPX_vec (i));
    set_param_value('RHORPE',RHORPE_vec (i));
    set_param_value('RHORPK',RHORPK_vec (i));
    set_param_value('RHOUCAP0',RHOUCAP0_vec (i));
    set_param_value('RPREME',RPREME_vec (i));
    set_param_value('RPREMK',RPREMK_vec (i));
    set_param_value('SE',SE_vec (i));
    set_param_value('SFPE',SFPE_vec (i));
    set_param_value('SFPME',SFPME_vec (i));
    set_param_value('SFPXE',SFPXE_vec (i));
    set_param_value('SFWE',SFWE_vec (i));
    set_param_value('SIGC',SIGC_vec (i));
    set_param_value('SIGEXE',SIGEXE_vec (i));
    set_param_value('SIGIME',SIGIME_vec (i));
    set_param_value('SLC',SLC_vec (i));
    set_param_value('RHOTR',RHOTR_vec (i));
    set_param_value('TYE1',TYE1_vec (i));
    set_param_value('TYE2',TYE2_vec (i));
    set_param_value('WRLAG',WRLAG_vec (i));
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG'});
    res_y(:,i)=cumsum(E_GY_E_EPS_G);
    res_gg(:,i)=cumsum(E_GG_E_EPS_G);
    res_phi(:,i)=cumsum(E_PHI_E_EPS_G);
    res_fm(:,i)=res_y(:,i)./res_gg(:,i)/GSN;
end

y_mean=mean(res_y,2);
y_CI_lower=quantile(res_y,0.05,2); 
y_CI_upper=quantile(res_y,0.95,2);
t = 1:1:n;
figure
plot(t,y_mean,'LineWidth',2); hold on
plot(t,y_CI_lower,'LineWidth',2); hold on
plot(t,y_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Y response','FontSize',10);


gg_mean=mean(res_gg,2);
gg_CI_lower=quantile(res_gg,0.05,2); 
gg_CI_upper=quantile(res_gg,0.95,2);
t = 1:1:n;
figure
plot(t,gg_mean,'LineWidth',2); hold on
plot(t,gg_CI_lower,'LineWidth',2); hold on
plot(t,gg_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('GG response','FontSize',10);


horizon=30;
t_small=t(1:horizon);
fm_mean=mean(res_fm,2);
fm_CI_lower=quantile(res_fm,0.05,2); 
fm_CI_upper=quantile(res_fm,0.95,2);
t = 1:1:n;
figure
plot(t_small,fm_mean(1:horizon),'LineWidth',2); hold on
plot(t_small,fm_CI_lower(1:horizon),'LineWidth',2); hold on
plot(t_small,fm_CI_upper(1:horizon),'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Government spending multipliers','FontSize',10);

figure

plot(t_small,fm_baseline(1:horizon),'LineWidth',2); hold on
plot(t_small,fm_mean(1:horizon),'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('Compare baseline and mean','FontSize',10);


