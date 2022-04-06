% Run and plot IRF of replication file 


% EA_QUEST3 model

s = 30; % enter a seed (for random draws)
clc;
close all;


% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);


%% Replication of the paper results

% Run the .mod file for the first time
    % Parameters value = posterior mean
    % We replicate what is done in the paper
dynare EA_Quest3_rep.mod

% Store reponses to a government consumption shock (E_EPS_G)
    % We call the response baseline
    % Therefore baseline responses = paper responses
    
y_baseline = cumsum(E_GY_E_EPS_G); % Response of output growth
phi_baseline = E_PHI_E_EPS_G; % Response of inflation
gg_baseline = cumsum(E_GG_E_EPS_G); % Response of government consumption growth

% Store the fiscal multipliers
    % We should recover the 0.73 for Q1 and 0.45 for Q4
fm_baseline=y_baseline./gg_baseline/GSN; % 1st formula
GY0_vec = ones(1, length(y_baseline)) * GY0;
fm2_baseline=(cumprod(1+E_GY_E_EPS_G+GY0)-cumprod(1+GY0_vec)')./(cumprod(1+E_GG_E_EPS_G+GY0)-cumprod(1+GY0_vec)')/GSN; % 2 nd formula
    % Check
    fm_baseline(1);
    fm_baseline(4);
    fm2_baseline(1);
    fm2_baseline(4);

% Plot the IRF of fiscal multipliers

t = 1:1:length(y_baseline); % x-axis

figure
subplot(2,1,1) % 1st formula
plot(t,fm_baseline,'LineWidth',2);
xlabel('quarters','FontSize',8);
title('Government consumption multipliers 1','FontSize',10);
subplot(2,1,2) % 2nd formula
plot(t,fm2_baseline,'LineWidth',2);
xlabel('quarters','FontSize',8);
title('Government spending multipliers 2','FontSize',10);

%% Reponse with uncertainty bounds
    % For each parameter, we draw a value
    % Either draw from prior distribution
    % Either draw from a normal distribution
    % Paper Table 1, posterior mean and st.dev of parameters are reported
    % From mean and variance we deduce the distribution inputs
    % See Excel file for computation of the beta, gamma, k and theta
    % We are now able to draw a sequence of value for each parameter
    % We can simulate the .mod file for each vector of parameters value
    % We obtain several responses to shocks
    % We can a Bayesian confidence interval around the mean response

H=100; % Set the number of the simulations
rng(s); % Specify the seed

% Draws
    % For each parameter, draw H values from the distribution
    % See the Excel file for the distribution inputs
    % Store the H values in a vector called VarName_vec
        
% Draw vectors from normal law
G1E_vec  =   normrnd(-0.0754,0.1066,1,H);
GSLAG_vec  =   normrnd(-0.4227,0.1041,1,H);
GVECM_vec  =   normrnd(-0.1567,0.0442,1,H);
IGVECM_vec  =   normrnd(-0.1222,0.0461,1,H);
A2E_vec  =   normrnd(0.0453,0.0128,1,H);
GAMIE_vec  =   normrnd(76.0366,20.5526,1,H);
GAMI2E_vec  =   normrnd(1.1216,0.5185,1,H);
GAMLE_vec  =   normrnd(58.2083,12.2636,1,H);
GAMPE_vec  =   normrnd(61.4415,10.4208,1,H);
GAMPME_vec  =   normrnd(1.6782,0.9092,1,H);
GAMPXE_vec  =   normrnd(26.1294,16.8398,1,H);
GAMWE_vec  =   normrnd(1.2919,0.8261,1,H);
HABE_vec  =   normrnd(0.5634,0.0412,1,H);
HABLE_vec  =   normrnd(0.8089,0.0778,1,H);
IG1E_vec  =   normrnd(0.1497,0.0996,1,H);
IGSLAG_vec  =   normrnd(0.4475,0.0895,1,H);
ILAGE_vec  =   normrnd(0.9009,0.0155,1,H);
KAPPAE_vec  =   normrnd(1.9224,0.4438,1,H);
RHOCE_vec  =   normrnd(0.9144,0.0295,1,H);
RHOETA_vec  =   normrnd(0.1095,0.0771,1,H);
RHOETAM_vec  =   normrnd(0.9557,0.0164,1,H);
RHOETAX_vec  =   normrnd(0.8109,0.0668,1,H);
RHOGE_vec  =   normrnd(0.2983,0.1,1,H);
RHOIG_vec  =   normrnd(0.853,0.0797,1,H);
RHOLE_vec  =   normrnd(0.975,0.0159,1,H);
RHOL0_vec  =   normrnd(0.9334,0.0188,1,H);
RHOPCPM_vec  =   normrnd(0.6652,0.1409,1,H);
RHOPWPX_vec  =   normrnd(0.2159,0.0686,1,H);
RHORPE_vec  =   normrnd(0.9842,0.0103,1,H);
RHORPK_vec  =   normrnd(0.9148,0.0233,1,H);
RHOUCAP0_vec  =   normrnd(0.9517,0.0187,1,H);
RPREME_vec  =   normrnd(0.02,0.0074,1,H);
RPREMK_vec  =   normrnd(0.0245,0.0026,1,H);
SE_vec  =   normrnd(0.8588,0.0196,1,H);
SFPE_vec  =   normrnd(0.8714,0.0567,1,H);
SFPME_vec  =   normrnd(0.7361,0.1227,1,H);
SFPXE_vec  =   normrnd(0.918,0.0473,1,H);
SFWE_vec  =   normrnd(0.7736,0.1565,1,H);
SIGC_vec  =   normrnd(4.0962,0.813,1,H);
SIGEXE_vec  =   normrnd(2.5358,0.32,1,H);
SIGIME_vec  =   normrnd(1.1724,0.2136,1,H);
SLC_vec  =   normrnd(0.3507,0.0754,1,H);
RHOTR_vec  =   normrnd(0.8636,0.0428,1,H);
TYE1_vec  =   normrnd(0.4274,0.1141,1,H);
TYE2_vec  =   normrnd(0.0783,0.0277,1,H);
WRLAG_vec  =   normrnd(0.2653,0.1315,1,H);


% Draw from prior law
% G1E_vec  =   -betarnd(0.3871750148721,4.74777213197273,1,H);
% GSLAG_vec  =   -betarnd(9.09571408283434,12.4224171753496,1,H);
% GVECM_vec  =   -betarnd(10.4425602204091,56.1979006628652,1,H);
% IGVECM_vec  =   -betarnd(6.04568315131211,43.4279923913401,1,H);
% A2E_vec  =   betarnd(11.9122825378418,251.052011895752,1,H);
% GAMIE_vec  =   gamrnd(13.6871125370781,5.5553426476197,1,H);
% GAMI2E_vec  =   gamrnd(4.67927403523814,0.239695301355207,1,H);
% GAMLE_vec  =   gamrnd(22.5285830778624,2.58375326130466,1,H);
% GAMPE_vec  =   gamrnd(34.7633401512157,1.76742222504333,1,H);
% GAMPME_vec  =   gamrnd(3.40697211803127,0.492578143248719,1,H);
% GAMPXE_vec  =   gamrnd(2.40760377777554,10.8528655093496,1,H);
% GAMWE_vec  =   gamrnd(2.44564013067148,0.528246156823283,1,H);
% HABE_vec  =   betarnd(81.080357597323,62.8322401970025,1,H);
% HABLE_vec  =   betarnd(19.8492375075171,4.68931794744286,1,H);
% IG1E_vec  =   betarnd(1.77116614802584,10.0602710465355,1,H);
% IGSLAG_vec  =   betarnd(13.365,16.5009217877095,1,H);
% ILAGE_vec  =   betarnd(333.88212714256,36.7274045952133,1,H);
% KAPPAE_vec  =   gamrnd(18.7634597430808,0.102454452767374,1,H);
% RHOCE_vec  =   betarnd(81.3292104751508,7.61349564378053,1,H);
% RHOETA_vec  =   betarnd(1.68669468122152,13.7169097134955,1,H);
% RHOETAM_vec  =   betarnd(149.482871932629,6.92904805547293,1,H);
% RHOETAX_vec  =   betarnd(27.0550083713561,6.30916522755388,1,H);
% RHOGE_vec  =   betarnd(5.9456293913,13.9860816087,1,H);
% RHOIG_vec  =   betarnd(15.9853198285919,2.75479720375499,1,H);
% RHOLE_vec  =   betarnd(93.0308739765042,2.38540702503857,1,H);
% RHOL0_vec  =   betarnd(163.236723064735,11.6472742191037,1,H);
% RHOPCPM_vec  =   betarnd(6.79701820447584,3.42098871746619,1,H);
% RHOPWPX_vec  =   betarnd(7.55065652002992,27.4222778015538,1,H);
% RHORPE_vec  =   betarnd(143.276939711565,2.30011750400604,1,H);
% RHORPK_vec  =   betarnd(130.420074851259,12.1466882130819,1,H);
% RHOUCAP0_vec  =   betarnd(124.150243398439,6.30078465498012,1,H);
% RPREME_vec  =   betarnd(7.13850986121256,349.786983199416,1,H);
% RPREMK_vec  =   betarnd(86.5944164201183,3447.87155991124,1,H);
% SE_vec  =   betarnd(270.226910453978,44.4294827155352,1,H);
% SFPE_vec  =   betarnd(29.5031576539166,4.35403497164756,1,H);
% SFPME_vec  =   betarnd(8.76171924619187,3.14117335833451,1,H);
% SFPXE_vec  =   betarnd(29.9690857153074,2.67697715539782,1,H);
% SFWE_vec  =   betarnd(4.75838075897478,1.39257678882095,1,H);
% SIGC_vec  =   gamrnd(25.3852365844692,0.16136150578585,1,H);
% SIGEXE_vec  =   gamrnd(62.795719140625,0.0403817335752031,1,H);
% SIGIME_vec  =   gamrnd(30.126530741068,0.0389158648925282,1,H);
% SLC_vec  =   betarnd(13.6959979217823,25.3573180798781,1,H);
% RHOTR_vec  =   betarnd(54.6695232771421,8.63469543191545,1,H);
% TYE1_vec  =   betarnd(7.60691935607018,10.1912073544356,1,H);
% TYE2_vec  =   betarnd(7.28637478137341,85.7707744060264,1,H);
% WRLAG_vec  =   betarnd(2.72512582214576,7.54673931975307,1,H);



% Run the simulations

n=length(y_baseline); % Lenght of the response vector (i.e # of quarters)

% Empty response matrix
    % Results will be store as follow:
    % h-th column of res_VarName = response vector of the h-th simulation
    % i-th row of res_VarName = response at the i-th quarter 
res_y=nan(n,H);
res_phi=nan(n,H);
res_gg=nan(n,H);
res_fm=nan(n,H);
res_fm2=nan(n,H);
for h=1:H
    % We modify the parameters values with the i-th draw
    % Parameters instable (see appendix 1) -> to be updated
    %set_param_value('G1E',G1E_vec(h));
    %set_param_value('GSLAG',GSLAG_vec(h)); % Remove it for small H
    %set_param_value('GVECM',GVECM_vec(h));
    %set_param_value('IGVECM',IGVECM_vec(h));

    % Positive posterior mean
    set_param_value('A2E',A2E_vec(h));
    set_param_value('GAMIE',GAMIE_vec(h));
    set_param_value('GAMI2E',GAMI2E_vec (h));
    set_param_value('GAMLE',GAMLE_vec (h));
    set_param_value('GAMPE',GAMPE_vec (h));
    set_param_value('GAMPME',GAMPME_vec (h));
    set_param_value('GAMPXE',GAMPXE_vec (h));
    set_param_value('GAMWE',GAMWE_vec (h));
    set_param_value('HABE',HABE_vec (h));
    set_param_value('HABLE',HABLE_vec (h));
    set_param_value('IG1E',IG1E_vec (h));
    set_param_value('IGSLAG',IGSLAG_vec (h));
    set_param_value('ILAGE',ILAGE_vec (h));
    set_param_value('KAPPAE',KAPPAE_vec (h));
    set_param_value('RHOCE',RHOCE_vec (h));
    set_param_value('RHOETA',RHOETA_vec (h));
    set_param_value('RHOETAM',RHOETAM_vec (h));
    set_param_value('RHOETAX',RHOETAX_vec (h));
    set_param_value('RHOGE',RHOGE_vec (h));
    set_param_value('RHOIG',RHOIG_vec (h));
    set_param_value('RHOLE',RHOLE_vec (h));
    set_param_value('RHOL0',RHOL0_vec (h));
    set_param_value('RHOPCPM',RHOPCPM_vec (h));
    set_param_value('RHOPWPX',RHOPWPX_vec (h));
    set_param_value('RHORPE',RHORPE_vec (h));
    set_param_value('RHORPK',RHORPK_vec (h));
    set_param_value('RHOUCAP0',RHOUCAP0_vec (h));
    set_param_value('RPREME',RPREME_vec (h));
    set_param_value('RPREMK',RPREMK_vec (h));
    set_param_value('SE',SE_vec (h));
    set_param_value('SFPE',SFPE_vec (h));
    set_param_value('SFPME',SFPME_vec (h));
    set_param_value('SFPXE',SFPXE_vec (h));
    set_param_value('SFWE',SFWE_vec (h));
    set_param_value('SIGC',SIGC_vec (h));
    set_param_value('SIGEXE',SIGEXE_vec (h));
    set_param_value('SIGIME',SIGIME_vec (h));
    set_param_value('SLC',SLC_vec (h));
    set_param_value('RHOTR',RHOTR_vec (h));
    set_param_value('TYE1',TYE1_vec (h));
    set_param_value('TYE2',TYE2_vec (h));
    set_param_value('WRLAG',WRLAG_vec (h));
    
    % Run the simulation for the new values
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG'});

    % Store the response of the h-th simulation
    res_y(:,h)=cumsum(E_GY_E_EPS_G);
    res_phi(:,h)=E_PHI_E_EPS_G;
    res_gg(:,h)=cumsum(E_GG_E_EPS_G);
    res_fm(:,h)=res_y(:,h)./res_gg(:,h)/GSN;
    res_fm2(:,h)=(cumprod(1+E_GY_E_EPS_G)-1)./(cumprod(1+E_GG_E_EPS_G)-1)/GSN;
end

%% Compute mean and confidence intervals
    % Note:
    % mean(matrix,2) -> sum over the raws
    % So, i-th raw mean(matrix,2) = mean response value at the i-th quarter
    % Idem for quantiles(matrix,p,2)

% output growth responses
y_mean=mean(res_y,2);
y_CI_lower=quantile(res_y,0.05,2); 
y_CI_upper=quantile(res_y,0.95,2);

% inflation responses
phi_mean=mean(res_phi,2);
phi_CI_lower=quantile(res_phi,0.05,2); 
phi_CI_upper=quantile(res_phi,0.95,2);

% government consumption responses
gg_mean=mean(res_gg,2);
gg_CI_lower=quantile(res_gg,0.05,2); 
gg_CI_upper=quantile(res_gg,0.95,2);

% fiscal multipliers (1st formula)
fm_mean=mean(res_fm,2);
fm_CI_lower=quantile(res_fm,0.05,2); 
fm_CI_upper=quantile(res_fm,0.95,2);

% fiscal multipliers (2nd formula)
fm2_mean=mean(res_fm2,2);
fm2_CI_lower=quantile(res_fm2,0.05,2); 
fm2_CI_upper=quantile(res_fm2,0.95,2);


%% Plot output, inflation and government consumption responses


t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval
    % It should look like the paper IRF

figure
subplot(3,1,1) % Output growth
plot(t,y_mean,'LineWidth',2); hold on
plot(t,y_CI_lower,'LineWidth',2); hold on
plot(t,y_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('$$Y_t$$','interpreter','latex','FontSize',10);

subplot(3,1,2) % Inflation
plot(t,phi_mean,'LineWidth',2); hold on
plot(t,phi_CI_lower,'LineWidth',2); hold on
plot(t,phi_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('$$\pi_t$$','interpreter','latex','FontSize',10);

subplot(3,1,3) % Government consumption growth
plot(t,gg_mean,'LineWidth',2); hold on
plot(t,gg_CI_lower,'LineWidth',2); hold on
plot(t,gg_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('$$G_t$$','interpreter','latex','FontSize',10);


% Compare the mean of the simulations with the baseline response

figure
subplot(3,1,1) % Output growth (mean vs. baseline)
plot(t,y_baseline,'LineWidth',2); hold on
plot(t,y_mean,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('$$Y_t$$ (mean vs. baseline)','interpreter','latex','FontSize',10);

subplot(3,1,2) % Inflation (mean vs. baseline)
plot(t,phi_baseline,'LineWidth',2); hold on
plot(t,phi_mean,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('$$\pi_t$$ (mean vs. baseline)','interpreter','latex','FontSize',10);

subplot(3,1,3) % Government consumption growth (mean vs. baseline)
plot(t,phi_baseline,'LineWidth',2); hold on
plot(t,phi_mean,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('$$G_t$$ (mean vs. baseline)','interpreter','latex','FontSize',10);

%% Plot fiscal multiplier responses


t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval

figure

subplot(2,2,1) % fiscal multiplier 1st formula
plot(t,fm_mean,'LineWidth',2); hold on
plot(t,fm_CI_lower,'LineWidth',2); hold on
plot(t,fm_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('FM 1st formula','FontSize',10);

subplot(2,2,2) % fiscal multiplier 2nd formula
plot(t,fm2_mean,'LineWidth',2); hold on
plot(t,fm2_CI_lower,'LineWidth',2); hold on
plot(t,fm2_CI_upper,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('FM 2nd formuma','FontSize',10);

% Add comparison with baseline

subplot(2,2,3) % fiscal multiplier 1st formula (mean vs. baseline)
plot(t,fm_baseline,'LineWidth',2); hold on
plot(t,fm_mean,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('FM 1st formula (mean vs. baseline)','FontSize',10);

subplot(2,2,4) % fiscal multiplier 2nd formula (mean vs. baseline)
plot(t,fm2_baseline,'LineWidth',2); hold on
plot(t,fm2_mean,'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('FM 2nd formula (mean vs. baseline)','FontSize',10);

% Same plots with a shorter time-horizon

horizon=30;
t_small=t(1:horizon);

figure
subplot(2,2,1)
plot(t_small,fm_mean(1:horizon),'LineWidth',2); hold on
plot(t_small,fm_CI_lower(1:horizon),'LineWidth',2); hold on
plot(t_small,fm_CI_upper(1:horizon),'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('FM 1st formula','FontSize',10);

subplot(2,2,2)
plot(t_small,fm2_mean(1:horizon),'LineWidth',2); hold on
plot(t_small,fm2_CI_lower(1:horizon),'LineWidth',2); hold on
plot(t_small,fm2_CI_upper(1:horizon),'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('FM 2nd formula','FontSize',10);

subplot(2,2,3)
plot(t_small,fm_baseline(1:horizon),'LineWidth',2); hold on
plot(t_small,fm_mean(1:horizon),'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('FM 1st formula (mean vs. baseline)','FontSize',10);

subplot(2,2,4)
plot(t_small,fm2_baseline(1:horizon),'LineWidth',2); hold on
plot(t_small,fm2_mean(1:horizon),'LineWidth',2); hold on
xlabel('quarters','FontSize',8);
title('FM 2nd formula (mean vs. baseline)','FontSize',10);


% RIGHT WAY TO PLOT CONFIDENCE INTERVALS
% figure
% 
% subplot(2,2,1) % fiscal multiplier 1st formula
% plot(t,fm_mean,'LineWidth',2); hold on
% fill([t fliplr(t)], [fm_CI_lower' flipud(fm_CI_upper)'], 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.4); hold on
% xlabel('quarters','FontSize',8);
% title('FM 1st formula','FontSize',10);
% 
% subplot(2,2,2) % fiscal multiplier 2nd formula
% plot(t,fm2_mean,'LineWidth',2); hold on
% fill([t fliplr(t)], [fm2_CI_lower' flipud(fm2_CI_upper)'], 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.4); hold on
% xlabel('quarters','FontSize',8);
% title('FM 2nd formuma','FontSize',10);

% figure 
% subplot(2,1,1)
% f_plot_response(t,fm_mean,fm_CI_lower,fm_CI_upper,"FM 1st formula")
% subplot(2,1,2)
% f_plot_response(t,fm2_mean,fm2_CI_lower,fm2_CI_upper,"FM 2nd formula")
% 
