% Run and plot IRF of replication file 


% EA_QUEST3 model

s = 10; % enter a seed (for random draws)
clc;
close all;


% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);


%% Replication of the paper results

% Run the .mod file for the first time
    % Parameters value = posterior mean
    % We replicate what is done in the paper
dynare EA_Quest3_rep.mod

% % Store reponses to a government consumption shock (E_EPS_G)
%     % We call the response baseline
%     % Therefore baseline responses = paper responses
%    
% GY0_vec = ones(1, length(E_GY_E_EPS_G))' * GY0;
% 
% Y_baseline=cumprod(1+E_GY_E_EPS_G+GY0)-cumprod(1+GY0_vec);% Response Output
% Infl_baseline = E_PHI_E_EPS_G; % Response Inflation
% G_baseline =cumprod(1+E_GG_E_EPS_G+GY0)-cumprod(1+GY0_vec); % Response Gov. consumption
% 
% % Store the fiscal multipliers
%     % We should recover the 0.73 for Q1 and 0.45 for Q4 (period multiplier)
%     % We should recover 0.58 for Q4 (cumulative multiplier)
% 
% fm_baseline=Y_baseline./G_baseline/GSN; % period multiplier
% FM_baseline=cumsum(Y_baseline)./cumsum(G_baseline)/GSN; % cumulative multiplier
%     % Check
%     fm_baseline(1);
%     fm_baseline(4);
%     FM_baseline(1);
%     FM_baseline(4);

% % Plot the fiscal multipliers
% 
% t = 1:1:length(Y_baseline); % x-axis
% 
% figure('Name','Period and cumulative government consumption multipliers')
% subplot(2,1,1) % Period multiplier
% plot(t,fm_baseline,'LineWidth',2);
% xlabel('quarters','FontSize',8);
% title('Government consumption period multiplier','FontSize',10);
% subplot(2,1,2) % 2nd formula
% plot(t,FM_baseline,'LineWidth',2);
% xlabel('quarters','FontSize',8);
% title('Government consumption cumulative multiplier','FontSize',10);

%% Reponse with uncertainty bounds
    % For each parameter, we draw a value
    % We draw from a normal distribution
    % Based on posterior mean and st.dev
    % See Excel File : randomize paramaters
    % We can simulate the .mod file for each vector of parameters value
    % We obtain several responses to shocks
    % We can a Bayesian confidence interval around the mean response

H=400; % Set the number of the simulations
rng(s); % Specify the seed

% Draws
    % For each parameter, draw H values from the distribution
    % See the Excel file for the distribution inputs
    % Store the H values in a vector called VarName_vec
        
% Draw vectors from normal law
% Instable parameters
GSLAG_vec  =   normrnd(-0.4227,0.1041,1,H);
SFPE_vec  =   normrnd(0.8714,0.0567,1,H);
SFPXE_vec  =   normrnd(0.918,0.0473,1,H);
SFWE_vec  =   normrnd(0.7736,0.1565,1,H);
SIGEXE_vec  =   normrnd(2.5358,0.32,1,H);

% Other parameters
G1E_vec  =   normrnd(-0.0754,0.1066,1,H);
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
SFPME_vec  =   normrnd(0.7361,0.1227,1,H);
SIGC_vec  =   normrnd(4.0962,0.813,1,H);
SIGIME_vec  =   normrnd(1.1724,0.2136,1,H);
SLC_vec  =   normrnd(0.3507,0.0754,1,H);
RHOTR_vec  =   normrnd(0.8636,0.0428,1,H);
TYE1_vec  =   normrnd(0.4274,0.1141,1,H);
TYE2_vec  =   normrnd(0.0783,0.0277,1,H);
WRLAG_vec  =   normrnd(0.2653,0.1315,1,H);

% Run the simulations

n=length(E_GY_E_EPS_G); % Lenght of the response vector (i.e # of quarters)
GY0_vec = ones(1, n)' * GY0; % to be used in the fiscal multiplier formulas 

% Empty results matrix
    % Results will be store as follow:
    % h-th column of res_VarName = response vector of the h-th simulation
    % i-th row of res_VarName = response at the i-th quarter 
res_Y=nan(n,H);
res_Infl=nan(n,H);
res_G=nan(n,H);
res_fm=nan(n,H);
res_FM=nan(n,H);
for h=1:H
    % We modify the parameters values with the i-th draw
    % Instable parameters (see appendix 2) (remove for small H)
%     set_param_value('GSLAG',GSLAG_vec(h));
%     set_param_value('SFPE',SFPE_vec(h));
%     set_param_value('SFPXE',SFPXE_vec(h));
%     set_param_value('SFWE',SFWE_vec(h));
%     set_param_value('SIGEXE',SIGEXE_vec(h));
  
    % Other parameters
    set_param_value('G1E',G1E_vec(h));
    set_param_value('GVECM',GVECM_vec(h));
    set_param_value('IGVECM',IGVECM_vec(h));
    set_param_value('A2E',A2E_vec(h));
    set_param_value('GAMIE',GAMIE_vec(h));
    set_param_value('GAMI2E',GAMI2E_vec(h));
    set_param_value('GAMLE',GAMLE_vec(h));
    set_param_value('GAMPE',GAMPE_vec(h));
    set_param_value('GAMPME',GAMPME_vec(h));
    set_param_value('GAMPXE',GAMPXE_vec(h));
    set_param_value('GAMWE',GAMWE_vec(h));
    set_param_value('HABE',HABE_vec(h));
    set_param_value('HABLE',HABLE_vec(h));
    set_param_value('IG1E',IG1E_vec(h));
    set_param_value('IGSLAG',IGSLAG_vec(h));
    set_param_value('ILAGE',ILAGE_vec(h));
    set_param_value('KAPPAE',KAPPAE_vec(h));
    set_param_value('RHOCE',RHOCE_vec(h));
    set_param_value('RHOETA',RHOETA_vec(h));
    set_param_value('RHOETAM',RHOETAM_vec(h));
    set_param_value('RHOETAX',RHOETAX_vec(h));
    set_param_value('RHOGE',RHOGE_vec(h));
    set_param_value('RHOIG',RHOIG_vec(h));
    set_param_value('RHOLE',RHOLE_vec(h));
    set_param_value('RHOL0',RHOL0_vec(h));
    set_param_value('RHOPCPM',RHOPCPM_vec(h));
    set_param_value('RHOPWPX',RHOPWPX_vec(h));
    set_param_value('RHORPE',RHORPE_vec(h));
    set_param_value('RHORPK',RHORPK_vec(h));
    set_param_value('RHOUCAP0',RHOUCAP0_vec(h));
    set_param_value('RPREME',RPREME_vec(h));
    set_param_value('RPREMK',RPREMK_vec(h));
    set_param_value('SE',SE_vec(h));
    set_param_value('SFPME',SFPME_vec(h));
    set_param_value('SIGC',SIGC_vec(h));
    set_param_value('SIGIME',SIGIME_vec(h));
    set_param_value('SLC',SLC_vec(h));
    set_param_value('RHOTR',RHOTR_vec(h));
    set_param_value('TYE1',TYE1_vec(h));
    set_param_value('TYE2',TYE2_vec(h));
    set_param_value('WRLAG',WRLAG_vec(h));

    
    % Run the simulation for the new values
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG'});

    % Store the response of the h-th simulation
    res_Y(:,h)=cumprod(1+E_GY_E_EPS_G+GY0)-cumprod(1+GY0_vec); % Response Output
    res_Infl(:,h)=E_PHI_E_EPS_G; % Response Inflation
    res_G(:,h)=cumprod(1+E_GG_E_EPS_G+GY0)-cumprod(1+GY0_vec); % Response Gov. consumption
    res_fm(:,h)=res_Y(:,h)./res_G(:,h)/GSN; % Response period multiplier
    res_FM(:,h)=cumsum(res_Y(:,h))./cumsum(res_G(:,h))/GSN; % Response cumulative multiplier
end

%% Compute mean and confidence intervals based on the H simulations
    % Note:
    % mean(matrix,2) -> sum over the raws
    % So, i-th raw mean(matrix,2) = mean response value at the i-th quarter
    % Idem for quantiles(matrix,p,2)

% Output response (simulation output)
Y_mean=mean(res_Y,2);
Y_CI_lower=quantile(res_Y,0.05,2); 
Y_CI_upper=quantile(res_Y,0.95,2);

% Inflation response (simulation output)
Infl_mean=mean(res_Infl,2);
Infl_CI_lower=quantile(res_Infl,0.05,2); 
Infl_CI_upper=quantile(res_Infl,0.95,2);

% government consumption responses (simulation output)
G_mean=mean(res_G,2);
G_CI_lower=quantile(res_G,0.05,2); 
G_CI_upper=quantile(res_G,0.95,2);

% Period multiplier (simulation output)
fm_mean=mean(res_fm,2);
fm_CI_lower=quantile(res_fm,0.05,2); 
fm_CI_upper=quantile(res_fm,0.95,2);

% Cumulative multiplier (simulation output)
FM_mean=mean(res_FM,2);
FM_CI_lower=quantile(res_FM,0.05,2); 
FM_CI_upper=quantile(res_FM,0.95,2);


%% Plot output, inflation and government consumption responses


t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval
    % It should look like the paper IRF

figure ('Name','Response to a transitory gov. consumption shock')
subplot(3,1,1) % Output growth
f_plot_with_CI(t,Y_mean,Y_CI_lower,Y_CI_upper,"$$(Y_t-Y_t^{SS})/Y_t^{SS}$$","Latex")

subplot(3,1,2) % Inflation
f_plot_with_CI(t,Infl_mean,Infl_CI_lower,Infl_CI_upper,"$$\pi_t- \pi_t^{SS}$$","Latex")

subplot(3,1,3) % Government consumption growth
f_plot_with_CI(t,G_mean,G_CI_lower,G_CI_upper,"$$(G_t-G_t^{SS})/G_t^{SS}$$","Latex")


%% Plot fiscal multiplier responses


t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval

figure ('Name','Fiscal multipliers to a transitory gov. consumption shock')

subplot(2,1,1) % period multiplier
f_plot_with_CI(t,fm_mean,fm_CI_lower,fm_CI_upper,"Period gov. consumption multiplier","none")


subplot(2,1,2) % cumulative multiplier
f_plot_with_CI(t,FM_mean,FM_CI_lower,FM_CI_upper,"Cumulative gov. consumption multiplier","none")

% Same plots with a shorter time-horizon

horizon=30;
t_small=t(1:horizon);

figure
subplot(2,1,1)
f_plot_with_CI(t_small,fm_mean(1:horizon),fm_CI_lower(1:horizon),fm_CI_upper(1:horizon),"Period gov. consumption multiplier","none")

subplot(2,1,2)
f_plot_with_CI(t_small,FM_mean(1:horizon),FM_CI_lower(1:horizon),FM_CI_upper(1:horizon),"Cumulative gov. consumption multiplier","none")
