%-------------------------------------------------------------------------%
%-------------------- GOVERNMENT CONSUMPTION SHOCK -----------------------% 
%-------------------------------------------------------------------------%

%% Initialize
s = 0; % enter a seed (for random draws)
clc;
close all;

% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run the .mod file for the first time
    % Parameters value = posterior mean
    % We replicate what is done in the paper
dynare EA_Quest3_rep.mod

%% Reponse with uncertainty bounds
    % For each parameter, we draw a value
    % We draw from a normal distribution
    % Based on posterior mean and st.dev
    % See Excel File : randomize paramaters
    % We can simulate the .mod file for each vector of parameters value
    % We obtain several responses to shocks
    % We can a Bayesian confidence interval around the mean response

H=50; % Set the number of the simulations
rng(s); % Specify the seed

% Draws
    % For each parameter, draw H values from the distribution
    % See the Excel file for the distribution inputs
    % Store the H values in a vector called VarName_vec
        
% Draw vectors from normal law
% Instable parameters
% GSLAG_vec  =   normrnd(-0.4227,0.1041,1,H);
% SFPE_vec  =   normrnd(0.8714,0.0567,1,H);
% SFPXE_vec  =   normrnd(0.918,0.0473,1,H);
% SFWE_vec  =   normrnd(0.7736,0.1565,1,H);
% SIGEXE_vec  =   normrnd(2.5358,0.32,1,H);

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
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG'}); % Store output, inflation and gov. consumption responses

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
subplot(3,1,1) % Output
f_plot_with_CI(t,Y_mean,Y_CI_lower,Y_CI_upper,"Output response : $$(Y_t-Y_t^{SS})/Y_0^{SS}$$","Latex",1)

subplot(3,1,2) % Inflation
f_plot_with_CI(t,Infl_mean,Infl_CI_lower,Infl_CI_upper,"Inflation response : $$\pi_t- \pi_t^{SS}$$","Latex",1)

subplot(3,1,3) % Government consumption
f_plot_with_CI(t,G_mean,G_CI_lower,G_CI_upper,"Gov. consumption response : $$(G_t-G_t^{SS})/G_0^{SS}$$","Latex",1)
saveas(gcf,'Figure_1_IRF_cons','png')

%% Plot fiscal multiplier responses

horizon=30; % use shorter time horizon to plot fiscal multipliers
t_small=t(1:horizon);

figure ('Name','Gov. consumption fiscal multipliers')
subplot(2,1,1)
f_plot_with_CI(t_small,fm_mean(1:horizon),fm_CI_lower(1:horizon),fm_CI_upper(1:horizon),"Period gov. consumption multiplier : $$fm_t$$","Latex",0)

subplot(2,1,2)
f_plot_with_CI(t_small,FM_mean(1:horizon),FM_CI_lower(1:horizon),FM_CI_upper(1:horizon),"Cumulative gov. consumption multiplier : $$FM_t$$","Latex",0)
saveas(gcf,'Figure_2_multipliers_cons','png')

%% Fiscal multipliers in a table

% Short run
varname={'5%','Mean','95%'};
Q1 = table([fm_CI_lower(1);FM_CI_lower(1)],[fm_mean(1);FM_mean(1)],[fm_CI_upper(1);FM_CI_upper(1)], 'VariableNames',varname);
Q2 = table([fm_CI_lower(2);FM_CI_lower(2)],[fm_mean(2);FM_mean(2)],[fm_CI_upper(2);FM_CI_upper(2)],'VariableNames',varname);
Q3 = table([fm_CI_lower(3);FM_CI_lower(3)],[fm_mean(3);FM_mean(3)],[fm_CI_upper(3);FM_CI_upper(3)],'VariableNames',varname);
Q4 = table([fm_CI_lower(4);FM_CI_lower(4)],[fm_mean(4);FM_mean(4)],[fm_CI_upper(4);FM_CI_upper(4)],'VariableNames',varname);

Table_sr=table(Q1,Q2,Q3,Q4,'VariableNames',{'Quarter 1','Quarter 2','Quarter 3','Quarter 4'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

%long-run 
Q8 = table([fm_CI_lower(8);FM_CI_lower(8)],[fm_mean(8);FM_mean(8)],[fm_CI_upper(8);FM_CI_upper(8)], 'VariableNames',varname);
Q12 = table([fm_CI_lower(12);FM_CI_lower(12)],[fm_mean(12);FM_mean(12)],[fm_CI_upper(12);FM_CI_upper(12)],'VariableNames',varname);
Q16 = table([fm_CI_lower(16);FM_CI_lower(16)],[fm_mean(16);FM_mean(16)],[fm_CI_upper(16);FM_CI_upper(16)],'VariableNames',varname);
Q20 = table([fm_CI_lower(20);FM_CI_lower(20)],[fm_mean(20);FM_mean(20)],[fm_CI_upper(20);FM_CI_upper(20)],'VariableNames',varname);

Table_lr=table(Q8,Q12,Q16,Q20,'VariableNames',{'Year 2','Year 3','Year 4','Year 5'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

% Store the values in Tables.xlsx

writetable(Q1,'Tables_12_multipliers.xlsx','Sheet','SR','Range','C3')
writetable(Q2,'Tables_12_multipliers.xlsx','Sheet','SR','Range','F3')
writetable(Q3,'Tables_12_multipliers.xlsx','Sheet','SR','Range','I3')
writetable(Q4,'Tables_12_multipliers.xlsx','Sheet','SR','Range','L3')

writetable(Q8,'Tables_12_multipliers.xlsx','Sheet','LR','Range','C3')
writetable(Q12,'Tables_12_multipliers.xlsx','Sheet','LR','Range','F3')
writetable(Q16,'Tables_12_multipliers.xlsx','Sheet','LR','Range','I3')
writetable(Q20,'Tables_12_multipliers.xlsx','Sheet','LR','Range','L3')

% Go back to original path
cd('..');

% ------------------------------------------------------------------------%
%-------------------- GOVERNMENT CONSUMPTION SHOCK -----------------------% 
%--------- no labor wage rigidities and no labor adjustment costs --------%
%-------------------------------------------------------------------------%

%% Initialize
s = 0; % enter a seed (for random draws)


% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run the .mod file for the first time
    % Parameters value = posterior mean
    % We replicate what is done in the paper
dynare EA_Quest3_rep.mod

%% Reponse with uncertainty bounds
    % For each parameter, we draw a value
    % We draw from a normal distribution
    % Based on posterior mean and st.dev
    % See Excel File : randomize paramaters
    % We can simulate the .mod file for each vector of parameters value
    % We obtain several responses to shocks
    % We can a Bayesian confidence interval around the mean response

H=50; % Set the number of the simulations
rng(s); % Specify the seed

% Draws
    % For each parameter, draw H values from the distribution
    % See the Excel file for the distribution inputs
    % Store the H values in a vector called VarName_vec
        
% Draw vectors from normal law
% Instable parameters
% GSLAG_vec  =   normrnd(-0.4227,0.1041,1,H);
% SFPE_vec  =   normrnd(0.8714,0.0567,1,H);
% SFPXE_vec  =   normrnd(0.918,0.0473,1,H);
% SFWE_vec  =   normrnd(0.7736,0.1565,1,H);
% SIGEXE_vec  =   normrnd(2.5358,0.32,1,H);

% Other parameters
G1E_vec  =   normrnd(-0.0754,0.1066,1,H);
GVECM_vec  =   normrnd(-0.1567,0.0442,1,H);
IGVECM_vec  =   normrnd(-0.1222,0.0461,1,H);
A2E_vec  =   normrnd(0.0453,0.0128,1,H);
GAMIE_vec  =   normrnd(76.0366,20.5526,1,H);
GAMI2E_vec  =   normrnd(1.1216,0.5185,1,H);
%GAMLE_vec  =   normrnd(58.2083,12.2636,1,H); not drawn
GAMPE_vec  =   normrnd(61.4415,10.4208,1,H);
GAMPME_vec  =   normrnd(1.6782,0.9092,1,H);
GAMPXE_vec  =   normrnd(26.1294,16.8398,1,H);
%GAMWE_vec  =   normrnd(1.2919,0.8261,1,H); not drawn
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
    set_param_value('GAMLE',0); % No labor adjustment cost
    set_param_value('GAMPE',GAMPE_vec(h));
    set_param_value('GAMPME',GAMPME_vec(h));
    set_param_value('GAMPXE',GAMPXE_vec(h));
    set_param_value('GAMWE',0); % No wage rigidities
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
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG'}); % Store output, inflation and gov. consumption responses

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

figure ('Name','Response to a transitory gov. consumption shock (no nominal wage rigidities, no labor adjustment cost)')
subplot(3,1,1) % Output
f_plot_with_CI(t,Y_mean,Y_CI_lower,Y_CI_upper,"Output response : $$(Y_t-Y_t^{SS})/Y_0^{SS}$$","Latex",1)

subplot(3,1,2) % Inflation
f_plot_with_CI(t,Infl_mean,Infl_CI_lower,Infl_CI_upper,"Inflation response : $$\pi_t- \pi_t^{SS}$$","Latex",1)

subplot(3,1,3) % Government consumption
f_plot_with_CI(t,G_mean,G_CI_lower,G_CI_upper,"Gov. consumption response : $$(G_t-G_t^{SS})/G_0^{SS}$$","Latex",1)
saveas(gcf,'Figure_3_IRF_cons','png')


%% Plot fiscal multiplier responses

horizon=30; % use shorter time horizon to plot fiscal multipliers
t_small=t(1:horizon);

figure ('Name','Gov. consumption fiscal multipliers (no nominal wage rigidities, no labor adjustment cost)')
subplot(2,1,1)
f_plot_with_CI(t_small,fm_mean(1:horizon),fm_CI_lower(1:horizon),fm_CI_upper(1:horizon),"Period gov. consumption multiplier (noNWRnoLAC) : $$fm_t$$","Latex",0)

subplot(2,1,2)
f_plot_with_CI(t_small,FM_mean(1:horizon),FM_CI_lower(1:horizon),FM_CI_upper(1:horizon),"Cumulative gov. consumption multiplier (noNWRnoLAC) : $$FM_t$$","Latex",0)
saveas(gcf,'Figure_4_multipliers_cons_noNWRnoLAC','png')

%% Fiscal multipliers in a table

% Short run
varname={'5%','Mean','95%'};
Q1 = table([fm_CI_lower(1);FM_CI_lower(1)],[fm_mean(1);FM_mean(1)],[fm_CI_upper(1);FM_CI_upper(1)], 'VariableNames',varname);
Q2 = table([fm_CI_lower(2);FM_CI_lower(2)],[fm_mean(2);FM_mean(2)],[fm_CI_upper(2);FM_CI_upper(2)],'VariableNames',varname);
Q3 = table([fm_CI_lower(3);FM_CI_lower(3)],[fm_mean(3);FM_mean(3)],[fm_CI_upper(3);FM_CI_upper(3)],'VariableNames',varname);
Q4 = table([fm_CI_lower(4);FM_CI_lower(4)],[fm_mean(4);FM_mean(4)],[fm_CI_upper(4);FM_CI_upper(4)],'VariableNames',varname);

Table_sr=table(Q1,Q2,Q3,Q4,'VariableNames',{'Quarter 1','Quarter 2','Quarter 3','Quarter 4'},'RowNames',{'period fiscal multiplier (noNWRnoLAC) ','cumulative fiscal multiplier (noNWRnoLAC)'});

%long-run 
Q8 = table([fm_CI_lower(8);FM_CI_lower(8)],[fm_mean(8);FM_mean(8)],[fm_CI_upper(8);FM_CI_upper(8)], 'VariableNames',varname);
Q12 = table([fm_CI_lower(12);FM_CI_lower(12)],[fm_mean(12);FM_mean(12)],[fm_CI_upper(12);FM_CI_upper(12)],'VariableNames',varname);
Q16 = table([fm_CI_lower(16);FM_CI_lower(16)],[fm_mean(16);FM_mean(16)],[fm_CI_upper(16);FM_CI_upper(16)],'VariableNames',varname);
Q20 = table([fm_CI_lower(20);FM_CI_lower(20)],[fm_mean(20);FM_mean(20)],[fm_CI_upper(20);FM_CI_upper(20)],'VariableNames',varname);

Table_lr=table(Q8,Q12,Q16,Q20,'VariableNames',{'Year 2','Year 3','Year 4','Year 5'},'RowNames',{'period fiscal multiplier (noNWRnoLAC)','cumulative fiscal multiplier (noNWRnoLAC)'});

% Store the values in Tables.xlsx

writetable(Q1,'Tables_34_multipliers_noNWRnoLAC.xlsx','Sheet','SR','Range','C3')
writetable(Q2,'Tables_34_multipliers_noNWRnoLAC.xlsx','Sheet','SR','Range','F3')
writetable(Q3,'Tables_34_multipliers_noNWRnoLAC.xlsx','Sheet','SR','Range','I3')
writetable(Q4,'Tables_34_multipliers_noNWRnoLAC.xlsx','Sheet','SR','Range','L3')

writetable(Q8,'Tables_34_multipliers_noNWRnoLAC.xlsx','Sheet','LR','Range','C3')
writetable(Q12,'Tables_34_multipliers_noNWRnoLAC.xlsx','Sheet','LR','Range','F3')
writetable(Q16,'Tables_34_multipliers_noNWRnoLAC.xlsx','Sheet','LR','Range','I3')
writetable(Q20,'Tables_34_multipliers_noNWRnoLAC.xlsx','Sheet','LR','Range','L3')

% Go back to original path
cd('..');

%-------------------------------------------------------------------------%
%---------------------- GOVERNMENT INVESTMENT SHOCK ----------------------% 
%-------------------------------------------------------------------------%

%% Initialize
s = 0; % enter a seed (for random draws)


% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run the .mod file for the first time
    % Parameters value = posterior mean
    % We replicate what is done in the paper
dynare EA_Quest3_rep.mod

%% Reponse with uncertainty bounds
    % For each parameter, we draw a value
    % We draw from a normal distribution
    % Based on posterior mean and st.dev
    % See Excel File : randomize paramaters
    % We can simulate the .mod file for each vector of parameters value
    % We obtain several responses to shocks
    % We can a Bayesian confidence interval around the mean response

H=50; % Set the number of the simulations
rng(s); % Specify the seed

% Draws
    % For each parameter, draw H values from the distribution
    % See the Excel file for the distribution inputs
    % Store the H values in a vector called VarName_vec
        
% Draw vectors from normal law
% Instable parameters
% GSLAG_vec  =   normrnd(-0.4227,0.1041,1,H);
% SFPE_vec  =   normrnd(0.8714,0.0567,1,H);
% SFPXE_vec  =   normrnd(0.918,0.0473,1,H);
% SFWE_vec  =   normrnd(0.7736,0.1565,1,H);
% SIGEXE_vec  =   normrnd(2.5358,0.32,1,H);

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

n=length(E_GY_E_EPS_IG); % Lenght of the response vector (i.e # of quarters)
GY0_vec = ones(1, n)' * GY0; % to be used in the fiscal multiplier formulas 

% Empty results matrix
    % Results will be store as follow:
    % h-th column of res_VarName = response vector of the h-th simulation
    % i-th row of res_VarName = response at the i-th quarter 
res_Y=nan(n,H);
res_Infl=nan(n,H);
res_IG=nan(n,H);
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
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GIG'}); % Store output, inflation and gov. investment responses

    % Store the response of the h-th simulation
    res_Y(:,h)=cumprod(1+E_GY_E_EPS_IG+GY0)-cumprod(1+GY0_vec); % Response Output
    res_Infl(:,h)=E_PHI_E_EPS_IG; % Response Inflation
    res_IG(:,h)=cumprod(1+E_GIG_E_EPS_IG+GY0)-cumprod(1+GY0_vec); % Response Gov. investment
    res_fm(:,h)=res_Y(:,h)./res_IG(:,h)/IGSN; % Response period multiplier
    res_FM(:,h)=cumsum(res_Y(:,h))./cumsum(res_IG(:,h))/IGSN; % Response cumulative multiplier
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
IG_mean=mean(res_IG,2);
IG_CI_lower=quantile(res_IG,0.05,2); 
IG_CI_upper=quantile(res_IG,0.95,2);

% Period multiplier (simulation output)
fm_mean=mean(res_fm,2);
fm_CI_lower=quantile(res_fm,0.05,2); 
fm_CI_upper=quantile(res_fm,0.95,2);

% Cumulative multiplier (simulation output)
FM_mean=mean(res_FM,2);
FM_CI_lower=quantile(res_FM,0.05,2); 
FM_CI_upper=quantile(res_FM,0.95,2);


%% Plot output, inflation


t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval
    % It should look like the paper IRF

figure ('Name','Response to a transitory gov. investment shock')
subplot(3,1,1) % Output
f_plot_with_CI(t,Y_mean,Y_CI_lower,Y_CI_upper,"Output response : $$(Y_t-Y_t^{SS})/Y_0^{SS}$$","Latex",1)

subplot(3,1,2) % Inflation
f_plot_with_CI(t,Infl_mean,Infl_CI_lower,Infl_CI_upper,"Inflation response : $$\pi_t- \pi_t^{SS}$$","Latex",1)

subplot(3,1,3) % Government investment
f_plot_with_CI(t,IG_mean,IG_CI_lower,IG_CI_upper,"Gov. investment response : $$(IG_t-IG_t^{SS})/IG_0^{SS}$$","Latex",1)
saveas(gcf,'Figure_5_IRF_cons','png')


%% Plot fiscal multiplier responses

horizon=30; % use shorter time horizon to plot fiscal multipliers
t_small=t(1:horizon);

figure ('Name','Gov. investment fiscal multipliers')
subplot(2,1,1)
f_plot_with_CI(t_small,fm_mean(1:horizon),fm_CI_lower(1:horizon),fm_CI_upper(1:horizon),"Period gov. investment multiplier : $$fm_t$$","Latex",0)

subplot(2,1,2)
f_plot_with_CI(t_small,FM_mean(1:horizon),FM_CI_lower(1:horizon),FM_CI_upper(1:horizon),"Cumulative gov. investment multiplier : $$FM_t$$","Latex",0)
saveas(gcf,'Figure_6_multipliers_inv','png')

%% Fiscal multipliers in a table

% Short run
varname={'5%','Mean','95%'};
Q1 = table([fm_CI_lower(1);FM_CI_lower(1)],[fm_mean(1);FM_mean(1)],[fm_CI_upper(1);FM_CI_upper(1)], 'VariableNames',varname);
Q2 = table([fm_CI_lower(2);FM_CI_lower(2)],[fm_mean(2);FM_mean(2)],[fm_CI_upper(2);FM_CI_upper(2)],'VariableNames',varname);
Q3 = table([fm_CI_lower(3);FM_CI_lower(3)],[fm_mean(3);FM_mean(3)],[fm_CI_upper(3);FM_CI_upper(3)],'VariableNames',varname);
Q4 = table([fm_CI_lower(4);FM_CI_lower(4)],[fm_mean(4);FM_mean(4)],[fm_CI_upper(4);FM_CI_upper(4)],'VariableNames',varname);

Table_sr=table(Q1,Q2,Q3,Q4,'VariableNames',{'Quarter 1','Quarter 2','Quarter 3','Quarter 4'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

%long-run 
Q8 = table([fm_CI_lower(8);FM_CI_lower(8)],[fm_mean(8);FM_mean(8)],[fm_CI_upper(8);FM_CI_upper(8)], 'VariableNames',varname);
Q12 = table([fm_CI_lower(12);FM_CI_lower(12)],[fm_mean(12);FM_mean(12)],[fm_CI_upper(12);FM_CI_upper(12)],'VariableNames',varname);
Q16 = table([fm_CI_lower(16);FM_CI_lower(16)],[fm_mean(16);FM_mean(16)],[fm_CI_upper(16);FM_CI_upper(16)],'VariableNames',varname);
Q20 = table([fm_CI_lower(20);FM_CI_lower(20)],[fm_mean(20);FM_mean(20)],[fm_CI_upper(20);FM_CI_upper(20)],'VariableNames',varname);

Table_lr=table(Q8,Q12,Q16,Q20,'VariableNames',{'Year 2','Year 3','Year 4','Year 5'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

% Store the values in Tables.xlsx

writetable(Q1,'Tables_56_multipliers_investment.xlsx','Sheet','SR','Range','C3')
writetable(Q2,'Tables_56_multipliers_investment.xlsx','Sheet','SR','Range','F3')
writetable(Q3,'Tables_56_multipliers_investment.xlsx','Sheet','SR','Range','I3')
writetable(Q4,'Tables_56_multipliers_investment.xlsx','Sheet','SR','Range','L3')

writetable(Q8,'Tables_56_multipliers_investment.xlsx','Sheet','LR','Range','C3')
writetable(Q12,'Tables_56_multipliers_investment.xlsx','Sheet','LR','Range','F3')
writetable(Q16,'Tables_56_multipliers_investment.xlsx','Sheet','LR','Range','I3')
writetable(Q20,'Tables_56_multipliers_investment.xlsx','Sheet','LR','Range','L3')

% Go back to original path
cd('..');

%-------------------------------------------------------------------------%
%---------(------------- GOVERNMENT TRANSFER SHOCK -----------------------% 
%-------------------------------------------------------------------------%

%% Initialize
s = 0; % enter a seed (for random draws)


% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run the .mod file for the first time
    % Parameters value = posterior mean
    % We replicate what is done in the paper
dynare EA_Quest3_rep.mod

%% Reponse with uncertainty bounds
    % For each parameter, we draw a value
    % We draw from a normal distribution
    % Based on posterior mean and st.dev
    % See Excel File : randomize paramaters
    % We can simulate the .mod file for each vector of parameters value
    % We obtain several responses to shocks
    % We can a Bayesian confidence interval around the mean response

H=50; % Set the number of the simulations
rng(s); % Specify the seed

% Draws
    % For each parameter, draw H values from the distribution
    % See the Excel file for the distribution inputs
    % Store the H values in a vector called VarName_vec
        
% Draw vectors from normal law
% Instable parameters
% GSLAG_vec  =   normrnd(-0.4227,0.1041,1,H);
% SFPE_vec  =   normrnd(0.8714,0.0567,1,H);
% SFPXE_vec  =   normrnd(0.918,0.0473,1,H);
% SFWE_vec  =   normrnd(0.7736,0.1565,1,H);
% SIGEXE_vec  =   normrnd(2.5358,0.32,1,H);

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

n=length(E_GY_E_EPS_IG); % Lenght of the response vector (i.e # of quarters)
GY0_vec = ones(1, n)' * GY0; % to be used in the fiscal multiplier formulas 

% Empty results matrix
    % Results will be store as follow:
    % h-th column of res_VarName = response vector of the h-th simulation
    % i-th row of res_VarName = response at the i-th quarter 
res_Y=nan(n,H);
res_Infl=nan(n,H);
res_TR=nan(n,H);
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
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GTR'}); % Store output, inflation and gov. transfer responses

    % Store the response of the h-th simulation
    res_Y(:,h)=cumprod(1+E_GY_E_EPS_TR+GY0)-cumprod(1+GY0_vec); % Response Output
    res_Infl(:,h)=E_PHI_E_EPS_TR; % Response Inflation
    res_TR(:,h)=cumprod(1+E_GTR_E_EPS_TR+GY0)-cumprod(1+GY0_vec); % Response Gov. investment
    res_fm(:,h)=res_Y(:,h)./res_TR(:,h)/0.1685; % Response period multiplier
    res_FM(:,h)=cumsum(res_Y(:,h))./cumsum(res_TR(:,h))/0.1685; % Response cumulative multiplier
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
TR_mean=mean(res_TR,2);
TR_CI_lower=quantile(res_TR,0.05,2); 
TR_CI_upper=quantile(res_TR,0.95,2);

% Period multiplier (simulation output)
fm_mean=mean(res_fm,2);
fm_CI_lower=quantile(res_fm,0.05,2); 
fm_CI_upper=quantile(res_fm,0.95,2);

% Cumulative multiplier (simulation output)
FM_mean=mean(res_FM,2);
FM_CI_lower=quantile(res_FM,0.05,2); 
FM_CI_upper=quantile(res_FM,0.95,2);


%% Plot output, inflation


t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval
    % It should look like the paper IRF

figure ('Name','Response to a transitory gov. transfer shock')
subplot(3,1,1) % Output
f_plot_with_CI(t,Y_mean,Y_CI_lower,Y_CI_upper,"Output response : $$(Y_t-Y_t^{SS})/Y_0^{SS}$$","Latex",1)

subplot(3,1,2) % Inflation
f_plot_with_CI(t,Infl_mean,Infl_CI_lower,Infl_CI_upper,"Inflation response : $$\pi_t- \pi_t^{SS}$$","Latex",1)

subplot(3,1,3) % Government investment
f_plot_with_CI(t,TR_mean,TR_CI_lower,TR_CI_upper,"Gov. transfer response : $$(IG_t-IG_t^{SS})/IG_0^{SS}$$","Latex",1)
saveas(gcf,'Figure_5_IRF_cons','png')


%% Plot fiscal multiplier responses

horizon=30; % use shorter time horizon to plot fiscal multipliers
t_small=t(1:horizon);

figure ('Name','Gov. investment fiscal multipliers')
subplot(2,1,1)
f_plot_with_CI(t_small,fm_mean(1:horizon),fm_CI_lower(1:horizon),fm_CI_upper(1:horizon),"Period gov. investment multiplier : $$fm_t$$","Latex",0)

subplot(2,1,2)
f_plot_with_CI(t_small,FM_mean(1:horizon),FM_CI_lower(1:horizon),FM_CI_upper(1:horizon),"Cumulative gov. investment multiplier : $$FM_t$$","Latex",0)
saveas(gcf,'Figure_6_multipliers_inv','png')

%% Fiscal multipliers in a table

% Short run
varname={'5%','Mean','95%'};
Q1 = table([fm_CI_lower(1);FM_CI_lower(1)],[fm_mean(1);FM_mean(1)],[fm_CI_upper(1);FM_CI_upper(1)], 'VariableNames',varname);
Q2 = table([fm_CI_lower(2);FM_CI_lower(2)],[fm_mean(2);FM_mean(2)],[fm_CI_upper(2);FM_CI_upper(2)],'VariableNames',varname);
Q3 = table([fm_CI_lower(3);FM_CI_lower(3)],[fm_mean(3);FM_mean(3)],[fm_CI_upper(3);FM_CI_upper(3)],'VariableNames',varname);
Q4 = table([fm_CI_lower(4);FM_CI_lower(4)],[fm_mean(4);FM_mean(4)],[fm_CI_upper(4);FM_CI_upper(4)],'VariableNames',varname);

Table_sr=table(Q1,Q2,Q3,Q4,'VariableNames',{'Quarter 1','Quarter 2','Quarter 3','Quarter 4'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

%long-run 
Q8 = table([fm_CI_lower(8);FM_CI_lower(8)],[fm_mean(8);FM_mean(8)],[fm_CI_upper(8);FM_CI_upper(8)], 'VariableNames',varname);
Q12 = table([fm_CI_lower(12);FM_CI_lower(12)],[fm_mean(12);FM_mean(12)],[fm_CI_upper(12);FM_CI_upper(12)],'VariableNames',varname);
Q16 = table([fm_CI_lower(16);FM_CI_lower(16)],[fm_mean(16);FM_mean(16)],[fm_CI_upper(16);FM_CI_upper(16)],'VariableNames',varname);
Q20 = table([fm_CI_lower(20);FM_CI_lower(20)],[fm_mean(20);FM_mean(20)],[fm_CI_upper(20);FM_CI_upper(20)],'VariableNames',varname);

Table_lr=table(Q8,Q12,Q16,Q20,'VariableNames',{'Year 2','Year 3','Year 4','Year 5'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

% Store the values in Tables.xlsx

writetable(Q1,'Tables_78_multipliers_transfer.xlsx','Sheet','SR','Range','C3')
writetable(Q2,'Tables_78_multipliers_transfer.xlsx','Sheet','SR','Range','F3')
writetable(Q3,'Tables_78_multipliers_transfer.xlsx','Sheet','SR','Range','I3')
writetable(Q4,'Tables_78_multipliers_transfer.xlsx','Sheet','SR','Range','L3')

writetable(Q8,'Tables_78_multipliers_transfer.xlsx','Sheet','LR','Range','C3')
writetable(Q12,'Tables_78_multipliers_transfer.xlsx','Sheet','LR','Range','F3')
writetable(Q16,'Tables_78_multipliers_transfer.xlsx','Sheet','LR','Range','I3')
writetable(Q20,'Tables_78_multipliers_transfer.xlsx','Sheet','LR','Range','L3')

% Go back to original path
cd('..');