%-------------------------------------------------------------------------%
%---------------------------- GENERAL CASE -------------------------------% 
%-------------------------------------------------------------------------%

%% Initialize
s = 50; % enter a seed (for random draws)
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

H=500; % Set the number of the simulations
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
% Responses to government consumption shock
res_Y_G=nan(n,H);
res_Infl_G=nan(n,H);
res_G_G=nan(n,H);
res_fm_G=nan(n,H);
res_FM_G=nan(n,H);

% Responses to government investment shock
res_Y_GI=nan(n,H);
res_Infl_GI=nan(n,H);
res_IG_GI=nan(n,H);
res_fm_GI=nan(n,H);
res_FM_GI=nan(n,H);

% Responses to government transfer shock
res_Y_TR=nan(n,H);
res_Infl_TR=nan(n,H);
res_TR_TR=nan(n,H);
res_fm_TR=nan(n,H);
res_FM_TR=nan(n,H);

for h=1:H
    % We modify the parameters values with the i-th draw
    % Instable parameters (see appendix 2) (remove for small H)
    set_param_value('GSLAG',GSLAG_vec(h));
    set_param_value('SFPE',SFPE_vec(h));
    set_param_value('SFPXE',SFPXE_vec(h));
    set_param_value('SFWE',SFWE_vec(h));
    set_param_value('SIGEXE',SIGEXE_vec(h));
  
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
    stoch_simul(M_,options_,oo_,{'E_GY','E_PHI','E_GG','E_GIG','E_GTR'}); % Store output, inflation and government variable responses

    % Store the response of the h-th simulation
    % IRF to gov. consumption shock
    res_Y_G(:,h)=cumprod(1+E_GY_E_EPS_G+GY0)-cumprod(1+GY0_vec); % Response Output
    res_Infl_G(:,h)=E_PHI_E_EPS_G; % Response Inflation
    res_G_G(:,h)=cumprod(1+E_GG_E_EPS_G+GY0)-cumprod(1+GY0_vec); % Response Gov. consumption
    res_fm_G(:,h)=res_Y_G(:,h)./res_G_G(:,h)/GSN; % Response period multiplier
    res_FM_G(:,h)=cumsum(res_Y_G(:,h))./cumsum(res_G_G(:,h))/GSN; % Response cumulative multiplier
    
    % IRF to gov. consumption shock
    res_Y_GI(:,h)=cumprod(1+E_GY_E_EPS_IG+GY0)-cumprod(1+GY0_vec); % Response Output
    res_Infl_GI(:,h)=E_PHI_E_EPS_IG; % Response Inflation
    res_IG_GI(:,h)=cumprod(1+E_GIG_E_EPS_IG+GY0)-cumprod(1+GY0_vec); % Response Gov. investment
    res_fm_GI(:,h)=res_Y_GI(:,h)./res_IG_GI(:,h)/IGSN; % Response period multiplier
    res_FM_GI(:,h)=cumsum(res_Y_GI(:,h))./cumsum(res_IG_GI(:,h))/IGSN; % Response cumulative multiplier


    % IRF to gov. transfer shock
    res_Y_TR(:,h)=cumprod(1+E_GY_E_EPS_TR+GY0)-cumprod(1+GY0_vec); % Response Output
    res_Infl_TR(:,h)=E_PHI_E_EPS_TR; % Response Inflation
    %res_TR_TR(:,h)=cumprod(1+E_GTR_E_EPS_TR)-cumprod(1); % Response Gov. transfer
    res_TR_TR(:,h)=cumsum(E_GTR_E_EPS_TR); % Response Gov. transfer
    res_fm_TR(:,h)=res_Y_TR(:,h)./res_TR_TR(:,h)/0.1685; % Response period multiplier
    res_FM_TR(:,h)=cumsum(res_Y_TR(:,h))./cumsum(res_TR_TR(:,h))/0.1685; % Response cumulative multiplier
end

%% Compute mean and confidence intervals based on the H simulations
    % Note:
    % mean(matrix,2) -> sum over the raws
    % So, i-th raw mean(matrix,2) = mean response value at the i-th quarter
    % Idem for quantiles(matrix,p,2)

% Get confidence interval --> IRF to governement consumption shock
[Y_G_mean,Y_G_CI_lower,Y_G_CI_upper]=get_CI(res_Y_G); % Output IRFs
[Infl_G_mean,Infl_G_CI_lower,Infl_G_CI_upper]=get_CI(res_Infl_G); % Inflation IRF
[G_G_mean,G_G_CI_lower,G_G_CI_upper]=get_CI(res_G_G); % Gov. consumption IRF
[fm_G_mean,fm_G_CI_lower,fm_G_CI_upper]=get_CI(res_fm_G); % Period multiplier
[FM_G_mean,FM_G_CI_lower,FM_G_CI_upper]=get_CI(res_FM_G); % Cumulative multiplier

% Get confidence interval --> IRF to governement investment shock
[Y_GI_mean,Y_GI_CI_lower,Y_GI_CI_upper]=get_CI(res_Y_GI); % Output IRFs
[Infl_GI_mean,Infl_GI_CI_lower,Infl_GI_CI_upper]=get_CI(res_Infl_GI); % Inflation IRF
[IG_GI_mean,IG_GI_CI_lower,IG_GI_CI_upper]=get_CI(res_IG_GI); % Gov. investment IRF
[fm_GI_mean,fm_GI_CI_lower,fm_GI_CI_upper]=get_CI(res_fm_GI); % Period multiplier
[FM_GI_mean,FM_GI_CI_lower,FM_GI_CI_upper]=get_CI(res_FM_GI); % Cumulative multiplier

% Get confidence interval --> IRF to governement transfer shock
[Y_TR_mean,Y_TR_CI_lower,Y_TR_CI_upper]=get_CI(res_Y_TR); % Output IRFs
[Infl_TR_mean,Infl_TR_CI_lower,Infl_TR_CI_upper]=get_CI(res_Infl_TR); % Inflation IRF
[TR_TR_mean,TR_TR_CI_lower,TR_TR_CI_upper]=get_CI(res_TR_TR); % Gov. investment IRF
[fm_TR_mean,fm_TR_CI_lower,fm_TR_CI_upper]=get_CI(res_fm_TR); % Period multiplier
[FM_TR_mean,FM_TR_CI_lower,FM_TR_CI_upper]=get_CI(res_FM_TR); % Cumulative multiplier


%% Build plots (IRF and multipliers) and tables--> Gov. consumption shock

t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval
    % It should look like the paper IRF

figure ('Name','Response to a transitory gov. consumption shock')
subplot(3,1,1) % Output
f_plot_with_CI(t,Y_G_mean,Y_G_CI_lower,Y_G_CI_upper,"Output response : $$(Y_t-Y_t^{SS})/Y_0^{SS}$$","Latex",1)

subplot(3,1,2) % Inflation
f_plot_with_CI(t,Infl_G_mean,Infl_G_CI_lower,Infl_G_CI_upper,"Inflation response : $$\pi_t- \pi_t^{SS}$$","Latex",1)

subplot(3,1,3) % Government consumption
f_plot_with_CI(t,G_G_mean,G_G_CI_lower,G_G_CI_upper,"Gov. consumption response : $$(G_t-G_t^{SS})/G_0^{SS}$$","Latex",1)
saveas(gcf,'Figure_1_IRF_cons','png')

% Plot fiscal multipliers

horizon=30; % use shorter time horizon to plot fiscal multipliers
t_small=t(1:horizon);

figure ('Name','Gov. consumption fiscal multipliers')
subplot(2,1,1)
f_plot_with_CI(t_small,fm_G_mean(1:horizon),fm_G_CI_lower(1:horizon),fm_G_CI_upper(1:horizon),"Period gov. consumption multiplier : $$fm_t$$","Latex",0)

subplot(2,1,2)
f_plot_with_CI(t_small,FM_G_mean(1:horizon),FM_G_CI_lower(1:horizon),FM_G_CI_upper(1:horizon),"Cumulative gov. consumption multiplier : $$FM_t$$","Latex",0)
saveas(gcf,'Figure_2_multipliers_cons','png')

% Fiscal multipliers in a table

% Short run
varname={'5%','Mean','95%'};
Q1_G = table([fm_G_CI_lower(1);FM_G_CI_lower(1)],[fm_G_mean(1);FM_G_mean(1)],[fm_G_CI_upper(1);FM_G_CI_upper(1)], 'VariableNames',varname);
Q2_G = table([fm_G_CI_lower(2);FM_G_CI_lower(2)],[fm_G_mean(2);FM_G_mean(2)],[fm_G_CI_upper(2);FM_G_CI_upper(2)],'VariableNames',varname);
Q3_G = table([fm_G_CI_lower(3);FM_G_CI_lower(3)],[fm_G_mean(3);FM_G_mean(3)],[fm_G_CI_upper(3);FM_G_CI_upper(3)],'VariableNames',varname);
Q4_G = table([fm_G_CI_lower(4);FM_G_CI_lower(4)],[fm_G_mean(4);FM_G_mean(4)],[fm_G_CI_upper(4);FM_G_CI_upper(4)],'VariableNames',varname);

Table_sr_G=table(Q1_G,Q2_G,Q3_G,Q4_G,'VariableNames',{'Quarter 1','Quarter 2','Quarter 3','Quarter 4'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

%long-run 
Q8_G = table([fm_G_CI_lower(8);FM_G_CI_lower(8)],[fm_G_mean(8);FM_G_mean(8)],[fm_G_CI_upper(8);FM_G_CI_upper(8)], 'VariableNames',varname);
Q12_G = table([fm_G_CI_lower(12);FM_G_CI_lower(12)],[fm_G_mean(12);FM_G_mean(12)],[fm_G_CI_upper(12);FM_G_CI_upper(12)],'VariableNames',varname);
Q16_G = table([fm_G_CI_lower(16);FM_G_CI_lower(16)],[fm_G_mean(16);FM_G_mean(16)],[fm_G_CI_upper(16);FM_G_CI_upper(16)],'VariableNames',varname);
Q20_G = table([fm_G_CI_lower(20);FM_G_CI_lower(20)],[fm_G_mean(20);FM_G_mean(20)],[fm_G_CI_upper(20);FM_G_CI_upper(20)],'VariableNames',varname);

Table_lr_G=table(Q8,Q12,Q16,Q20,'VariableNames',{'Year 2','Year 3','Year 4','Year 5'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

% Store the values in Tables.xlsx

writetable(Q1_G,'Tables_12_multipliers.xlsx','Sheet','SR','Range','C3')
writetable(Q2_G,'Tables_12_multipliers.xlsx','Sheet','SR','Range','F3')
writetable(Q3_G,'Tables_12_multipliers.xlsx','Sheet','SR','Range','I3')
writetable(Q4_G,'Tables_12_multipliers.xlsx','Sheet','SR','Range','L3')

writetable(Q8_G,'Tables_12_multipliers.xlsx','Sheet','LR','Range','C3')
writetable(Q12_G,'Tables_12_multipliers.xlsx','Sheet','LR','Range','F3')
writetable(Q16_G,'Tables_12_multipliers.xlsx','Sheet','LR','Range','I3')
writetable(Q20_G,'Tables_12_multipliers.xlsx','Sheet','LR','Range','L3')

%% Build plots (IRF and multipliers) and tables--> Gov. investment shock

t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval
    % It should look like the paper IRF

figure ('Name','Response to a transitory gov. investment shock')
subplot(3,1,1) % Output
f_plot_with_CI(t,Y_GI_mean,Y_GI_CI_lower,Y_GI_CI_upper,"Output response : $$(Y_t-Y_t^{SS})/Y_0^{SS}$$","Latex",1)

subplot(3,1,2) % Inflation
f_plot_with_CI(t,Infl_GI_mean,Infl_GI_CI_lower,Infl_GI_CI_upper,"Inflation response : $$\pi_t- \pi_t^{SS}$$","Latex",1)

subplot(3,1,3) % Government investment
f_plot_with_CI(t,IG_GI_mean,IG_GI_CI_lower,IG_GI_CI_upper,"Gov. investment response : $$(GI_t-GI_t^{SS})/GI_0^{SS}$$","Latex",1)
saveas(gcf,'Figure_3_IRF_inv','png')

% Plot fiscal multipliers

horizon=30; % use shorter time horizon to plot fiscal multipliers
t_small=t(1:horizon);

figure ('Name','Gov. investment fiscal multipliers')
subplot(2,1,1)
f_plot_with_CI(t_small,fm_GI_mean(1:horizon),fm_GI_CI_lower(1:horizon),fm_GI_CI_upper(1:horizon),"Period gov. investment multiplier : $$fm_t$$","Latex",0)

subplot(2,1,2)
f_plot_with_CI(t_small,FM_GI_mean(1:horizon),FM_GI_CI_lower(1:horizon),FM_GI_CI_upper(1:horizon),"Cumulative gov. investment multiplier : $$FM_t$$","Latex",0)
saveas(gcf,'Figure_4_multipliers_inv','png')

% Fiscal multipliers in a table

% Short run
varname={'5%','Mean','95%'};
Q1_GI = table([fm_GI_CI_lower(1);FM_GI_CI_lower(1)],[fm_GI_mean(1);FM_GI_mean(1)],[fm_GI_CI_upper(1);FM_GI_CI_upper(1)], 'VariableNames',varname);
Q2_GI = table([fm_GI_CI_lower(2);FM_GI_CI_lower(2)],[fm_GI_mean(2);FM_GI_mean(2)],[fm_GI_CI_upper(2);FM_GI_CI_upper(2)],'VariableNames',varname);
Q3_GI = table([fm_GI_CI_lower(3);FM_GI_CI_lower(3)],[fm_GI_mean(3);FM_GI_mean(3)],[fm_GI_CI_upper(3);FM_GI_CI_upper(3)],'VariableNames',varname);
Q4_GI = table([fm_GI_CI_lower(4);FM_GI_CI_lower(4)],[fm_GI_mean(4);FM_GI_mean(4)],[fm_GI_CI_upper(4);FM_GI_CI_upper(4)],'VariableNames',varname);

Table_sr_GI=table(Q1_GI,Q2_GI,Q3_GI,Q4_GI,'VariableNames',{'Quarter 1','Quarter 2','Quarter 3','Quarter 4'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

%long-run 
Q8_GI = table([fm_GI_CI_lower(8);FM_GI_CI_lower(8)],[fm_GI_mean(8);FM_GI_mean(8)],[fm_GI_CI_upper(8);FM_GI_CI_upper(8)], 'VariableNames',varname);
Q12_GI = table([fm_GI_CI_lower(12);FM_GI_CI_lower(12)],[fm_GI_mean(12);FM_GI_mean(12)],[fm_GI_CI_upper(12);FM_GI_CI_upper(12)],'VariableNames',varname);
Q16_GI = table([fm_GI_CI_lower(16);FM_GI_CI_lower(16)],[fm_GI_mean(16);FM_GI_mean(16)],[fm_GI_CI_upper(16);FM_GI_CI_upper(16)],'VariableNames',varname);
Q20_GI = table([fm_GI_CI_lower(20);FM_GI_CI_lower(20)],[fm_GI_mean(20);FM_GI_mean(20)],[fm_GI_CI_upper(20);FM_GI_CI_upper(20)],'VariableNames',varname);

Table_lr_GI=table(Q8_GI,Q12_GI,Q16_GI,Q20_GI,'VariableNames',{'Year 2','Year 3','Year 4','Year 5'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

% Store the values in Tables.xlsx

writetable(Q1_GI,'Tables_34_multipliers_investment.xlsx','Sheet','SR','Range','C3')
writetable(Q2_GI,'Tables_34_multipliers_investment.xlsx','Sheet','SR','Range','F3')
writetable(Q3_GI,'Tables_34_multipliers_investment.xlsx','Sheet','SR','Range','I3')
writetable(Q4_GI,'Tables_34_multipliers_investment.xlsx','Sheet','SR','Range','L3')

writetable(Q8_GI,'Tables_34_multipliers_investment.xlsx','Sheet','LR','Range','C3')
writetable(Q12_GI,'Tables_34_multipliers_investment.xlsx','Sheet','LR','Range','F3')
writetable(Q16_GI,'Tables_34_multipliers_investment.xlsx','Sheet','LR','Range','I3')
writetable(Q20_GI,'Tables_34_multipliers_investment.xlsx','Sheet','LR','Range','L3')


%% Build plots (IRF and multipliers) and tables--> Gov. investment shock

t = 1:1:n; % x-axis

% Plot simulation mean and confidence interval
    % It should look like the paper IRF

figure ('Name','Response to a transitory gov. transfer shock')
subplot(3,1,1) % Output
f_plot_with_CI(t,Y_TR_mean,Y_TR_CI_lower,Y_TR_CI_upper,"Output response : $$(Y_t-Y_t^{SS})/Y_0^{SS}$$","Latex",1)

subplot(3,1,2) % Inflation
f_plot_with_CI(t,Infl_TR_mean,Infl_TR_CI_lower,Infl_TR_CI_upper,"Inflation response : $$\pi_t- \pi_t^{SS}$$","Latex",1)

subplot(3,1,3) % Government investment
f_plot_with_CI(t,TR_TR_mean,TR_TR_CI_lower,TR_TR_CI_upper,"Gov. transfer response : $$(TR_t-TR_t^{SS})/TR_0^{SS}$$","Latex",1)
saveas(gcf,'Figure_5_IRF_transfer','png')

% Plot fiscal multipliers

horizon=30; % use shorter time horizon to plot fiscal multipliers
t_small=t(1:horizon);

figure ('Name','Gov. transfer fiscal multipliers')
subplot(2,1,1)
f_plot_with_CI(t_small,fm_TR_mean(1:horizon),fm_TR_CI_lower(1:horizon),fm_TR_CI_upper(1:horizon),"Period gov. transfer multiplier : $$fm_t$$","Latex",0)

subplot(2,1,2)
f_plot_with_CI(t_small,FM_TR_mean(1:horizon),FM_TR_CI_lower(1:horizon),FM_TR_CI_upper(1:horizon),"Cumulative gov. transfer multiplier : $$FM_t$$","Latex",0)
saveas(gcf,'Figure_6_multipliers_transfer','png')

% Fiscal multipliers in a table

% Short run
varname={'5%','Mean','95%'};
Q1_TR = table([fm_TR_CI_lower(1);FM_TR_CI_lower(1)],[fm_TR_mean(1);FM_TR_mean(1)],[fm_TR_CI_upper(1);FM_TR_CI_upper(1)], 'VariableNames',varname);
Q2_TR = table([fm_TR_CI_lower(2);FM_TR_CI_lower(2)],[fm_TR_mean(2);FM_TR_mean(2)],[fm_TR_CI_upper(2);FM_TR_CI_upper(2)],'VariableNames',varname);
Q3_TR = table([fm_TR_CI_lower(3);FM_TR_CI_lower(3)],[fm_TR_mean(3);FM_TR_mean(3)],[fm_TR_CI_upper(3);FM_TR_CI_upper(3)],'VariableNames',varname);
Q4_TR = table([fm_TR_CI_lower(4);FM_TR_CI_lower(4)],[fm_TR_mean(4);FM_TR_mean(4)],[fm_TR_CI_upper(4);FM_TR_CI_upper(4)],'VariableNames',varname);

Table_sr_TR=table(Q1_TR,Q2_TR,Q3_TR,Q4_TR,'VariableNames',{'Quarter 1','Quarter 2','Quarter 3','Quarter 4'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

%long-run 
Q8_TR = table([fm_TR_CI_lower(8);FM_TR_CI_lower(8)],[fm_TR_mean(8);FM_TR_mean(8)],[fm_TR_CI_upper(8);FM_TR_CI_upper(8)], 'VariableNames',varname);
Q12_TR = table([fm_TR_CI_lower(12);FM_TR_CI_lower(12)],[fm_TR_mean(12);FM_TR_mean(12)],[fm_TR_CI_upper(12);FM_TR_CI_upper(12)],'VariableNames',varname);
Q16_TR = table([fm_TR_CI_lower(16);FM_TR_CI_lower(16)],[fm_TR_mean(16);FM_TR_mean(16)],[fm_TR_CI_upper(16);FM_TR_CI_upper(16)],'VariableNames',varname);
Q20_TR = table([fm_TR_CI_lower(20);FM_TR_CI_lower(20)],[fm_TR_mean(20);FM_TR_mean(20)],[fm_TR_CI_upper(20);FM_TR_CI_upper(20)],'VariableNames',varname);

Table_lr_TR=table(Q8_TR,Q12_TR,Q16_TR,Q20_TR,'VariableNames',{'Year 2','Year 3','Year 4','Year 5'},'RowNames',{'period fiscal multiplier','cumulative fiscal multiplier'});

% Store the values in Tables.xlsx

writetable(Q1_TR,'Tables_56_multipliers_transfer.xlsx','Sheet','SR','Range','C3')
writetable(Q2_TR,'Tables_56_multipliers_transfer.xlsx','Sheet','SR','Range','F3')
writetable(Q3_TR,'Tables_56_multipliers_transfer.xlsx','Sheet','SR','Range','I3')
writetable(Q4_TR,'Tables_56_multipliers_transfer.xlsx','Sheet','SR','Range','L3')

writetable(Q8_TR,'Tables_56_multipliers_transfer.xlsx','Sheet','LR','Range','C3')
writetable(Q12_TR,'Tables_56_multipliers_transfer.xlsx','Sheet','LR','Range','F3')
writetable(Q16_TR,'Tables_56_multipliers_transfer.xlsx','Sheet','LR','Range','I3')
writetable(Q20_TR,'Tables_56_multipliers_transfer.xlsx','Sheet','LR','Range','L3')

% Go back to original path
cd('..');