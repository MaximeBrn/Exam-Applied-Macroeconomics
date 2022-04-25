% Run and plot IRF of replication file 

% EA_QUEST3 model

clc;
close all;

% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run replication dynare file


dynare EA_Quest3_rep_perfect_foresight.mod

GY0_vec = ones(1, length(oo_.endo_simul(1,:)))' * GY0;

Y_baseline=cumprod(1+oo_.endo_simul(find(strcmp(M_.endo_names,'E_GY')),:))'-cumprod(1+GY0_vec);
G_baseline=cumprod(1+oo_.endo_simul(find(strcmp(M_.endo_names,'E_GG')),:))'-cumprod(1+GY0_vec);
fm_baseline=Y_baseline./G_baseline/GSN;
FM_baseline=cumsum(Y_baseline)./cumsum(G_baseline)/GSN;

figure
subplot(2,1,1)
t=2:1:41;
plot(t,Y_baseline(2:41),'LineWidth',2);

subplot(2,1,2)
t=2:1:41;
plot(t,G_baseline(2:41),'LineWidth',2);

figure
subplot(2,1,1)
t=2:1:41;
plot(t,fm_baseline(2:41),'LineWidth',2);

subplot(2,1,2)
t=2:1:41;
plot(t,FM_baseline(2:41),'LineWidth',2);




