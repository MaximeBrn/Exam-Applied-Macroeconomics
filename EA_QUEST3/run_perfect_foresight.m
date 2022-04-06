% Run and plot IRF of replication file 

% EA_QUEST3 model

clc;
close all;

% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run replication dynare file


dynare EA_Quest3_rep_perfect_foresight.mod

gy_simul=oo_.endo_simul(find(strcmp(M_.endo_names,'E_GY')),:)';
gg_simul=oo_.endo_simul(find(strcmp(M_.endo_names,'E_GG')),:)';
figure
t=1:1:42;
plot(t,gy_simul,'LineWidth',2);

figure
t=1:1:42;
plot(t,cumsum(gg_simul-GY0),'LineWidth',2);




fm=cumsum(gy_simul-GY0)./cumsum(gg_simul-GY0)/GSN;



