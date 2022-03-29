# Exam-Applied-Macroeconomics
DSGE 3 Mutliplicateur Budgétaire / Fiscal Multiplier

# The topic

A l’aide d’un modèle DSGE (déjà) estimé par techniques bayesiennes sur données zone euro, déterminez la réponse du PIB et de l’inflation à des chocs transitoires et permanents de dépense publique. Quantifier les multiplicateurs budgétaires de court et long termes. Caractériser l’intervalle de confiance bayesien entourant les multiplicateurs.

Using a DSGE model (already) estimated with bayesian methods for the euro area, compute the reaction of GDP and inflation to transitory and permanent shocks to public spending. Estimate the short and long term fiscal multipliers. Identify the Bayesian confidence intervals around those multipliers.

# The choice of the DSGE model

The model should be such that:
1. it is estimated with bayesian methods;
2. it is estimated for the euro area;
3. it is tailored to study the impact of government spending on output and inflation.

In class talked about the [Macroeconomic Model dataBase (MBB)](https://www.macromodelbase.com/) which is an archive of macroeconomic models. On the website is uploaded a [short list](https://www.macromodelbase.com/files/documentation_source/mmb-model-list.pdf) and a [detailed list](https://www.macromodelbase.com/files/documentation_source/mmb-model-description.pdf) of the models available on the platform. We note that there are the main models we discussed in class. 

It seems that there are two models that match our criteria:
- EA_SW03: Smets and Wouters (2003);
- EA_QUEST3: QUEST III Euro Area Model of the DG-ECFIN EU, Ratto et al. (2009)

QUEST III seems to be a better choice because:
- it includes a share of non-Ricardian households (hand-to-mouth) which amplify the role of fiscal policy;
- it included 3 fiscal policy variables (government consumption, governement investment, and transfers) which is good to study fiscal multipliers;
- it has been used in a comparative analysis to assess the impact of the 2008-2009 fiscal packages (see Cwik and Wieland et al., 2011, Keynesian government spending multipliers and spillovers in the euro area).

# Replicating the paper results

On the github of the MMB, we can find a [replication of the model QUEST III](https://github.com/IMFS-MMB/mmb-rep/tree/master/EA_QUEST3). There is also the list of the variables used in the model [here](https://github.com/IMFS-MMB/mmb-gui-mlab/blob/master/MODELS/EA_QUEST3/list_of_variables.xls).

Before answering the question, we need to make sure that the code works and replicates the results of the authors.

## How to run the code?

In the folder `EA_QUEST3`, run the file `run.m`. This file plot the response function. It is based on the file `EA_QUEST3_rep.mod` which is located in the folder `EA_QUEST3_rep`. 

Note: make sure you have configured Dynare on Matlab permanently (i.e. using `Set Path` instead of `Addpath`).

## Compare results

Slightly modifying and adapting the code we can replicate the results of Ratto et al. (2008) on the response to government consumption.

The code replication       |  Ratto et al. paper
:-------------------------:|:-------------------------:
![image](https://user-images.githubusercontent.com/37322244/160674359-fc7b354d-7c41-45b3-ba9a-45b94932fe6c.png)  |  ![image](https://user-images.githubusercontent.com/37322244/160674054-fbb50d1c-ef76-4439-902e-a2ddf45c6905.png)
![image](https://user-images.githubusercontent.com/37322244/160674415-07a0a1a3-9df0-4973-91b5-8c08e2ddaeee.png) | ![image](https://user-images.githubusercontent.com/37322244/160674211-defb202d-d210-4073-9bc6-f67537f207d4.png)

Regarding fiscal multipliers, in Ratto et al. (2008) it written : 

>To assess the impact of the government spending shocks on output in terms of traditional “multipliers”, the impact effect for a 1% of government spending shock on GDP is 0.73 in the first quarter, falling to 0.45 in the fourth.
>
In our code we find  0.7365 in the first quarter and 0.4595 in the fourth. It seems consistent with the authors' results. 

We have the following plot for other time horizons.

![image](https://user-images.githubusercontent.com/37322244/160675225-b61cf418-ac1b-4db8-b837-893945bc9f6f.png)

For robustness, we also checked that we find the same response for a government investment shock.

# Awnsering the research question

##  Compute the reaction of GDP and inflation...

- [ ] Remove the other variable from the analysis (easy)

##  ...to transitory and permanent shocks to public spending.

- [X] Response to a transitory public consumption shock
- [ ] Response to a permanent public consumption shock 

Question to ask:
- How to make a shock permanent when using stoch_simul? Element of answer [here](https://forum.dynare.org/t/permanent-shock-in-stochastic-simulation/2399/3)
- On which parameters are we allowed to play?
- Should we only focus on government consumption? Or look also at investment and transfers?

## Estimate the short and long term fiscal multipliers

- [X] Find the formula used in Ratto et al. (2008)

## Identify the Bayesian confidence intervals around those multipliers

- [ ] Find the approach to compute the Bayesian confidence intervals

Question to ask:
- How to compute the Bayesian confidence intervals? Look at Appendix 2 of Ratto et al. (2008)?
- Should we take into account other source of uncertainty? Element of answer in Smets and Wouter, 2004, Forecasting with a bayesian DSGE model (page 10/32).

An idea to compute the Bayesian uncertainty bounds:
- Draw a list of parameters estimate from their posterior distribution
- Run for each draw a simulation and store the response
- Identify the 90% bound of responses
- Plot with confidence interval
- Deduce the confidence intervals of the fiscal multiplier

https://forum.dynare.org/t/dsge-irfs-from-bayesian-estimation/3477/3
