%%%% Nominal parameters of the AD model

P_env=95; % Environmental stress load 
gamma_B=1; % Barrier-mediated inhibition of pathogen infiltration 
kappa_P=0.6; % Nominal skin permeability 
alpha_I=0.25; %  Rate of pathogen eradication by innate immune responses 
delta_P= 1; % Basal pathogen death rate 
kappa_B=0.5; %Barrier production rate 
gamma_R= 10; %Innate immunity-mediated inhibition of barrier production 
delta_B= 0.1; % Rate of Kallikrein-dependent barrier degradation 
P_minus=26.6; % Receptor inactivation threshold 
P_plus=40; % Receptor activation threshold
R_off=0; %Receptor-off level 
R_on=16.7; % Receptor-on level 
K_off=0; % Kallikrein-off level 
m_on= 0.45; % Slope of the linear relation between $P(t)$ and $K_{\rm on}$ 
beta_on= 6.71; % Y-intercept of the linear relation between $P(t)$ and $K_{\rm on}$

par_vector=[P_env, gamma_B, kappa_P, delta_P, kappa_B,alpha_I,  gamma_R,delta_B, P_minus, P_plus, R_off,R_on, K_off,m_on, beta_on];

par_vector_off=[P_env, gamma_B, kappa_P, delta_P, kappa_B];

par_vector_on=[P_env, gamma_B, kappa_P, delta_P, kappa_B,alpha_I, gamma_R,delta_B, R_on, m_on, beta_on];
