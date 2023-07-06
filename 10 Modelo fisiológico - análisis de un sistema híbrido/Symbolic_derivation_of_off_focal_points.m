close all
clear all
clc


syms P_env gamma_B kappa_P delta_P kappa_B%%
syms P B
%%

dPdt_off=  (P_env*kappa_P/(1+gamma_B*B))-P*delta_P; % Pathogen load ODE
dBdt_off=  (kappa_B*(1-B));


%%
system_off=[dPdt_off dBdt_off];
sol=solve(system_off, [P B]);


%%
Pss_off=(P_env*kappa_P)/(delta_P + delta_P*gamma_B);
Bss_off=1;

%%
jacobian(system_off, [P B])

%%
Jac=[[ -delta_P, -(P_env*gamma_B*kappa_P)/(B*gamma_B + 1)^2]
     [        0,                                   -kappa_B]];

%%
Jac_off=subs(subs(Jac, B, Bss_off), P, Pss_off);
%%


eig(Jac_off);

