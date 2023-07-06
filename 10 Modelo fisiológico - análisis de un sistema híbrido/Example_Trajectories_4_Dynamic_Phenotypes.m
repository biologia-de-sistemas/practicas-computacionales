close all
clear all
clc

% load the parameter values
run Nominal_parameters.m

t_ODE=[0 10];
IC_States=[50 1];
IC_Switches=1;


jj=0;
%%
for ii=1:1:4
   
    %% select the values of the bifurcation parameters
    if ii==1
        % healthy control
        kappa_P=.6; % healthy control
        alpha_I=.3;
        Col='k';
    elseif ii==2
        % Damage
        kappa_P=1;
        alpha_I=.05;
        Col='r';
     elseif ii==3
        % Oscillations
        kappa_P=1;
        alpha_I=.3;
        Col='b';
    elseif ii==4
        % Bistability
        kappa_P=.6;
        alpha_I=.05;
        Col='y';
    end


%% Numerical integration
[T,Y,Switch_state] = Integrate_Full_Hybrid_Model(t_ODE,IC_States,IC_Switches,P_env, gamma_B, kappa_P, delta_P, kappa_B,alpha_I,  gamma_R,delta_B, P_minus, P_plus, R_off,R_on, K_off,m_on, beta_on);


%% Focal points
[Pss_off, Bss_off]=StableFocalPoints_off(P_env ,gamma_B, kappa_P, delta_P);
[Pss_on, Bss_on, stable_state]=StableFocalPoints_on(P_env, gamma_B, kappa_P, delta_P, kappa_B,alpha_I, gamma_R,delta_B, R_on, m_on, beta_on);



%%
%figure;
subplot(3,4, 1+jj)
plot( T, Y(1,:),'LineWidth',2, 'Color', Col)
hold on
line([0, T(end)], [P_minus, P_minus],'Color', 'm');
axis square
line([0, T(end)], [P_plus, P_plus],'Color', 'c');
line([0, T(end)], [Pss_on, Pss_on],'Color', 'm', 'LineStyle', '--');
axis square
line([0, T(end)], [Pss_off, Pss_off],'Color', 'c', 'LineStyle', '--');

ylim([0, 60])
ylabel('P(t)')
title(['\kappa_P=' num2str(kappa_P) '\alpha_I=' num2str(alpha_I) ])


subplot(3,4,5+jj)
plot( T, Y(2,:),'LineWidth',2,'Color',Col)
axis square
ylim([0, 1.1])
ylabel('B(t)')

subplot(3,4, 9+jj)
plot( T, Switch_state, 'LineWidth',2,'Color', Col)
axis square
ylim([0, 1.1])
xlabel('Time [h]')
ylabel('R switch state')


jj=jj+1;
end



