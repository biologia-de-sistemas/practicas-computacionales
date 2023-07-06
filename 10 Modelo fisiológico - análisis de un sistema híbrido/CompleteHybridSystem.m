function dydt = CompleteHybridSystem(~,y,R_t_state,P_env, gamma_B, kappa_P, delta_P, kappa_B,alpha_I, gamma_R,delta_B, R_on, m_on, beta_on)

%R_t_state is the state of the switch, either 0 or 1
P=y(1);
B=y(2);

dydt=[  (P_env*kappa_P/(1+gamma_B*B))-P*(alpha_I*R_t_state*R_on+delta_P); % Pathogen load ODE
        (kappa_B/((1+R_t_state*gamma_R*R_on)))*(1-B)-R_t_state*delta_B*(m_on*P-beta_on)*B]; % Barrier function ODE
                
         
   
end