
% Integration of complete Hybrid system: Includes Pathogen and Barrier
% Operates using the Events location mechanism in MATLAB and recursively
% integrates the system.

%%% Inputs %%%
% t_ODE = integration time
% IC_states = initial conditions for states
% IC_Switches = initial conditions for switches
% par = model parameters
% Pminus, Pplus = thresholds

%%% Outputs %%%
% T = final time vector, size(T) = [1 X], X = number of points in the
% integration
% Y = State vectors, size(Y) = [1 X]
% Switch_state = shows if switches are ON or OFF




function [T,Y,Switch_state] = Integrate_Full_Hybrid_Model(t_ODE,IC_States,IC_Switches,P_env, gamma_B, kappa_P, delta_P, kappa_B,alpha_I,  gamma_R,delta_B, P_minus, P_plus, R_off,R_on, K_off,m_on, beta_on);

T = [];
Y = [];
Switch_state = [];

IntegrationTime = t_ODE(2);

time = 0; % base case
    while time ~= IntegrationTime
        
        % Hybrid system approximations
        R_t_state = HybridParameters(IC_Switches); % this is a function declared below
        
        % Integration  & uses hybrid switching, declared below
        options = odeset('Events',@(t,y)HybridSwitchingEvent(t,y,IC_Switches,P_minus,P_plus),'RelTol',1e-10,'AbsTol',1e-10,'Refine',5);
                
        solution = ode23(@(t,y)CompleteHybridSystem(t,y,R_t_state,P_env, gamma_B, kappa_P, delta_P, kappa_B,alpha_I, gamma_R,delta_B, R_on, m_on, beta_on),t_ODE,IC_States,options);
        
        t_ODE_extended = linspace(solution.x(1),solution.x(end),5000); % NOTE: 10000 is arbitrary increase or decrease depending on requirements
        
        Solution = deval(solution,t_ODE_extended);
        
        
        % Update recursive variable
        time = solution.x(end); 
        
        % Concatinate solutions into a vector
        T = [T,t_ODE_extended]; 
        Y = [Y,Solution];
        
        % Concatinate switch state for the integration
        Switch_state = [Switch_state,repmat(IC_Switches',1,numel(t_ODE_extended))];
        
        if solution.x(end) ~= t_ODE(2) % check if while conditions have been achieved
            % Depending on what event happend change the governing system by manipuating the switches
            IC_Switches = ChangeSwitchState(solution.ie,IC_Switches);
                        
            % Re-initialize data to re-start integration
            t_ODE = [t_ODE_extended(end),IntegrationTime];
            IC_States = Solution(:,end);      
        end
    end
    
end


function Hybrid_par = HybridParameters(IC_Switches)
    % R-switch
    if IC_Switches== 0
         Hybrid_par = 0; % this defines an anonymous function. You can then write PAR_ss(1.2). The result will be 0 since x does not appear in the equality
       
    elseif IC_Switches == 1
        Hybrid_par = 1;
    end
        
   
end    


function IC_Switches = ChangeSwitchState(Switch_identity,IC_Switches)
    % Reversible switch
    if Switch_identity == 1 
        if IC_Switches(1) == 0 % if initially switched off then switch on
            IC_Switches(1) = 1;
        elseif IC_Switches(1) == 1 % if initially switched on then switch off
            IC_Switches(1) = 0;
        end
    % Irriversible switch
    elseif Switch_identity == 2
        IC_Switches(2) = 1; % switch on Gata3 switch
    end   
end


function [value,isterminal,direction] = HybridSwitchingEvent(~,y,IC_Switches,Pminus,Pplus)
    if IC_Switches(1) == 0 
        value = y(1) - Pplus;
        isterminal = 1; % 0 = don't stop
        direction = 0;  % 0 = both negative and positive direction
    elseif IC_Switches(1) == 1 
        value = y(1) - Pminus;
        isterminal = 1; % 0 = don't stop
        direction = 0;  % 0 = both negative and positive direction 
     end
end

