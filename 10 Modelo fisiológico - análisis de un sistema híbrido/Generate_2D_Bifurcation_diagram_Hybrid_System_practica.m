function Generate_2D_Bifurcation_diagram_Hybrid_System_practica

close all
clear 
clc

tic

%%%  Generate the bifurcation diagram with 4 qualitative behaviours.
%%% Elisa Domínguez Hüttinger, re-using code from Panayiotis Christodoulides.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Los parámetros se encuentran guardados en la estructura:
%load Nominal_parameters_AD_model.mat
P_env=95;
gamma_B=1;
delta_P=1;
kappa_B= 0.5;
gamma_lip=10;
delta_B=0.1; 

Ron=16.7000;
Pplus=40;
Pminus= 26.6000;

% Initialize
column=1;
row=1;

for kappa_P=0:0.03:1;
    for alpha_I=0.3:-0.01:0
         
               
        %%% Determine the qualitative behaviour, depending on the position
        %%% of the focal points
        par=[P_env, kappa_P, gamma_B, alpha_I, delta_P, kappa_B, gamma_lip, 1, delta_B]; % parameters in Pana notation     
        Behaviour=determine_qualitative_behaviour(par, Ron, Pplus, Pminus);
        % 1: Monostable healhty, 2: bistable (will be refined for specific
        % initial conditions to: 21: bistable healhty branch, 24: bistable
        % unhealhty branch), 3: oscillations, and 4: monostable unhealhty
        
               
        
        %% Fill in the matrices
        Parameter_kappaBmatrix(column, row)=kappa_P;
        Parameter_alphaImatrix(column, row)=alpha_I;
        
        % Qualitative behaviour
        Behaviour_matrix(column,row)= Behaviour;
        
               
        
        column=column+1;
    end;
    column=1;
    row=row+1;
    
end;

%%%% Behaviour matrix
AD_figure=figure;
set(AD_figure, 'Position', [100 100 500 500]);

imagesc(Parameter_kappaBmatrix(1,:),Parameter_alphaImatrix(:,1),Behaviour_matrix);
ylabel('Strength of immune response (alpha I) ','FontSize',12,'FontName','Arial');%
xlabel('FLG defficiency (kappa P) ','FontSize',12,'FontName','Arial');%
hold on

%%%% Find the lines separating these behaviours:
% limits 
kappa_P_crit_index_o=find(Behaviour_matrix(1,:)==3,1)-1;
kappa_P_crit_o=Parameter_kappaBmatrix(end,kappa_P_crit_index_o);
alpha_I_crit_o=Parameter_alphaImatrix(1,1);


kappa_P_crit_index_B=find(Behaviour_matrix(end,:)>1,1);
kappa_P_crit_B=(Parameter_kappaBmatrix(end,kappa_P_crit_index_B));
alpha_I_crit_B=0;

kappa_P_crit_index_U=find(Behaviour_matrix(1,:)==4,1);
kappa_P_crit_U=(Parameter_kappaBmatrix(end,kappa_P_crit_index_U));


alpha_I_crit_index_U=find(Behaviour_matrix(:,end)==4,1);
alpha_I_crit_U=(Parameter_alphaImatrix(alpha_I_crit_index_U,end));

% Use the obtained coordinates to plot those lines on the image sc plot
line([kappa_P_crit_o,kappa_P_crit_o],[0, 0.3],'Color','k', 'Linewidth', 3);
line([kappa_P_crit_B,1],[alpha_I_crit_B, alpha_I_crit_U], 'Color','k', 'Linewidth', 3)
title('Qualitative behaviour','FontSize',12,'FontName','Arial');%
axis square 
set(gca, 'YDir', 'normal');
cb = colorbar;
set(cb, 'YTick',[1,2,3,4],'YTickLabel',{'Mono-H';'Bi';'Osc';'Mono-U'} );%,'YTickLabel',{'Mono-H','Bi','Osc','Mono-U'})
set(gca,'FontSize',10,'FontName','Arial');%,'FontWeight','Bold');





toc

end


function Behaviour=determine_qualitative_behaviour(par, PAR_ss, Splus, Sminus)


% Order of Parameters par is important:
%par_name = {'Sout','tildeP','epsilonB','fIP','I_0','b_{pre}','kL','tildeB','dK'};
%%%% Given parameter values, returns the behaviour of the system:
% possibilities: 
%1  monostable healhty
%2 bistable healhty
%3 oscillations
% monostable unhealthy;
%Behaviour

%%%%% Get the steady states of the low and the high subsystems
[SS_Low,SS_High] = AnalyticalSteadyStates(par,PAR_ss);


%%%%% Long-term behaviour depending on location of steady-states
            
            if isreal(SS_High) == 0 
                if SS_Low(1,1) < Splus % monostable healthy
                    Behaviour=1;
                elseif SS_Low(1,1) >= Splus % oscillatory behaviour
                   Behaviour=3;
                end
            elseif isempty(SS_High) == 1  % monostable healthy  
                if SS_Low(1,1) < Splus % monostable healthy
                     Behaviour=1;
                elseif SS_Low(1,1) >= Splus % oscillatory behaviour
                      Behaviour=3;
                end                            
            elseif SS_High(1,1) > Splus && SS_Low(1,1) < Sminus % bistable
                    Behaviour=2;
            elseif (SS_High(1,1) <= Splus && SS_High(1,1) >= Sminus) && (SS_Low(1,1) <= Splus && SS_Low(1,1) >= Sminus)     % bistable
                    Behaviour=2;
            elseif (SS_High(1,1) <= Splus && SS_High(1,1) >= Sminus) && SS_Low(1,1) > Splus % monostable unhealthy
                    Behaviour=4;
            elseif (SS_High(1,1) <= Splus && SS_High(1,1) >= Sminus) && SS_Low(1,1) < Sminus % bistable
                   Behaviour=2;
            elseif (SS_Low(1,1) <= Splus && SS_Low(1,1) >= Sminus) && SS_High(1,1) > Splus % bistable
                    Behaviour=2;
            elseif (SS_Low(1,1) <= Splus && SS_Low(1,1) >= Sminus) && SS_High(1,1) < Sminus % monostable healthy
                    Behaviour=1;
            elseif SS_Low(1,1) > Splus && SS_High(1,1) < Sminus % oscillatory behaviour
                    Behaviour=3;
            elseif SS_Low(1,1) < Sminus && SS_High(1,1) < Sminus % monostable healthy
                   Behaviour=1;
            elseif SS_High(1,1) > Splus && SS_Low(1,1) > Splus % monostable unhealthy
                    Behaviour=4;
            end 
     
   
end

 function [SS_Low,SS_High] = AnalyticalSteadyStates(final_par,PAR_ss) 
 
   InflammationP_ss=PAR_ss;
    alpha = final_par(6)/(1+final_par(7)*PAR_ss);
    beta = final_par(4)*InflammationP_ss+final_par(5);
    gamma = final_par(1)*final_par(2);
    delta = alpha/final_par(8)-final_par(9)*6.71;
    epsilon = final_par(9)*0.45;

    %Low

    Bss_Low = final_par(8);
    Sss_Low = final_par(1)*final_par(2)/((final_par(8)+final_par(3))*final_par(5));
    SS_Low = [Sss_Low,Bss_Low];

    %High

    B_ss_High = [(-(delta*beta*final_par(3)+epsilon*gamma-alpha*beta)+sqrt((delta*beta*final_par(3)+epsilon*gamma-alpha*beta)^2+4*delta*beta^2*alpha*final_par(3)))/(2*delta*beta);
                (-(delta*beta*final_par(3)+epsilon*gamma-alpha*beta)-sqrt((delta*beta*final_par(3)+epsilon*gamma-alpha*beta)^2+4*delta*beta^2*alpha*final_par(3)))/(2*delta*beta)];

    B_ss_High = B_ss_High(B_ss_High >= 0); 
    S_ss_High = gamma./((B_ss_High+final_par(3))*beta);
    SS_High = [S_ss_High,B_ss_High];
end
