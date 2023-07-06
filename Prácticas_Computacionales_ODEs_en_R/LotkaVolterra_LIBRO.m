function LotkaVolterra_LIBRO
%  Como siempre, buena costumbre limpiar nuestro espacio de trabajo antes
%  de comenzar
close all; clear all; clc

% condicion inicial:
x0=1; y0=1; IC=[x0 y0]; 

%Parametros
alpha = 2/3; beta= 4/3; gamma=1; delta=1;

%tiempo de integracion
tspan = [0 50];

%integramos  numericamente:
[t,y] = ode45(@(t,y)Lotka_Volterra(t,y, alpha, beta, gamma, delta),tspan,IC);

% pintamos la funcion 
figure
plot(t, y(:,1),'color','b', 'LineWidth',2)
xlabel('Tiempo');
ylabel('Lotka Volterra');
hold on
plot(t, y(:,2),'color','r', 'LineWidth',2)
legend('Presa', 'Depredador')
axis square
set(gcf, 'Position', [100 100 300 300]); 
title(['x_0=' num2str(x0) ', y_0=' num2str(y0), ', \alpha=' num2str(alpha)   ', \beta=' num2str(beta) ', \gamma=' num2str(gamma) ', \delta=' num2str(delta)]);

%
figure
plot(y(:,1), y(:,2),'color','k', 'LineWidth',2)
xlabel('Presa');
ylabel('Depredador');
axis square
set(gcf, 'Position', [100 100 300 300]); 

end

function dydt = Lotka_Volterra(t,y, alpha, beta, gamma, delta)

dydt =[alpha*y(1)-beta*y(1)*y(2); % linear growth of pray, death by the predator 
       delta*y(1)*y(2)-gamma*y(2)]; % growth of predator depends on prey, linear death
   
end