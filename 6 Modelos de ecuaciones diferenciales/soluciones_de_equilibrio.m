clear all
close all
clc



%% Solución de equilibrios con Matlab


%%  ejemplos de ustedes
% 1
syms r1 s In r2 r
%ds suceptibles, di infectados, dr recuperados 

ds = -r1*s*In;
di =  r1*s*In - In*r2;
dr =            In*r2;

S=solve(ds == 0, di==0, dr==0,  s, In, r)


% 2
syms x y alfa gamma beta delta epsilon
% y poblacion sumidero, x poblacion fuente
% Jesus
%epsilon=0

dx = epsilon+x*(alfa - (gamma+ beta));
dy =         x*gamma - y*delta;

Sol=solve(dx == 0, dy==0, x, y)

% gamma 0 o 1
Sol.x
Sol.y


%% Ejemplo 1: degradación lineal
% sistema de dimensión 1, lineal
syms k1 x 
dx=-x*k1
solve(dx==0, x)

%% Ejemplo 2: Producción de orden 0, degradación constante
% sistema de dimensión 1, lineal
syms k1 k2 x
dx=k1-x*k2
solve(dx==0, x)


%% Ejemplo 3: circuito de adaptación perfecta 
% este es nuestro primer ejemplo de 2D
syms k1 k2 k3 k4 S R X
solve(k1*S-k2*X*R==0, 3*S-k4*X==0, R, X)     % %dX/dt=f2(R,X)
%%
% poner slución en una función 
y1ss=@(k1, k2, k3, k4,S)k1*k4/(k2*k3);
y2ss=@(k1, k2, k3, k4,S)k3*S/k4;
y1ss(k1, k2, k3, k4,S)
y2ss(k1, k2, k3, k4,S);

%% Ejemplo 4: Switch con histéresis
syms alpha IL4 k_g Gata3 k
dGata3dt = alpha.*IL4 + (k_g.*Gata3.^2/(1+Gata3.^2)) - k*Gata3;
 solve( dGata3dt ==0, Gata3)
 % no hay solución general - necesitamos argumentos geométricos o numéricos
 %% si damos 
 
%Parametros
alpha = 0.02; k_g = 5; k = 1; IL4=1;
dGata3dt = alpha.*IL4 + (k_g.*Gata3.^2/(1+Gata3.^2)) - k*Gata3;

dGata3dt = (1+Gata3.^2)*alpha.*IL4 + (k_g.*Gata3.^2) - (1+Gata3.^2)*k*Gata3;

solve( dGata3dt ==0, Gata3)



%% Ejemplo 5: Angeli
%% encontrar las ceroclinas de la función de angeli


%dx/dt = alpha1*(1-x) -beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
% dx/dt=fx(x,y) --> fx(x, y) =0 --> dx/dt=0 --> x está en equilibrio.
% xss=g(y)

%fy= alpha1*(1-x) - beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
% análogo
% ...
% yss=h(x)

% g(y) y h(x) se llaman ceroclinas (nullclines) de x y de y,
% respectivamente.

alpha1=1; alpha2=1; beta1=200; beta2=10; gamma1=4; gamma2=4; K1=30; K2=1; %v=3;



for  v=[0.1, 1, 3]
    
%%


%%
%%% ceroclina de x: dx=0
%0 = alpha1*(1-x)-beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
%alpha1*(1-x)=beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
%(1-x)/x=beta1*(v*y)^gamma1/(K1+(v*y)^gamma1)/alpha1
%1/x-1=beta1*(v*y)^gamma1/(K1+(v*y)^gamma1)/alpha1
%1/x=(beta1*(v*y)^gamma1/(K1+(v*y)^gamma1)/alpha1)+1
%x=((beta1*(v*y)^gamma1/(K1+(v*y)^gamma1)/alpha1)+1)^(-1)

xss=@(y)1./((beta1*(v.*y).^gamma1./(K1+(v.*y).^gamma1)./alpha1)+1);
%%
  
%%% encontrar la ceroclina de y
% dy = alpha2*(1-y)-beta2*y*x^gamma2/(K2+x^gamma2)=0
% alpha2*(1-y)=beta2*y*x^gamma2/(K2+x^gamma2)
% (1-y)/y=beta2*x^gamma2/(K2+x^gamma2)/alpha2
% 1/y-1=beta2*x^gamma2/(K2+x^gamma2)/alpha2
% 1/y=beta2*x^gamma2/(K2+x^gamma2)/alpha2+1
yss=@(x)1./(beta2*x.^gamma2./(K2+x.^gamma2)./alpha2+1);


%%
ww=0:0.01:1;
   
%%
figure;
subplot(1,3,1)
plot(ww, xss(ww),'LineWidth', 2,  'Color', 'r');
xlabel('y');
ylabel('xss(y): ceroclina de x')
axis square
xlim([0, 1.1]);
ylim([0, 1.1]);

%%
subplot(1,3,2)
plot(ww, yss(ww),'LineWidth', 2,  'Color', 'b');
xlabel('x');
ylabel('yss(x): ceroclina de y')
axis square
xlim([0, 1.1]);
ylim([0, 1.1]);

%%
subplot(1,3,3)
plot(ww, yss(ww),'LineWidth', 2,  'Color', 'b');
hold on
plot(xss(ww), ww, 'LineWidth', 2,  'Color', 'r');
xlabel('x');
ylabel('y');
legend('yss(y): ceroclina de y', 'xss(y): ceroclina de x')
axis square
xlim([0, 1.1]);
ylim([0, 1.1]);
%%

end;


%%

J=Jacobian_matrix_Angeli(0.5, 0.62)
%% evaluar la estabilidad de las raíces

%%%%% ahora vamos a ver la estabiliad de los puntos de equilibrio, que
%%%%% obtuvimos de manera numérica en R

%stable
J=Jacobian_matrix_Angeli(0.9946980, 0.1681565);

% 
J=Jacobian_matrix_Angeli(0.5057743, 0.6195077)






