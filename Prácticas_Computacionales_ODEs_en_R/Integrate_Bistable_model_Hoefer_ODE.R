## Integrate the model of Hoefer

library(deSolve)
# Parámeteros
alpha<-0.02;
k_g<-5;
k<- 1;
    
# Parámetro de bifurcaci?n:
IL4<-1; 

# Ecuación diferencial: 
dGata3dt<- function(t, y, parms) 
list(alpha*IL4+k_g*y^2/(( 1+y)^2)-k*y)

# Primero, un análisis dinámico - integración numérica del sistema, 


times <- seq(from = 0, to = 35, by = 0.1) 
yini<-0.25;

out <- ode(y = yini, times = times, func =dGata3dt, parms = NULL)


yini <- 1; 
out2 <- ode(y = yini, times = times, func = dGata3dt, parms = NULL)

plot(out, out2, main = "Hoefer model:Th2 diff", lwd = 2)

