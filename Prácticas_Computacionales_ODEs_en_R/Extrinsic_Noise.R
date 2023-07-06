library(deSolve)

#  Declara la función 
cocos <- function(t, X, parms){
  kP=parms[1]; kD=parms[2];
  dX= kP-X*kD
  list(dX)
}

# declara los valores nominales (ej: output de optimización paramétrica)
kP_nominal=5;
kD_nominal=7;

#.. y las condiciones iniciales
x0_nominal=1;

#  Define el intervalo de integración
tspan =seq(from = 0, to = 10, by = 0.01)

parms=c(kP_nominal, kD_nominal)
x_n= ode(y = x0_nominal, times = tspan, func = cocos, parms=parms)

par(pty="s")
plot(x_n[,1], x_n[,2], col="red", type="l", lwd=5) 

# define el número de iteraciones
iterations=100

# prealoca la matrix para guardar los resulados de las simulaciones
x_t_family=matrix(data=NA, nrow=length(tspan), ncol = iterations)

for (i in 1:iterations){
  # muestrea parámetros:
  kP=rnorm(1, mean=kP_nominal, sd=1)
  kD=rnorm(1, mean=kD_nominal, sd=1)
  x0=rnorm(1, mean=x0_nominal, sd=1)
  
  parms=c(kP, kD)
   # Invoca al integrador 
  x= ode(y = x0, times = tspan, func = cocos, parms = parms)
  x_t_family[,i]=x[,2]
  
  lines(x[,1], x[,2], col="grey", type="l", lwd=1) 
  
}# end for iteraciones

lines(x_n[,1], x_n[,2], col="red", type="l", lwd=5, ylab="cocos", xlab="time" )



# Graficar la distribución de valores finales
par(pty="s")
hist(x_t_family[1001,], freq=FALSE, col="blue", main=" ", xlab="x(t=10)")
points(kP_nominal/kD_nominal, 1, col="red")




