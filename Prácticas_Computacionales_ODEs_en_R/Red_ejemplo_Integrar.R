#### Preámbulo
# Borramos todo
rm(list=ls())

# Nos ubicamos donde queremos

# Cargar (antes, si es necesario, instalar con install.packages("deSolve")) la paquetería que necesitamos
library(deSolve)


#######
# (1) Declara la función $\frac{d\bar{y}(t)}{dt}=\bar{F}(\bar{y}(t), \bar{P})$

RedEjemplo <- function(t, y, parms){  
  # parámetros
        #  1        2           3       4         5        6   
  #parms=(Xpre, kprodEprod, kdegEdeg, kdimEdim, kdisEdis, Ytot)
  # Variables de estado - FORMAN UN VECTOR; \vec{y}=[y(1), y(2),...; y(n)]  ; en este sistema n=3
  # x--> y(1) y --> y(2)   xy--> y(3)
  # \bar{P}=[p1, p2, p3, ..., pm]
  X_t=y[1]; 
  Y_t=y[2];
  XY_t=(parms[6]-Y_t); # 
  dX <- parms[1]*parms[2]-X_t*parms[3]-X_t*Y_t*parms[4]+XY_t*parms[5];
  dY <- -X_t*Y_t*parms[4]+XY_t*parms[5];
  list(c(dX,dY))
}

# (2) Declara los valores de parámetro $\bar{P}$
Xpre=10; kprodEprod=.5; kdegEdeg=.5; kdimEdim=10; kdisEdis=1;
Ytot=6;
parms=c(Xpre, kprodEprod, kdegEdeg, kdimEdim, kdisEdis, Ytot);

# (3) Define las condiciones iniciales $\bar{x}(t_0)=\bar{x}_{\rm 0}$
y0=c(0, Ytot);

# (4) Define el intervalo de integración $t=[0:\Delta t:t_{end}]$.
tspan =seq(from = 0, to = 5, by = 0.01)

# (5) Invoca al integrador - para obtener $\bar{x}(t)$.  
# ¡A integrar!
out <- ode(y = y0, times = tspan, func = RedEjemplo, parms = parms)

# Grafiquemos los resultados
par(pty="s")
plot(out[,1], out[,2],type = "l", col="red", xlab = "Time", ylab ="[X(t), Y(t), XY(t)]")
lines(out[,1], out[,3],type = "l", col="blue")
lines(out[,1], Ytot-out[,3],type = "l", col="black")
legend=c("X(t)", "Y(t)", "XY(t)")
