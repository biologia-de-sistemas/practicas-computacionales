## parameter estimation in R
# http://tbb.bio.uu.nl/rdb/grindR/
#rm(list=ls())

# (1) Run or Source Grind.R

# (2) Define the function (i.e., the model)
######### PONER ACÁ LAS ECUACIONES
model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
   dy1=k1*S-k2*y2*y1;
   dy2=k3*y1-k4*y2; 
return(list(c(dy1, dy2)))
  })
}

# (3) Declare first guess of parameter values: 
# PARAMETROS RAZONABLES
p=c(k1=5, k2=13, k3=1, k4=.3, S=2.5);
# ...and initial conditions (which, effectively, are also parameters)

# CONDICIONES INICIALES QUE SEAN CONGRUENTES CON LOS DATOS EXPERIMENTALES
# (SI ES POSIBLE, SE LEEN DE LOS DATOS EXPERIMENTALES)
s <- c(y1=0,y2=0) # s is the state

# (4) Input the experimental data
#t_exp=[   0    0.5000    1.0000    2.0000    4.0000    6.0000 ];%  24.0000]; % Time of measurement, in hours
#y_exp=[0.0126    1.6059    0.8196    0.7323    0.5339    0.5613 ];%   0.1898]; % pStatStat3

# make a dataframe 
data <- data.frame("time"= c( 0, 0.5, 1, 2, 4,6) , "y1"=c(0.0126,1.6059,0.8196, 0.7323, 0.5339,0.5613) )
# Names of variables same as in the model

# (5) Run the optimization!
w <- c("k1", "k2", "k3", "k4") # w provides the names of free parameters
# to be optimized
f=fit(legend=FALSE, free=w, tstep=0.0001,method="BFGS") # see ?modFit for additional options
summary(f) # Check confidence ranges, etcetera
# we can check the results for example with f$par
f$par #just returns the estimated parameters 
f$ssr # summed squared residuals (error) 



## LOS VALORES DE PARAMETROS SON:
f$par #just returns the estimated parameters 