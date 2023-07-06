## parameter estimation in R
# http://tbb.bio.uu.nl/rdb/grindR/
#rm(list=ls())

# (1) Run or Source Grind.R

# (2) Define the function (i.e., the model)
model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
  dX_t=Xpre*kprodEprod-X_t*kdegEdeg-X_t*Y_t* kdimEdim+(Ytot-Y_t)*kdisEdis;
  dY_t=-X_t*Y_t*kdimEdim+(Ytot-Y_t)*kdisEdis;
  dXY_t=-dY_t; #-(-X_t*Y_t*kdimEdim+(Ytot-Y_t)*kdisEdis) #-dY_t;
  return(list(c(dX_t, dY_t, dXY_t)))
  #return(list(c(dX_t, dXY_t)))
  })
}

# (3) Declare first guess of parameter values:  
p=c(Xpre=10, kprodEprod=.5, kdegEdeg=.5, kdimEdim=10, kdisEdis=1, Ytot=6);
# ...and initial conditions (which, effectively, are also parameters)
s <- c(X_t=0,Y_t=5, XY_t=1) # s is the state

# (4) Input the experimental data
# ref. descrip exp...
t_exp= c(0, 1, 3, 6) # hours post-stimulation
XY_exp=c(.5, 2.5, 3, 1.8 )/.5; # Smad1- promoter complex

# make a dataframe 
data <- data.frame("time"= c(0, 1, 3,6) , "XY_t"=c(.5, 2.5, 3 , 1.8)/.5 )
# Names of variables same as in the model


# (5) Run the optimization!
w <- c("kdimEdim", "Ytot") # w provides the names of free parameters
# to be optimized
f=fit(legend=FALSE, free=w, tstep=0.01) #method="Nelder-Mead" # see ?modFit for additional options
summary(f) # Check confidence ranges, etcetera
# we can check the results for example with f$par
f$par #just returns the estimated parameters 
f$ssr # summed squared residuals (error) 

# (6) some options
# assume you have more datasets
#data_1 <- data; # the original data
# invent a new dataset (only for practice purposes of course)
#data_2 <- data; # the original data
#data_2$XY_t <- c(1.4, 4, 5) # empty to fill it in again
#f <- fit(list(data_1,data_2),  free=w, tstep=0.01, differ = "Ytot") # which gives the same result
# differ are the parameters that are allowed to be different

p_opt=p
p_opt[4]=f$par[1]
p_opt[6]=f$par[2]

p_opt
par(pty="s") # axis square
run(parms=p_opt, legend=FALSE, tstep=0.001, tmax=7, main=c(paste("cost=", toString(round(f$ssr, digits = 3) ), sep=" ")))
points(t_exp, XY_exp)