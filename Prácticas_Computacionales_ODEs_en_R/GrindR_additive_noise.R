## Additive noise

# http://tbb.bio.uu.nl/rdb/grindR/
# Run or Source Grind.R

# Define the function (i.e., the model)
# Note: the derivatives should be specified in the same order asthe variables in the state.
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N
    return(list(c(dR, dN)))
  })
}

p <- c(r=1,K=1,a=1,c=1,delta=0.5) # p is a named vector of parameters

s <- c(R=1,N=0.01) # s is the state


##### sample parameter from a random distribution AT EACH TIME STEP###

# r was our alpha parameter #
run(after="parms[\"r\"]<-rnorm(1,1,0.1)", legend=FALSE) # sets the parameter r to a random value,
#drawn from a normal distribution with a mean of one and a standard deviation of 0.1

#Since r is the first parameter, one can also just write
run(after="parms[1]<-rnorm(1,1,0.1)", legend=FALSE) 
# This random resetting of r is done every timestep 
# (as defined by the parameter tstep=1 in run()).


###### ADD GAUSSIAN NOISE TO BOTH VARIABLES ##############
### ADDITIVE NOISE ######## 

# START FROM THE EQUILIBRIUM STATE
f <- newton(c(R=0.5,N=0.5))
X=run(state=f,after="state<-state+rnorm(2,mean=0,sd=0.01)", legend=FALSE) #,ymax=1)
#Note that rnorm(2,0,0.01) provides two
#random values, that are added to the two variables, respectively.


########### ahora, repitamos el experimento un par de veces; vamos a recolectar sólo el valor final

# define el número de iteraciones
iterations=500

# prealoca la matrix para guardar los resulados de las simulaciones
x_t_family=matrix(data=NA, nrow=iterations, ncol = 2)

for (i in 1:iterations){
  
  X=run(state=f,after="state<-state+rnorm(2,mean=0,sd=0.01)", legend=FALSE)
  x_t_family[i, ]=X
}# end for iteraciones


lines(x_t_family[,1], x_t_family[,2], col="red", type="l", lwd=5, ylab="cocos", xlab="time" )



# Graficar la distribución de valores finales
par(pty="s")
hist(x_t_family[,1], freq=TRUE, col="blue", main=" ", xlab="Presa steady state")
points(f[1], 10, col="red", lwd=5)

par(pty="s")
hist(x_t_family[,2], probability = TRUE, col="red", main=" ", xlab="Depredador steady state")
points(f[2], 10, col="blue", lwd=5)




