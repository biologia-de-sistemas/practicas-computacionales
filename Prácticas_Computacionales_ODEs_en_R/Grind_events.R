# http://tbb.bio.uu.nl/rdb/grindR/
rm(list=ls())



# (1) Run or Source Grind.R
# Define the function (i.e., the model)
# Note: the derivatives should be specified in the same order asthe variables in the state.
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dP <- alpha*P - beta*P*D
    dD <- gamma*P*D - delta*D
    return(list(c(dP, dD)))
  })
}
p <- c(alpha=2/3,beta=4/3,gamma=1,delta=1) # p is a named vector of parameters
s <- c(P=1,D=1) # s is the state

#run(after="a condition") sets a condition after each time step 
run(xlab="Time", ylab="Predator / Prey", ymax=2, tmax=25, tstep=0.01)
run(after="if(t==20)state[\"D\"]<-0", legend=FALSE, tstep=0.01, xlab="Time", ylab="Predator / Prey", ymax=5, tmax=25) # sets the predators N = 0 when
#time t = 20.  Note again that s is called state in run() (see the Manual), and the
#backslashes in \"N\". Again, the more simple state[2]<-0 would achieve the same eect, because
#N is the second variable.

#For example, 
run(after="parms[\"r\"]<-rnorm(1,1,0.1)", legend=FALSE) # sets the parameter r to a random value,
#drawn from a normal distribution with a mean of one and a standard deviation of 0.1
#Note that p is called parms within run() (see the Manual), and the back-
#  slashes in \"r\" before the quotes around the parameter name, because these quotes would otherwise
#mark the beginning or ending of a text. 
#(Since r is the first parameter, one can also just write

run(after="parms[1]<-rnorm(1,1,0.1)") # to achieve the same effect). This random resetting of r is done every
#timestep (as defined by the parameter tstep=1 in run()).

#The second example, 
run(after="if(t==20)state[\"N\"]<-0") # sets the predators N = 0 when
#time t = 20.  Note again that s is called state in run() (see the Manual), and the
#backslashes in \"N\". Again, the more simple state[2]<-0 would achieve the same eect, because
#N is the second variable.

# Use arrest to handle events at time points within time steps:
run(50,arrest=33.14,after="if(t==33.14)state[\"N\"]<-0",table=T)

#Finally, the third example adds Gaussian noise to both variables.
# START FROM THE EQUILIBRIUM STATE
f <- newton(c(R=0.5,N=0.5))
run(state=f,after="state<-state+rnorm(2,mean=0,sd=0.01)") #,ymax=1)
#Note that rnorm(2,0,0.01) provides two
#random values, that are added to the two variables, respectively.


# Here is an example of similar event handling in deSolve:
fun<-function(t, y, parms){y["N"]<-0;return(y)}
run(events=list(func=fun,time=20))



##ADDITional \event" to tweak the numerical data. To modify the numerical solution data com-
#puted within run() one can pass on any R-command as a text with the tweak option. 
#For instance,
run(tweak="nsol$T<-nsol[,2]+nsol[,3]",timeplot=TRUE )
#adds a fourth column to the solution by summing the
#first and second variable, and calls this column \T" (for total).
#Note that the numerical solution is called nsol, and that the first column 
#contains the time. This manipulated nsol table is subsequently
#passed on to timePlot, or printed to screen (with the table=TRUE option in run()). 


#This tweak option can be very helpful when fitting data containing columns 
# representing (transformed) combina  tions of the variables of the model. 
#Finally, note that this same tweak can be done by using apply to
#call sum: 
run(tweak="nsol$T<-apply(nsol[,2:3],1,sum)") #, which shows how one can generalize this to sum over a subset of columns.

############ MAPS (DISCRETE DYNAMICAL SYSTEM)

#One can study maps by simply switching to Euler integration:
  model <- function(t, state, parms) {
    with(as.list(c(state,parms)), {
      dN <- r*N*(1 - N) - N
      return(list(c(dN)))
    })
  }
p <- c(r=3.75)
s <- c(N=0.01)
data <- run(1000,method="euler",table=TRUE)
plot(data$N[1:999],data$N[2:1000],pch=".")

#Be careful with using the steady state functions (newton(), continue(), and plane()) because they
#are not aware of the fact that this model defines a map. Steady state values of this model will be
#correct because of the -N in the equation (i.e., dN=0 has the same solution as Nt+1 = rNt(1 ???? Nt)),
#but the reported eigenvalues no longer deFIne the stability of steady states.

################### delay differential equations
#One can study delay dierential equations (DDEs) because the deSolve package implements a general
#solver (dede()) using the same syntax as the general ODE solver (ode()). For instance, the Lotka
#Volterra model with delayed growth of the predator,


model <- function(t, state, parms) {
with(as.list(c(state,parms)), {
  tlag <- t - Delta
  if (tlag < 0) lags <- c(0,0) # no initial predation
  else lags <- lagvalue(tlag) # returns lags of R and N
  dR <- r*R*(1 - R/K) - a*R*N
  dN <- a*lags[1]*lags[2] - d*N
  return(list(c(dR, dN)))
})
}
p <- c(r=1,K=1,a=1,c=1,d=0.5,Delta=10)
s <- c(R=1,N=0.1)
run(delay=TRUE)

#One can study delay diferential equations (DDEs) because the deSolve package implements a general
#solver (dede()) using the same syntax as the general ODE solver (ode()). For instance, the Lotka
#Volterra model with delayed growth of the predator.


################Vectors of equations
#Here is an example of a model with n = 3 prey populations, Ri, that are competing with each other
#via a logistic term (see the left panel in Fig. 4). Each prey has its own predator Ni

model <- function(t, state, parms){
  with(as.list(c(state,parms)),{
    R <- state[1:n]
    N <- state[(n+1):(2*n)]
    S <- sum(R)
    dR <- b*R*(1-S) - d1*R - a*R*N
    dN <- a*R*N - d2*N
    return(list(c(dR,dN)))
  })
}
n <- 3 # number of species
b <- rnorm(n,mean=1,sd=0.1) # b is a global parameter
p <- c(d1=0.1,d2=0.2,a=1) # other parameters
R <- rep(0.1/n,n) # initial condition of R
names(R) <- paste("R",seq(1,n),sep="") # Name them R1, R2, ... Rn
N <- rep(0.01/n,n) # initial condition of N
names(N) <- paste("N",seq(1,n),sep="") # Name them N1, N2, ... Nn
s <- c(R,N) # combine R and N into s
run()

###### FITTING PARAMETERS

# the ode model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N
    return(list(c(dR, dN)))
  })
}

# nominal parameters
p=c(r=1,K=1,a=1,c=1,delta=0.5)

# initial conditions
s=c(R=1,N=0.01)

data <- run(20,table=T) # Make a FAKE data set  to fit the model to
s <- s*abs(rnorm(2,1,0.1));s # Random guess for initial condition
p <- p*abs(rnorm(5,1,0.1));p # Random guess for parameters
f <- fit() # Fit data with all 7 parameters free
summary(f) # Check confidence ranges, etcetera
# we can check the results for example with f$par


p <- f$par[3:7];p # Store estimates in p

p <- p*abs(rnorm(5,1,0.1));p # Another random guess for the parameters
w <- c(names(s),names(p)) # w provides the names of free parameters
f <- fit(data,free=w) # Fit the data again
f <- fit(data,initial=T,free=names(p)) # Take initial condition from data

dataR <- data; dataR$N <- NULL # Make two data sets one with R,
dataN <- data; dataN$R <- NULL # and the other with N,

f <- fit(list(dataR,dataN)) # which gives the same result
p <- c(r=1,K=1,a=1,c=1,delta=0.5) # Start again with same parameters
p["K"] <- 0.75 # Change K,
s <- c(R=0.5,N=0.05) # and the initial condition,
data2 <- run(25,table=T) # and make a new data set
s <- c(R=0.75,N=0.02) # An "average" guess for the 2 initial conditions
par(mfrow=c(1,2)) # Show two panels next to each other
f <- fit(list(data,data2),differ=c("R","N","K"),main=c("A","B"))
f$par # Show parameters only
# Provide individual initial guesses as a list:
differ <- list(R=c(0.9,0.55),N=c(0.02,0.04),K=c(1.1,0.7))
f <- fit(list(data,data2),differ=differ,main=c("A","B"))
# Provide fixed parameters as a list:
fixed <- list(R=c(1,0.5),N=c(0.01,0.05))
differ <- "K" # one unknown parameter (K)
free <- names(p)[-2];free # and four shared unknown parameters
f <- fit(list(data,data2),free=free,differ=differ,fixed=fixed,main=c("A","B"))
# The latter is identical to taking the initial condition from the data:
f <- fit(list(data,data2),free=free,differ=differ,initial=T,main=c("A","B"))
# Finally perform a 100 bootstrap simulations:
fit(list(data,data2),free=free,differ=differ,fixed=fixed,main=c("A","B"),bootstrap=100)$par
par(mfrow=c(1,1))