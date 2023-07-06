## parameter estimation in R

# http://tbb.bio.uu.nl/rdb/grindR/
#rm(list=ls())

# (1) Run or Source Grind.R

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

run()

## Basic Phase plane stuff


##### Integrating the ODE: run()command
run(legend=FALSE) # run the model and make a time plot 

#Solve model nummerically.
# Output= timeplot, trajectory or tabl; + final state.
# Full definition:
#run <- function(tmax=100, tstep=1, state=s, parms=p, odes=model, ymin=0, ymax=NULL,
 #                 log="", x=1, y=2, xlab="Time", ylab="Density", tmin=0, draw=lines, times=NULL,
  #                show=NULL, arrest=NULL, after=NULL, tweak=NULL, timeplot=TRUE, traject=FALSE,
   #               table=FALSE, add=FALSE, legend=TRUE, solution=FALSE, delay=FALSE, lwd=2,
    #              col="black", pch=20, ...)
# uses:    run() calls the ode() function from the deSolve package. 
  
############ PHASE PLANE
plane() # make a phase plane with nullclines
plane(xmin=-0.001,ymin=-0.001) # include the full axis in the phase plane 
plane(tstep=0.5,portrait=T) # make a phase portrait 
plane() # make a clean phase plain again (Fig 1d)
p["K"] <- 0.75 # change the parameter K from 1 to 0.75
plane(add=T) # add the new nullclines
s["R"] <- 0.1 # change the initial state to (R=0.1,N=0.01)


#plane() sets up a space with the First variable on the horizontal axis, and
#the second on the vertical axis. 
#  plane <- function(xmin=0, xmax=1.1, ymin=0, ymax=1.1, log="", npixels=500, state=s,
 #                   parms=p, odes=model, x=1, y=2, time=0, grid=5, eps=0, show=NULL, portrait=FALSE,
  #                  vector=FALSE, add=FALSE, legend=TRUE, zero=TRUE, lwd=2, col="black", pch=20, ...)

# Additional arguments (...) are passed on to run() (for the phase portrait) and to plot(). Note that
#plane() calls the \vectorized" R-function outer(), which implies that if one calls functions in the
#ODEs they should also be vectorized, e.g., one should use pmax() instead of max(). Finally, note that
#there is an extension, cube.R, for 3-dimensional nullclines.


run(traject=T) # run the model and plot a trajectory


#### find stady states
newton(c(R=0.5,N=0.5),plot=T) # find a steady state around (R=0.5,N=0.5)a 
#The function newton() finds a steady state from a nearby initial state, and can report the Jacobi
#matrix with its eigenvalues and eigenvectors. The full definition of newton() is:
 # newton <- function(state=s, parms=p, odes=model, time=0, x=1, y=2, positive=FALSE,
  #                   jacobian=FALSE, vector=FALSE, plot=FALSE, ...)
   # newton() calls the function steady() from the rootSolve package (which calls stode()). Additional
#arguments (...) are passed on to both of them. newton() needs an initial state close to an equilibrium
#point.

f <- newton(c(R=0.5,N=0.5)) # store this steady state in f

#### Continuation
continue(f,x="K",xmax=2,y="N") # continue this steady state while varying K (Fig 1e)
continue(f,x="K",xmax=2,y="N",step=0.001) # get a better value with a smaller step size

#The function continue() continues a steady state by changing a \bifurcation" parameter dened by
#the horizontal axis of the bifurcation diagram. The full de nition of continue() is:
#continue <- function(state=s, parms=p, odes=model, x=1, step=0.01, xmin=0, xmax=1,
#y=2, ymin=0, ymax=1.1, log="", time=0, positive=FALSE, add=FALSE, ...)
#continue() calls the function steady()from the rootSolve package (additional arguments (...) are
#passed on), and needs an initial state close to an equilibrium point. Note that there is much more
#proper software for bifurcation analysis like XPPAUT or MatCont, which reports the type of bifurca-
#tions encountered, and automatically continues all branches of branch points.

p["K"] <- 0.5 # set K to the value at which N goes extinct
plane(vector=T) # make a phase plane for this value of K (Fig 1f)

######### Events location function


p <- c(r=1,K=1,a=1,c=1,delta=0.5) # p is a named vector of parameters
s <- c(R=1,N=0.01) # s is the state

#run(after="a condition") sets a condition after each time step 
#ej
#run(after="state<-ifelse(state<1e-9,0,state)") will set small variables
#to zero after each time step. One can also stop when any of the variables becomes negative by after
#<- "if(min(state)<0)break".

#run(after="state<-ifelse(state<1e-9,0,state)") # will set small variables to zero after each time step. 

#One can also stop when any of the variables becomes negative by after
#run(after= "if(min(state)<0)break")

#This new option is illustrated by the following three examples each
#providing an R-command as a text

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
run(after="if(t==20)state[\"R\"]<-0") # sets the predators N = 0 when
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
    dP <- alpha*P - beta*P*D
    dD <- gamma*P*lags[1] - delta*D
    return(list(c(dP, dD)))
  })
}
p <- c(alpha=2/3,beta=4/3,gamma=1,delta=1,Delta=2)
s <- c(P=1,D=0)
run(delay=TRUE, ylab="Predator / Prey", ymax=5, tmax=50, tstep=0.1)


#One can study delay diferential equations (DDEs) 
#because the deSolve package implements a general
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






model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    R = 1/(1+A^n) # Repressor
    dA = M*L - delta*A - v*M*A # Allolactose
    dM = c0 + c*(1-R) - d*M # mRNA
    return(list(c(dA, dM)))
  })
}
p <- c(L=1,c=1,c0=0.05,d=1,delta=0.2,n=5,v=0.25)
s <- c(A=0,M=0)
plane(xmax=4)
low <- newton(s,plot=T)
mid <- newton(c(A=0.8,M=0.2),plot=T)
hig <- newton(c(A=2,M=1),plot=T)
continue(mid,x="L",y="A",xmax=2,ymax=4)