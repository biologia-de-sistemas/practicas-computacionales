

IC=1; 


c1=100; 
c2=5;


th=c(c1,c2); 
n=50; # Parameters / settings


h=function(y, th){return(c(th[1], th[2]*y))} # Hazard functions
# Gillespie algorithm
gillespieGFp <- function(IC, th, h, n)
  
{v=length(th) # number of reactions (IN THIS CASE: same as number of kinetic parameters)
tt=0; x=IC; tvec=vector("numeric",n); xvec=vector("numeric", n) # matrix for the u variables of the system
  for (i in 1:n)    {
    totH= sum(h(x, th))# total hazard
    tt=tt+rexp(1,totH) # sample time to next event from exp probability distribution
    j=sample(v,1, prob=h(x, th)) # choose a reaction to occur
    # update state vector
    if ( j==1) {x=x+1} # production
    else # j==2
      x=x-1 # degradation
        # done, just fill in the vectors:
    tvec[i]=tt; xvec[i]=x  }
return(list(t=tvec, x=xvec))
} # end function
# one realization
out=  gillespieGFp(IC, th, h, n)
plot(out$t, out$x, type="l",ylim=range(c(out$t,out$x)),  main=" ", xlab="time ", ylab="GFP" )
# second plot  EDIT: needs to have same ylim
iterations=100
# Now we will do several iterations to generate the distributions
# For this we have to disretize the output and put on a regular grit:
discretize <- function(out)
{
  # out: el output de una realización del algoritmo de gillespie, i.e. tu GFP(t)
  events=length(out$t) #número de eventos - igual para todas las realizaciones del algoritmo
    start=0
  end=out$t[events] # tiempo final
  dt=0.01 # choose this delta t wisely
  len=(end-start)%/%dt # gradilla regular
  x=vector("numeric", len) # prealocar vector x
  target=0;   t=0;   j=1
  for(i in 1:events){
    while (out$t[i]>=target){
      x[j]=out$x[i] # fill in the vector until "the time comes"
      j=j+1
      target=target+dt
      t[j]=target
      } # end while
     #else {x[j]=NaN}
    } # end for
  #xdisc=ts(x, start=0, deltat=dt) #create a time-series object with the output
  return(list(tdisc=t, xdisc=x))
  } # end function


xENDvec=vector("numeric", iterations)
numPoints=50

regularOutMatrix=matrix(0, nrow=numPoints, ncol=iterations)
for (ii in 1:iterations)
{
par(new = TRUE)
out2=  gillespieGFp(IC, th, h, n)

regularOutMatrix[,ii]=discretize(out2)$xdisc[1:numPoints]
plot(out2$t, out2$x,type="l", ylim=range(c(out$t,out$x)), axes = FALSE, xlab = "", ylab = "")
xENDvec[ii]=out2$x[n]
}
par(new = TRUE)
plot(discretize(out2)$tdisc[1:numPoints], rowMeans(regularOutMatrix,1),type="l",col="blue", ylim=range(c(out$t,out$x)), axes = FALSE, xlab = "", ylab = "")


library(deSolve)
## From: Soetaert, K., Cash, J. & Mazzia, F. Use R! Solving Differential equations in R. (2012). doi:10.1007/978-0-387-78171-6

#simple differential equation is implemented in an R function called derivs that takes as arguments the current time (t), the value of the dependent variable (y) and a parameter vector (parms), and that returns the derivative, as a list.The parameters r and K, although defined outside of function derivs, are also known within the derivative function.1

detGFP <- function(t, y, parms) {
  list(c1-c2*y)}


times <- discretize(out2)$tdisc[1:numPoints]#seq(from = 0, to = out$t[n], by = 0.2) 
outD <- ode(y = IC, times = times, func = detGFP, parms = NULL)
par(new = TRUE)
#par(new = FALSE)
plot(outD,type="l",col="red", ylim=range(c(out$t,out$x)), axes = FALSE, xlab = "", ylab = "", main= "")

par(new = FALSE)
hist(xENDvec,prob = TRUE)
par(new = TRUE)
lines(density(xENDvec), col="red")
abline(v=c1/c2, col="blue")


par(new = FALSE)
plot(discretize(out2)$tdisc[1:numPoints], apply(regularOutMatrix[1:numPoints,], 1, sd),type="l",col="blue", ylim=range(c(out$t,out$x)), xlab = "time", ylab = "standard deviation")

