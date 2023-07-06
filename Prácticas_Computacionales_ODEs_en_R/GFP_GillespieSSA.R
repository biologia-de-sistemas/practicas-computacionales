library(GillespieSSA)

  

##### Declarar el modelo Determinista
model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
    dGFP = c1-c2*GFP
    return(list(c(dGFP)))
  })
}

c1=10
c2=.5


# Constant parameters
p <-  c(c1=10, c2=.5)

# queremos ver la distibución estacionaria, por ello empezamos en el punto de equilibrio
x0=c(GFP=c1/c2);

# state-change matrix (nu), 

#defines the change in the number of individuals in each state (rows) 
#as caused by one reaction of a given type (columns). 
#For example, the state-change matrix for system with the species S1 and S2
#with two reactions
nu= matrix(c(+1,-1),nrow=1,byrow=TRUE) 


#The propensity vector, a, 
#defines the probabilities that a particular
#reaction will occur over the next infinitesimal time interval [t; t + dt]. 
#a= c("c1*X1","c2*X2")
#where c1 and c2 are the per capita reaction probabilities.
h=c("c1", "c2*GFP")

# final time
tf=100


out <- ssa(x0=x0,a=h,nu=nu,parms=p,tf=tf)


plot(out$data[,1], out$data[,2])


### integra numérciamente la EDO:
out1=run(state=x0,table=TRUE, timeplot=FALSE, tstep=0.001)


lines(out1[,1], out1[,2],type = "l", col="red")

# saco la varianza 
var(out$data[,2])


hist(out$data[,2])
points(c1/c2, 100, col="red")



