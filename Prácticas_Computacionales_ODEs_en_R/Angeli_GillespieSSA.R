library(GillespieSSA)


# Constant parameters
v_value=1.8


p <-  c(alpha1=1, alpha2=1, beta1=200, beta2=10, gamma1=4, gamma2=4, K1=30, K2=1, v=v_value)


##### initial state vector (x0),
#The elements of the vector have to be labelled using
#the same notation as the state variables used in the propensity functions

#mid <- newton(c(x=1,y=30),plot=F)
### UTILIZAR GRIND.R (CORRELO UNA VEZ PARA USAR SUS FUNCIONES)
## Y DECLARAR EL MODELO DETERMINISTA:

##### Declarar el modelo Determinista
model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
    dx = alpha1*(100-x)-beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
    dy = alpha2*(100-y)-beta2*y*x^gamma2/(K2+x^gamma2)
    return(list(c(dx, dy)))
  })
}



low <- newton(c(x=100,y=0),plot=F, jacobian=FALSE, vector=FALSE)
#hig <- newton(c(x=0,y=100),plot=F)

xT=100; yT=100;
x0=c(y1=low[[1]],y2=low[[2]],y3=xT-low[[1]],y4=yT-low[[2]]);
x0=round(x0)

# state-change matrix (nu), 

#defines the change in the number of individuals in each state (rows) 
#as caused by one reaction of a given type (columns). 
#For example, the state-change matrix for system with the species S1 and S2
#with two reactions
nu= matrix(c(+1,-1,0,0,
             0,0,+1,-1,
             -1,+1,0,0,
             0,0,-1,+1),nrow=4,byrow=TRUE) 


#The propensity vector, a, 
#defines the probabilities that a particular
#reaction will occur over the next infinitesimal time interval [t; t + dt]. 
#a= c("c1*X1","c2*X2")
#where c1 and c2 are the per capita reaction probabilities.
h=c("alpha1*y3", "beta1*y1*(v*y2)^gamma1/(K1+(v*y2)^gamma1)", "alpha2*y4",  "beta2*y2*y1^gamma2/(K2+y1^gamma2)")

# final time
tf=4


out <- ssa(x0=x0,a=h,nu=nu,parms=p,tf=tf)


plot(out$data[,2])



par(pty="s")
plot(out$data[,1], out$data[,2], type = "l" ,  col="blue", xlab = "Time", ylim=c(0, 100), ylab = "X(t)", main = paste("v=", toString(v_value), sep=" "))

par(pty="s")
plot(out$data[,1], out$data[,3], type = "l" ,  col="blue", xlab = "Time", ylim=c(0, 100), ylab = "Y(t)", main = paste("v=", toString(v_value), sep=" "))


### integra numérciamente la EDO:

out1=run(state=low,table=TRUE, timeplot=FALSE, tstep=0.001)


lines(out1[,1], out1[,3],type = "l", col="red")

# saco la varianza 
var(out$data[,1])
var(out$data[,2])

hist(out$data[,1])

hist(out$data[,2])


