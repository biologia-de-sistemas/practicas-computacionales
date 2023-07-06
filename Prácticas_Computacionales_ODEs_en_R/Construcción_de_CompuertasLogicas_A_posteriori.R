# (1) Run or Source Grind.R
library(deSolve)     # run() calls the ode() function
library(rootSolve)   # newton() and continue() call steady()
library(FME)         # fit() calls modFit() and modCost()

# Define the function (i.e., the model)
# Note: the derivatives should be specified in the same order as the variables in the state.

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    A1<-(A1T-(A1A2+A1D+A1A2D))
    A2<-(A2T-(A1A2+A1D+A1A2D))#Algebraicas
    D <-(DT-(A2D+A1D+A1A2D))
    dA1D<- -(K4dim*A1D*A2)+(K4dis*A1A2D) #Dinamicas
    dA2D<- -(K5dim*A2D*A1)+(K5dis*A1A2D)
    dA1A2<- (K3dim*A1*A2)+ (K6dis*A1A2D)-(K3dis*A1A2)-(K6dim*A1A2*D)
    dA1A2D<- (K5dim*A2D*A1)-(K5dis*A1A2D)+(K4dim*A1D*A2)-(K4dis*A1A2D)+(K6dim*D*A1A2)-(K6dis*A1A2D)
    return(list(c(dA1D,dA2D,dA1A2,dA1A2D)))
  })
}
#

p <- c(A1T=1,A2T=1,DT=1,K3dim=1,K3dis=1,K4dim=1,K4dis=1,K5dim=1,K5dis=1,K6dim=1,K6dis=1) # p is a named vector of parameters
s <- c(A1D=0.5,A2D=0.5,A1A2=0.5,A1A2D=0.5) # s is the state
valor_final_t=run(tmax=10, tstep=0.01, state=s, parms=p, odes=model, xlab="tiempo", ylab="Heteromeros", add=FALSE, legend=FALSE)
#Algoritmo de newton para puntos estables 
ss=newton(valor_final_t)

CRIF<- function(A1c, A2c){
  alpha=1
  beta=1
  gamma=1
  p["A1T"]<-A2c
  p["A2T"]<-A2c
  temp<-run(tmax=10, tstep=0.01, state=s, parms=p, odes=model, xlab="tiempo", ylab="Heteromeros", add=FALSE, legend=FALSE)
  ss<-newton(temp)
  criff<-alpha*ss[[1]]+beta*ss[[2]]+gamma*ss[[3]]
  f_est<-c("F",criff)
  return(criff)
}

#A1T y a2T SIEMPRE DEBEN SER MAYORES A LOS COMPLEJOA QUE FORMAN
A1T=0.1
A2T=0.1
p["A1T"]<-0.1
p["A2T"]<-0.1
s <- c(A1D=0,A2D=0,A1A2=0,A1A2D=0) # 
valor_final_t=run(tmax=10, tstep=0.01, state=s, parms=p, odes=model, xlab="tiempo", ylab="Heteromeros", add=FALSE, legend=FALSE)
sscc<-newton(valor_final_t)
#Agregar valores de alpha beta y gamma
alpha<-1
beta<-1
gamma<-1

F1=alpha*sscc[[1]]+beta*sscc[[2]]+gamma*sscc[[3]]
#falta cambiar los otros valores de A1T Y A2T
#Así podemos graficar la sabana
A1T<-0.9
A2T<-0.9
p["A1T"]<-0.9
p["A2T"]<-0.9
s <- c(A1D=0,A2D=0,A1A2=0,A1A2D=0) # 
valor_final_g=run(tmax=10, tstep=0.01, state=s, parms=p, odes=model, xlab="tiempo", ylab="Heteromeros", add=FALSE, legend=FALSE)
ssgg<-newton(valor_final_g)
F2=alpha*ssgg[[1]]+beta*ssgg[[2]]+gamma*ssgg[[3]]
#Ahorauno grande y uno chico
A1T<-0.9
A2T<-0.1
p["A1T"]<-0.9
p["A2T"]<-0.1
s <- c(A1D=0,A2D=0,A1A2=0,A1A2D=0) # 
valor_final_gc=run(tmax=10, tstep=0.01, state=s, parms=p, odes=model, xlab="tiempo", ylab="Heteromeros", add=FALSE, legend=FALSE)
ssgc<-newton(valor_final_gc)
F3=alpha*ssgc[[1]]+beta*ssgc[[2]]+gamma*ssgc[[3]]
#Ahora chico contra grande
A1T<-0.1
A2T<-0.9
p["A1T"]<-0.1
p["A2T"]<-0.9
s <- c(A1D=0,A2D=0,A1A2=0,A1A2D=0) # 
valor_final_cg=run(tmax=10, tstep=0.01, state=s, parms=p, odes=model, xlab="tiempo", ylab="Heteromeros", add=FALSE, legend=FALSE)
sscg<-newton(valor_final_cg)
F4=alpha*sscg[[1]]+beta*sscg[[2]]+gamma*sscg[[3]]
#POner uno mediano contra mediano
p_eq<-matrix(nrow=2,ncol = 2,c(F1,F5,F2,F3,F6,F4))
p_eq<-matrix(nrow=2,ncol = 2,c(F1,F4,F3,F2))

valor_final_t
valor_final_g
valor_final_cg
valor_final_gc

#Llamamos a la paqueteria para poder imprimir con la matriz
library(plot3D)

vector<-c(0.1,0.9)
image2D(p_eq,vector,vector,border="black")
persp3D(x = vector, y = vector, z = p_eq, xlab = "S3", ylab = "S2", zlab = "S1_crit", xlim = c(0,0.8), ylim = c(0,0.8))#, zlim = c(0,2.5))


#persp3D(z = z, x = x, y = y)
#Probar la función
x<-CRIF(0.1,0.2)
x[1]
vector<-seq(0,2.8,0.3)
#Matriz para guardar los valores
stable_points<-matrix(ncol = 10,nrow = 10)
for (i in 1:10){
  #print(i)
  for (z in 1:10){
    #print(z)
    x<-CRIF(vector[i],vector[z]) 
    print(x)
    stable_points[i,z]<-x
  }
}
image2D(stable_points, vector, vector, border = "black")
persp3D(x = vector, y = vector, z = stable_points, xlab = "S3", ylab = "S2", zlab = "S1_crit", xlim = c(0,3), ylim = c(0,3))#, zlim = c(0,2.5))
