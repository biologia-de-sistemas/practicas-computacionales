
# ejemplo de EDO

# A --> 0; dA/dt=-A*k1

library(deSolve)
## From: Soetaert, K., Cash, J. & Mazzia, F. Use R! Solving Differential equations in R. (2012). doi:10.1007/978-0-387-78171-6

#simple differential equation is implemented in an R function called derivs that takes as arguments the current time (t), the value of the dependent variable (y) and a parameter vector (parms), and that returns the derivative, as a list.The parameters r and K, although defined outside of function derivs, are also known within the derivative function.1

# (1) darle valor al parámetro
k1=1;

# (2) declarar nuestra función
dA <- function(t, y, parms) {
  list(-y*k1)}

# (3) definir el intervalo de integración, 
times <- seq(from = 0, to = 5, by = 0.2) 

# (4) definir nuestra condición inicial 
IC= 10; 


# (5) integrar numéricamente
At <- ode(y= IC, times = times, func = dA, parms = NULL)


#(6) graficar
plot(At,type="l",col="red")


