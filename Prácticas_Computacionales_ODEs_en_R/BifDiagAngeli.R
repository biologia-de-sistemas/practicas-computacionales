getwd()

setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/Posgrado sostenibilidad 2020/Modulo EDOs Elisa Domínguez Hüttinger/Prácticas_Computacionales_ODEs_en_R")

source('Grind.r') # descargado de {http://theory.bio.uu.nl/rdb/grind.html}

# definir la función
model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
    dx = alpha1*(1-x)-beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
    dy = alpha2*(1-y)-beta2*y*x^gamma2/(K2+x^gamma2)
    return(list(c(dx, dy)))
  })
}


# Declarar los valores de parámetros que permanecen constantes
p <-  c(alpha1=1, alpha2=1, beta1=200, beta2=10, gamma1=4, gamma2=4, K1=30, K2=1, v=1)

# "Initial guess" para buscar raíces con el algoritmo de Newton Raphson
s <- c(x=0,y=0)

# Grafica las ceroclinas para este valor de v
par(pty="s") #ejes cuadrados
plane(xmin=0, xmax=1, ymin=0, ymax=1)

# ahora obtengamos estos tres puntos de equilibrio 
mid <- newton(s,plot=F)


low <- newton(c(x=1,y=0),plot=F)

hig <- newton(c(x=0,y=1),plot=F)



#####

par(pty="s")

continue(state=hig, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, add=F) # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)




continue(state=hig, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="x", ymin=0, ymax=1.1, add=F) # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="x", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="x", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
