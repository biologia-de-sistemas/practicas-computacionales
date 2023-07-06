library(BoolNet)
setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/EXAMENES Y RESPUESTAS A PRÁCTICAS PCBIOL 2021 1/ensamblaje de comunidades")

net<-loadNetwork("Red_planta_polinizador.txt")
attr<-getAttractors(net)
plotAttractors(attr)
plotStateGraph(attr) 
path <- getPathToAttractor(net, c(0,1,1,1,1))
plotSequence(sequence=path)
