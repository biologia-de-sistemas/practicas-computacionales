
### SIEMPRE COMENTAR SU CODIGO
#### #COMENTARIO



# check in which directory you are:
getwd() 

# Ir a donde queremos estar (ojo con las diagonales / )
setwd("C:/Users/Elisa/Dropbox/Docencia adscripci�n biom�dicas/Bio Mates PCBIOL 2021 1/Pr�cticas computacionales/1 Intro a BoolNet")


# cargar las librer�as
# insall.packages()
library(BoolNet)

# cargar la red
net<- loadNetwork("red_ejemplo.txt")

net  # look at the tables

plotNetworkWiring(net)


# get attractors
attr <- getAttractors(net)

attr  #% look at the attractors
plotAttractors(attr)

# perturb the network....
#Assuming that the c-gene is always off
mut=getAttractors(net, genesOFF=c(0,0,1))
mut

plotAttractors(mut)


getPathToAttractor(net, c(1,0,0)) 


plotStateGraph(attr, piecewise=TRUE)



# nota: actualizaci�n as�ncrona; fuente de ruido. repetir varias veces y ver si resultados
# coinciden con s�ncrono, para descartar que resultados sean artefactos del m�todo computacional

att_asymchron=getAttractors(net, type="asynchronous") # asynchronament

plotAttractors(att_asymchron)

# %ginsim.org/
# % www.colomoto.org/



## Ahora veremos qu� pasa con un atractor c�clico
cell_cyle <- loadNetwork("cellcycle.txt")

cell_cyle


plotNetworkWiring(cell_cyle)


# get attractors
attrCC <- getAttractors(cell_cyle)

attrCC  #% look at the attractors
plotAttractors(attrCC)



  
 