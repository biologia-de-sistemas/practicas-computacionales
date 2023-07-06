rm(list=ls())


setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/Bio Mates PCBIOL 2021 1/Prácticas computacionales/Prácticas_Computacionales_ODEs_en_R")

library(deSolve)
library(sensitivity)
library(checkmate)

library(ODEnetwork)
library(ODEsensitivity)

##### Lotka-Volterra equations #####
# The model function:
LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ingestion    <- rIng  * Prey * Predator
    GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
    MortPredator <- rMort * Predator
    
    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion * assEff - MortPredator
    
    return(list(c(dPrey, dPredator)))
  })
}

#The k = 5 parameters rG, rI, rM, kAE and K are considered as input variables for the sensitivity
#analysis. Hence, we will analyze the sensitivity of the prey's and the predator's population with
#regard to changes in these 5 parameters.


# The parameters to be included in the sensitivity analysis and their lower 
# and upper boundaries:
LVpars  <- c("rIng", "rGrow", "rMort", "assEff", "K")
# normalmente alrededor del valor nominal (output de la optimizacion)
LVbinf <- c(0.05, 0.05, 0.05, 0.05, 1)
LVbsup <- c(1.00, 3.00, 0.95, 0.95, 20)

# The initial values of the state variables:
LVinit  <- c(Prey = 1, Predator = 2)
# The timepoints of interest:
LVtimes <- c(0.01, seq(1, 5, by = 0.1))
set.seed(59281)
# Sobol' sensitivity analysis (here only with n = 500, but n = 1000 is
# recommended):
# Warning: The following code might take very long!

LVres_sobol <- ODEsobol(mod = LVmod,
                        pars = LVpars,
                        state_init = LVinit,
                        times = LVtimes,
                        n = 500,
                        rfuncs = "runif",
                        rargs = paste0("min = ", LVbinf,
                                       ", max = ", LVbsup),
                        sobol_method = "Martinez",
                        ode_method = "lsoda",
                        parallel_eval = TRUE,
                        parallel_eval_ncores = 2)

str(LVres_sobol, vec.len = 3, give.attr = FALSE)

# it is a list of class "ODEmorris" with one element for each state
#variable (here, Prey and Predator). Those elements are matrices of 3*length(LVpars) + 1  rows and length(LVtimes) columns. 
#The first row contains a copy of all timepoints. 
#The other rows contain the 3 Sobol sensitivity indices 
#  for all parameters at all 51 timepoints.

plot(LVres_sobol,  state_plot = "Predator")
plot(LVres_sobol,  state_plot = "Prey")


# now let's plot the distribution of indexes per time
library(vioplot)

x1 <- LVres_sobol$Prey$S[2, ]
x2 <- LVres_sobol$Prey$S[3, ]
x3 <- LVres_sobol$Prey$S[4, ]
x4 <- LVres_sobol$Prey$S[5, ]
x5 <- LVres_sobol$Prey$S[6, ]

par(pty="s") 
vioplot(x1, x2, x3,x4, x5, names=LVpars,  col="black")
title("Prey sensitivity")
LVres_sobol$Prey$S

y1 <- LVres_sobol$Predator$S[2, ]
y2 <- LVres_sobol$Predator$S[3, ]
y3 <- LVres_sobol$Predator$S[4, ]
y4 <- LVres_sobol$Predator$S[5, ]
y5 <- LVres_sobol$Predator$S[6, ]

par(pty="s") 
vioplot(y1, y2, y3,y4, y5, names=LVpars,  col="red")
title("Predator sensitivity")

