# definir el modelo

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    tlag <- t - Delta
    if (tlag < 0) lags <- 0 # no initial predation
    else lags <- lagvalue(tlag) # returns lags of P
    dP <- alpha*P - beta*P*lags[1]
    return(list(dP))
  })
}


par(pty="s")

p <- c(alpha=2/3,beta=4/3,Delta=3)

s <- c(P=1)

run(delay=TRUE, ylab="Prey", ymax=5, tmax=100, tstep=0.1, legend=FALSE)


p <- c(alpha=2/3,beta=4/3,Delta=0)
s <- c(P=1)
par(new=TRUE)
run(delay=TRUE, ylab="Prey", ymax=5, tmax=100, tstep=0.1, legend=FALSE)
