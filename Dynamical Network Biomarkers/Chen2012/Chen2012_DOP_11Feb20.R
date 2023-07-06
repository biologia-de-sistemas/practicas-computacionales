########################Dinamical network biomarkers##########################
#---------------------------Chen et al model-------------------------------
rm(list = ls()); graphics.off()#clear all; close all

#First, we write the deterministic model 
model <- function(t,state,parms){
  with(as.list(c(state,parms)),{
    dA = (90*abs(P)-1236) + (240-120*abs(P))/(1+C) + 1488*D/(1+D) - 30*abs(P)*A
    dB = (75*abs(P)-150) + (60-30*abs(P))/(1+A) + ((240-120*abs(P))*C)/(1+C) - 60*B
    dC = -1056 + (1488*D)/(1+D) - 60*C
    dD = -600 + (1350*E)/(1+E) - 100*D
    dE = 108 + 160/(1+A) + 40/(1+B) + 1488/(1+D) -300*E
    return(list(c(dA,dB,dC,dD,dE)))
  })
} #The variables are named as follows: A=z1, B=z2, C=z3, D=z4, E=z5.

#Set the bifurcation parameter
p <- c(P=0.4) 
#Look for the steady state founded in the paper 
low <- newton(state=c(A=0,B=0,C=0,D=2,E=1), parms = p, odes = model, plot=F)
  #SS=(1,0,1,3,2) is the steady state indeed

#----------------Deterministic model using GrindR

s<-c(A=0,B=0,C=0,D=3,E=2) #Initial condition
run() #Deterministic trajectories
 
#----------------Stochastic model using GrindR

s<-c(A=1,B=0,C=1,D=3,E=2) #initial condition in the steady state

p <- c(P=0.4) # Set P in the "normal" region
out<-run(tmax=1000, tstep=0.1, after="state<-state+rnorm(5,mean=0,sd=0.1)",timeplot=FALSE,table = T)
    #we use the event location function to create the stochastic version

#Plotting trajectories
plot(out$time,out$A,type='l',col=1,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "",ylab = "")
par(new=T)
plot(out$time,out$B+1,type='l',col=2,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "",ylab = "")
par(new=T)
plot(out$time,out$C,type='l',col=3,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "",ylab = "")
par(new=T)
plot(out$time,out$D/3,type='l',col=4,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "",ylab = "")
par(new=T)
plot(out$time,out$E/2,type='l',col=5,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "Time",ylab = "Z",
     main = "P=0.4")
abline(h = 1, col = "black",lwd=3)
legend(915,3.8, legend=c("Z1", "Z2", "Z3", "Z4", "Z5"),
       col=c(1,2,3,4,5), lty=1,cex=0.65)

#Now we are near the bifurcation
p <- c(P=0.01) # Set P in the "pre-critical" state
s<-c(A=1,B=0,C=1,D=3,E=2)
outNearBif<-run(tmax=1000, tstep=0.1, after="state<-state+rnorm(5,mean=0,sd=0.1)",timeplot=FALSE,table = T)
plot(outNearBif$time,outNearBif$A,type='l',col=1,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "",ylab = "")
par(new=T)
plot(outNearBif$time,outNearBif$B+1,type='l',col=2,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "",ylab = "")
par(new=T)
plot(outNearBif$time,outNearBif$C,type='l',col=3,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "",ylab = "")
par(new=T)
plot(outNearBif$time,outNearBif$D/3,type='l',col=4,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "",ylab = "")
par(new=T)
plot(outNearBif$time,outNearBif$E/2,type='l',col=5,ylim = c(min(outNearBif$A),max(outNearBif$A)),xlab = "Time",ylab = "Z",
     main = "P=0.01")
abline(h = 1, col = "black",lwd=3)
legend(915,3.8, legend=c("Z1", "Z2", "Z3", "Z4", "Z5"),
       col=c(1,2,3,4,5), lty=1,cex=0.5)

#------------------------Early warning signals

#NOTE: we start at the same initial condition as before
Parameters<-seq(0.4,0,by=-0.01) # Parameters
variables=5
auxSD<-matrix(0,length(Parameters),variables+1) # To store the standard deviation 
auxCOR<-matrix(0,length(Parameters),10)# To store the correlation
                                        # 10 = nCm (Binomial Coefficient)
auxMarkers<-matrix(0,length(Parameters),4)
#auxComIndex <- matrix(0,length(Parameters),1) # To store the composite index

cont=1;
for (i in 1:length(Parameters)) {
  p <- c(P=Parameters[i]) 
  outSD<-run(tmax=5000, tstep=0.1, after="state<-state+rnorm(5,mean=0,sd=0.1)",timeplot=FALSE,table = T)
    
    #Computing standard deviation
    dev<-apply(outSD,2,sd) # 2 to calculate sd over the colums 
    for (j in 1:variables+1) {
      auxSD[i,j]<-dev[j]
    }
  
  #Computing Pearson correlation
  corAB<-cor(outSD$A,outSD$B); corCD<-cor(outSD$C,outSD$D)
  corAC<-cor(outSD$A,outSD$C); corCE<-cor(outSD$C,outSD$E)
  corAD<-cor(outSD$A,outSD$D); corDE<-cor(outSD$D,outSD$E)
  corAE<-cor(outSD$A,outSD$E); corBD<-cor(outSD$B,outSD$D)
  corBC<-cor(outSD$B,outSD$C); corBE<-cor(outSD$B,outSD$E)
  
  auxCOR[i,cont]<-corAB
  auxCOR[i,cont+1]<-corAC
  auxCOR[i,cont+2]<-corAD
  auxCOR[i,cont+3]<-corAE
  auxCOR[i,cont+4]<-corBC
  auxCOR[i,cont+5]<-corBD
  auxCOR[i,cont+6]<-corBE
  auxCOR[i,cont+7]<-corCE
  auxCOR[i,cont+8]<-corCD
  auxCOR[i,cont+9]<-corDE
  
  #Computing composite index
  SDd <- mean(c(auxSD[i,2],auxSD[i,3])) #This should increase
  PCCd <- abs(mean(c(corAB))) #This should increase
  PCCo <- abs(mean(c(corAC,corAD,corAE,corBC,corBD,corBE))) #This should decrease
  CompIndex <- (SDd*PCCd)/(PCCo)
  
  auxMarkers[i,1] <- SDd
  auxMarkers[i,2] <- PCCd
  auxMarkers[i,3] <- PCCo
  auxMarkers[i,4] <- CompIndex
}

#Store the info collected
write.csv(auxSD,"SD_Chen2012_12Feb20.csv")
write.csv(auxCOR,"PCC_Chen2012_12Feb20.csv")
write.csv(auxMarkers,"CompIndex_Chen2012__12Feb20.csv")

SD<-read.csv("SD_Chen2012_12Feb20.csv",sep = ",")
PCC<-read.csv("PCC_Chen2012_12Feb20.csv",sep = ",")
CompIndex<-read.csv("CompIndex_Chen2012__12Feb20.csv",sep = ",")

plot(Parameters,SD$V2,type="o",pch=16,col=1,ylim = c(min(SD$V2[1:40]),max(SD$V2[1:40])),xlab = "",ylab = "")
par(new=T)
plot(Parameters,SD$V3,type="o",pch=16,col=2,ylim = c(min(SD$V2[1:40]),max(SD$V2[1:40])),xlab = "",ylab = "")
par(new=T)
plot(Parameters,SD$V4,type="o",pch=16,col=3,ylim = c(min(SD$V2[1:40]),max(SD$V2[1:40])),xlab = "",ylab = "")
par(new=T)
plot(Parameters,SD$V5,type="o",pch=16,col=4,ylim = c(min(SD$V2[1:40]),max(SD$V2[1:40])),xlab = "",ylab = "")
par(new=T)
plot(Parameters,SD$V6,type="o",pch=16,col=4,ylim = c(min(SD$V2[1:40]),max(SD$V2[1:40])),xlab = "Parameter P",ylab = "Standard Deviation")
legend(.35,.65, legend=c("Z1", "Z2", "Z3", "Z4","Z5"),
       col=c(1,2,3,4,5), lty=1,cex=0.8)

plot(Parameters,(PCC$V1),type="o",pch=16,col=1,ylim = c(min((PCC$V1[1:40])),max((PCC$V5[1:40]))),xlab = "",ylab = "")
par(new=T)
plot(Parameters,(PCC$V2),type="o",pch=16,col=2,ylim = c(min((PCC$V1[1:40])),max((PCC$V5[1:40]))),xlab = "",ylab = "")
par(new=T)
plot(Parameters,(PCC$V3),type="o",pch=16,col=3,ylim = c(min((PCC$V1[1:40])),max((PCC$V5[1:40]))),xlab = "",ylab = "")
par(new=T)
plot(Parameters,(PCC$V4),type="o",pch=16,col=4,ylim = c(min((PCC$V1[1:40])),max((PCC$V5[1:40]))),xlab = "",ylab = "")
par(new=T)
plot(Parameters,(PCC$V5),type="o",pch=16,col=5,ylim = c(min((PCC$V1[1:40])),max((PCC$V5[1:40]))),xlab = "",ylab = "")
par(new=T)
plot(Parameters,(PCC$V6),type="o",pch=16,col=6,ylim = c(min((PCC$V1[1:40])),max((PCC$V5[1:40]))),xlab = "",ylab = "")
par(new=T)
plot(Parameters,(PCC$V7),type="o",pch=16,col=7,ylim = c(min((PCC$V1[1:40])),max((PCC$V5[1:40]))),xlab = "Parameter P",ylab = "Pearson's correlation coeficients")
# par(new=T)
# plot(Parameters,abs(PCC$V8),type="o",pch=16,col=8,ylim = c(0,max(abs(PCC$V1[1:38]))),xlab = "",ylab = "")
# par(new=T)
# plot(Parameters,abs(PCC$V9),type="o",pch=16,col=9,ylim = c(0,max(abs(PCC$V1[1:38]))),xlab = "",ylab = "")
# par(new=T)
# plot(Parameters,abs(PCC$V10),type="o",pch=16,col=10,ylim = c(0,max(abs(PCC$V1[1:38]))),xlab = "Parameters P",ylab = "Pearson's correlation coeficients")
legend(.3,-0.38, legend=c("PCC(Z1,Z2)", "PCC(Z1,Z3)", "PCC(Z1,Z4)", "PCC(Z1,Z5)",
                         "PCC(Z2,Z3)","PCC(Z2,Z4)","PCC(Z2,Z5)"),col=c(1,2,3,4,5,6,7), lty=1,cex=0.9)

#Plot the composite index
plot(Parameters[1:40],CompIndex$V4[1:40],type="o",pch=16,col=1,xlab = "Parameter P",ylab = "Composite index")
legend(.3,35, legend=c("Composite index"),cex=0.8,pch=16,lty=1)

#SDd
plot(Parameters[1:40],CompIndex$V1[1:40],type="o",pch=16,col=1,xlab = "Parameter P",ylab = "SDd")

#PCCo
plot(Parameters[1:40],CompIndex$V3[1:40],type="o",pch=16,col=1,xlab = "Parameter P",ylab = "PCCo")

#PCCd
plot(Parameters[1:40],CompIndex$V2[1:40],type="o",pch=16,col=1,xlab = "Parameter P",ylab = "PCCd")
                        
