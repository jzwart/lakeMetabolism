#Return DOHat, residuals, atmFlux, and NLL given par estimates
#CTS 18 Aug 2009

metabPredix_v4 <- function (parsIn,dataIn) {

##
#Inputs

#Unpack parameters and exponentiate to force positive,
iota <- exp(parsIn[1])
rho <- exp(parsIn[2])
DOInit <- exp(parsIn[3])

#Label useful things
nObs <- dim(dataIn)[1]
kO2 <- dataIn$kO2
DOSat <- dataIn$DOSat
zMix <- dataIn$zMix
irr <- dataIn$irr
DOObs <- dataIn$DOObs
fluxDummy <- dataIn$fluxDummy

##
#Calculate predictions and residuals

#Set up output
DOHat <- rep(NA,nObs)
atmFlux <- rep(NA,nObs)  
#Initialize DOHat
DOHat[1] <- DOInit

#Calculate atmFlux and predicted DO for each time point
#Fluxes out of lake have negative sign
for (i in 1:(nObs-1)) {
  atmFlux[i] <- fluxDummy[i] * -kO2[i] * (DOHat[i] - DOSat[i]) / zMix[i]  
  DOHat[i+1] <- DOHat[i] + iota*irr[i] - rho + atmFlux[i]
  }

#Compare observed and predicted DO; calculate residuals and NLL
#Exclude from calculation any cases where DOObs=NA
if (any(is.na(DOObs)))
  {
  NAObs <- which(is.na(DOObs))
  res <- DOObs[-NAObs] - DOHat[-NAObs]
  } else
  
  {
  res <- DOObs - DOHat
  }

nRes <- length(res)
SSE <- sum(res^2)
sigma2 <- SSE/nRes
NLL <- 0.5*((SSE/sigma2) + nRes*log(2*pi*sigma2))

#Set up output structure
dataOut <- list(atmFlux=atmFlux,DOHat=DOHat,res=res,NLL=NLL)

#Return dataOut
return(dataOut)

}
