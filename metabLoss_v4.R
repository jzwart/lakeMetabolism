#Estimating metabolism from DO data
#CTS 29 Jan 2009
#Version 4 Aug 2009

#This function calculates the likelihoods; use in conjunction
#with an optimization routine to minimize NLL

#Use this for one day's data at a time

#dataIn needs to include the following columns
#  DOObs:     DO (mg/L)
#  DOSat:     DO at saturation (mg/L)
#  irr:       irradiance (PAR) (mmol m-2 s-1) >>>>>note that GLEON data may actually be in umol; check
#  kO2:       piston velocity for O2, m (timeStep)-1
#  zMix:      depth of the mixed layer (m)
#  fluxDummy: 0 if zMix is above DO sensor depth (prevents atm flux); 1 otherwise

#Description of parameters
#Note timeStep is number of minutes between DO readings
#Note units given here are those that the model works with; parsIn is given in log units
#  iota   - primary productivity per unit of PAR
#           units are (mg L-1 timeStep-1) / (mmol m-2 s-1)
#  rho    - nighttime respiration
#           units are (mg L-1 (timeStep)-1)
#  DOInit - initial DOHat value
#           units are mg L-1


metabLoss_v4 <- function (parsIn,dataIn) {

##
#Inputs

#Unpack parameters and exponentiate to force positive
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

#Return NLL
return(NLL)

}
