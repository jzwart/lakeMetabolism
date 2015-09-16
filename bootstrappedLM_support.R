### supporting functions for parallelized bootstrapping of lake metabolism model
### 11/18/2014
### SEJ



######## setup dataframes for model fitting 
metabDataSetup<-function(rootDir,year,lake,lat=46.15,elev=518,windHeight=2,timeStep=10,sensorDepth=0.5){
	### setup directories
	outName=paste(lake,year,sep='_')
	# dirData       Directory where data file are located, e.g. 'C:/GLEON/Acton'. Character.
	dirData=file.path(rootDir,'Data/Pelagic',year,lake)
	# dataIn        Names of files that contain the data. This should be a character vector of
	dataIn=c("DO.txt","PAR.txt","WS.txt","TEMP.txt","TEMP_PROF.txt")
	
	### load data
	setwd(dirData)
	#DO
	dataDO <- read.table(dataIn[1],header=T,sep='\t',colClasses=c(datetime="POSIXct"))
	#PAR
	#Divide PAR by 1000 to convert from measured units (umol m-2 s-1) to model units (mmol m-2 s-1)
	dataPAR <- read.table(dataIn[2],header=T,sep='\t',colClasses=c(datetime="POSIXct"))
	dataPAR$PAR <- dataPAR$PAR/1000
	#Wind speed
	dataWind <- read.table(dataIn[3],header=T,sep='\t',colClasses=c(datetime="POSIXct"))
	#Water temp at depth of DO sensor
	dataSensorTemp <- read.table(dataIn[4],header=T,sep='\t',colClasses=c(datetime="POSIXct"))
	#Temp profile
	dataTempProfile <- read.table(dataIn[5],header=T,sep='\t',colClasses=c(datetime="POSIXct"))

	#### Remove duplicate rows
	dataDO=dataDO[!duplicated(dataDO$datetime),]
	dataPAR=dataPAR[!duplicated(dataPAR$datetime),]
	dataWind=dataWind[!duplicated(dataWind$datetime),]
	dataSensorTemp=dataSensorTemp[!duplicated(dataSensorTemp$datetime),]
	dataTempProfile=dataTempProfile[!duplicated(dataTempProfile$datetime),]
					
	#Make all data sets extend from startTime to endTime by timeStep
	dataDO$datetime <- floorMins(dataDO)
	dataPAR$datetime <- floorMins(dataPAR)
	dataWind$datetime <- floorMins(dataWind)
	dataSensorTemp$datetime <- floorMins(dataSensorTemp)
	dataTempProfile$datetime <- floorMins(dataTempProfile)

	#Repeat check for dup rows in case any introduced during floorMins
	dataDO=dataDO[!duplicated(dataDO$datetime),]
	dataPAR=dataPAR[!duplicated(dataPAR$datetime),]
	dataWind=dataWind[!duplicated(dataWind$datetime),]
	dataSensorTemp=dataSensorTemp[!duplicated(dataSensorTemp$datetime),]
	dataTempProfile=dataTempProfile[!duplicated(dataTempProfile$datetime),]

	#Find the latest first time point and the earliest last time point of all the data
	startTime <- max(min(dataDO$datetime),min(dataPAR$datetime),min(dataWind$datetime),min(dataSensorTemp$datetime),min(dataTempProfile$datetime))
	endTime <- min(max(dataDO$datetime),max(dataPAR$datetime),max(dataWind$datetime),max(dataSensorTemp$datetime),max(dataTempProfile$datetime))

	#Data.frame with one column "datetime" which is sequence of times at time interval of timeStep, from startTime to endTime
	completeTimes <- data.frame(datetime=seq(startTime,endTime,paste(timeStep,"mins")))

	#Merge all of input data.frames with completeTimes, so that they now all extend from startTime to endTime by timeStep
	dataDO <- merge(completeTimes,dataDO,by="datetime",all.x=T)
	dataPAR <- merge(completeTimes,dataPAR,by="datetime",all.x=T)
	dataWind <- merge(completeTimes,dataWind,by="datetime",all.x=T)
	dataSensorTemp <- merge(completeTimes,dataSensorTemp,by="datetime",all.x=T)
	dataTempProfile <- merge(completeTimes,dataTempProfile,by="datetime",all.x=T)

	#### Calculate sunrise, sunset
	#Days of year for which to calculate sunrise and sunset
	daysVec <- seq.POSIXt(trunc(startTime,"day"),trunc(endTime,"day"),"1 day")
	
	
	
	sun=data.frame(day=daysVec,sunrise=NA,sunset=NA)
	for(i in 1:length(daysVec)){
		sun[i,2:3]=findSunriseSunset(dataPAR,daysVec[i])
	}
	sun$sunrise=as.POSIXct(sun$sunrise,origin="1970-01-01")
	sun$sunset=as.POSIXct(sun$sunset,origin="1970-01-01")
	
	###Re-trim data sets so that they start at or after first sunrise, and end at last time before last sunrise
	# i.e. lop off partial day at end
	startTrim <- min(sun$sunrise)
	endTrim <- max(sun$sunrise)
	dataDO <- dataDO[dataDO$datetime >= startTrim & dataDO$datetime < endTrim,]
	dataPAR <- dataPAR[dataPAR$datetime >= startTrim & dataPAR$datetime < endTrim,]
	dataWind<- dataWind[dataWind$datetime >= startTrim & dataWind$datetime < endTrim,]
	dataSensorTemp <- dataSensorTemp[dataSensorTemp$datetime >= startTrim & dataSensorTemp$datetime < endTrim,]
	dataTempProfile <- dataTempProfile[dataTempProfile$datetime >= startTrim & dataTempProfile$datetime < endTrim,]
	completeTimes <- data.frame(datetime=completeTimes[completeTimes$datetime >= startTrim & completeTimes$datetime < endTrim,])

	#(Useful later) Vector giving which solar day each time in completeTimes belongs to
	solarDaysBreaks <- sun$sunrise[sun$sunrise <= endTrim]
	solarDaysVec <- cut.POSIXt(completeTimes$datetime,breaks=solarDaysBreaks)

	### Fill gaps in data
	#DO - do not fill gaps
	#PAR - linearly interpolate gaps up to 60 min long
	dataPAR <- fillHoles(dataPAR,maxLength=60,timeStep=timeStep)
	#sensorTemp - linearly interpolate gaps up to 60 min long
	dataSensorTemp <- fillHoles(dataSensorTemp,maxLength=60,timeStep=timeStep)
	#windSpeed - fill with daily average as long as at least 80% of data are available
	#Loop over days
	for (i in 1:length(unique(solarDaysVec))){
  		#Extract data between sunrise on day i and sunrise on day i+1
  		timeSlice <- c(sun$sunrise[i], sun$sunrise[i+1])
  		dataTemp <- dataWind[dataWind$datetime>=timeSlice[1] & dataWind$datetime<timeSlice[2],]
  
  		#Determine total number of observations, and number that are NA
  		nTot <- length(dataTemp$WS)
  		nNA <- length(which(is.na(dataTemp$WS)))
  
  		#If >20% of obs are NA, skip to next i
  		if (nNA/nTot > 0.20) next else{
 			#Calculate mean windSpeed and sub in for NA values
 			meanSpeed <- mean(dataTemp$WS,na.rm=T)
  			naRows <- as.numeric(row.names(dataTemp[is.na(dataTemp$WS),]))
  			dataWind$WS[naRows] <- meanSpeed
  		}
  	}
	#tempProfile - linearly interpolate gaps up to 60 min long 
	nCols <- dim(dataTempProfile)[2]
	#Loop over the columns of dataTempProfile
	for (i in 2:nCols){
		dataTemp <- dataTempProfile[,c(1,i)]
  		dataTempFilled <- fillHoles(dataTemp,maxLength=1000000,timeStep=timeStep)	#***** set maxLength to enable long interpolations or not
  		dataTempProfile[,i] <- dataTempFilled[,2]
	}
  

	####Calculate zMix and fluxDummy
	#If temperature measured at only one depth, use maxZMix as zMix at every time
	if (ncol(dataTempProfile) <= 2){
  		dataZMix <- data.frame(datetime=dataTempProfile$datetime,zMix=rep(maxZMix,length(dataTempProfile$datetime)))
  	}else{
	#Otherwise calculate zMix from data
  		#Convert tempProfile data to density
  		#Density of water (kg m-3) as function of temp from McCutcheon (1999)
  		#Note there is a different formula if salinity is appreciable; formula below ignores that
  		dataDensProfile <- dataTempProfile 
  		dataDensProfile[,-1] <- 1000*(1-((dataDensProfile[,-1]+288.9414)/(508929.2*(dataDensProfile[,-1]+68.12963)))*(dataDensProfile[,-1]-3.9863)^2)

	  	#Calc zMix
  		dataZMix <- calcZMixDens(dataDensProfile)
  	}

	####Identify when to shut off atmospheric flux
	# If zMix > sensorDepth, then sensor is in mixed layer and fluxDummy = 1
	# If zMix <= sensorDepth, then there is stratification at or above sensor and fluxDummy = 0 -- shut off atmosphere at this time step
	fluxDummy <- as.numeric(dataZMix$zMix>sensorDepth)
	
	###Merge data for convenience
	data1 <- merge(dataDO,dataPAR,by="datetime",all=T)
	data1 <- merge(data1,dataWind,by="datetime",all=T)
	data1 <- merge(data1,dataSensorTemp,by="datetime",all=T)
	data1 <- merge(data1,dataZMix,by="datetime",all=T)
	
	####Calculate DOSat and kO2 at each time step

	#Calculate average atmospheric pressure at elevation of lake
	#Using the 'barometric formula' - see e.g. U.S. Standard Atmosphere 1976 or
	# Jacobs 1999 Atmospheric Chemistry, Eqn 2.9
	#Values of Rstar, g0, M are according to US Standard Atmosphere 1976; could use SI instead

	#Constants
	Pb <- 101325        #static pressure, pascals
	Tb <- 288.15        #standard temp, K
	Lb <- -0.0065       #standard temp lapse rate, K m-1
	h <- elev           #elevation above sea level, m
	hb <- 0             #elevation at bottom of atmospheric layer 0, m (note layer 0 extends to 11000 masl)
	Rstar <-  8.31432   #universal gas constant, N m mol-1 K-1 (equiv to J K-1 mol-1)  SI: 8.314472
	g0 <- 9.80665       #acceleration of gravity, m s-1
	M <- 0.0289644      #molar mass of Earth's air, kg mol-1

	#Pressure, in Pa (pascals)
	P <- Pb * (Tb/(Tb+Lb*(h-hb)))^(g0*M/(Rstar*Lb))
	# In mmHg
	atmPres <- P*0.00750061683
	
	###Calculate DO saturation
	#Use eqn from Weiss 1970 Deep Sea Res. 17:721-735; simplified since salinity=0
	# ln DO = A1 + A2 100/T + A3 ln T/100 + A4 T/100
	attach(data1)

	#Convert sensorTemp to Kelvin
	sensorTempK <- dataSensorTemp[,2] + 273.15

	#Weiss equation
	A1 <- -173.4292;  A2 <- 249.6339;  A3 <- 143.3483;  A4 <- -21.8492
	DOSat <- exp(((A1 + (A2*100/sensorTempK) + A3*log(sensorTempK/100) + A4*(sensorTempK/100))))

	#Correction for local average atmospheric pressure
	u <- 10^(8.10765 - (1750.286/(235+dataSensorTemp[,2])))
	DOSat <- (DOSat*((atmPres-u)/(760-u)))   #ml/L
	DOSat <- DOSat/1000                      #L/L

	#Convert using standard temperature and pressure. 
	#Similar to calculating saturation DO at STP in ml/L, converting to mg?L (at STP),
	#and then doing the above temperature and pressure conversions.
	R <- 0.082057  #L atm deg-1 mol-1
	O2molWt <- 15.999*2
	convFactor <- O2molWt*(1/R)*(1/273.15)*(760/760) #g/L
	DOSat <- DOSat*convFactor*1000                   #mg/L

	###Calculate kO2
	wp <- 0.15                       #exponent of wind profile power relationship, Smith 1985 Plant, Cell & Environment 8:387-398
	wind10 <- (10/windHeight)^wp * dataWind[,2]
	k600 <- 2.07 + 0.215*wind10^1.7  #k600 in cm hr-1 per Cole and Caraco 1998;
	k600 <- k600*24/100              #k600 in m day-1
	schmidt <- 1800.6 - 120.1*dataSensorTemp[,2] + 3.7818*dataSensorTemp[,2]^2 - 0.047608*dataSensorTemp[,2]^3
	kO2 <- k600*(schmidt/600)^-0.5   #Jahne et al. 87. exp could also be -.67
	kO2 <- kO2*(timeStep/1440)       #change kO2 to units of m/(timeStep*min)

	detach(data1)
	
	data2=data.frame(datetime=data1$datetime,DOObs=data1$DO,DOSat=DOSat,irr=data1$PAR,kO2=kO2,zMix=data1$zMix,fluxDummy=fluxDummy)
	
	return(list(data2=data2,sun=sun))
}


###### fit metabolism model to data using output of metabDataSetup() for a given day
		# this is used by mclapply()
fitMetabMLEparallel<-function(date,data2,sun,iotaGuess,rhoGuess,bootstrap=FALSE,ar1.resids=TRUE,nBoot=1000){
	# pull day's data
	timeSlice<-c(sun$sunrise[sun$day==date],sun$sunrise[which(sun$day==date)+1])
	dataTemp <- data2[data2$datetime>=timeSlice[1] & data2$datetime<timeSlice[2],]	

	#Run optimization
  	#If more than 20% of DOObs is missing, or if any NA in DOSat, irr, kO2, or zMix, 
  	# return NA for optimization results and plot blank plots
 	nTot <- length(dataTemp$DOObs)
  	nNA <- length(which(is.na(dataTemp$DOObs)))
  	if ((nNA/nTot > 0.20) |  any(is.na(dataTemp[,3:6]))){
    	stop("too many NAs")
    }else{
	#Otherwise, fit model
    
    	#For guess of initial DOHat, use first obs unless that is NA, in which case use min obs
    	if (is.na(dataTemp$DOObs[1])==F) (DOInit <- dataTemp$DOObs[1]) else (DOInit <- min(dataTemp$DOObs,na.rm=T))

		#Find parameter values by minimizing nll
    	parGuess <- log(c(iotaGuess,rhoGuess,DOInit))
    	optimTemp <- optim(parGuess,metabLoss_v4,dataIn=dataTemp)
    
    	curFit=c(optimTemp$value,exp(optimTemp$par[1:2])*(1440/timeStep),exp(optimTemp$par[1])*(1440/timeStep)/(24*60*60)*sum(dataTemp$irr*timeStep*60),exp(optimTemp$par[3]),optimTemp$convergence)
    	#  Multiply by 1440/timeStep to get from units of timeStep^-1 to units of day^-1
   	}

	# if bootstrapping do so...
	if(bootstrap){
		bootResult=bootstrap.metab(parGuess=parGuess,dataTemp=dataTemp,optimTemp=optimTemp,n=nBoot,ar1.resids=ar1.resids)
		curFitBoot=c(curFit,sd(bootResult$iota),sd(bootResult$rho),sd(bootResult$GPP))
		return(curFitBoot)
	}else{
		return(curFit)
	}
}

# Bootstrapping function for metabolism model; JAZ 2014-11-3; modified from Winslow and GLEON fellows 
bootstrap.metab <- function(parGuess, dataTemp, optimTemp, n=1000, ar1.resids=FALSE){
	#Calculate atmFlux and DOHat given max likelihood parameter estimates
  	predix <- metabPredix_v4(optimTemp$par,dataTemp)
    DOHat <- predix$DOHat
    atmFlux <- predix$atmFlux
    res <- predix$res
    
	n.obs = length(dataTemp$DOObs)
  
	doHat  = metabPredix_v4(optimTemp$par,dataTemp)$DOHat
	resids = dataTemp$DOObs - doHat
  
	#If we are maintaining the ar1 component of the residuals, 
	# we must estimate ar1 coeff and the ar1 residual standard deviation
	if(ar1.resids){
    	ar1.lm    = lm(resids[1:n.obs-1] ~ resids[2:n.obs]-1)
    	ar1.coeff = ar1.lm$coefficients
    	ar1.sd    = sd(ar1.lm$residuals)
  	}
  
  	#Pre-allocate the result data frame
  	result <- data.frame(boot.iter = 1:n,iota = rep(NA,n),rho = rep(NA,n),DOInit = rep(NA,n),covergence = rep(NA,n),nll = rep(NA,n),GPP = rep(NA,n))
  
	for(i in 1:n){
    	#Randomize the residuals using one of two methods
    	if(ar1.resids){ #residual randomization keeping the ar1 data structure
      		simRes = rep(NA, n.obs)
      		simRes[1] = sample(resids[!is.na(resids)],1)
      		for(j in 2:n.obs){
        		simRes[j] = ar1.coeff*simRes[j-1] + rnorm(n=1, sd=ar1.sd)
      		}      
    	}else{ #Raw residual randomization
      		#Randomize residuals without replacement
      		simRes = sample(resids[!is.na(resids)], length(resids), replace=FALSE) 
    	}
    
    	doSim = doHat + simRes
    
    	#Run optim again with new simulated DO signal
    	dataBoot<-dataTemp
    	dataBoot$DOObs<-doSim
    	optimTemp <- optim(parGuess,metabLoss_v4,dataIn=dataBoot)
    
    	result[i,2:3] <- exp(optimTemp$par[1:2])*(1440/timeStep) #iota and rho to units of mg O2 L-1 day-1 
    	result[i,4] <- exp(optimTemp$par[3]) #initial DO estimate 
    	result[i,5] <- optimTemp$convergence  #did model converge or not (0=yes, 1=no)
    	result[i,6] <- optimTemp$value #value of nll 
    	result[i,7] <- result[i,2]/(24*60*60)*sum(dataTemp$irr*60*timeStep) #GPP in units of mg O2 L-1 d-1  
  	}
  return(result)
}



####### Other support functions

#Round all time down to nearest timeStep (e.g. if timeStep is 5, round 00:07 to 00:05)
floorMins <- function(dataIn){
  #Pull out datetime column and name it x
  x <- dataIn$datetime
  nRows <- length(x)
  #Truncate each datetime to hour; convert to class numeric
  floorHour <- as.POSIXct(trunc(x[1:nRows],"hour"))
  floorNumeric <- as.numeric(floorHour)
  #Create sequence from floorNumeric to next hour by timeStep (in seconds)
  seqSec <- seq(0,3600,60*timeStep)
  #Create matrix where each row is floorNumeric + the elements of seqSec
  matSec <- matrix(rep(seqSec,nRows),nrow=nRows,byrow=T)
  matSec <- floorNumeric + matSec
  #Calculate abs(time difference) between each element of x and the timeStep intervals
  difs <- abs(as.numeric(x) - matSec)
  #Find the minimum absolute difference in each row and select the corresponding time from matSec
  whichMin <- apply(difs,1,which.min)
  rowNames <- as.numeric(rownames(data.frame(whichMin)))
  matIndex <- (whichMin-1)*nRows + rowNames
  matSecFlat <- matrix(matSec,ncol=1)
  outTime <- as.POSIXct(matSecFlat[matIndex],origin="1970-01-01")
  #Return outTime
  return(outTime)
}



#MCV and CTS
#Version 27 Aug 2009

#Using density instead of temp to indicate mixed layer
#Though note that objects in this script are still labeled e.g. 'temps', not 'dens'

#Inputs:
# dataIn - A density profile.
#          -First column is datetime, POSIXct
#          -Second and subsequent columns are density (kg m-3) at depth. Column headers
#           for these columns are e.g. temp0, temp0.5, temp1.0, temp2, temp2.3
# thresh - Threshold density change to indicate end of mixed layer. Units (kg m-3) m-1.
#          Defaults to 0.075 (value Jordan Read is using in their project)

#Outputs:
# metaDepthOut - A data.frame, first column is datetime, second column is zMix calculated at that datetime

calcZMixDens <- function(dataIn,thresh=0.075)
{

#Useful things
nObs <- dim(dataIn)[1]
nCols <- dim(dataIn)[2]

#Extract vector of depths from column names of dataIn
depths <- as.numeric(substr(colnames(dataIn)[2:nCols],5,10))
# depths=c(0.1,2.0,3.0,4.0,10.0)

#Set up structure to hold calculated metaDepth for each time point
metaDepthOut <- data.frame(datetime=dataIn$datetime,zMix=rep(NA,nObs))

#Calculate zMix for each time in dataIn
for (i in 1:nObs)       #nObs
  {
  
  #If all temp data at this time point are NA, return NA for zMix at this time
  if (all(is.na(dataIn[i,2:nCols])))
    {
    metaDepthOut[i,2] <- NA
    next
    }
    
  #Get vectors of depths and temps, excluding cols where temp is NA
  temps <- unlist(dataIn[i,2:nCols])
  DEPTHS <- depths[!(is.na(temps))]
  TEMPS <- temps[!(is.na(temps))]

  #Find change in depth, change in temp, and change in temp per depth
  dz <- diff(DEPTHS)
  dT <- diff(TEMPS)
  dTdz <- dT/dz
 
  #If no rate of change exceeds thresh, return max depth for zMix at this time
  if (all(dTdz<thresh))
    {
    metaDepthOut[i,2] <- max(depths)
    next
    }
  
  #Identify first depth where rate of change exceeds thresh
  whichFirstBigChange <- min(which(dTdz >= thresh))
  metaDepth <- mean(c(DEPTHS[whichFirstBigChange],DEPTHS[(whichFirstBigChange+1)]))
  
  #Save result to metaDepthOut
  metaDepthOut[i,2] <- metaDepth
  }
  
#Return metaDepthOut
return(metaDepthOut)

}


#Linearly interpolate values for strings (up to specified length) of missing data
#CTS 31 Jul 2009

#Inputs
# dataIn:     A data.frame with two columns, first is "datetime", second is data
# maxLength:  Maximum length of NA string that you are willing to interpolate across.
#             NOTE that this is given in minutes
# timeStep:   The time step of the data


fillHoles <- function(dataIn,maxLength,timeStep)
{

#Number of rows in dataIn
nObs <- dim(dataIn)[1]

#Express maxLength as number of time steps instead of number of minutes
maxLength <- maxLength/timeStep

#Temporarily replace NAs in data with -99999
whichNA <- which(is.na(dataIn[,2]))
if(length(whichNA)){
dataIn[whichNA,2] <- -99999
}
#Identify strings of NA (-99999) values
rleOut <- rle(dataIn[,2])
which9 <- which(rleOut$values==-99999)

#If no NA strings in data, return dataIn
if (length(which9)==0)
  {
  return(dataIn)
  } else

#Otherwise, continue
  {

  #Interpolate valus for each string of NA values
  for (i in 1:length(which9))
    {
    
    #Determine start and end index of string i, and calculate length of string
    if (which9[i]==1)
      {
      stringStart <- 1
      } else
      
      {
      stringStart <- 1 + sum(rleOut$lengths[1:(which9[i]-1)])
      }
      
    stringEnd <- sum(rleOut$lengths[1:which9[i]])
    stringLength <- stringEnd-stringStart+1
    
    #Skip to next i if:
    #  -length of string exceeds maxLength,
    #  -stringStart is the first obs in dataIn
    #  -stringEnd is the last obs in dataIn
    if (stringLength > maxLength | stringStart==1 | stringEnd==nObs) next else
    {
    
    #Linearly interpolate missing values
    interp <- approx(x=c(dataIn[stringStart-1,"datetime"],dataIn[stringEnd+1,"datetime"]),
                     y=c(dataIn[stringStart-1,2],dataIn[stringEnd+1,2]),
                     xout=dataIn[stringStart:stringEnd,"datetime"],
                     method="linear")
    dataIn[stringStart:stringEnd,2] <- interp$y
    }
    }

  #Replace any remaining -99999 with NA
  dataIn[which(dataIn[,2]==-99999),2] <- NA

  #Return result
  return(dataIn)
  }

}

findSunriseSunset<-function(par,day){
	#pull day of interest
	cur=par[as.POSIXct(trunc(par[,1],"day"))==day,]
	cur=cur[!is.na(cur$PAR),]
	
	# lowest PAR measurement Long Lake met station returns is 0.0012 mM
	return(range(cur[cur[,2]>0.0012,1]))
}


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
