#metabFunc_v6
#Script for running metabolism model for a single day
#CTS 9 Apr 2010 based on earlier versions


#Inputs are
# lat           Latitude of lake, decimal degrees, north positive
# elev          Elevation of lake surface, m above sea level
# windHeight    Height above lake surface at which wind speed is measured, m
# timeStep      Time interval between DO measurements, minutes
# sensorDepth   Depth of DO sensor, m
# maxZMix       Maximum depth to use for zMix, if temp only measured at one depth
#               This input is not required unless temp is only measured at one depth
# outName       Text to use in labeling outputs, e.g. 'Acton2008'. Character.
# dirData       Directory where data file are located, e.g. 'C:/GLEON/Acton'. Character.
# dataIn        Names of files that contain the data. This should be a character vector of
#                 length 5, e.g. c('Acton_2008_DO.txt','Acton_2008_PAR_5min.txt','Acton_2008_windSpeed.txt',
#                 'Acton_2008_sensorTemp.txt','Acton_2008_tempProfile.txt')
# dirFxns       Directory where functions are located, e.g. 'C:/GLEON/functions'
# dirDump       Directory where outputs should be dumped, e.g. 'C:/GLEON/Acton/Results'

#All of these should be loaded to the workspace before running this script



########################################
#Set up

#Set console output to sink to file 'consoleOutput.txt'
setwd(dirDump)
sink(file='consoleOutput.txt')

#Load in required functions
setwd(dirFxns)
source('metabLoss_v4.R')
source('metabPredix_v4.R')
source('fillHoles.R')
source('calcZMixDens.R')

#Set environment tz variable to GMT
Sys.setenv(tz="GMT")


########################################
#Read and organize data

##
#Read data
setwd(dirData)

#DO
dataDO <- read.table(dataIn[1],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))

#PAR
#Divide PAR by 1000 to convert from measured units (umol m-2 s-1) to model units (mmol m-2 s-1)
dataPAR <- read.table(dataIn[2],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))
dataPAR$PAR <- dataPAR$PAR/1000

#Wind speed
dataWind <- read.table(dataIn[3],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))

#Water temp at depth of DO sensor
dataSensorTemp <- read.table(dataIn[4],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))

#Temp profile
dataTempProfile <- read.table(dataIn[5],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))


##
#Display some info about time grain of measurements

#Print first five time readings for each variable
firstFiveTimes <- data.frame(DO=dataDO$dateTime[1:5],PAR=dataPAR$dateTime[1:5],windSpeed=dataWind$dateTime[1:5],sensorTemp=dataSensorTemp$dateTime[1:5],tempProfile=dataTempProfile$dateTime[1:5])
print('First five time readings for each variable'); print(firstFiveTimes)

#Calculate first differences of time readings, display unique values
cat('\n','First differences of dateTime','\n')
difTimesDO <- diff(dataDO$dateTime); print(table(difTimesDO))
difTimesPAR <- diff(dataPAR$dateTime); print(table(difTimesPAR))
difTimesWindSpeed <- diff(dataWind$dateTime); print(table(difTimesWindSpeed))
difTimesSensorTemp <- diff(dataSensorTemp$dateTime); print(table(difTimesSensorTemp))
difTimesTempProfile <- diff(dataTempProfile$dateTime); print(table(difTimesTempProfile))


##
#Remove rows with duplicate dateTime stamps (and warn)

#Function to find duplicate dateTime stamps
findNotDupRows <- function(dataInName)
  {
  #This function returns the indexes of rows where the dateTime is NOT a duplicate
  #of the dateTime in a previous row
  #dataInName is character, e.g. "dataPAR"
  dataIn <- eval(parse(text=dataInName))
  #Find duplicated time stamps
  dups <- duplicated(dataIn$dateTime)
  #If no duplicated time stamps, notDupRows=all rows in dataIn
  if (all(dups==FALSE))
    {
    notDupRows <- c(1:dim(dataIn)[1])
    } else
  #If at least one time stamp is duplicated, warn, and notDupRows is indexes
  #of rows where dateTime is not duplicate of the dateTime in a previous row
    {
    notDupRows <- which(dups==FALSE)
    nDups <- dim(dataIn)[1]-length(notDupRows)
    print(paste("Warning:",nDups,"rows with duplicate time stamps in",dataInName,"will be removed"))
    }
  #Return dupRows
  return(notDupRows)
  }

notDupRows <- findNotDupRows("dataDO")
dataDO <- dataDO[notDupRows,]

notDupRows <- findNotDupRows("dataPAR")
dataPAR <- dataPAR[notDupRows,]

notDupRows <- findNotDupRows("dataWind")
dataWind <- dataWind[notDupRows,]

notDupRows <- findNotDupRows("dataSensorTemp")
dataSensorTemp <- dataSensorTemp[notDupRows,]

notDupRows <- findNotDupRows("dataTempProfile")
dataTempProfile <- dataTempProfile[notDupRows,]


##
#Make all data sets extend from startTime to endTime by timeStep
#Note that for some lakes it may be necessary to aggregate some variables to coarser time scale to get match up

#Round all time down to nearest timeStep (e.g. if timeStep is 5, round 00:07 to 00:05)
floorMins <- function(dataIn)
  {
  #Pull out dateTime column and name it x
  x <- dataIn$dateTime
  nRows <- length(x)
  #Truncate each dateTime to hour; convert to class numeric
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

dataDO$dateTime <- floorMins(dataDO)
dataPAR$dateTime <- floorMins(dataPAR)
dataWind$dateTime <- floorMins(dataWind)
dataSensorTemp$dateTime <- floorMins(dataSensorTemp)
dataTempProfile$dateTime <- floorMins(dataTempProfile)

#Repeat check for dup rows in case any introduced during floorMins
notDupRows <- findNotDupRows("dataDO")
dataDO <- dataDO[notDupRows,]
notDupRows <- findNotDupRows("dataPAR")
dataPAR <- dataPAR[notDupRows,]
notDupRows <- findNotDupRows("dataWind")
dataWind <- dataWind[notDupRows,]
notDupRows <- findNotDupRows("dataSensorTemp")
dataSensorTemp <- dataSensorTemp[notDupRows,]
notDupRows <- findNotDupRows("dataTempProfile")
dataTempProfile <- dataTempProfile[notDupRows,]

#Find the latest first time point and the earliest last time point of all the data
startTime <- max(min(dataDO$dateTime),min(dataPAR$dateTime),min(dataWind$dateTime),min(dataSensorTemp$dateTime),min(dataTempProfile$dateTime))
endTime <- min(max(dataDO$dateTime),max(dataPAR$dateTime),max(dataWind$dateTime),max(dataSensorTemp$dateTime),max(dataTempProfile$dateTime))

#Data.frame with one column "dateTime" which is sequence of times at time interval of timeStep, from startTime to endTime
completeTimes <- data.frame(dateTime=seq(startTime,endTime,paste(timeStep,"mins")))

#Merge all of input data.frames with completeTimes, so that they now all extend from startTime to endTime by timeStep
dataDO <- merge(completeTimes,dataDO,by="dateTime",all.x=T)
dataPAR <- merge(completeTimes,dataPAR,by="dateTime",all.x=T)
dataWind <- merge(completeTimes,dataWind,by="dateTime",all.x=T)
dataSensorTemp <- merge(completeTimes,dataSensorTemp,by="dateTime",all.x=T)
dataTempProfile <- merge(completeTimes,dataTempProfile,by="dateTime",all.x=T)


########################################
#Calculate sunrise, sunset

#Days of year for which to calculate sunrise and sunset
daysVec <- seq.POSIXt(trunc(startTime,"day"),trunc(endTime,"day"),"1 day")
day <- as.numeric(format(daysVec,format="%j"))

#Factors to convert degrees to radians and vice versa
degToRad <- 2*pi/360
radToDeg <- 180/pi

#Day angle "gamma" (radians). Iqbal 1983 Eq. 1.2.2
dayAngle <- 2*pi*(day-1)/365

#Declination of the sun "delta" (radians). Iqbal 1983 Eq. 1.3.1
dec <- 0.006918 - 0.399912*cos(dayAngle) + 0.070257*sin(dayAngle) - 0.006758*cos(2*dayAngle) +  0.000907*sin(2*dayAngle) - 0.002697*cos(3*dayAngle) + 0.00148*sin(3*dayAngle)

#Sunrise hour angle "omega" (degrees). Iqbal 1983 Eq. 1.5.4
latRad <- lat*degToRad
sunriseHourAngle <- acos(-tan(latRad)*tan(dec))*radToDeg

#Sunrise and sunset times (decimal hours, relative to solar time) Iqbal 1983 Ex. 1.5.1
sunrise <- 12 - sunriseHourAngle/15
sunset <- 12 + sunriseHourAngle/15
# As number of seconds from midnight
sunrise <- sunrise/24*86400
sunset <- sunset/24*86400
# As number of seconds from beginning of year
sunrise <- sunrise+(day-1)*86400
sunset <- sunset+(day-1)*86400
# Convert to POSIXct and round to nearest minute
yr <- format(daysVec,format="%Y")
origin <- paste(yr,"01","01",sep="-")
sunrise <- round(as.POSIXct(sunrise,origin=origin),"mins")
sunset <- round(as.POSIXct(sunset,origin=origin),"mins")

#Create data.frame with sunrise, sunset times for each day
sun <- data.frame(day=daysVec,sunrise,sunset)


########################################
#Re-trim data sets so that they start at or after first sunrise, and end at last time before last sunrise
# i.e. lop off partial day at end

#Trim
startTrim <- min(sun$sunrise)
endTrim <- max(sun$sunrise)
dataDO <- dataDO[dataDO$dateTime >= startTrim & dataDO$dateTime < endTrim,]
dataPAR <- dataPAR[dataPAR$dateTime >= startTrim & dataPAR$dateTime < endTrim,]
dataWind<- dataWind[dataWind$dateTime >= startTrim & dataWind$dateTime < endTrim,]
dataSensorTemp <- dataSensorTemp[dataSensorTemp$dateTime >= startTrim & dataSensorTemp$dateTime < endTrim,]
dataTempProfile <- dataTempProfile[dataTempProfile$dateTime >= startTrim & dataTempProfile$dateTime < endTrim,]
completeTimes <- data.frame(dateTime=completeTimes[completeTimes$dateTime >= startTrim & completeTimes$dateTime < endTrim,])

#(Useful later) Vector giving which solar day each time in completeTimes belongs to
solarDaysBreaks <- sun$sunrise[sun$sunrise <= endTrim]
solarDaysVec <- cut.POSIXt(completeTimes$dateTime,breaks=solarDaysBreaks)


########################################
#Fill gaps in data

##
#DO - do not fill gaps

##
#PAR - linearly interpolate gaps up to 60 min long
dataPAR <- fillHoles(dataPAR,maxLength=60,timeStep=timeStep)

##
#sensorTemp - linearly interpolate gaps up to 60 min long
dataSensorTemp <- fillHoles(dataSensorTemp,maxLength=60,timeStep=timeStep)

##
#windSpeed - fill with daily average as long as at least 80% of data are available

#Loop over days
for (i in 1:length(unique(solarDaysVec)))
  {
  
  #Extract data between sunrise on day i and sunrise on day i+1
  timeSlice <- c(sun$sunrise[i], sun$sunrise[i+1])
  dataTemp <- dataWind[dataWind$dateTime>=timeSlice[1] & dataWind$dateTime<timeSlice[2],]
  
  #Determine total number of observations, and number that are NA
  nTot <- length(dataTemp$windSpeed)
  nNA <- length(which(is.na(dataTemp$windSpeed)))
  
  #If >20% of obs are NA, skip to next i
  if (nNA/nTot > 0.20) next else
  
  {
  #Calculate mean windSpeed and sub in for NA values
  meanSpeed <- mean(dataTemp$windSpeed,na.rm=T)
  naRows <- as.numeric(row.names(dataTemp[is.na(dataTemp$windSpeed),]))
  dataWind$windSpeed[naRows] <- meanSpeed
  }
  }

##
#tempProfile - linearly interpolate gaps up to 60 min long 

nCols <- dim(dataTempProfile)[2]

#Loop over the columns of dataTempProfile
for (i in 2:nCols)
  {
  dataTemp <- dataTempProfile[,c(1,i)]
  dataTempFilled <- fillHoles(dataTemp,maxLength=60,timeStep=timeStep)
  dataTempProfile[,i] <- dataTempFilled[,2]
  }
  

########################################
#Calculate zMix and fluxDummy

#If temperature measured at only one depth, use maxZMix as zMix at every time
if (ncol(dataTempProfile) <= 2)
  {
  dataZMix <- data.frame(dateTime=dataTempProfile$dateTime,zMix=rep(maxZMix,length(dataTempProfile$dateTime)))
  } else

#Otherwise calculate zMix from data
  {
  #Convert tempProfile data to density
  #Density of water (kg m-3) as function of temp from McCutcheon (1999)
  #Note there is a different formula if salinity is appreciable; formula below ignores that
  dataDensProfile <- dataTempProfile 
  dataDensProfile[,-1] <- 1000*(1-((dataDensProfile[,-1]+288.9414)/(508929.2*(dataDensProfile[,-1]+68.12963)))*(dataDensProfile[,-1]-3.9863)^2)

  #Calc zMix
  dataZMix <- calcZMixDens(dataDensProfile)
  }

#Plot zMix
setwd(dirDump)
maxDepth <- max(as.numeric(substr(colnames(dataTempProfile)[2:nCols],5,10)))
pdf(file=paste(outName,'zMix.pdf'))
plot(zMix~dateTime,data=dataZMix,ylim=c(maxDepth,0))
dev.off()

#Identify when to shut off atmospheric flux
#If zMix > sensorDepth, then sensor is in mixed layer and fluxDummy = 1
#If zMix <= sensorDepth, then there is stratification at or above sensor and fluxDummy = 0 -- shut off atmosphere at this time step
fluxDummy <- as.numeric(dataZMix$zMix>sensorDepth)


########################################
#Merge data for convenience

#Merge
data1 <- merge(dataDO,dataPAR,by="dateTime",all=T)
data1 <- merge(data1,dataWind,by="dateTime",all=T)
data1 <- merge(data1,dataSensorTemp,by="dateTime",all=T)
data1 <- merge(data1,dataZMix,by="dateTime",all=T)


########################################
#Report on lengths of NA strings in data

#For each variable in data1, find the length of each run of NA in the data

#Have to temporarily sub in -99999 for NA to get rle() to work as desired
data1Temp <- as.matrix(data1[,2:6])
whichNA <- which(is.na(data1Temp))
data1Temp[whichNA] <- -99999

cat('\n','Lengths of NA strings after fillHoles etc.','\n')

#DO
rleOut <- rle(data1Temp[,"DO"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of DO NA strings')
print(sort(rleOut$lengths[whichNA]))
#PAR
rleOut <- rle(data1Temp[,"PAR"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of PAR NA strings')
print(sort(rleOut$lengths[whichNA]))
#windSpeed
rleOut <- rle(data1Temp[,"windSpeed"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of windSpeed NA strings')
print(sort(rleOut$lengths[whichNA]))
#sensorTemp
rleOut <- rle(data1Temp[,"sensorTemp"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of sensorTemp NA strings')
print(sort(rleOut$lengths[whichNA]))
#zMix
rleOut <- rle(data1Temp[,"zMix"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of zMix NA strings')
print(sort(rleOut$lengths[whichNA]))

rm(data1Temp)


########################################
#Calculate DOSat and kO2 at each time step

##
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


##
#Calculate DO saturation
#Use eqn from Weiss 1970 Deep Sea Res. 17:721-735; simplified since salinity=0
# ln DO = A1 + A2 100/T + A3 ln T/100 + A4 T/100

attach(data1)

#Convert sensorTemp to Kelvin
sensorTempK <- sensorTemp + 273.15

#Weiss equation
A1 <- -173.4292;  A2 <- 249.6339;  A3 <- 143.3483;  A4 <- -21.8492
DOSat <- exp(((A1 + (A2*100/sensorTempK) + A3*log(sensorTempK/100) + A4*(sensorTempK/100))))

#Correction for local average atmospheric pressure
u <- 10^(8.10765 - (1750.286/(235+sensorTemp)))
DOSat <- (DOSat*((atmPres-u)/(760-u)))   #ml/L
DOSat <- DOSat/1000                      #L/L

#Convert using standard temperature and pressure. 
#Similar to calculating saturation DO at STP in ml/L, converting to mg?L (at STP),
#and then doing the above temperature and pressure conversions.
R <- 0.082057  #L atm deg-1 mol-1
O2molWt <- 15.999*2
convFactor <- O2molWt*(1/R)*(1/273.15)*(760/760) #g/L
DOSat <- DOSat*convFactor*1000                   #mg/L

##
#Calculate kO2
wp <- 0.15                       #exponent of wind profile power relationship, Smith 1985 Plant, Cell & Environment 8:387-398
wind10 <- (10/windHeight)^wp * windSpeed
k600 <- 2.07 + 0.215*wind10^1.7  #k600 in cm hr-1 per Cole and Caraco 1998;
k600 <- k600*24/100              #k600 in m day-1
schmidt <- 1800.6 - 120.1*sensorTemp + 3.7818*sensorTemp^2 - 0.047608*sensorTemp^3
kO2 <- k600*(schmidt/600)^-0.5   #Jahne et al. 87. exp could also be -.67
kO2 <- kO2*(timeStep/1440)       #change kO2 to units of m/(timeStep*min)

detach(data1)


########################################
#Fit model to find parameter estimates for each day

##
#Set up

#Organize input data
data2 <- data.frame(dateTime=data1$dateTime, DOObs=data1$DO, DOSat=DOSat, irr=data1$PAR, kO2=kO2, zMix=data1$zMix, fluxDummy=fluxDummy)

#Set up data.frame to store output of optimizations
nDays <- dim(sun)[1] - 1  #Not sure if this indexing works appropriately for all lakes
dateRange <- c(sun$day[1],sun$day[nDays])
outDays <- seq(dateRange[1],dateRange[2],"1 day")
optimOut <- data.frame(solarDay=outDays,nll=rep(NA,nDays), iotaEst=rep(NA,nDays), rhoEst=rep(NA,nDays), DOInitEst=rep(NA,nDays), optimCode=rep(NA,nDays),R2=rep(NA,nDays))

#Calculate appropriate starting guesses for iota and rho parameters
#  This takes as reasonable guesses
#    iota = 1 (mg L-1 d-1)/(mmol m-2 s-1)
#    rho  = 0.5 mg L-1 d-1
#  And converts them to appropriate units given time step of model
iotaGuess <- 1*timeStep/1440
rhoGuess <- 0.5*timeStep/1440

#Remove some stuff
rm(DOSat, kO2, fluxDummy)

#Change directory in prep for dumping results
setwd(dirDump)

##
#Useful stuff for plotting

#Limits for y-axis for drivers
irrLims <- range(data2$irr,na.rm=T)
zMixLims <- c(max(data2$zMix,na.rm=T),0)
atmFluxLims <- c(-10*max(data2$kO2,na.rm=T)*min(data2$zMix,na.rm=T),10*max(data2$kO2,na.rm=T)*min(data2$zMix,na.rm=T))

#Set up pdf device
pdf(file=paste(outName,'daily fits.pdf'),width=11,height=8.5)
layout(rbind(matrix(c(1:14),nrow=2,byrow=F),matrix(c(15:28),nrow=2,byrow=F)),heights=c(1,1.5,1,1.5))
par(mar=c(1,2,0,0)+0.1)

#Save workspace
save(list=ls(),file=paste(outName,'.RData',sep=""))

##
#Run optimization for each day

for (i in 1:nDays)
  {
  
  #Print occasional progress report on i
  if (i %in% seq(1,nDays,10)) (print(paste("Starting day",i)))
  
  #Extract data between sunrise on day i and sunrise on day i+1
  timeSlice <- c(sun$sunrise[i], sun$sunrise[i+1])
  dataTemp <- data2[data2$dateTime>=timeSlice[1] & data2$dateTime<timeSlice[2],]
  
  #If more than 20% of DOObs is missing, or if any NA in DOSat, irr, kO2, or zMix, 
  # return NA for optimization results and plot blank plots
  nTot <- length(dataTemp$DOObs)
  nNA <- length(which(is.na(dataTemp$DOObs)))
  if ((nNA/nTot > 0.20) |  any(is.na(dataTemp[,3:6])))
    {
    optimOut[i,2:6] <- NA
    frame(); frame()
    next
    } else
  
  #Otherwise, fit model and make plots
    {
    
    #For guess of initial DOHat, use first obs unless that is NA, in which case use min obs
    if (is.na(dataTemp$DOObs[1])==F) (DOInit <- dataTemp$DOObs[1]) else (DOInit <- min(dataTemp$DOObs,na.rm=T))

    #Find parameter values by minimizing nll
    parGuess <- log(c(iotaGuess,rhoGuess,DOInit))
    optimTemp <- optim(parGuess,metabLoss_v4,dataIn=dataTemp)
    
    #Save min nll
    optimOut[i,2] <- optimTemp$value
    #Save parameter estimates
    #  Multiply by 1440/timeStep to get from units of timeStep^-1 to units of day^-1
    optimOut[i,3:4] <- exp(optimTemp$par[1:2])*(1440/timeStep)
    optimOut[i,5] <- exp(optimTemp$par[3]) 
    #Save code indicating whether nlm completed successfully
    optimOut[i,6] <- optimTemp$convergence
    
    #Calculate atmFlux and DOHat given max likelihood parameter estimates
    predix <- metabPredix_v4(optimTemp$par,dataTemp)
    DOHat <- predix$DOHat
    atmFlux <- predix$atmFlux
    res <- predix$res
    
    #Calculate SST, SSE, and R2, and save R2
    SST <- sum(dataTemp$DOObs^2,na.rm=T)
    SSE <- sum(res^2,na.rm=T)
    R2 <- (SST-SSE)/SST
    optimOut[i,"R2"] <- R2
    
    #Plot irradiance (orange points), zMix (dashed line), atmFlux (hollow black points)
    #y-axis tick labels are for atmFlux; positive values are flux into lake and negative values are flux out of lake
    par(mar=c(1,2,0,0)+0.1)
    plot(dataTemp$irr~dataTemp$dateTime, ylim=irrLims, axes=F, xlab="", ylab="", pch=18, col="dark orange")
    axis.POSIXct(1,dataTemp$dateTime,labels=F); box()
    text(x=min(dataTemp$dateTime),y=irrLims[2],labels=format(dataTemp$dateTime[1],format="%d-%b"),adj=c(0,1))
    par(new=T); plot(dataTemp$zMix~dataTemp$dateTime, ylim=zMixLims, type="l", lty=2, axes=F, xlab="", ylab="")
    par(new=T); plot(atmFlux~dataTemp$dateTime, ylim=atmFluxLims, axes=F, xlab="", ylab=""); axis(2)
    
    #Plot observed and predicted DO
    yLims <- range(c(DOHat,dataTemp$DOObs),na.rm=T)
    par(mar=c(2,2,0,0)+0.1)
    plot(DOHat ~ dataTemp$dateTime, ylim=yLims, type="l", axes=F, xlab="", ylab="")
    axis.POSIXct(1,dataTemp$dateTime,format="%H:%M")
    axis(2)
    box()
    points(dataTemp$DOObs ~ dataTemp$dateTime)
    meanDOSat <- round(mean(dataTemp$DOSat,na.rm=T),1)
    text(x=min(dataTemp$dateTime),y=yLims[2],labels=paste('DOSat',meanDOSat),adj=c(0,1))
    
    }

  }  #end loop over nDays


#Close pdf graphics device
dev.off()

#Dump optimOut to dirDump
write.table(optimOut,paste(outName,'optimOut.txt'))


#######################################################
# Bookkeeping Calculations below, MCV 6/11/2009
attach(data2)

# calculate gasflux for each datapoint
gasflux <- kO2/timeStep *(DOSat - DOObs) / zMix  # grams O2 m^-3 minute^-1

# Find delta O2 / delta time
dDO<-diff(DOObs)
dT=diff(dateTime)
dDOdT <- dDO / as.numeric(dT,units="mins") # grams O2 m^-3 minute^-1

# calcuate "metabolism" value for each timestep, where m = R in dark, NEP in light
# since metabolism is calculated for a time interval (opposed to a point), use 
#   mean of zMix, and gasflux from start and end of each interval
m <- dDOdT - gasflux[1:length(gasflux)-1]+diff(gasflux)/2
# grams O2 m^-3 minute^-1

# For each night, find areal rate of respiration
R_nightly <- rep(NA,nDays)
for (i in 1:nDays)
  {
  r <- which(dateTime>sun$sunset[i] & dateTime<sun$sunrise[i+1])
  R_nightly[i] <- mean(m[r],na.rm=T)*1440    # grams O2 m^-3 day^-1
  rm(r)
  }

# For each daylight period, average R from night before and after
Rmean<-rep(NA,length(R_nightly)-1)
for (i in 1:length(R_nightly)-1)
  {
  Rmean[i]<-(R_nightly[i]+R_nightly[i+1])/2 # grams O2 m^-3 day^-1
  }

# Fill in R values for first and last days of deployment based on single nights
R <- c(R_nightly[1], Rmean)  # grams O2 m^-3 day^-1

# Calculate daylight NEP (this is NOT 24hr NEP)
Nd <- rep(NA,nDays)
for (i in 1:length(Nd))
  {
   r <- which(dateTime>sun$sunrise[i] & dateTime<sun$sunset[i])
   Nd[i] <- mean(m[r],na.rm=T) * (as.double(sun$sunset[i]-sun$sunrise[i],units="mins"))
   # grams O2 m^-3 daylight period^-1
   # line above multiplies the rate per minute by the number of daylight minutes
   rm(r)
   }
   
# Calculate GPP & true NEP (24hr NEP)  # grams O2 m^-3 day^-1    
GPP <- Nd + (-R * as.double(sun$sunset[1:nDays]-sun$sunrise[1:nDays],units="days") )
# line above multiplies R (per day) by fraction of day that is daylight.
NEP <- GPP + R

#Redefine R so that positive respiration rates are displayed as positive
R <- -R

#Detach data2
detach(data2)

#Write out results of bookkeeping estimates
bookOut <- data.frame(solarDay=optimOut$solarDay,GPP,R)
write.table(bookOut,paste(outName,'bookOut.txt'))


########################################
#Plots

##
#Plots to compare bookkeeping and fitting estimates

#Calculate model-fitting version of GPP from iotaEst and total daily PAR
#First calculate total mmol photons m-2 for each time step
solarFlux <- data2$irr*timeStep*60  #mmol m-2 timeStep-1
#Convert iota units
iotaNewUnits <- optimOut$iotaEst/(60*60*24) #(mg L-1 s-1) / (mmol m-2 s-1)
#Aggregate solarFlux - sum by day
aggSolarFlux <- aggregate(solarFlux,by=list(solarDay=solarDaysVec),sum,na.rm=T)
aggSolarFlux$solarDay <- as.POSIXct(trunc.POSIXt(as.POSIXct(aggSolarFlux$solarDay),"days"))
aggSolarFlux <- merge(aggSolarFlux,data.frame(solarDay=optimOut$solarDay),all.y=T)
totalSolarFlux <- aggSolarFlux[,2] #mmol m-2 d-1
#Calc GPP from fitted iota
GPPFit <- iotaNewUnits*totalSolarFlux #mg L-1 d-1

#Write out GPPFit
GPPFitOut <- data.frame(solarDay=optimOut$solarDay,GPPFit)
write.table(GPPFitOut,paste(outName,'GPPFitOut.txt'))

#Plot R vs rhoEst and GPP vs. GPPFit
pdf(file=paste(outName,'model comparison.pdf'),width=10,height=5)
par(mfrow=c(1,2))
plot(R ~ optimOut$rhoEst, xlab="R (max likelihood est)", ylab="R (bookkeeping est)"); abline(0,1)
plot(GPP ~ GPPFit, xlab="GPP (max likelihood est)", ylab="GPP (bookkeeping est)"); abline(0,1)
dev.off()


##
#Plot estimates of iota, GPP, and rho
pdf(file=paste(outName,'time series of iota GPP rho.pdf'),width=8.5,height=11)
par(mfrow=c(3,1),mar=c(3,4,1,1))
plot(iotaEst ~ solarDay, data=optimOut, xlab="", ylab="iotaEst, (mg L-1 d-1) / (mmol m-2 s-1)")
plot(GPPFit ~ optimOut$solarDay, xlab="", ylab="GPP, (mg L-1 d-1)")
plot(rhoEst ~ solarDay, data=optimOut, xlab="", ylab="rhoEst, (mg L-1 d-1)")
dev.off()


##
#Turn off console sink
sink()


##END
