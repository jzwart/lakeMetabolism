#MCV and CTS
#Version 27 Aug 2009

#Using density instead of temp to indicate mixed layer
#Though note that objects in this script are still labeled e.g. 'temps', not 'dens'

#Inputs:
# dataIn - A density profile.
#          -First column is dateTime, POSIXct
#          -Second and subsequent columns are density (kg m-3) at depth. Column headers
#           for these columns are e.g. temp0, temp0.5, temp1.0, temp2, temp2.3
# thresh - Threshold density change to indicate end of mixed layer. Units (kg m-3) m-1.
#          Defaults to 0.075 (value Jordan Read is using in their project)

#Outputs:
# metaDepthOut - A data.frame, first column is dateTime, second column is zMix calculated at that dateTime

calcZMixDens <- function(dataIn,thresh=0.075)
{

#Useful things
nObs <- dim(dataIn)[1]
nCols <- dim(dataIn)[2]

#Extract vector of depths from column names of dataIn
depths <- as.numeric(substr(colnames(dataIn)[2:nCols],5,10))

#Set up structure to hold calculated metaDepth for each time point
metaDepthOut <- data.frame(dateTime=dataIn$dateTime,zMix=rep(NA,nObs))

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
