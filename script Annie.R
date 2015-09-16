#Annie script
#9 Apr 2010

#Directory where results of model fitting should get dumped
dirDump <- 'J:/Research projects/GLEON/Respiration/Metabolism model/Model v6/Results/Annie'

#Lake name and year, to be used in labeling outputs
outName <- 'Annie2008'

#Directory where data files are located
dirData <- 'J:/Research projects/GLEON/Respiration/Core data - QAd FINAL/Annie'

#Names of files to import
dataIn <- c('Annie_2008_DO.txt','Annie_2008_PAR.txt','Annie_2008_windSpeed.txt',
            'Annie_2008_sensorTemp.txt','Annie_2008_tempProfile.txt')

#Set pars
lat <- 27.207       #Latitude of lake, decimal degrees. N lats are positive, S lats are negative
elev <- 3.7         #Elevation above sea level at surface of lake, m
windHeight <- 10   #height above lake at which wind speed is meaured, m
timeStep <- 15       #number of minutes between DO measurements
sensorDepth <- 1.35  #depth of DO sensor, m 

#Directory where functions are located - this can be the same for any lake you run
dirFxns <- 'J:/Research projects/GLEON/Respiration/Metabolism model/Model v6'

#Run main metab script
setwd(dirFxns)
source('metabFunc_v6.R')
