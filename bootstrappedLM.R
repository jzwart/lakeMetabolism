### code for parallelized bootstrapping of lake metabolism model
###		-almost all the code taken from Solomon's metabolism scripts from his 25 lakes R paper
### 11/18/2014
### SEJ

rm(list=ls())

library('parallel')

# arguments
lat=	46.15	#latitutde in decimal degrees
elev=518	#elevatin in m
windHeight=2		#height of wind sensor in m
timeStep=10		#time between observations in minutes
sensorDepth=0.73	#depth sensor is located at in m
rootDir="/Users/JonesLab/Documents/zwart/Long Lake Metabolism/"

doBootstrapping=TRUE

setwd(rootDir)

# set up dirDump       Directory where outputs should be dumped, e.g. 'C:/GLEON/Acton/Results'
if(!file.exists(file.path(rootDir,"Results",'Pelagic'))){
	dir.create(file.path(rootDir,'Results','Pelagic'))
}
dirDump=file.path(rootDir,'Results','Pelagic')

lakes=c('EL','WL')
year=c(2011,2012,2013,2014)

source(file.path(rootDir,"R Code/bootLMsupport.R"))

#Set environment tz variable to GMT
Sys.setenv(tz="GMT")

for(i in 1:length(year)){
	for(j in 1:length(lakes)){
		# generate input forcing data and observations for metabolism fits
	cleanForcings=metabDataSetup(rootDir,year[i],lakes[j],lat,elev,windHeight,timeStep,sensorDepth)

	# fit metabolism model to observations
	# Set up
	#Set up data.frame to store output of optimizations
	nDays <- dim(cleanForcings$sun)[1] - 1  #Not sure if this indexing works appropriately for all lakes
	dateRange <- c(cleanForcings$sun$day[1],cleanForcings$sun$day[nDays])
	outDays <- seq(dateRange[1],dateRange[2],"1 day")
	#optimOut <- data.frame(solarDay=outDays,nll=rep(NA,nDays), iotaEst=rep(NA,nDays),rhoEst=rep(NA,nDays), DOInitEst=rep(NA,nDays), optimCode=rep(NA,nDays))

	
	#Calculate appropriate starting guesses for iota and rho parameters
	#  This takes as reasonable guesses
	#    iota = 1 (mg L-1 d-1)/(mmol m-2 s-1)
	#    rho  = 0.5 mg L-1 d-1
	#  And converts them to appropriate units given time step of model
	iotaGuess <- 3*timeStep/1440 
	rhoGuess <- 0.5*timeStep/1440
	
	#Change directory in prep for dumping results
	setwd(file.path(dirDump,year[i],lakes[j]))

	out=mclapply(X=as.list(outDays),FUN=fitMetabMLEparallel,mc.preschedule=FALSE,mc.cores=3,data2=cleanForcings$data2,sun=cleanForcings$sun,iotaGuess=iotaGuess,rhoGuess=rhoGuess,bootstrap=doBootstrapping,nBoot=1000,ar1.resids=TRUE)


	if(doBootstrapping){
	
		outMat=matrix(NA,length(out),9)
	
		for(q in 1:length(out)){
			outMat[q,]=out[[q]]
		}
	
		colnames(outMat)=c('NLL','iota','R','GPP','DO0','convergence','sd_iota','sd_R','sd_GPP')
		rownames(outMat)=strftime(outDays,format="%Y-%m-%d")
		write.table(outMat,file.path(dirDump,year[i],lakes[j],'results.txt'),sep='\t',quote=F)
	}else{
		outMat=matrix(NA,length(out),6)
	
		for(q in 1:length(out)){
			outMat[q,]=out[[q]]
		}
	
		colnames(outMat)=c('NLL','iota','R','GPP','DO0','convergence')
		rownames(outMat)=strftime(outDays,format="%Y-%m-%d")
		write.table(outMat,file.path(dirDump,year[i],lakes[j],'results.txt'),sep='\t',quote=F)
	}
	}
}
