# Background respriation estimation from  bootstrapped metabolism data; original code from CTS; modified by JAZ, 2015-04-14

# loading in bootstrapping data for background R and priming calculation 
lakes<-c('EL','WL') 
years<-c(2011,2012,2013,2014)
trt<-c('Pre','Post') 

for(i in 1:length(lakes)){
  curBootOut<-data.frame()
  for(j in 1:length(trt)){
    curFiles<-list.files(file.path('/Users/Jake/Documents/Jake/MyPapers/Long Lake Metabolism/Results/PelagicBootstrapped/',trt[j],lakes[i]))
    curFiles<-curFiles[grep('bootOut.txt',curFiles)]
    if(j==1){ #pre temp data 
      curTemp<-read.table(file.path('/Users/Jake/Documents/Jake/MyPapers/Long Lake Metabolism/Data/Pelagic/','2011',lakes[i],'TEMP.txt'),
                          sep='\t',header=T,stringsAsFactor=F)
      curTemp2<-read.table(file.path('/Users/Jake/Documents/Jake/MyPapers/Long Lake Metabolism/Data/Pelagic/','2012',lakes[i],'TEMP.txt'),
                           sep='\t',header=T,stringsAsFactor=F)
      curTemp<-rbind(curTemp,curTemp2)
      curTemp<-subset(curTemp,!duplicated(curTemp$datetime))
      curTemp$datetime<-strftime(curTemp$datetime,'%Y-%m-%d')
      curTemp<-aggregate(curTemp$TEMP,by=list(curTemp$datetime),FUN=mean,na.rm=T) # mean daily temperature 
      colnames(curTemp)<-c('solarDay','temp')
    }else{ #post temp data 
      curTemp<-read.table(file.path('/Users/Jake/Documents/Jake/MyPapers/Long Lake Metabolism/Data/Pelagic/','2013',lakes[i],'TEMP.txt'),
                          sep='\t',header=T,stringsAsFactor=F)
      curTemp2<-read.table(file.path('/Users/Jake/Documents/Jake/MyPapers/Long Lake Metabolism/Data/Pelagic/','2014',lakes[i],'TEMP.txt'),
                           sep='\t',header=T,stringsAsFactor=F)
      curTemp<-rbind(curTemp,curTemp2)
      curTemp<-subset(curTemp,!duplicated(curTemp$datetime))
      curTemp$datetime<-strftime(curTemp$datetime,'%Y-%m-%d')
      curTemp<-aggregate(curTemp$TEMP,by=list(curTemp$datetime),FUN=mean,na.rm=T) # mean daily temperature 
      colnames(curTemp)<-c('solarDay','temp')
      
    }

    dataBoot<-data.frame()
    for(q in 1:length(curFiles)){
      cur<-read.table(file.path('/Users/Jake/Documents/Jake/MyPapers/Long Lake Metabolism/Results/PelagicBootstrapped/',trt[j],lakes[i],
                               curFiles[q]),sep='\t',header=T,stringsAsFactor=F)
      dataBoot<-rbind(dataBoot,cur)
    }
    dataBoot<-merge(dataBoot,curTemp,by='solarDay',all.x=T)
    
    regrCoefs <- data.frame(b0=rep(NA,1000),b1=rep(NA,1000))
    for (b in 1:1000) {
      if (b %in% c(100,200,300,400,500,600,700,800,900)) (print(b))
      #Pull out rhoEst and iotaEst for each day from the bth bootstrapped dataset
      dataB <- dataBoot[dataBoot$boot.iter==b,c("rho","iota",'GPP','temp')]
      #Temperature-correct R and GPP to 20C - Holtgrieve et al. 2011, Venkiteswaran et al. 2007, etc
      R20 <- dataB$rho*1.047^(20-dataB$temp)
      P20 <- dataB$GPP*1.047^(20-dataB$temp)
      #Fit R~P regression to the bth boostrapped data set
      lm1 <- lm(R20~P20)
      #Save results from the bth regression
      regrCoefs[b,] <- coef(lm1)
    } #end loop over bootstrap samples
    assign(paste(lakes[i],trt[j],'bootR',sep='_'),regrCoefs)
    curBootOut<-rbind(curBootOut,dataBoot)
  }
  assign(paste(lakes[i],'dataBoot',sep='_'),curBootOut)
}
