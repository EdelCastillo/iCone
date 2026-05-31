

# statisticalQuality()
# Determines the deviation of the mzTest vector from the mzRef vector
# Returns statistical values of the deviations
# Resolution is used to identify repeated masses (sharing a reference mass)
# mzRef and mzTest are derived from getPeakMatrix()
#' @export
statisticalQuality<-function(mzRef, mzTest, resolution)
{
  txt1=sprintf("    refSize=%d  testSize=%d", length(mzRef), length(mzTest))
  
  Mx=fitQuality(mzRef, mzTest, resolution)
  m=mean(Mx[,3])
  sigma=sd(Mx[,3])
  md=median(Mx[,3])
  logic=Mx[,5]==1
  repes=length((Mx[,5])[logic])
  
  txt2=sprintf("    mean=%.3f  sigma=%.3f  median=%.3f  repes:%.0f (%.1f%%)", m, sigma, md, repes, 100*repes/length(mzTest));
  cat(txt1, txt2)
}

# fitQuality()
fitQuality<-function(refCentroids, testCentroids, referenceMassResolution, histo=FALSE)
{
  fail<-matrix(nrow=2, ncol=5);
  testLength=length(testCentroids);
  refLength=length(refCentroids);

  deviation <-matrix(nrow = testLength, ncol = 5);
  deviation[,1]=rep(0, times=testLength);
  deviation[,2]=rep(0, times=testLength);
  deviation[,3]=rep(0, times=testLength);
  deviation[,4]=rep(0, times=testLength);
  deviation[,5]=rep(0, times=testLength);
  
  for(iPk in 1:testLength) #para cada pico del test
  {
    offset=0;
    testPPM=referenceMassResolution;
     testMass=testCentroids[iPk];

    retMass<-nearestValue(testMass, refCentroids);
    if(retMass==-1) {offset=-1;}
    else {offset<-abs(testMass-retMass)}
    ppm<-1e6*offset/testMass;
    
    deviation[iPk, 1]=testMass;
    deviation[iPk, 2]=retMass;
    deviation[iPk, 3]=ppm;
    
    if(ppm>1.5*testPPM) #desviation > 1.5*mass_resolution (low resolution)
    {deviation[iPk, 4]=2;}
    else if(ppm>testPPM) #desviation >1 && <= 1.5*mass_resolución 
    {deviation[iPk, 4]=1;}
    else                #desviation <=1*mass_resolución 
    {deviation[iPk, 4]=0;}
    
    deviation[iPk, 5]=0;
    if(iPk>1)
      if(deviation[iPk, 2]==deviation[iPk-1, 2]) #share the same reference mass.
      {deviation[iPk, 5]=1;}
    
  }
  if(histo==TRUE)
  {
    hist(deviation[, 3], main="Histogram of rSirem deviations", xlab="ppm", ylab="Frequency");
    #  legend("topright", legend="SNR=1");
  }
  colnames(deviation)<-c("mzTest", "mzRef", "ppm", "maxDev", "repe")
  return(deviation);
}


#' nearestValue()
#' Return the nearest value in data
#' Successive approximation algorithm.
#' 
#' @param value -> reference value
#' @param data  -> array of sort values
#'
#' @return nearest value in data to value; -1 if value es out of range
#' export
#' 
nearestValue<-function(value, data)
{
  indexLow<-1;
  indexHigh<-length(data);
  
  if(indexHigh==indexLow) return(data[1]);
  if(indexHigh==indexLow+1)
  {
    if(indexLow==-1)indexLow=0;
    if(value-data[indexLow] <= data[indexHigh]-value) {return(data[indexLow]);}
    else {return(data[indexHigh]);}
  }
  
  if(indexLow!=-1 & value==data[indexLow])       return(value);
  if(indexLow!=-1 & value<data[indexLow])       {return(data[indexLow]) ;}
  else if(value>data[indexHigh]) {return(data[indexHigh]);}
  else if(value==data[indexHigh]) return(value);
  
  while(1)
  {
    indexCenter<-round((indexHigh+indexLow)/2);
    if(value==data[indexCenter]) return(value);
    if(value<data[indexCenter]) {indexHigh<-indexCenter; }
    else {indexLow <-indexCenter;}
    if(indexHigh==indexLow+1)
    {
      if(indexLow!=-1 & value-data[indexLow] <= data[indexHigh]-value) {return(data[indexLow]);}
      else {return(data[indexHigh]);}
    }
  }
}

