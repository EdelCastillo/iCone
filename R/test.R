

# statisticalQuality() 
# It reports statistical values regarding the deviation between the two m/z vectors being compared. 
# This deviation is determined from the nearest centroids in the high-resolution sample to those in the low-resolution sample.
# arguments: 
# mzRef      is the high-resolution centroid vector.
# mzTest     is the low-resolution centroid vector.
# tolerance  indicates the minimum distance for the result to be considered a false positive (expressed in ppm).

# report: 
# refSize    size of the high-resolution centroid vector.
# testSize   size of the low-resolution centroid vector.
# mean       mean of the deviations between centroids.
# sigma      standard deviation of the deviations between centroids.
# median     median of the deviations between centroids.
# FP         False Positive: number of low-resolution centroids with a deviation exceeding 1.5 times the tolerance relative to the nearest high-resolution centroid.
#' @export
statisticalQuality<-function(mzRef, mzTest, tolerance)
{
  txt1=sprintf("refSize=%d  testSize=%d\n", length(mzRef), length(mzTest))
  cat(txt1)
  
  Mx=fitQuality(mzRef, mzTest, tolerance)
  v1=Mx[,3]
  logic=Mx[,4]<2 
  v2=v1[logic] #good data
  
  m=mean(v1)
  sigma=sd(v1)
  md=median(v1)
  logic=Mx[,4]==2
  FP=length((Mx[,4])[logic])
  
  txt2=sprintf("mean=%9.4f  sigma=%9.4f  median=%9.4f  FP:%.0f (%.1f%%)\n", m, sigma, md, FP, 100*FP/length(mzTest));
  cat("        All data:", txt2)

  m=mean(v2)
  sigma=sd(v2)
  md=median(v2)
  logic=Mx[,4]==2
  FP=length((Mx[,4])[logic])
  
  txt2=sprintf("mean=%9.4f  sigma=%9.4f  median=%9.4f  FP:%.0f (%.1f%%)\n", m, sigma, md, FP, 100*FP/length(mzTest));
  cat("Without bad data:", txt2)
  return(Mx)
}

# fitQuality()
fitQuality<-function(refCentroids, testCentroids, tolerance, histo=FALSE)
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
    testPPM=tolerance;
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
