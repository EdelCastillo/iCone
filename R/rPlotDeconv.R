#########################################################################
#     iCone - R package for MSI data deconvolution
#     Copyright (C) november 2025, Esteban del Castillo Pérez
#     esteban.delcastillo@urv.cat
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

# functions for graphical presentation.
#---------------------------------------------------------------------------


#' @name rPlotGaussianSpectrum
#' @title Presents an image with the Gaussians from a spectra fragment,
#' along with the curve resulting from the sum of all of them (in red).
#' The blue circles indicate the original intensity data to be adjusted.
#' The drawInfo parameters come from getPixelGaussians() or getGaussiansFromSpectrum()
#'
#' @param drawInfo$gaussians: matrix
#'  column 1 -> mean values of each Gaussian.
#'  column 2 -> standard deviation.
#'  column 3 -> factor associated with each Gaussian.
#' @param drawInfo$mass: vector (in Daltons).
#' @param drawInfo$intensity: vector of intensities. It can be null or non-existent.
#' @param sum: if TRUE, the sum of Gaussian is displayed.
#' @return -1 if KO; 0 if OK
#' @export
#' 
rPlotGaussianSpectrum<-function(drawInfo, minMass=0, maxMass=0, sum=TRUE)
{
  if(length(drawInfo$mass)<2) {print("no data to draw"); return(-1);}
  
  nPoints<-length(drawInfo$mass); #number of original data on each axis.
  if(minMass==0)
    {minX<-drawInfo$mass[1];}         #minimum value on the X axis.
  else
    {minX=minMass}
  if(maxMass==0)
    {maxX<-drawInfo$mass[nPoints];}   #maximum value on the X axis.
  else
  {maxX=maxMass}
  
  logic=drawInfo$mass>=minX & drawInfo$mass<=maxX
  mass=drawInfo$mass[logic]
  intensity=c(0)
  if(exists("intensity", where=drawInfo) && length(intensity)>1)
    intensity=drawInfo$intensity[logic]
  nPoints=length(mass)
  
  if(nPoints<2) {print("no data to draw."); return(-1);}
  
  logicG=drawInfo$gaussians[,1]>=minX & drawInfo$gaussians[,1]<=maxX
  gMean =drawInfo$gaussians[,1][logicG]
  if(is.null(gMean)) {print("There are no Gaussian within the interval."); return(-1);}
  
  gSigma=drawInfo$gaussians[,2][logicG]
  gInt  =drawInfo$gaussians[,3][logicG]
  nGauss=length(gMean)
  maxIntG=max(gInt)
  maxIntP=max(intensity)
  k=maxIntG/maxIntP;
  
#  nGauss<-nrow(drawInfo$gaussians); #number of Gaussians to represent.
  maxGaussValue=0;
  maxGaussIndex=0;
  minSigma<-min(gSigma); #required to determine the increments in X.
  
  Y=c(0, max(gInt));
  X=c(min(mass), max(mass))
  
  #Blue circles are drawn corresponding to the original data.
  plot(X, Y, type="p",col="white",lwd=.01, main="spectrum", xlab="mz(Daltons)", 
       ylab="Intensity", las=1, col.axis="black");  
  totalY=0; #used to represent the sum curve.
  if(length(intensity)>1)
    lines(mass, intensity, type="p",col="blue",lwd=1)
  
  deltaX<-minSigma/10; #5 points fit within the lowest sigma (resolution in X).
  if(deltaX<1e-6) deltaX=1e-6;
  
  X<-seq(minX, maxX, deltaX); # X axis.
  for(i in 1 : nGauss)
    {
    mean =gMean[i]; #gaussians parameters
    sigma=gSigma[i];
    value=gInt[i];
      
    Y=value*exp(-((X-mean)^2)/(2*(sigma^2))); #magnitude of the current Gaussian.
    totalY=totalY+Y; #the sum is updated.
    color=((i-1)%%8)+1;
    lines(X, Y, type="l", col=color, lwd=1); #the current Gaussian is drawn.
  }
  if(sum==TRUE)
    lines(X, totalY, type="l", col="red", lwd=2); #the sum curve is drawn.
  return(0);
}


#' @name rPlotGaussianSpectrum2
#' @title Presents an image with the Gaussians of two spectrum fragments,
#' along with the curve resulting from the sum of them (in red).
#' The blue circles indicate the original concentration data to be adjusted.
#' The drawInfo parameters come from getPixelGaussians() or getGaussiansFromSpectrum()
#' 
#' @param drawInfo$gaussians: matrix
#'  column 1 -> mean values of each Gaussian.
#'  column 2 -> standard deviation.
#'  column 3 -> factor associated with each Gaussian.
#' @param drawInfo$mass: vector (in Daltons).
#' @param drawInfo$intensity: vector of intensities
#' @param minMass: minimun mz
#' @param maxMass: maximun mz
#'
#' @return -1 if KO; 0 if OK
#' @export
#' 
rPlotGaussianSpectrum2<-function(drawInfo1, drawInfo2, minMass=0, maxMass=0)
{
  #Los tres ejes X se unifican 
  minMass1=min(drawInfo1$mass);
  minMass2=min(drawInfo2$mass);
  maxMass1=max(drawInfo1$mass);
  maxMass2=max(drawInfo2$mass);
  
  if(minMass==0) {minMass=min(c(minMass1, minMass2));} #eje completo
  else 
  {
    minMass=minMass; 
    if(minMass<minMass1 & minMass<minMass2)
    {minMass=min(c(minMass1, minMass2));
    print("warning: the low mass was update");}
  }
  if(maxMass==0) {maxMass=max(c(maxMass1, maxMass2));}
  else 
  {
    maxMass=maxMass; 
    if(maxMass>maxMass1 & maxMass>maxMass2)
    {maxMass=max(c(maxMass1, maxMass2));
    print("warning: the high mass was update");}
  }
  
  if(maxMass<=minMass) {print("mass range is wrong"); return(-1);}
  
  #mz extremas de cada eje
  x1Range=c(minMass1, maxMass1);
  x2Range=c(minMass2, maxMass2);
  
  #se extrae un subconjunto de datos delimitado por minMz y maxMz
  #para arg1
  gauss1_logic=drawInfo1$gaussians[,1]>=minMass & drawInfo1$gaussians[,1]<=maxMass;
  mean1_sub=(drawInfo1$gaussians[,1])[gauss1_logic];
  sigma1_sub=(drawInfo1$gaussians[,2])[gauss1_logic];
  value1_sub=(drawInfo1$gaussians[,3])[gauss1_logic];
  axis1_logic=drawInfo1$mass>=minMass & drawInfo1$mass<=maxMass;
  X1_axis=drawInfo1$mass[axis1_logic];
  Y1_axis=drawInfo1$intensity[axis1_logic];
  if(length(mean1_sub)==0) {print("warning: no peaks in arg1 for this mass range"); return(-1);}
  maxRelVal1=100*max(value1_sub)/max(drawInfo1$gaussians[,3]);
  
  #para arg2
  gauss2_logic=drawInfo2$gaussians[,1]>=minMass & drawInfo2$gaussians[,1]<=maxMass;
  mean2_sub=(drawInfo2$gaussians[,1])[gauss2_logic];
  sigma2_sub=(drawInfo2$gaussians[,2])[gauss2_logic];
  value2_sub=(drawInfo2$gaussians[,3])[gauss2_logic];
  axis2_logic=drawInfo2$mass>=minMass & drawInfo2$mass<=maxMass;
  X2_axis=drawInfo2$mass[axis2_logic];
  Y2_axis=drawInfo2$intensity[axis2_logic];
  if(length(mean2_sub)==0) {print("warning: no peaks in arg2 for this mass range"); return(-1);}
  maxRelVal2=100*max(value2_sub)/max(drawInfo2$gaussians[,3]);
  
  #mensaje informativo
  msg1<-sprintf("arg1 max relative value in range:%.1f %%", maxRelVal1)  
  msg2<-sprintf("arg2 max relative value in range:%.1f %%", maxRelVal2) 
  print(msg1);
  print(msg2);
  
  #normalización de los datos en rango de 0..100
  value1_sub=value1_sub*100/max(Y1_axis);
  value2_sub=value2_sub*100/max(Y2_axis);
  Y1Factor=100/max(Y1_axis);
  Y2Factor=100/max(Y2_axis);
  Y1_axis=Y1_axis*Y1Factor; #100/max(Y1_axis);
  Y2_axis=Y2_axis*Y2Factor; #100/max(Y2_axis);
  
  X=c(minMass, maxMass);  #extremos del eje X
  Y=c(0, 225);      #extremos del eje Y
  
  #marco de imagen: sólo con eje X
  plot(X, Y, type="p", col="white", cex=0.01, axes = FALSE, main="spectra", xlab="mz(Da)", ylab="relative intensity", las=1, col.axis="black");
  axis(1); ##visualiza el eje X
  legendY1=sprintf("top factor=%.4f", Y1Factor);
  legendY2=sprintf("bottom factor=%.4f", Y2Factor);
  #  print(legendY1, legendY2);
  #  legend("topleft", legend=c(legendY1, legendY2));
  mtext(legendY1, side=3, adj=1);
  mtext(legendY2, side=3, adj=0);
  
  #resolución en X ajustada a 1/5 de la sigma inferior
  minSigma1<-min(sigma1_sub); #required to determine the increments in X.
  minSigma2<-min(sigma2_sub); 
  minSigma=min(c(minSigma1, minSigma2));
  deltaX<-minSigma/5; #5 points fit within the lowest sigma (resolution in X).
  if(deltaX<1e-6) deltaX=1e-6;
  
  #para cada muestra
  for(res in 1:2)
  {
    if(res==1) #low resolution
    {
      mean_t=mean1_sub; #info de gausianas
      sigma_t=sigma1_sub;
      value_t=value1_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y1_axis;        #eje Y
      X=X1_axis;        #eje X
    }
    if(res==2) #median resolution
    {
      mean_t=mean2_sub; #info de gausianas
      sigma_t=sigma2_sub;
      value_t=value2_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y2_axis;        #eje Y
      X=X2_axis;
    }
    #Blue circles are drawn corresponding to the original data.
    points(X, Y+115*(res-1), pch=21, col="blue",lwd=0.75, cex=0.75);  
    
    X<-seq(minMass, maxMass, deltaX); # X axis común
    totalY=0; #used to represent the sum curve.
    
    #para cada gausiana de cada pico
    for(i in 1 : nGauss) 
    {
      mean =mean_t[i]; #gaussians parameters
      sigma=sigma_t[i];
      value=value_t[i];
      
      Y=value*exp(-((X-mean)^2)/(2*(sigma^2))); #magnitude of the current Gaussian.
      totalY=totalY+Y; #the sum is updated.
      color=((i-1)%%8)+1;
      lines(X, Y+115*(res-1), type="l", col=color, lwd=1); #the current Gaussian is drawn.
    }
    lines(X, totalY+115*(res-1), type="l", col="red", lwd=2); #the sum curve is drawn.
  }
  
  #Lineas discontinuas verticales: mean of gaussians from mean spectrum 3
  for(i in 1:length(mean2_sub))
  {
    yHigh=value2_sub[i]+115
    Y<-1:yHigh;
    X<-rep(mean2_sub[i], times=yHigh)
    lines(X, Y, type="l", col="green", lwd=0.75, lty=2); 
  }
  return(0);
}



#' @name rPlotGaussianSpectrum3
#' @title Presents an image with the Gaussians of three spectrum fragments,
#' along with the curve resulting from the sum of them (in red).
#' The blue circles indicate the original intensity data to be adjusted.
#' The drawInfo parameters come from getPixelGaussians() or getGaussiansFromSpectrum()
#' 
#' @param drawInfo$gaussians: matrix
#'  column 1 -> mean values of each Gaussian.
#'  column 2 -> standard deviation.
#'  column 3 -> factor associated with each Gaussian.
#' @param drawInfo$mass: vector (in Daltons).
#' @param drawInfo$intensity: vector of intensities
#' @param minMass: minimun mz
#' @param maxMass: maximun mz
#' @param rMSI2_peaks1: low    resolution peaks from rMSI peak matrix
#' @param rMSI2_peaks1: median resolution peaks from rMSI peak matrix
#' @param rMSI3_peaks1: high   resolution peaks from rMSI peak matrix
#'
#' @return -1 if KO; 0 if OK
#' @export
#' 
rPlotGaussianSpectrum3<-function(drawInfo1, drawInfo2, drawInfo3, minMass=0, maxMass=0, rMSI2_peaks1=-1, rMSI2_peaks2=-1, rMSI2_peaks3=-1)
{
  #Los tres ejes X se unifican 
  mMass=c(0,0,0)
  mMass[1]=min(drawInfo1$mass);
  mMass[2]=min(drawInfo2$mass);
  mMass[3]=min(drawInfo3$mass);
  
  if(minMass==0) {minMass=min(mMass);} #eje completo
  else 
  {
    minMass=minMass; 
    if(minMass<min(mMass))
    {minMass=min(mMass);
    print("warning: the low mass was update");}
  }
  
  mMass[1]=max(drawInfo1$mass);
  mMass[2]=max(drawInfo2$mass);
  mMass[3]=max(drawInfo3$mass);

  if(maxMass==0) {maxMass=max(mMass);}
  else 
  {
    maxMass=maxMass; 
    if(maxMass>max(mMass))
    {maxMass=max(mMass);
    print("warning: the high mass was update");}
  }
  
  if(maxMass<=minMass) {print("mass range is wrong"); return(-1);}
  
 
  #se extrae un subconjunto de datos delimitado por minMz y maxMz
  #para arg1
  gauss1_logic=drawInfo1$gaussians[,1]>=minMass & drawInfo1$gaussians[,1]<=maxMass;
  mean1_sub=(drawInfo1$gaussians[,1])[gauss1_logic];
  sigma1_sub=(drawInfo1$gaussians[,2])[gauss1_logic];
  value1_sub=(drawInfo1$gaussians[,3])[gauss1_logic];
  axis1_logic=drawInfo1$mass>=minMass & drawInfo1$mass<=maxMass;
  X1_axis=drawInfo1$mass[axis1_logic];
  Y1_axis=drawInfo1$intensity[axis1_logic];
  if(length(mean1_sub)==0) {print("warning: no peaks in arg1 for this mass range"); return(-1);}
  maxRelVal1=100*max(value1_sub)/max(drawInfo1$gaussians[,3]);
  
  #para arg2
  gauss2_logic=drawInfo2$gaussians[,1]>=minMass & drawInfo2$gaussians[,1]<=maxMass;
  mean2_sub=(drawInfo2$gaussians[,1])[gauss2_logic];
  sigma2_sub=(drawInfo2$gaussians[,2])[gauss2_logic];
  value2_sub=(drawInfo2$gaussians[,3])[gauss2_logic];
  axis2_logic=drawInfo2$mass>=minMass & drawInfo2$mass<=maxMass;
  X2_axis=drawInfo2$mass[axis2_logic];
  Y2_axis=drawInfo2$intensity[axis2_logic];
  if(length(mean2_sub)==0) {print("warning: no peaks in arg2 for this mass range"); return(-1);}
  maxRelVal2=100*max(value2_sub)/max(drawInfo2$gaussians[,3]);
  
  #para arg3
  gauss3_logic=drawInfo3$gaussians[,1]>=minMass & drawInfo3$gaussians[,1]<=maxMass;
  mean3_sub=(drawInfo3$gaussians[,1])[gauss3_logic];
  sigma3_sub=(drawInfo3$gaussians[,2])[gauss3_logic];
  value3_sub=(drawInfo3$gaussians[,3])[gauss3_logic];
  axis3_logic=drawInfo3$mass>=minMass & drawInfo3$mass<=maxMass;
  X3_axis=drawInfo3$mass[axis3_logic];
  Y3_axis=drawInfo3$intensity[axis3_logic];
  if(length(mean3_sub)==0) {print("warning: no peaks in arg3 for this mass range"); return(-1);}
  maxRelVal3=100*max(value3_sub)/max(drawInfo3$gaussians[,3]);
  
  #mensaje informativo
  msg1<-sprintf("arg1 max relative value in range:%.1f %%", maxRelVal1)  
  msg2<-sprintf("arg2 max relative value in range:%.1f %%", maxRelVal2) 
  msg3<-sprintf("arg3 max relative value in range:%.1f %%", maxRelVal3) 
  print(msg1);
  print(msg2);
  print(msg3);
  
  #normalización de los datos en rango de 0..100
  value1_sub=value1_sub*100/max(Y1_axis);
  value2_sub=value2_sub*100/max(Y2_axis);
  value3_sub=value3_sub*100/max(Y3_axis);
  Y1Factor=100/max(Y1_axis);
  Y2Factor=100/max(Y2_axis);
  Y3Factor=100/max(Y3_axis);
  Y1_axis=Y1_axis*Y1Factor; #100/max(Y1_axis);
  Y2_axis=Y2_axis*Y2Factor; #100/max(Y2_axis);
  Y3_axis=Y3_axis*Y3Factor; #100/max(Y2_axis);
  
  X=c(minMass, maxMass);  #extremos del eje X
#  Y=c(0, 325);      #extremos del eje Y
  Y=c(-30, 325);      #extremos del eje Y
  
  #marco de imagen: sólo con eje X
  plot(X, Y, type="p", col="white", cex=0.01, axes = FALSE, main="Flowchart", 
       xlab="mz(Da)", ylab="relative intensity", las=1, col.axis="black");
  axis(1); ##visualiza el eje X
#  legendY1=sprintf("top factor=%.4f", Y1Factor);
#  legendY2=sprintf("central factor=%.4f", Y2Factor);
#  legendY3=sprintf("bottom factor=%.4f", Y3Factor);
#  mtext(legendY1, side=3, adj=1);
#  mtext(legendY2, side=3);
#  mtext(legendY3, side=3, adj=0);
  
  #resolución en X ajustada a 1/5 de la sigma inferior
  minSigma1<-min(sigma1_sub); #required to determine the increments in X.
  minSigma2<-min(sigma2_sub); 
  minSigma3<-min(sigma3_sub); 
  minSigma=min(c(minSigma1, minSigma2, minSigma3));
  deltaX<-minSigma/5; #5 points fit within the lowest sigma (resolution in X).
  if(deltaX<1e-6) deltaX=1e-6;
  
  colors=c("red", "green", "blue", "cyan")
  
  #para cada muestra
  for(res in 1:3)
  {
    if(res==1) #low resolution
    {
      mean_t=mean1_sub; #info de gausianas
      sigma_t=sigma1_sub;
      value_t=value1_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y1_axis;        #eje Y
      X=X1_axis;        #eje X
    }
    if(res==2) #median resolution
    {
      mean_t=mean2_sub; #info de gausianas
      sigma_t=sigma2_sub;
      value_t=value2_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y2_axis;        #eje Y
      X=X2_axis;
    }
    if(res==3) #high resolution
    {
      mean_t=mean3_sub; #info de gausianas
      sigma_t=sigma3_sub;
      value_t=value3_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y3_axis;        #eje Y
      X=X3_axis;
    }
    #Blue circles are drawn corresponding to the original data.
    points(X, Y+115*(res-1), pch=21, col="bisque4", lwd=0.75, cex=0.75);  
    
    X<-seq(minMass, maxMass, deltaX); # X axis común
    totalY=0; #used to represent the sum curve.
    
    #para cada gausiana de cada pico
    for(i in 1 : nGauss) 
    {
      mean =mean_t[i]; #gaussians parameters
      sigma=sigma_t[i];
      value=value_t[i];
      
      Y=value*exp(-((X-mean)^2)/(2*(sigma^2))); #magnitude of the current Gaussian.
      totalY=totalY+Y; #the sum is updated.
      color=((i-1)%%8)+1;
      lines(X, Y+115*(res-1), type="l", col=color, lwd=1); #the current Gaussian is drawn.
    }
    lines(X, totalY+115*(res-1), type="l", col=colors[res], lwd=2); #the sum curve is drawn.
  ceros=rep(0, times=length(X))
  lines(X, ceros+115*(res-1), type="l", col="black", lwd=1); #the sum curve is drawn.
  }
#return()
  
  #centroides from rMSI2
  centros=rep(0, times=5);
  if(rMSI2_peaks1!=-1 && length(rMSI2_peaks1)>=1)
  {
    logicMass<-rMSI2_peaks1>=minMass & rMSI2_peaks1<=maxMass;
    rMSI2_mass<-rMSI2_peaks1[logicMass];
    lineLength=5;
    Y<-1:lineLength;
    for(i in 1:length(rMSI2_mass))
    {
      X<-rep(rMSI2_mass[i], times=lineLength)
      centros[i]=centros[i]+rMSI2_mass[i];
      lines(X, Y, type="l", col="red", lwd=3)
      lines(X, Y-20, type="l", col="red", lwd=3)
    }
  }
  if(rMSI2_peaks2!=-1 && length(rMSI2_peaks2)>=1)
  {
    logicMass<-rMSI2_peaks2>=minMass & rMSI2_peaks2<=maxMass;
    rMSI2_mass<-rMSI2_peaks2[logicMass];
    for(i in 1:length(rMSI2_mass))
    {
      X<-rep(rMSI2_mass[i], times=lineLength)
      centros[i]=centros[i]+rMSI2_mass[i];
      lines(X, Y+115, type="l", col="green", lwd=3)
      lines(X, Y-20, type="l", col="green", lwd=3)
    }
  }
  if(rMSI2_peaks3!=-1 && length(rMSI2_peaks3)>=1)
  {
    logicMass<-rMSI2_peaks3>=minMass & rMSI2_peaks3<=maxMass;
    rMSI2_mass<-rMSI2_peaks3[logicMass];
    for(i in 1:length(rMSI2_mass))
    {
      X<-rep(rMSI2_mass[i], times=lineLength)
      centros[i]=centros[i]+rMSI2_mass[i];
      lines(X, Y+230, type="l", col="blue", lwd=3)
      lines(X, Y-20, type="l", col="blue", lwd=3)
    }
  }
#  return()
  #Lineas discontinuas verticales: mean of gaussians from mean spectrum 3
  for(i in 1:length(mean3_sub))
  {
    yHigh=value3_sub[i]+230
    Y<-1:yHigh;
    X<-rep(mean3_sub[i], times=yHigh)
    lines(X, Y, type="l", col="bisque4", lwd=0.75, lty=2); 
  }
  if(0)#for(i in 1:length(mean2_sub))
  {
    yHigh=value2_sub[i]+115
    Y<-1:yHigh;
    X<-rep(mean2_sub[i], times=yHigh)
    lines(X, Y, type="l", col="bisque4", lwd=0.75, lty=2); 
  }
  if(0)#for(i in 1:length(mean1_sub))
  {
    yHigh=value1_sub[i]+0
    Y<-1:yHigh;
    X<-rep(mean1_sub[i], times=yHigh)
    lines(X, Y, type="l", col="bisque4", lwd=0.75, lty=2); 
  }
#  return()
  if(0)#if(length(rMSI2_peaks1)>=1 & length(rMSI2_peaks2)>=1 & length(rMSI2_peaks3)>=1)
  {  
    centros=centros/3
    Y<-1:(lineLength+3);
    for(i in 1:length(centros))
    {
      X<-rep(centros[i], times=lineLength+3)
      lines(X, Y-40, type="l", col="cyan", lwd=3)
    }
  }
  
  return(0);
}

#' @name rPlotSpectrum
#' @title Presents an image with the info of two spectrum fragments,
#' @param List: values come from getAverageSpectrum() or getAverageGaussianSpectrum()
#'  values$averageMz: vector of mases(in Daltons).
#'  values$averageIntensity: vector of intesities
#' @param mzLow : minimun mz
#' @param mzHigh: maximun mz
#'
#' @export
#' 
rPlotSpectrum<-function(values, mzLow, mzHigh)
{
  logic=values$averageMz>=mzLow & values$averageMz<=mzHigh
  X=values$averageMz[logic];
  Y=values$averageIntensity[logic]
  #marco de imagen: sólo con eje X
  plot(X, Y, type="l", col="green", cex=0.01, axes = FALSE, main="Mean spectrum", 
       xlab="mz(Da)", ylab="Intensity", las=1, col.axis="black");
  axis(1); ##visualiza el eje X
  axis(2)
#  lines(mass, intensity, type="p",col="blue",lwd=1)
}

#' @name rPlotIon
#' @title visualize an ion.
#' @param pMatrix: peak matrix from getPeakMatrix().
#' @param mz: desired mass (Da).
#' @comment very poor solution.
#' @export
rPlotIon<-function(pMatrix, mz)
{
  colors=c("red3", "red",  "orange1", "yellow", "forestgreen", "olivedrab", "yellowgreen", "cyan", "blue")
  ionIndex=which.min(abs(pMatrix$mass - mz))

  intensity=pMatrix$peakMatrix[,ionIndex]
  coord=pMatrix$coordinates
  #maxInt=max(intensity)
  maxInt=mean(intensity)+3*sd(intensity)
  co=(maxInt/9)
  X=coord[,1]
  Y=coord[,2]
  data=abs(intensity/co)
  plot(X, Y,  type="p", col="white")
  for (i in 1:length(pMatrix$coordinates[,1]))
  {
    A=intensity[i]/co;
    A=9-(A%%9)
    if(A>9) A=9;
    lines(X[i], Y[i], type="p", cex=0.1, col=colors[A])
  }
}


