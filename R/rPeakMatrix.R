#*********************************************************************************
#     iCone - R package for MSI data processing
#     Copyright (C) 2026 Esteban del Castillo Pérez (esteban.delcastillo@urv.cat)
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
#********************************************************************************/
 

#' @name getPeakMatrix
#' @title obtains the peak matrix from imzML files.
#' @param data_file: A list of the absolute paths of the files to be processed is required. 
#'                   Alternatively, a single filename with the 'txt' extension containing the absolute paths of all files to be processed is acceptable. 
#'                   In file.txt, lines are considered comments if they begin with the '#' character or the blank character.
#'                   Only files in 'imzML' format are recognized.
#'                   The attached binary file, with the 'ibd' extension, must be in the same directory that 'imzML' files.
#' @param params  
#'      "massResolution": desired mass resolution.
#'                 "SNR": signal-to-noise ratio (by defect=1)
#'         "noiseMethod": method for estimating noise (by defect="estnoise_mad").
#'    "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered (by defect=1).
#'         "linkedPeaks": two peaks are considered linked if they are closer than the given standard deviation (by defect=3).   
#' @param initMass:  Initial  mass to consider. By default, the minimum value from the entire range of masses is used.
#' @param finalMass: Final    mass to consider. By default, the maximum value from the entire range of masses is used.
#' @param pxList:    List of pixels. First pixel=1. By default everyone.
#' @param nThreads:  Number of threads for parallel processing (by default maxCores-1)
#' @param imzMLChecksum: if the binary file checksum must be verified, it can be disabled for convenice with really big files.
#' @param fixBrokenUUID: set to FALSE by default to automatically fix an uuid mismatch between the ibd and the imzML files (a warning message will be raised).
#'
#' @return a list: 
#'     peakMatrix: Matrix of peak (centroids) rows = pixels, columns = intensity of each pixel.
#'           mass: Vector with the masses associated with each column of peakMatrix.
#' massResolution: mass resolution associated with each mass.
#'  pixelsSupport: Vector with the number of pixels in each column with non-zero intensity. 
#'    coordinates: Matrix with pixel coordinates (X/Y).
#'   pixelsSample: Vector with the number of pixels in each of the samples. In the peakMatrix and coordinates they appear in the same order.
#'   pixelSize_um: Vector with the size of the pixels (spatial resolution) in micrometers.
#'
#' @export
getPeakMatrix<-function(data_file,
                   params,
                   initMass=0,
                   finalMass=0,
                   pxList=c(0),
                   nThreads=0,
                   imzMLChecksum = F, 
                   fixBrokenUUID = F)
{
  if(length(data_file)==0)
  {
    stop("No information is provided about the files to be processed.\n")
  }

  #parameters control 
  if(!(exists("SNR", where=params)))
  {
    cat("warning: by default, SNR parameter will be 1.\n")
    params=c(params, "SNR"=1)
  }
  if(!(exists("minPixelsSupport", where=params)))
  {
    cat("warning: by default, minPixelsSupport parameter will be 1%.\n")
    params=c(params, "minPixelsSupport"=1)
  }
  if(!(exists("noiseMethod", where=params)))
  {
    cat("warning: by default, noiseMethod parameter will be MAD type.\n")
    params=c(params, "noiseMethod"="estnoise_mad")
  }
  if(!(exists("linkedPeaks", where=params)))
  {
    cat("warning: by default, linkedPeaks parameter will be 3 sd.\n") 
    params=c(params, "linkedPeaks"=3)
  }
  if(!(exists("massResolution", where=params)))
  {
    cat("ERROR: massResolution parameter is required.\n") 
    return (0)
  }
  if(params$massResolution<=0)
  {
    cat("ERROR: massResolution must be greater than zero.\n") 
    return (0)
  }

  imgData <- NULL
  pt<-proc.time()
  
  samples=1; #default number of samples to be analyzed.
  totalPixels=0;
  samples=c();
  tmpFile=c();
  
  if(length(data_file)==1)
  {
    fileExtension <- unlist(strsplit(basename(data_file[1]), split = "\\."))
    fileExtension <- as.character(fileExtension[length(fileExtension)])
    if( fileExtension == "txt")
      {  
      fileCon <- file(data_file[1], open="r") # We opened the connection.
      
      tmpFile <- readLines(fileCon) #The file names are loaded line by line.
      for(i in 1:length(tmpFile))
      {
        if(substr(tmpFile[i], 1, 1)=='#' || substr(tmpFile[i], 1, 1)==' ' || nchar(tmpFile[i])<=2) {next;}
        if(!file.exists(tmpFile[i]))
        {
          cat(sprintf("Warning: %s file not found\n", tmpFile[i]))
          next
        }
        else
        {
          samples=c(samples, tmpFile[i])
        }
      }
      close(fileCon);
      }
    else
    {
      samples=data_file[1]; #unique sample.
    }
  }
  else if(length(data_file)>1)
  {
    for(i in 1:length(data_file))
    {
    samples=c(samples, data_file[i]);
    }
  }
  
  
  nSamples=length(samples);
  if(nSamples<=0)
  {
    cat("warning: there is no valid file.\n")
    return(0);
  }
  else if(nSamples>150)
  {
    cat("Warning: the maximum number of samples allowed has been reached. It is limited to 150.\n)");
    nSamples=150;
  }
  
  baseDir=rGetDirectory(data_file[1])
  
  fileExtension <- unlist(strsplit(basename(data_file[1]), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension != "bin")
  {
    #Deletion of temporary files, if any. 
    #These files should not contain any content since they are cumulative.
    fileA=paste0(baseDir, "tmpGaussians.bin")
    if(file.exists(fileA)) file.remove(fileA)
    fileB=paste0(baseDir, "tmpPixelsCoordinates.bin")
    if(file.exists(fileB)) file.remove(fileB)
  }
  
  minMass=1e32
  maxMass=0
  lowMz=0
  highMz=0
  totalSamples=0
  pixelSize=c()
  hit=FALSE;
  
  for( sample in 1 : nSamples)
    {
    fileExtension <- unlist(strsplit(basename(samples[sample]), split = "\\."))
    fileExtension <- as.character(fileExtension[length(fileExtension)])
    cat(sprintf("\nSample file name: %s\n", samples[sample]))
    if( fileExtension == "imzML")
    {
      if(!file.exists(samples[sample]))
      {
        cat("Error: file not found\n")
        return (0)
      }
      
      #capturing information from an imzML file.
      in_img <- import_imzML(path.expand(samples[sample]),  fun_progress = NULL, fun_text = NULL, 
                             close_signal = NULL, verifyChecksum = imzMLChecksum, subImg_rename = NULL, 
                             subImg_Coords = NULL, fixBrokenUUID = fixBrokenUUID)
      
      #binary data file name.
      file=path.expand(file.path(in_img$data$path, paste0(in_img$data$imzML$file, ".ibd")));
      if(!file.exists(file))
      {
        cat(sprintf("Error: %s file not found\n", file))
        return (0)
      }
      pixelSize=c(pixelSize, in_img$pixel_size_um)
      
      #dimensioned (minPixel, maxPixel; minMz, maxMz)
      basicInfo=rGetBasicInfo(file, in_img$data$imzML, pxList-1)
      if(initMass ==0 || initMass<basicInfo$minMz  || initMass>basicInfo$maxMz)   {lowMz =basicInfo$minMz;}
      else {lowMz=initMass}
      if(finalMass==0 || finalMass>basicInfo$maxMz || finalMass<basicInfo$minMz)  {highMz=basicInfo$maxMz;}
      else {highMz=finalMass}
      if(lowMz <minMass) minMass=lowMz
      if(highMz>maxMass) maxMass=highMz
      
      if(lowMz>highMz)
      {
        cat("warning: the lower mass exceeds the upper mass.\n")
        return(0)
      }
      
     #imzML file to Gaussians (for each sample).
     nPixels=rawToGaussiansR(baseDir, file, in_img$data$imzML, params, lowMz, highMz, pxList-1, nThreads);
     totalPixels=totalPixels+nPixels;
     totalSamples=totalSamples+1;
     hit=TRUE;
     
     fileA=paste0(baseDir, "tmpMassRange.bin")
     rSaveMassRange(fileA, minMass, maxMass); #save the temporary file with the peaks converted to Gaussians.
     
     gc();
     rm(in_img, basicInfo)
    }
  else if (fileExtension == "bin") #allows execution to continue at step 2 (Gaussians to centroids)
    {
    fileA=paste0(baseDir, "tmpGaussians.bin")
    cat(sprintf("The information contained in the %s file will be used.\nThis eliminates phase 1 of the processing.\n", fileA))
    hit=TRUE;
    totalSamples=1;
    break;
    }
  else 
   {
    cat(sprintf("The file extension must be imzML or bin.\n"))
     hit=FALSE;
    }
}
  gc()  
  if(totalSamples==0) 
  {
    cat("warning: there is no valid file.\n")
    return(0);
  }
  else if(totalSamples>1) 
  {
    cat(sprintf("\t\t\tunified mass range(Da):%.4f to %.4f\n", minMass, maxMass))
  }
  if(!hit) return(0);
  
  #rm(list = ls())
  gc()
  fileA=paste0(baseDir, "tmpMassRange.bin")
  massRange=rLoadMassRange(fileA); #mass range to file
  
  #step 2 and 3 (gaussians to centroids). Peak matrix to file tmpPeakMatrix.bin
  peakMatrix=peakMatrixR(baseDir, params, massRange[1], massRange[2], nThreads);

  gc() 
  
  #The peak array is captured in pkMatrix object from the tmpPeakMatrix.bin file.
  fileA=paste0(baseDir, "tmpPeakMatrix.bin")
  metaInfo=rGetMetaDataFromFile(fileA) #generic information from peak matrix
  pkMatrix=matrix(ncol=metaInfo$nIons, nrow=metaInfo$totalPx);
  mass=vector(length = metaInfo$nIons)
  massResolution=vector(length = metaInfo$nIons)
  pxSupport=vector(length = metaInfo$nIons)
  
  #load ion by ion (This reduces the memory required)
  for(ion in 1:metaInfo$nIons)
  {
    tmpIon=rGetIntensityFromFile(fileA, ion)
    mass[ion]=tmpIon[1];
    massResolution[ion]=tmpIon[2];
    pxSupport[ion]=tmpIon[3];
    pkMatrix[,ion]=tmpIon[4:length(tmpIon)]
  }
  
  fileA=paste0(baseDir, "tmpPixelsCoordinates.bin") 
  coordinates=rGetCoordinatesFromFile(fileA, 0); #load coordinates for all samples
  peakMatrix <- list("peakMatrix" = pkMatrix, "mass" = mass, "massResolution" = massResolution,
                   "pixelsSupport"=pxSupport, "coordinates"=coordinates, 
                   "pixelsSample"=metaInfo$pixelsSample, "pixelSize_um"=pixelSize)  
  pt<-proc.time() - pt
  display_processing_time(pt, "Data processing time")
  return(peakMatrix)
}

#//////////////////////////////////////////////////////////////

#' @name getPxGaussians
#' @title It retrieves information from a single pixel.
#' @param data_file: absolute reference to the file with the imzML extension.
#' @param params  
#'      "massResolution": mass resolution with which the spectra were acquired (mz/deltaMz).
#'                 "SNR": signal-to-noise ratio
#'         "noiseMethod": method for estimating noise.
#' @param initMass:  Initial  mass to consider. By default, the minimum value from the entire range of masses is used.
#' @param finalMass: Final    mass to consider. By default, the maximum value from the entire range of masses is used.
#' @param pixel:    spectrum to evaluate. First pixel=1
#' @param imzMLChecksum: if the binary file checksum must be verified, it can be disabled for convenice with really big files.
#' @param fixBrokenUUID: set to FALSE by default to automatically fix an uuid mismatch between the ibd and the imzML files (a warning message will be raised).
#'
#' @return a list: 
#'    gaussians: parameters of each Gaussian.
#'         mass: Vector with raw masses.
#'    intensity: vector with raw intensities;
#'          SNR: signal-to-noise ratio for each point of mass.
#'        noise: noise estimation.
#' @export
getPixelGaussians<-function(data_file,
                        params=c(),
                        initMass=0,
                        finalMass=0,
                        pixel=1,
                        imzMLChecksum = F, 
                        fixBrokenUUID = F)
{
  if(!file.exists(data_file))
  {
    stop("File not found\n")
  }
  if(length(params)==0)
  {
    params=list("SNR"=1, "noiseMethod"="estnoise_mad", "minPixelsSupport"=1, "linkedPeaks"=3, "massResolution"=30000)
    cat("defect parameters: SNR=1, noiseMethod=estnoise_mad, minPixelsSupport=1, linkedPeaks=3, massResolution=30000\n")
  }
  else
  {
    #control de parámetros
    if(!(exists("SNR", where=params)))
    {
      cat("warning: by default, SNR parameter will be 1.\n")
      params=c(params, "SNR"=1)
    }
    if(!(exists("minPixelsSupport", where=params)))
    {
      cat("warning: by default, minPixelsSupport parameter will be 1%.\n")
      params=c(params, "minPixelsSupport"=1)
    }
    if(!(exists("noiseMethod", where=params)))
    {
      cat("warning: by default, noiseMethod will be MAD type.\n")
      params=c(params, "noiseMethod"="estnoise_mad")
    }
    if(!(exists("massResolution", where=params)))
    {
      cat("ERROR: massResolution parameter is required.\n") 
      return (0)
    }
    if(params$massResolution<=0)
    {
      cat("ERROR: massResolution must be greater than zero.\n") 
      return (0)
    }
  }
  imgData <- NULL
  

  pt<-proc.time()
  
  fileExtension <- unlist(strsplit(basename(data_file), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension == "imzML")
  {
    #capturing information from an imzML file.
    in_img <- import_imzML(path.expand(data_file),  fun_progress = NULL, fun_text = NULL, 
                           close_signal = NULL, verifyChecksum = imzMLChecksum, subImg_rename = NULL, 
                           subImg_Coords = NULL, fixBrokenUUID = fixBrokenUUID)
    
    #binary data file name.
    file=path.expand(file.path(in_img$data$path, paste0(in_img$data$imzML$file, ".ibd")));
    
    #dimensioned (lowMz, highMz)
    basicInfo=rGetBasicInfo(file, in_img$data$imzML, pixel-1)
    if(initMass ==0 || initMass<basicInfo$minMz  || initMass>basicInfo$maxMz)   {lowMz =basicInfo$minMz;}
    else {lowMz=initMass}
    if(finalMass==0 || finalMass>basicInfo$maxMz || finalMass<basicInfo$minMz)  {highMz=basicInfo$maxMz;}
    else {highMz=finalMass}
    
    if(lowMz>highMz)
    {
      cat("warning: the lower mass exceeds the upper mass.\n")
      return(0)
    }
    
    if(pixel<1) 
    {
      cat("warning: pixel must be greater than zero. Thus, pixels will be 1.\n")
      pixel=1;
    }
    #info about pixel (the first pixel is zero)
    pxGauss<-rGetPixelGaussians(file, in_img$data$imzML, params, lowMz, highMz, pixel-1);

    pt<-proc.time() - pt
    display_processing_time(pt, "Data processing time")
    return(pxGauss)
  }
  cat("warning: imzML file type is required.\n" )
}

#//////////////////////////////////////////////////////////////

#' @name getAverageGaussianSpectrum
#' @title It obtains the average value of the Gaussian from all data into .imzML file.
#' @param data_file: absolute reference to the file with the imzML extension.
#' @param params  
#'      "massResolution": mass resolution with which the spectra were acquired (mz/deltaMz).
#'                 "SNR": signal-to-noise ratio
#'         "noiseMethod": method for estimating noise.
#' @param initMass:  Initial  mass to consider. By default, the minimum value from the entire range of masses is used.
#' @param finalMass: Final    mass to consider. By default, the maximum value from the entire range of masses is used.
#' @param pxList:    list of pixels. First pixel=1. By default everyone.
#' @param overSampling: interval between points on the mass axis = massResolution/overSampling.
#' @param nThreads:  number of threads
#' @param imzMLChecksum: if the binary file checksum must be verified, it can be disabled for convenice with really big files.
#' @param fixBrokenUUID: set to FALSE by default to automatically fix an uuid mismatch between the ibd and the imzML files (a warning message will be raised).
#'
#' @return a list: averageMz and averageIntensity
#'            averageMz: array of masses at intervals of 1/4 of the resolution
#'     averageIntensity: array of average Gaussians values 
#' @export
getAverageGaussianSpectrum<-function(data_file,
                         params=c(),
                         initMass=0,
                         finalMass=0,
                         pxList=c(0),
                         overSampling=4,
                         nThreads=0,
                         imzMLChecksum = F, 
                         fixBrokenUUID = F)
{
  if(!file.exists(data_file))
  {
    stop("File not found\n")
  }
  
  if(length(params)==0)
  {
    params=list("SNR"=1, "noiseMethod"="estnoise_mad", "minPixelsSupport"=1, "linkedPeaks"=3, "massResolution"=30000)
    cat("defect parameters: SNR=1, noiseMethod=estnoise_mad, minPixelsSupport=1, linkedPeaks=3, massResolution=30000\n")
  }
  else
  {
    #control de parámetros
    if(!(exists("SNR", where=params)))
    {
      cat("warning: by default, SNR parameter will be 1.\n")
      params=c(params, "SNR"=1)
    }
    if(!(exists("minPixelsSupport", where=params)))
    {
      cat("warning: by default, minPixelsSupport parameter will be 1%.\n")
      params=c(params, "minPixelsSupport"=1)
    }
    if(!(exists("noiseMethod", where=params)))
    {
      cat("warning: by default, noiseMethod will be MAD type.\n")
      params=c(params, "noiseMethod"="estnoise_mad")
    }
    if(!(exists("linkedPeaks", where=params)))
    {
      #cat("warning: by default, linkedPeaks parameter will be 3 sd.\n") 
      params=c(params, "linkedPeaks"=3)
    }
    if(!(exists("massResolution", where=params)))
    {
      cat("ERROR: massResolution parameter is required.\n") 
      return (0)
    }
  }
  if(params$massResolution<=0)
  {
    cat("ERROR: massResolution must be greater than zero.\n") 
    return (0)
  }

  imgData <- NULL
  
  pt<-proc.time()
  
  fileExtension <- unlist(strsplit(basename(data_file), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension == "imzML")
  {
    #capturing information from an imzML file.
    in_img <- import_imzML(path.expand(data_file),  fun_progress = NULL, fun_text = NULL, 
                           close_signal = NULL, verifyChecksum = imzMLChecksum, subImg_rename = NULL, 
                           subImg_Coords = NULL, fixBrokenUUID = fixBrokenUUID)
    
    #binary data file name.
    file=path.expand(file.path(in_img$data$path, paste0(in_img$data$imzML$file, ".ibd")));
    
    #dimensioned (lowMz, highMz)
    basicInfo=rGetBasicInfo(file, in_img$data$imzML, pxList-1)
    if(initMass ==0 || initMass<basicInfo$minMz  || initMass>basicInfo$maxMz)   {lowMz =basicInfo$minMz;}
    else {lowMz=initMass}
    if(finalMass==0 || finalMass>basicInfo$maxMz || finalMass<basicInfo$minMz)  {highMz=basicInfo$maxMz;}
    else {highMz=finalMass}
    
    if(lowMz>highMz)
    {
      cat("warning: the lower mass exceeds the upper mass.\n")
      return(0)
    }
    
    avSp=rGetAverageGaussianSpectrum(file, in_img$data$imzML, params, lowMz, highMz, pxList-1, overSampling, nThreads);

    pt<-proc.time() - pt
    display_processing_time(pt, "Data processing time")
    return(avSp)
  }
  cat("warning: imzML file type is required.\n" )
}

#//////////////////////////////////////////////////////////////

#' @name getAverageSpectrum
#' @title It obtains the average value of the intensities from all data into .imzML file.
#'        Noise is not taken into account.
#' @param data_file: absolute reference to the file with the imzML extension.
#' @param initMass:  Initial  mass to consider. By default, the minimum value from the entire range of masses is used.
#' @param finalMass: Final    mass to consider. By default, the maximum value from the entire range of masses is used.
#' @param pxList:    list of pixels. First pixel=1. By default everyone.
#' @param massResolution: mass resolution with which the spectra were acquired (mz/deltaMz).
#' @param overSampling:  interval between points on the mass axis = massResolution/overSampling.
#' @param imzMLChecksum: if the binary file checksum must be verified, it can be disabled for convenice with really big files.
#' @param fixBrokenUUID: set to FALSE by default to automatically fix an uuid mismatch between the ibd and the imzML files (a warning message will be raised).
#'
#' @return a list: averageMz and averageIntensity
#'            averageMz: array of masses at intervals of the resolution/overSampling
#'     averageIntensity: array of average Gaussians values 
#' @export
getAverageSpectrum<-function(data_file,
                             initMass=0,
                             finalMass=0,
                             pxList=c(0),
                             massResolution=30000,
                             overSampling=4,
                             imzMLChecksum = F, 
                             fixBrokenUUID = F)
{
  if(!file.exists(data_file))
  {
    stop("File not found\n")
  }
  
  if(!(exists("params"))) #parameters not considered but must exist.
  {
    params=list("SNR"=1, "noiseMethod"="estnoise_mad", "minPixelsSupport"=1, "linkedPeaks"=3)
    params=c(params, "massResolution"=massResolution)
  }
  
  imgData <- NULL
  
  pt<-proc.time()
  
  baseDir=rGetDirectory(data_file)
  
  fileExtension <- unlist(strsplit(basename(data_file), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension == "imzML")
  {
    #capturing information from an imzML file.
    in_img <- import_imzML(path.expand(data_file),  fun_progress = NULL, fun_text = NULL, 
                           close_signal = NULL, verifyChecksum = imzMLChecksum, subImg_rename = NULL, 
                           subImg_Coords = NULL, fixBrokenUUID = fixBrokenUUID)
    
    #binary data file name.
    file=path.expand(file.path(in_img$data$path, paste0(in_img$data$imzML$file, ".ibd")));
    
    #dimensioned (lowMz, highMz)
    basicInfo=rGetBasicInfo(file, in_img$data$imzML, pxList-1)
    if(initMass ==0 || initMass<basicInfo$minMz  || initMass>basicInfo$maxMz)   {lowMz =basicInfo$minMz;}
    else {lowMz=initMass}
    if(finalMass==0 || finalMass>basicInfo$maxMz || finalMass<basicInfo$minMz)  {highMz=basicInfo$maxMz;}
    else {highMz=finalMass}
    
    if(lowMz>highMz)
    {
      cat("warning: the lower mass exceeds the upper mass.\n")
      return(0)
    }
    
    avSp=rGetAverageSpectrum(file, in_img$data$imzML, params, lowMz, highMz, pxList-1, overSampling);
    
    pt<-proc.time() - pt
    display_processing_time(pt, "Data processing time")
    return(avSp)
  }
  cat("warning: imzML file type is required.\n" )
}

#//////////////////////////////////////////////////////////////

#' @name getGaussiansFromSpectrum
#' @title Form Gaussians over each peak of the given spectrum.
#' @param Intensity:   vector with value associated with each mass point.
#' @param mz:          vector with mass/charge information
#' @param initMass:  Initial  mass to consider. By default, the minimum value from the entire range of masses is used.
#' @param finalMass: Final    mass to consider. By default, the maximum value from the entire range of masses is used.
#' @param SNR:         signal-to-noise ratio
#' @param noiseMethod: method for estimating noise (estnoise_diff, estnoise_sd, estnoise_mad).
#'
#' @return List: gaussians, mass intensity, SNR, noise
#'       gaussians: matrix with parameters for each Gaussian
#'            mass: vector with the masses of the raw spectrum
#'       intensity: intensity associated with each mass of the raw spectrum
#'             SNR: signal-to-noise ratio associated with each mass of the raw spectrum.
#'           noise: noise estimation
#' @export
getGaussiansFromSpectrum<-function(
                              intensity,
                              mz,
                              initMass=0,
                              finalMass=0,
                              SNR=3,
                              noiseMethod="estnoise_mad")
{

  #control de parámetros
    params=list("SNR"=3, "noiseMethod"="estnoise_mad", "massResolution"=30000, "minPixelsSupport"=1, "linkedPeaks"=3)

    params$SNR=SNR
    params$noiseMethod=noiseMethod
  if(length(intensity) != length(mz))
  {
    cat("warning: vectors must have the same length.\n")
    return(0)
  }
    pt<-proc.time()
  
    if(initMass>finalMass)
    {
      cat("warning: the lower mass exceeds the upper mass.\n")
      initMass=finalMass
    }
    if(initMass==0)  initMass =min(mz);
    if(finalMass==0) finalMass=max(mz);
  
    Gauss=rGetGaussiansFromSpectrum(intensity, mz, params, initMass, finalMass);
    
    pt<-proc.time() - pt
    display_processing_time(pt, "Data processing time")
    return(Gauss)
}

#//////////////////////////////////////////////////////////////

#' @name getMetaDataFromFile
#' @title returns generic information about the peak matrix located in a file.
#' @param file -> file name with peak matrix  (tmpPeakMatrix.bin)
#' @return a list:
#'   nSamples     -> number of samples analyzed.
#'   totalPx      -> total number of pixels (cumulative of each sample).
#'   nIons        -> number of columns in the matrix.
#'   pixelsSample -> vector with the pixels in each sample.
#'
#' @export
getMetaDataFromFile<-function(file)
{
  return(rGetMetaDataFromFile(file));
}

#//////////////////////////////////////////////////////////////

#' @name getPeakMatrixFromFile
#' @title return the peak matrix into file
#' @param file -> file name with peak matrix (tmpPeakMatrix.bin)
#' @return a list:
#'     peakMatrix: Matrix of peak (centroids) rows = pixels, columns = intensity of each pixel.
#'           mass: Vector with the masses associated with each column of peakMatrix.
#' massResolution: mass resolution associated with each mass.
#'  pixelsSupport: Vector with the number of pixels in each column with non-zero intensity. 

#' @export
getPeakMatrixFromFile<-function(file)
{
metaInfo=rGetMetaDataFromFile(file) #generic information from peak matrix

pkMatrix=matrix(ncol=metaInfo$nIons, nrow=metaInfo$totalPx);
mass=vector(length = metaInfo$nIons)
massResolution=vector(length = metaInfo$nIons)
pxSupport=vector(length = metaInfo$nIons)

#load ion by ion (This reduces the memory required)
for(ion in 1:metaInfo$nIons)
{
  tmpIon=rGetIntensityFromFile(file, ion)
  mass[ion]=tmpIon[1];
  massResolution[ion]=tmpIon[2];
  pxSupport[ion]=tmpIon[3];
  pkMatrix[,ion]=tmpIon[4:length(tmpIon)]
}
peakMatrix <- list("peakMatrix" = pkMatrix, "mass" = mass, "massResolution" = massResolution,
                   "pixelsSupport"=pxSupport) 
return(peakMatrix)
}

#//////////////////////////////////////////////////////////////

#' @name getMassVectorFromFile()
#' @title returns a vector with all the masses in peak matrix 
#' @param file     -> file name with peak matrix (tmpPeakMatrix.bin)
#' @return mass vector

#' @export
getMassVectorFromFile<-function(file)
{
  return (rGetMassVectorFromFile(file))
}

#//////////////////////////////////////////////////////////////

#' @name getColumnFromFile()
#' @title returns a column information of the peak matrix: 
#' 
#' @param file     -> file name with peak matrix (tmpPeakMatrix.bin)
#' @param mass     -> reference to the desired initial column of the peak matrix (Da).
#' @param sample   -> just download the pixels from this sample.
#'                 if sample < 0, all sample coordinates are returned
#' @return a list:
#'     intensity: vector of intesities 
#'          mass: mass associated with the column of peakMatrix.
#'massResolution: final mass resolution at centroid.
#' pixelsSupport: number of pixels in column with non-zero intensity. 
#' 
#' @export
getColumnFromFile<-function(file, mass, sample)
{
  return (rGetColumnFromFile(file, mass, sample))
}

#//////////////////////////////////////////////////////////////

#' @name getIntensityFromFile()
#' @title returns a column whit pixels intensities from the peak matrix: 
#' 
#' @param file -> file name with peak matrix (tmpPeakMatrix.bin)
#' @param column -> desired column (first = 1)
#' @return a vector with intensity of each pixel.
#' 
#' @export
getIntensityFromFile<-function(file, column)
{
  return (rGetIntensityFromFile(file, column))
}

#//////////////////////////////////////////////////////////////

#' @name getCoordinatesFromFile()
#' @title returns a matrix with the coordinates of all pixels (X/Y).
#' If there are multiple samples, they appear sequentially; that is, the matrix has as many rows 
#' as the cumulative number of pixels in each sample and two columns.
#' @param file   -> file name with pixels coordinates (tmpPixelsCoordinates.bin)
#' @param sample -> just download the pixels from this sample.
#'                  if sample < 0, all sample coordinates are returned
#' @return a matrix with the coordinates (X/Y) of pixels.
#' 
#' @export
getCoordinatesFromFile<-function(file, sample)
{
return(rGetCoordinatesFromFile(file, sample))
}


#//////////////////////////////////////////////////////////////
#' uuid.
#' 
#' Generates a timecode-based 16-bytes UUID.
#' The generated bytes are generated using the following pattern:
#' bytes 0 and 1: Year
#' byte 2: Month
#' byte 3: Day
#' byte 4: Hour
#' byte 5: Minute
#' byte 6: Second
#' bytes 7 to 15: random.
#'
#' @return the generated UUID encoded in a text string.
#'
uuid_timebased <- function()
{
  currentTime <- Sys.time()
  sUUID <- sprintf( "%04X", as.integer( format(currentTime, "%Y")))
  sUUID <- paste0( sUUID, sprintf( "%02X", as.integer( format(currentTime, "%m"))))
  sUUID <- paste0( sUUID, sprintf( "%02X", as.integer( format(currentTime, "%d"))))
  sUUID <- paste0( sUUID, sprintf( "%02X", as.integer( format(currentTime, "%H"))))
  sUUID <- paste0( sUUID, sprintf( "%02X", as.integer( format(currentTime, "%M"))))
  sUUID <- paste0( sUUID, sprintf( "%02X", as.integer( format(currentTime, "%S"))))
  sUUID <-paste0(sUUID, paste0(sprintf("%02X",sample(0:255, 9)), collapse = ""))
  
  Sys.sleep(1) #force to sleep one seconf to ensure each execition provides a unique uuid
  
  return(sUUID)
}

#Function to display processing time properly
display_processing_time <- function(time_elapsed, message)
{
  if(class(time_elapsed) != "proc_time")
  {
    stop("time_elapsed invalid class time")  
  }
  
  dsec <- time_elapsed["elapsed"]
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- round(dsec - 3600*hours - 60*minutes, digits = 3)
  cat(paste0(message, ": "))
  if(hours > 0) cat(paste0(hours, " hours,\t"))
  if(minutes > 0 || hours > 0) cat(paste0(minutes, " minutes,\t"))
  cat(paste0(seconds, " seconds\n"))
}


#//////////////////////////////////////////////////////////////
#retorna el índice más próximo en massVector a la masa dada
getIndexFromMass <- function(massVector, mass)
{
  highLimit=length(massVector);
  lowLimit=0;
  index=trunc(highLimit/2);
  while(1)
  {
    mz=massVector[index];
    if(mass>mz)
      {lowLimit=index;}
    else
      {highLimit=index;}
    index=trunc((highLimit+lowLimit)/2);
    if(trunc(highLimit-lowLimit)<=1) break;
  }
  return(index+1); #+1 por R
}

#//////////////////////////////////////////////////////////////
#' Create an empty rMSI object with defined mass axis and size.
#'
#' @param x_size the number of pixel in X direction.
#' @param y_size the number of pixel in Y direction.
#' @param mass_axis the mass axis.
#' @param pixel_resolution defined pixel size in um.
#' @param img_name the name for the image.
#' @param rMSIXBin_path where the rMSI files will be stored.
#' @param uuid a string containing an universal unique identifier for the image. If it is not provided it will be created using a time code.
#'
#' Creates an empty rMSI object with the provided parameters. This method is usefull to implement importation of new data formats
#' and synthetic datasets to test and develop processing methods and tools.
#'
#' @return the created rMSI object
#' export
#'
CreateEmptyImage<-function(x_size,
                           y_size,
                           mass_axis,
                           pixel_resolution,
                           img_name = "New empty image",
                           uuid = NULL)
{
  
  #TODO: create empty image with given X Y sizes and mas axis, so this function should actually create the rMSIXBin and the imzML files!
  #Document this and thing about the other CreateEmptyImage() image function
  
  
  img <- CreateEmptyImage( num_of_pixels = x_size*y_size,
                           mass_axis = mass_axis, 
                           pixel_resolution = pixel_resolution, 
                           img_name = img_name, 
                           uuid = uuid)
  
  img$size <- c( x_size, y_size )
  i <- 1
  for( xi in 1:x_size)
  {
    for( yi in 1:y_size)
    {
      img$pos[i,]<- c(xi, yi)
      i<-i+1
    }
  }
  
  return(img)
}

#//////////////////////////////////////////////////////////////
#' Create an empty rMSI object with defined mass axis and total number of pixels.
#'
#' @param num_of_pixels Total number of spectrums/pixels.
#' @param mass_axis the mass axis.
#' @param pixel_resolution defined pixel size in um.
#' @param img_name the name for the image.
#' @param rMSIXBin_path where the rMSI files will be stored.
#' @param uuid a string containing an universal unique identifier for the image. If it is not provided it will be created using a time code.
#'
#' Creates an empty rMSI object with the provided parameters. This method is usefull to implement importation of new data formats
#' and synthetic datasets to test and develop processing methods and tools.
#' img$size is initialized with c(NA, NA) and the pos matrix with NA coords. Size and pos matrix must be filled by user.
#'
#' @return the created rMSI object
#' export
#'
CreateEmptyImage<-function(num_of_pixels,
                           mass_axis, 
                           pixel_resolution, 
                           img_name = "New empty image", 
                           img_path = getwd(), 
                           uuid = NULL)
{
  img<-list()
  class(img) <- "rMSIObj"
  img$rMSI_format_version <- 2
  img$name <- img_name
  img$mass <- mass_axis
  img$size <- c( NA, NA )
  names(img$size) <- c("x", "y")
  
  #Prepare the pos matrix
  img$pos <- matrix( NA, ncol = 2, nrow = num_of_pixels )
  img$posMotors <- matrix( NA, ncol = 2, nrow = num_of_pixels )
  colnames(img$pos)<- c("x", "y")
  colnames(img$posMotors)<- c("x", "y")
  
  img$pixel_size_um <-  pixel_resolution
  img$mean <- rep(0, length(mass_axis))
  img$base <- rep(0, length(mass_axis))
  
  img$data <- list()
  class(img$data) <- "rMSIData"
  img$data$path <- img_path
  
  img$data$imzML <- list()
  class(img$data$imzML) <- "imzMLData"
  img$data$imzML$file <- NULL
  if(is.null(uuid))
  {
    img$data$imzML$uuid <- uuid_timebased()
  }
  else
  {
    img$data$imzML$uuid <- uuid  
  }
  
  #Init not availble imZML data to NULL (will be set latter outside this function)
  img$data$imzML$SHA <- NULL
  img$data$imzML$MD5 <- NULL
  img$data$imzML$continuous_mode <- NULL
  img$data$imzML$mz_dataType <- NULL
  img$data$imzML$int_dataType <- NULL
  img$data$imzML$run <- NULL
  
  img$data$peaklist <- list()
  class(img$data$peaklist) <- "peakList"
  #The peaklist field are created outside this function.
  
  img$normalizations <- data.frame();
  
  return(img)
}

#//////////////////////////////////////////////////////////////
#' remap2ImageCoords.
#'
#' @param dataPos a pos matrix as it is in rMSI data object.
#'
#' This function should be only used to implement data importers from foreign formats.
#' This functions maps a MALDI motors coors space to a image coord space.
#' dataPos matrix is a two columns matrix where first column stores x positions and second y pixel positions.
#' a remapped dataPos matrix do not contain empty raster positions neighter offsets.
#'
#' @return the dataPos matrix remapped.
#' export
#'
remap2ImageCoords <- function(dataPos)
{
  #1- Calc offsets and subtract it
  x_offset<-min(dataPos[,"x"])
  y_offset<-min(dataPos[,"y"])
  for(i in 1:nrow(dataPos))
  {
    dataPos[i, "x"] <- dataPos[i, "x"] - x_offset + 1
    dataPos[i, "y"] <- dataPos[i, "y"] - y_offset + 1
  }
  
  #2- Compute Motor coords range
  x_size<-max(dataPos[,"x"])
  y_size<-max(dataPos[,"y"])
  
  #3- Map MALDI motor coords to image cords (1-pixels steps)
  #It is important to map MALDI motors coords to image coords.
  #Otherwise, null extra pixels may be added leading to bad reconstruction
  px_map <- matrix( 0, nrow = x_size, ncol = y_size)
  for(i in 1:nrow(dataPos))
  {
    xi <- dataPos[i, "x"]
    yi <- dataPos[i, "y"]
    px_map[xi, yi]<- i
  }
  
  colNull <- which( base::colSums(px_map) == 0)
  rowNull <- which( base::rowSums(px_map) == 0)
  remap<-FALSE
  if( length(colNull) > 0 && length(rowNull) > 0 )
  {
    px_map_ <- px_map[ -rowNull , -colNull ]
    remap<-TRUE
  }
  if( length(colNull) > 0 && length(rowNull) == 0 )
  {
    px_map_ <- px_map[ , -colNull ]
    remap<-TRUE
  }
  if( length(colNull) == 0 && length(rowNull) > 0 )
  {
    px_map_ <- px_map[ -rowNull , ]
    remap<-TRUE
  }
  
  if(remap)
  {
    for(ix in 1:nrow(px_map_))
    {
      for(iy in 1:ncol(px_map_))
      {
        if(px_map_[ix, iy] > 0)
        {
          dataPos[px_map_[ix, iy], "x"] <- ix
          dataPos[px_map_[ix, iy], "y"] <- iy
        }
      }
    }
  }
  
  return(dataPos)
}

massX<-function(mzVect, initMass, finalMass, sp)
{
  if(exists("sp")) rPlotSpectrum(sp, initMass, finalMass)  
  mzVect[mzVect>=initMass & mzVect<=finalMass]  
}
