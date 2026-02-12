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
#' @param data_file: absolute reference to the file with the imzML extension.
#' @param params  
#'      "massResolution": mass resolution with which the spectra were acquired.
#'                 "SNR": signal-to-noise ratio (by defect=1)
#'         "noiseMethod": method for estimating noise (by defect="estnoise_mad").
#'    "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered (by defect=1).
#'         "linkedPeaks": two peaks are considered linked if they are closer than the given standard deviation (by defect=3).   
#' @param lowMass:  lower mass to consider
#' @param highMass: higher mass to consider
#' @param pxList:   list of pixels. First pixel=1. By default everyone.
#' @param nThreads: number of threads for parallel processing (if zero, nThreads=maxCores-1)
#' @param imzMLChecksum: if the binary file checksum must be verified, it can be disabled for convenice with really big files.
#' @param fixBrokenUUID: set to FALSE by default to automatically fix an uuid mismatch between the ibd and the imzML files (a warning message will be raised).
#'
#' @return a list: 
#'   peakMatrix: Matrix of peak (centroids) rows = pixels, columns = intensity of each pixel.
#'         mass: Vector with the masses associated with each column of peakMatrix.
#'pixelsSupport: Vector with the number of pixels in each column with non-zero intensity. 
#'
#' @export
getPeakMatrix<-function(data_file,
                   params,
                   lowMass=0,
                   highMass=0,
                   pxList=c(0),
                   nThreads=0,
                   imzMLChecksum = F, 
                   fixBrokenUUID = F)
{
  if(!file.exists(data_file))
  {
    stop("File not found\n")
  }
  
  #parameters control 
  if(lowMass==0 && highMass==0)
    {print("ERROR: lowMass and highMass must be given."); return (0);}
  if(!(exists("SNR", where=params)))
  {
    print("warning: by default, SNR parameter will be 1.")
    params=c(params, "SNR"=1)
  }
  if(!(exists("minPixelsSupport", where=params)))
  {
    print("warning: by default, minPixelsSupport parameter will be 1%.")
    params=c(params, "minPixelsSupport"=1)
  }
  if(!(exists("noiseMethod", where=params)))
  {
    print("warning: by default, noiseMethod parameter will be MAD type.")
    params=c(params, "noiseMethod"="estnoise_mad")
  }
  if(!(exists("linkedPeaks", where=params)))
  {
    print("warning: by default, linkedPeaks parameter will be 3 sd.") 
    params=c(params, "linkedPeaks"=3)
  }
  if(!(exists("massResolution", where=params)))
  {
    print("ERROR: massResolution parameter is required.") 
    return (0)
  }
  if(params$massResolution<=0)
  {
    print("ERROR: massResolution must be greater than zero.") 
    return (0)
  }

  imgData <- NULL
  
  #Default  fun_label()
  #fun_label This is a callback function to update the progress bar dialog text.
  fun_label <- function(text)
    {
      cat(text)
      cat("\n")
    }
  #fun_progress This is a callback function to update the progress of loading data. See details for more information.
  fun_progress = NULL; 
  #close_signal function to be called if the loading process is aborted.
  close_signal = NULL;
  #imzMLRename the image name, if NULL a default name based on the file name will be used.
  imzMLRename  = NULL;
  pt<-proc.time()

  fileExtension <- unlist(strsplit(basename(data_file), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension == "imzML")
  {
    #capturing information from an imzML file.
    in_img <- import_imzML(path.expand(data_file),  fun_progress = fun_progress, fun_text = fun_label, 
                           close_signal = close_signal, verifyChecksum = imzMLChecksum, subImg_rename = imzMLRename, 
                           subImg_Coords = NULL, fixBrokenUUID = fixBrokenUUID)
    
    #binary data file name.
    file=path.expand(file.path(in_img$data$path, paste0(in_img$data$imzML$file, ".ibd")));
    if(lowMass>highMass)
    {
      print("warning: the lower mass exceeds the upper mass")
      lowMass=highMass
    }
    #peak matrix generation.
    peakMatrix<-peakMatrix(file, in_img$data$imzML, params, lowMass, highMass, pxList-1, nThreads);
  
  gc()  
  pt<-proc.time() - pt
  display_processing_time(pt, "Data processing time")
  return(peakMatrix)
  }
  print("warning: imzML file type is required." )
}


#' @name getPxGaussians
#' @title It retrieves information from a single pixel.
#' @param data_file: absolute reference to the file with the imzML extension.
#' @param params  
#'      "massResolution": mass resolution with which the spectra were acquired (mz/deltaMz).
#'                 "SNR": signal-to-noise ratio
#'         "noiseMethod": method for estimating noise.
#' @param lowMass:  lower mass to consider
#' @param highMass: higher mass to consider
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
                        params,
                        lowMass=0,
                        highMass=0,
                        pixel=1,
                        imzMLChecksum = F, 
                        fixBrokenUUID = F)
{
  if(!file.exists(data_file))
  {
    stop("File not found\n")
  }
  
  #control de parámetros
  if(!(exists("SNR", where=params)))
  {
    print("warning: by default, SNR parameter will be 1.")
    params=c(params, "SNR"=1)
  }
  if(!(exists("minPixelsSupport", where=params)))
  {
    print("warning: by default, minPixelsSupport parameter will be 1%.")
    params=c(params, "minPixelsSupport"=1)
  }
  if(!(exists("noiseMethod", where=params)))
  {
    print("warning: by default, noiseMethod will be MAD type.")
    params=c(params, "noiseMethod"="estnoise_mad")
  }
  if(!(exists("massResolution", where=params)))
  {
    print("ERROR: massResolution parameter is required.") 
    return (0)
  }
  if(params$massResolution<=0)
  {
    print("ERROR: massResolution must be greater than zero.") 
    return (0)
  }
  
  imgData <- NULL
  
  #Default  fun_label()
  #fun_label This is a callback function to update the progress bar dialog text.
  fun_label <- function(text)
  {
    cat(text)
    cat("\n")
  }
  #fun_progress This is a callback function to update the progress of loading data. See details for more information.
  fun_progress = NULL; 
  #close_signal function to be called if the loading process is aborted.
  close_signal = NULL;
  #imzMLRename the image name, if NULL a default name based on the file name will be used.
  imzMLRename  = NULL;
  
  pt<-proc.time()
  
  fileExtension <- unlist(strsplit(basename(data_file), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension == "imzML")
  {
    #capturing information from an imzML file.
    in_img <- import_imzML(path.expand(data_file),  fun_progress = fun_progress, fun_text = fun_label, 
                           close_signal = close_signal, verifyChecksum = imzMLChecksum, subImg_rename = imzMLRename, 
                           subImg_Coords = NULL, fixBrokenUUID = fixBrokenUUID)
    
    #binary data file name.
    file=path.expand(file.path(in_img$data$path, paste0(in_img$data$imzML$file, ".ibd")));
    if(lowMass>highMass)
    {
      print("warning: the lower mass exceeds the upper mass")
      lowMass=highMass
    }
    if(pixel<1) 
    {
      print("warning: pixel must be greater than zero. Thus, pixels will be 1.")
      pixel=1;
    }
    #info about pixel (the first pixel is zero)
    pxGauss<-rGetPixelGaussians(file, in_img$data$imzML, params, lowMass, highMass, pixel-1);

    pt<-proc.time() - pt
    display_processing_time(pt, "Data processing time")
    return(pxGauss)
  }
  print("warning: imzML file type is required." )
}

#' @name getAverageGaussianSpectrum
#' @title It obtains the average value of the Gaussian from all data into .imzML file.
#' @param data_file: absolute reference to the file with the imzML extension.
#' @param params  
#'      "massResolution": mass resolution with which the spectra were acquired (mz/deltaMz).
#'                 "SNR": signal-to-noise ratio
#'         "noiseMethod": method for estimating noise.
#' @param lowMass:   lower mass to consider
#' @param highMass:  higher mass to consider
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
                         params,
                         lowMass=0,
                         highMass=0,
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
  
  #control de parámetros
  if(!(exists("SNR", where=params)))
  {
    print("warning: by default, SNR parameter will be 1.")
    params=c(params, "SNR"=1)
  }
  if(!(exists("minPixelsSupport", where=params)))
  {
    print("warning: by default, minPixelsSupport parameter will be 1%.")
    params=c(params, "minPixelsSupport"=1)
  }
  if(!(exists("noiseMethod", where=params)))
  {
    print("warning: by default, noiseMethod will be MAD type.")
    params=c(params, "noiseMethod"="estnoise_mad")
  }
  if(!(exists("linkedPeaks", where=params)))
  {
    #print("warning: by default, linkedPeaks parameter will be 3 sd.") 
    params=c(params, "linkedPeaks"=3)
  }
  if(!(exists("massResolution", where=params)))
  {
    print("ERROR: massResolution parameter is required.") 
    return (0)
  }
  if(params$massResolution<=0)
  {
    print("ERROR: massResolution must be greater than zero.") 
    return (0)
  }

  imgData <- NULL
  
  #Default  fun_label()
  #fun_label This is a callback function to update the progress bar dialog text.
  fun_label <- function(text)
  {
    cat(text)
    cat("\n")
  }
  #fun_progress This is a callback function to update the progress of loading data. See details for more information.
  fun_progress = NULL; 
  #close_signal function to be called if the loading process is aborted.
  close_signal = NULL;
  #imzMLRename the image name, if NULL a default name based on the file name will be used.
  imzMLRename  = NULL;
  
  pt<-proc.time()
  
  fileExtension <- unlist(strsplit(basename(data_file), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension == "imzML")
  {
    #capturing information from an imzML file.
    in_img <- import_imzML(path.expand(data_file),  fun_progress = fun_progress, fun_text = fun_label, 
                           close_signal = close_signal, verifyChecksum = imzMLChecksum, subImg_rename = imzMLRename, 
                           subImg_Coords = NULL, fixBrokenUUID = fixBrokenUUID)
    
    #binary data file name.
    file=path.expand(file.path(in_img$data$path, paste0(in_img$data$imzML$file, ".ibd")));
    if(lowMass>highMass)
    {
      print("warning: the lower mass exceeds the upper mass")
      lowMass=highMass
    }

    avSp=rGetAverageGaussianSpectrum(file, in_img$data$imzML, params, lowMass, highMass, pxList-1, overSampling, nThreads);

    pt<-proc.time() - pt
    display_processing_time(pt, "Data processing time")
    return(avSp)
  }
  print("warning: imzML file type is required." )
}

#' @name getAverageSpectrum
#' @title It obtains the average value of the intensities from all data into .imzML file.
#'        Noise is not taken into account.
#' @param data_file: absolute reference to the file with the imzML extension.
#' @param lowMass:   lower mass to consider
#' @param highMass:  higher mass to consider
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
                             lowMass=0,
                             highMass=0,
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
  
  #Default  fun_label()
  #fun_label This is a callback function to update the progress bar dialog text.
  fun_label <- function(text)
  {
    cat(text)
    cat("\n")
  }
  #fun_progress This is a callback function to update the progress of loading data. See details for more information.
  fun_progress = NULL; 
  #close_signal function to be called if the loading process is aborted.
  close_signal = NULL;
  #imzMLRename the image name, if NULL a default name based on the file name will be used.
  imzMLRename  = NULL;
  
  pt<-proc.time()
  
  fileExtension <- unlist(strsplit(basename(data_file), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension == "imzML")
  {
    #capturing information from an imzML file.
    in_img <- import_imzML(path.expand(data_file),  fun_progress = fun_progress, fun_text = fun_label, 
                           close_signal = close_signal, verifyChecksum = imzMLChecksum, subImg_rename = imzMLRename, 
                           subImg_Coords = NULL, fixBrokenUUID = fixBrokenUUID)
    
    #binary data file name.
    file=path.expand(file.path(in_img$data$path, paste0(in_img$data$imzML$file, ".ibd")));
    if(lowMass>highMass)
    {
      print("warning: the lower mass exceeds the upper mass")
      lowMass=highMass
    }
    
    avSp=rGetAverageSpectrum(file, in_img$data$imzML, params, lowMass, highMass, pxList-1, overSampling);
    
    pt<-proc.time() - pt
    display_processing_time(pt, "Data processing time")
    return(avSp)
  }
  print("warning: imzML file type is required." )
}

#' @name getGaussiansFromSpectrum
#' @title Form Gaussians over each peak of the given spectrum.
#' @param Intensity:   vector with value associated with each mass point.
#' @param mz:          vector with mass/charge information
#' @param lowMass:     lower mass to consider
#' @param highMass:    higher mass to consider
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
                              lowMass=0,
                              highMass=0,
                              SNR=3,
                              noiseMethod="estnoise_mad")
{

  #control de parámetros
    params=list("SNR"=3, "noiseMethod"="estnoise_mad", "massResolution"=30000, "minPixelsSupport"=1, "linkedPeaks"=3)

    params$SNR=SNR
    params$noiseMethod=noiseMethod

  pt<-proc.time()
  
    if(lowMass>highMass)
    {
      print("warning: the lower mass exceeds the upper mass")
      lowMass=highMass
    }
    
    Gauss=rGetGaussiansFromSpectrum(intensity, mz, params, lowMass, highMass);
    
    pt<-proc.time() - pt
    display_processing_time(pt, "Data processing time")
    return(Gauss)
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

