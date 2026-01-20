---
title: "iCone: basic information"
author: "E del Castillo"
date: 'January 2026'
output: html_document
---

>iCone is an R package oriented to the field of mass spectrometry and focused on image processing (MSI).
Its purpose is to generate a centroid matrix from the information provided by a spectrometer in imzML format.
This is a new peak selection algorithm with a dual purpose: to obtain centroids with high accuracy and to deconvolve overlapping ions.

![](./centroids.png)

## Installation

> devtools::install_github("EdelCastillo/iCone")

> Please note that is a development version and no release has been made yet. So, keep looking at this page for future updates and releases.

## Use guide


### **Step 1**.- A list with the necessary parameters is generated.

> For **example**:  
>**params**\ <-list("massResolution"=30000, "minPixelsSupport"=5, "SNR"=3,  "SNRmethod"="estnoise_mad", "linkedPeaks"=3);

**Description of the parameters:**

```         
massResolution:     Desired mass resolution for the peak matrix (mz/delta_mz).
minPixelSupport:    Minimum percentage of pixels that must provide intensity information for each mass (mz).
                    That is, in the peak matrix, the relative number of pixels with non-zero intensity must be equal to or greater than this.
SNR:                Signal-to-noise ratio. Discards information based on its proximity to the estimated noise level.
noiseMethod:        Procedure used for noise estimation:
                      estnoise_diff:  Average of the absolute differences with respect to the mean value.
                      estnoise_sd:    Standard deviation of the absolute differences (Gaussian filter).
                      estnoise_mad:   Median of the absolute deviations (Gaussian filter).
linkedPeaks:        Two peaks are considered linked if they are closer than the given standard deviation.                         
```

### **Step 2**.- The peak selection algorithm is applied and its peak matrix is obtained:

> lst <-**getPeakMatrix**(imzML_file_path, params, initMass, finalMass, pxList, nThreads)

**Description of the parameters:**   
```
imzML_file_path: Absolute path to the filename with the imzML extension.
                 The attached binary file, with the ibd extension, must be in the same directory.
         params: List of parameters given in Step 1.
       initMass: Initial mass to consider.
      finalMass: Final   mass to consider.
         pxList: list of pixels. First pixel=1. By default everyone.
       nThreads: Number of threads for parallel processing (if zero, nThreads=maxCores-1)

  return a list: 
    peakMatrix: Matrix of peak (centroids) rows = pixels, columns = intensity of each pixel.
          mass: Vector with the masses associated with each column of peakMatrix.
massResolution: Vector with the mass resolution achieved in the peak matrix (mz/delta_mz).
 pixelsSupport: Vector with the number of pixels in each column with non-zero intensity.
```
### **complementary functions.**

> lst <-**getPixelGaussians**(imzML_file_path, params, initMass, finalMass, pixel)

Reports information about a single spectrum.

**Description of the parameters:**   
```
imzML_file_path: Absolute path to the filename with the imzML extension.
                 The attached binary file, with the ibd extension, must be in the same directory.
         params: List of parameters given in Step 1.
       initMass: Initial mass to consider.
      finalMass: Final   mass to consider.
          pixel: Reference to the spectrum (the first one is 1)
          
  return a list: 
      gaussians: Matrix with the Gaussian spectra.
           mass: Vector with raw masses.
      intensity: Vector with raw intensities.
            SNR: Vector with signal-to-noise ratio for each mass point.
          noise: Noise estimation.
```

> lst <-**getAverageGaussianSpectrum**(imzML_file_path, params, initMass, finalMass, pxList, overSampling, nThreads)

Report the average of all conformal Gaussians in the indicated spectra.

**Description of the parameters:**   
```
imzML_file_path: Absolute path to the filename with the imzML extension.
                 The attached binary file, with the ibd extension, must be in the same directory.
         params: List of parameters.
            "massResolution": mass resolution (mz/deltaMz).
                       "SNR": signal-to-noise ratio
               "noiseMethod": method for estimating noise.
       initMass: Initial mass to consider.
      finalMass: Final   mass to consider.
         pxList: list of pixels. First pixel=1. By default everyone.
   overSampling: interval between points on the mass axis = massResolution/overSampling.
          
  return a list: 
       averageMz: Vector with mass axis. A mass bin is 1/4 of the massResolution parameter
averageIntensity: Vector with intensity axis.
```

> lst <-**getAverageSpectrum**(imzML_file_path, params, initMass, finalMass, pxList, overSampling)

Report the average of all intensities in the indicated spectra.

**Description of the parameters:**   
```
imzML_file_path: Absolute path to the filename with the imzML extension.
                 The attached binary file, with the ibd extension, must be in the same directory.
         params: 
          "massResolution": mass resolution (mz/deltaMz).
       initMass: Initial mass to consider.
      finalMass: Final   mass to consider.
         pxList: list of pixels. First pixel=1. By default everyone.
   overSampling: interval between points on the mass axis = massResolution/overSampling.
          
  return a list: 
       averageMz: Vector with mass axis. A mass bin is 1/2 of the massResolution parameter
averageIntensity: Vector with intensity axis.
```
### **functions for data visualization.**

> **rPlotSpectrum**(values, initMass, finalMass)

View the information returned by the getAverageGaussianSpectrum() and getAverageSpectrum() functions.

**Description of the parameters:**   
```
         values: list returned by the getAverageGaussianSpectrum() and getAverageSpectrum() functions.
       initMass: Initial mass to consider.
      finalMass: Final   mass to consider.
```

> **rPlotGaussianSpectrum**(drawInfo, initMass, finalMass)

View the information returned by the getPixelGaussians() function.

**Description of the parameters:**   
```
       drawInfo: information returned by the getPixelGaussians() function.
       initMass: Initial mass to consider.
      finalMass: Final   mass to consider.
```
> **rPlotGaussianSpectrum2**(drawInfo1, drawInfo2, minMass, maxMass)

Display information from two spectra on the same plot.

**Description of the parameters:**   
```
      drawInfo1: information returned by the getPixelGaussians() function.
      drawInfo2: information returned by the getPixelGaussians() function.
       initMass: Initial mass to consider.
      finalMass: Final   mass to consider.
```
> **rPlotGaussianSpectrum3**(drawInfo1, drawInfo2, drawInfo3, minMass, maxMass)

Display information from tree spectra on the same plot.

**Description of the parameters:**   
```
      drawInfo1: information returned by the getPixelGaussians() function.
      drawInfo2: information returned by the getPixelGaussians() function.
      drawInfo3: information returned by the getPixelGaussians() function.
       initMass: Initial mass to consider.
      finalMass: Final   mass to consider.
```

