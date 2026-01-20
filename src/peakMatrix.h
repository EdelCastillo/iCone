/*********************************************************************************
 *     peakMatrix
 *     iCone - R package for MSI data processing
 *     Copyright (C) 2025 Esteban del Castillo Pérez (esteban.delcastillo@urv.cat)
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *********************************************************************************/
#ifndef R_PEAK_MATRIX
#define R_PEAK_MATRIX

#include <Rcpp.h>

#include "peakInfo.h"
#include "gmmPeak.h"
#include <time.h>
#include <sys/time.h>    
#include "rGetImzMLData.h"
#include "kmeansR.h"
#include <thread>
#include <mutex>
#include <chrono>
#include "common_methods.h"
#include <stdlib.h>
#include "noiseestimation.h"

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace Rcpp;
using namespace std; 


//class to obtain a peak matrix from the imzML file information.
class PeakMatrix
{
public:
  
  //Constructor
  //captures input information, allocates memory and initializes.
 PeakMatrix(const char* ibdFname, Rcpp::List imzML, Rcpp::List params,  Rcpp::NumericVector pxList, float mzLow=0, float mzHigh=0, int nThreads=0);
  
  //destructor
  //free reserved memory
  ~PeakMatrix();
  
  //freeing buffer.
  void freeMemoryPeak();
  
  //mtGetGaussians()
  //Parallel processing.
  //Peak are delimited and their Gaussians are formed.
  //This thread remains active, processing spectra until none remain.
  //Each spectrum is converted into Gaussians that can overlap (join).
  //spIndex: thread
  //Returns -1 on failure, 0 = OK.
  int  mtGetGaussians(int spIndex);
  
  //mtSegmentation()
  //Kmeans segmentation in parallel processing.
  //For each mass segment, centroids are generated.
  //The information is stored in the array pointed to by m_ionEntry_p.
  void mtSegmentation(int thrIndex);
  
  //rawToGaussians
  //gets the intensity peak and converts them into Gaussians.
  //The intensity and mass data adjusted to the range of interest are loaded from the imzML file.
  //the SNR info is established for each point of the spectra.
  //The generated information is stored in the m_spectro_p structure.
  int rawToGaussians();
  
  //getMassRanges()
  //establishes the mass segments where overlapping Gaussians exist.
  //the information is stored in the array pointed to by m_massRange_p.
  //returns the number of segments.
  //procedure:
  //The process is iterated up to 10 times, halving the m_linkedPeaks parameter until mass segments smaller than 3 Da 
  //are obtained. It ends when linkedPeaks < 0.001.
  //In each iteration, a vector is established with information on the number of Gaussian waves supporting a given bin, 
  //considering its sigma weighted by linkedPeaks. Very wide Gaussians (sigma > 5*deltaMass) are discarded.
  //Then, its baseband is estimated and considered as noise. Lower intensities are canceled, and the composite peaks are determined. 
  //These peaks determine the width of the mass segments.
  int getMassRanges(float *linkedPeaks);
  
  //getMassRanges2()
  //establishes the mass segments where overlapping Gaussians exist.
  //generates a Boolean mass axis with a resolution of 1/4 of the spectrometer's mass resolution.
  //if any Gaussian invades its space, sets the corresponding +/-3 m_linkedPeaks*sigma checkboxes to true.
  //Very wide Gaussians (sigma > 5*deltaMass) are discarded.
  //the information is stored in the array pointed to by m_massRange_p.
  //returns the number of segments.
  int getMassRanges2();
  
  //massRangeToCentroids()
  //Parallel processing.
  //1) Sets the number of independent mass ranges (joined peak).
  //2) Establishes clusters within those ranges (kmeans segmentation).
  //Requires preprocessing by getGaussians()
  //Receives the total mass range to consider.
  //Returns a list: peakMatrix, massVector, pixelsSupport.
  //peakMatrix: Matrix of centroids and the intensity associated with each pixel.
  //massVector: The mz associated with each column of the peakMatrix.
  //pixelsSupport: Number of pixels with intensity > 0.
  List massRangeToCentroids(MASS_RANGE massRange);
  
  //Obtains the average values of the Gaussian masses on an artificial mass axis.
  //The mass axis is formed from the extreme masses to be considered and the desired resolution.
  //Each mass resolution segment is divided into four parts.
  //pxLow and pxHigh delimit the spectra to be considered.
  //Returns a list with two arrays: averageMz and averageIntensity.
  List getMeanGaussianSpectrum(float resolution, int overSampling);
  
  //Obtains the average values of the intensities on an artificial mass axis.
  //The mass axis is formed from the extreme masses to be considered and the desired resolution.
  //Each mass resolution segment is divided into two parts.
  //pxLow and pxHigh delimit the spectra to be considered.
  //Returns a list with two arrays: averageMz and averageIntensity.
  List getMeanSpectrum(float resolution, int overSampling);
  
  GAUSS_SP *getGaussiansPointer();
  int getPixelsNumber();
  NumericMatrix getPixelGaussians(int px, float mzLow, float mzHigh);
    
private:  
  //getGaussians()
  //Called from a thread.
  //Sets the Gaussians on the peak.
  //Uses the peak separation information (m_peakFG_p).
  //px: pixel.
  //spectro: pointer to the spectrum.
  //gaussians_p: pointer to the structure containing the Gaussians' parameters.
  //Returns the number of Gaussians or a value < 0 on failure.
  int getGaussians(int px, SPECTRO *spectro_p, GAUSS_PARAMS *gaussians_p);
  
  //getCentroidsIntoRange()
  //Extracts the existing Gaussians within a mass range from the information in m_gaussians_p.
  //massRange: Mass range from which to extract the Gaussians.
  //gaussians_p: Requested Gaussians.
  //Returns the number of Gaussians.
  int getCentroidsIntoRange(MASS_RANGE massRange, float **gaussians_p, int size);
  
  //getRawInfo()
  //Loads the full spectrum information associated with a pixel from an imzML file.
  //The information is stored in the m_spectro structure, set to the range [m_mzLow, m_mzHigh].
  //px: Pixel whose spectrum should be loaded.
  //spIndex: Threads that manage it
  //Returns the size of the spectrum.
  int getRawInfo(int px, int spIndex);
  
  //setGaussiansNumberIntoSegments()
  //Determines the maximum number of Gaussians over the given mass intervals and all pixels.
  //aproximaciones sucesivas
  //massRange: Mass range to consider.
  //Returns the maximum value.
  int setGaussiansNumberIntoSegments (MASS_RANGE massRange); 
  
  //setGaussiansNumberIntoSegments()
  //Determines the maximum number of Gaussians over the given mass intervals and all pixels.
  //brute force.
  //massRange: Mass range to consider.
  //Returns the maximum value.
  int setGaussiansNumberIntoSegments2(MASS_RANGE massRange); 
  
  //getCentroidsNumberIntoRange()
  //Returns the number of Gaussians in a mass range from the information in m_gaussians_p.
  //massRange: Mass range from which to extract Gaussians.
  //Returns the number of Gaussians.
  int getCentroidsNumberIntoRange(MASS_RANGE massRange);
  
  int centers(float *mass_p, int size, float segmentSize, float* centers_p, int *centerSize_p);

public:   
  int     
          m_maxPxGaussians,
          m_massRangeSize;
  float   m_linkedPeaks; //Two peaks are considered linked if they are closer than the given standard deviation.
  double     
          m_massResolution;
 bool     m_hit;
    
private:  
  //input info to the constructor.
  bool m_continuous;
  int 
    m_nThreads,
    m_pxSupport,
    m_NPixels,
    m_maxMzLength,
    *m_pxList,
    m_pxMax,
    m_pxMin;
  float     
    m_mzLow,
    m_mzHigh,
    m_SNR;
  double     
    m_maxMassResolution;  
  
  //info generated in the class.
  GetImzMLData  *m_getImzMLData_p;
  PEAK_F_GROUP  *m_peakFG_p=0;
  GAUSS_SP      *m_gaussians_p;

  MASS_SEGMENT  m_massSegment;
  ION_ENTRY     **m_ionEntry_p;
  SPECTRO  m_spectro[MAX_THREADS];
  MASS_RANGE    *m_massRange_p;
  bool          m_enable;
  int     
                m_totalIons[MAX_THREADS],
                m_SNRmethod;
  NoiseEstimation *m_noiseEst_p; 
  int           *m_x, 
                *m_y,
                m_coordSize;
}; 

#endif
