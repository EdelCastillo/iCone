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
#include <fstream>

using namespace Rcpp;
using namespace std; 


//class: converts the peaks of each spectrum in the imzML file into Gaussian waves.
class RawToGaussians
{
public:
  
  //Constructor
  //captures input information, allocates memory and initializes.
  RawToGaussians(char *baseDir, const char* ibdFname, Rcpp::List imzML, Rcpp::List params,  Rcpp::NumericVector pxList, double mzLow=0, double mzHigh=0, int nThreads=0);
  
  // Saves the Gaussian data to the given file
  // Adds it to any existing data
  // Returns false if failed 
  bool saveGaussians(char *fileName, GAUSS_SP *gauss_p);
  
  //Save the coordinate (XY) information of each pixel of m_pxList[] to the fileName file.
  //It is added to any existing data.
  //returns false if it failed.
  bool savePixelsCoordinates(char *fileName, NumericVector X, NumericVector Y);
    
    //destructor
  //free reserved memory
  ~RawToGaussians();
  
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
  
  //rawToGaussians
  //gets the intensity peak and converts them into Gaussians.
  //The intensity and mass data adjusted to the range of interest are loaded from the imzML file.
  //the SNR info is established for each point of the spectra.
  //The generated information is stored in the m_spectro_p structure.
  int rawToGaussians();
  
  //Obtains the average values of the Gaussian masses on an artificial mass axis.
  //The mass axis is formed from the extreme masses to be considered and the desired resolution.
  //Each mass resolution segment is divided into four parts.
  //pxLow and pxHigh delimit the spectra to be considered.
  //Returns a list with two arrays: averageMz and averageIntensity.
  List getMeanGaussianSpectrum(double resolution, int overSampling);
  
  //Obtains the average values of the intensities on an artificial mass axis.
  //The mass axis is formed from the extreme masses to be considered and the desired resolution.
  //Each mass resolution segment is divided into two parts.
  //pxLow and pxHigh delimit the spectra to be considered.
  //Returns a list with two arrays: averageMz and averageIntensity.
  List getMeanSpectrum(double resolution, int overSampling);
  
  GAUSS_SP *getGaussiansPointer();
  int getPixelsNumber();
  NumericMatrix getPixelGaussians(int px, double mzLow, double mzHigh);
  
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
  
  //getRawInfo()
  //Loads the full spectrum information associated with a pixel from an imzML file.
  //The information is stored in the m_spectro structure, set to the range [m_mzLow, m_mzHigh].
  //px: Pixel whose spectrum should be loaded.
  //spIndex: Threads that manage it
  //Returns the size of the spectrum.
  int getRawInfo(int px, int spIndex);
  
public:   
  int     
  m_maxPxGaussians,
  m_massRangeSize,
  m_NPixels;
  double   m_linkedPeaks; //Two peaks are considered linked if they are closer than the given standard deviation.
  double     
    m_massResolution;
  bool     m_hit;
  GAUSS_SP      *m_gaussians_p;
  
private:  
  //input info to the constructor.
  bool m_continuous;
  int 
    m_nThreads,
    m_pxSupport,
    m_maxMzLength,
    *m_pxList,
    m_pxMax,
    m_pxMin;
    double     
    m_mzLow,
    m_mzHigh,
    m_SNR;
  double     
    m_maxMassResolution;  
  
  //info generated in the class.
  GetImzMLData  *m_getImzMLData_p;
  PEAK_F_GROUP  *m_peakFG_p=0;
  
  MASS_SEGMENT  m_massSegment;
  SPECTRO  m_spectro[MAX_THREADS];
  bool          m_enable;
  int           m_SNRmethod;
  NoiseEstimation *m_noiseEst_p; 
  double *Z;
}; 

#endif
