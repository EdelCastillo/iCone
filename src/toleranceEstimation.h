
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
#ifndef TOLERANCE_EST
#define TOLERANCE_EST

#include <Rcpp.h>

#include "peakInfo.h"
#include "rGetImzMLData.h"
#include <thread>
#include <mutex>
#include "common_methods.h"
#include <stdlib.h>
#include "noiseestimation.h"

using namespace Rcpp;
using namespace std; 


//class: converts the peaks of each spectrum in the imzML file into Gaussian waves.
class ToleranceEstimation
{
public:
  typedef struct
  {
    double lowMz, maxMz, highMz, maxIntensity;
    int lowIdx, maxIdx, highIdx;
  }PEAKS_F;
  
  typedef struct
  {
    PEAKS_F *peaks_p;
    int size;
  }PEAKS_FG;
  
  //Constructor
  //captures input information, allocates memory and initializes.
  ToleranceEstimation(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, int nThreads=0);
  
  //destructor
  //free reserved memory
  ~ToleranceEstimation();
  
  //parallel processing
  //Estimates tolerance by counting the mass points that make up each peak
  //Analyzes the 1,000 spectra with the highest number of mass points
  //Returns a matrix containing the estimate for 10 consecutive mass ranges.
  NumericMatrix getTolerance();

private:  

  //getRawInfo()
  //Loads the full spectrum information associated with a pixel from an imzML file.
  //The information is stored in the m_spectro structure, set to the range [m_mzLow, m_mzHigh].
  //px: Pixel whose spectrum should be loaded.
  //spIndex: Threads that manage it
  //Returns the size of the spectrum.
  int getRawInfo(int px, int spIndex);
  
  //mtGetGaussians()
  //Parallel processing.
  //Peak are delimited and their Gaussians are formed.
  //This thread remains active, processing spectra until none remain.
  //Each spectrum is converted into Gaussians that can overlap (join).
  //spIndex: thread
  //Returns -1 on failure, 0 = OK.
  int  mtToleranceEstimation(int spIndex);
  
  //Returns the index of data.maxMz closest to value
  //If nearest bits 1:0 == 00, returns the nearest
  //If nearest bits 1:0 == 01, returns the nearest above
  //If nearest bits 1:0 == 10, returns the nearest below
  //If nearest bit    2 == 0 & value is out of range returns the nearest
  //If nearest bit    2 == 1 & value is out of range returns -1
  int nearestIndexPeak(float value, PEAKS_F *data, int size, int nearest);
  
public:   
  int     
  m_NPixels;
  bool  m_hit;

  
private:  
  //input info to the constructor.
  bool  m_continuous;
  int 
      m_nThreads,
      m_maxMzLength,
      m_pxMax,
      m_pxMin;
  double     
      m_SNR;
  
  //info generated in the class.
  GetImzMLData  *m_getImzMLData_p;
  SPECTRO       m_spectro[MAX_THREADS];
  bool          m_enable;
  int           m_SNRmethod;
  NoiseEstimation *m_noiseEst_p; 
  PeakMethod    m_peakMethod;
  PEAKS_FG      *m_peaksFG_p;
  int           *m_pxLen_p,
                *m_pxLenIdx_p;
}; 

#endif
