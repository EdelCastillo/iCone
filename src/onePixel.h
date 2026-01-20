/*********************************************************************************
 *     onePixel
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
#ifndef ONE_PIXEL
#define ONE_PIXEL

#include <Rcpp.h>
#include "peakInfo.h"
#include "gmmPeak.h"
#include "rGetImzMLData.h"
#include "common_methods.h"
#include <stdlib.h>
#include "noiseestimation.h"

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace Rcpp;
using namespace std; 

class OnePixel
{
public:
  //Constructor
  //captures input information, allocates memory and initializes.
  OnePixel(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, int pixel);
  
  //destructor
  //free reserved memory
  ~OnePixel();
  
  //Loads the raw data from a file, separates the peaks, and establishes the Gaussian.
  //results in internal structure. 
  //return the number of Gaussians or -1 is KO
  int rawToGaussians();
  
  //getGaussians()
  //Sets the Gaussians on the peak.
  //Uses the peak separation information (m_peakFG_p).
  //spectro: pointer to the spectrum.
  //gaussians_p: pointer to the structure containing the Gaussians' parameters.
  //Returns the number of Gaussians or a value < 0 on failure.
  int getGaussians(SPECTRO *spectro_p, GAUSS_PARAMS *gaussians_p);
  
  //Returns a list with information about a spectrum
  // gaussians: matrix with parameters for each Gaussian
  // mass: vector with the masses of the raw spectrum
  // intensity: intensity associated with each mass of the raw spectrum
  // SNR: signal-to-noise ratio associated with each mass of the raw spectrum.
  // require of rGetPixelGaussians()
  List getPixelGaussians(int px, float mzLow, float mzHigh);
  
  int       m_NPixels;
  bool      m_exit;
    
private:  
  //input info to the constructor.
  bool m_continuous;
  int 
      m_nThreads,
      m_mzLength;
  float     
      m_mzLow,
      m_mzHigh,
      m_SNR,
      m_noise;
  double     
      m_massResolution,
      m_maxMassResolution;  
  
  //info generated in the class.
  GetImzMLData  *m_getImzMLData_p;
  PEAK_F_GROUP  m_peakFG;
  GAUSS_SP      m_gaussians;
  
  MASS_SEGMENT  m_massSegment;
  ION_ENTRY     **m_ionEntry_p;
  SPECTRO       m_spectro;
  MASS_RANGE    *m_massRange_p;
  bool          m_enable;
  int     
              m_pixel,
             m_SNRmethod;
  NoiseEstimation *m_noiseEst_p; 
  
};

#endif
