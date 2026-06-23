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
#ifndef SEGMENTS_H
#define SEGMENTS_H

#include <Rcpp.h>

#include "peakInfo.h"
#include "common_methods.h"
#include "noiseestimation.h"
#include <stdlib.h>
#include <fstream>

using namespace Rcpp;
using namespace std; 


//class: extracts mass segments where overlapping Gaussians exist. Over all pixels.
//The source information is located in a temporary file.
class Segments
{
public:
  Segments(float massResolution, float mzLow, float mzHigh, float linkedPeaks);
  ~Segments();
  
  //Load the file with information about the Gaussian curves associated with each pixel
  //If there is more than one sample, all its pixels are integrated
  //fileName: temporary file name
  //totalPixels: total pixels in the file
  //Returns false if the load failed
  int loadGaussians(char *fileName);
    
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
  

  int   m_totalPixels;
  float m_massResolution, 
        m_mzHigh, 
        m_mzLow, 
        m_linkedPeaks;
  GAUSS_SP      *m_gaussians_p;
  MASS_RANGE    *m_massRange_p;
};
#endif

  
