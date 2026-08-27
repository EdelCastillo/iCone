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
#ifndef PEAK_MATRIX_H
#define PEAK_MATRIX_H

#include <Rcpp.h>

#include "common_methods.h"
#include <stdlib.h>
#include <fstream>

using namespace Rcpp;
using namespace std; 


//class to obtain a peak matrix. It requires the prior contribution of classes RawToGaussians and Segments.
class PeakMatrix
{
public:
  //Constructor
  //Generates the peak matrix; final step. It requires the prior contribution of classes RawToGaussians and Segments.
  //   totalPixels: cumulative value of the pixels of each sample analyzed.
  //massResolution: desired mass resolution (mz/delta_mz).
  //         mzLow: lower  mass to consider
  //        mzHigh: higher mass to consider
  //     pxSupport: minimum percentage of pixels that must support an ion for it to be considered.
  PeakMatrix(int totalPixels, double massResolution, double mzLow, double mzHigh, double pxSupport, int nSamples, char *baseDir);
    
  //destructor
  ~PeakMatrix();
  
  //Loads the file with the coordinates of each pixel
  //If there is more than one sample, all its pixels are integrated
  //totalPixels: total pixels in the file
  //Returns false if the load failed
  int loadPixelsCoordinates(char *fileName, int totalPixels);
  
  //Load the file with information about the Gaussian curves associated with each pixel
  //If there is more than one sample, all its pixels are integrated
  //fileName: temporary file name (".../tmpGaussians.bin")
  //Returns false if the load failed
  int loadGaussians(char *fileName);
  
  //getCentroidsIntoRange()
  //Extracts the existing Gaussians within a mass range from the information in m_gaussians_p.
  //massRange: Mass range from which to extract the Gaussians.
  //mass_p: vector of internally generated masses.
  //px_p: vector of internally generated pixels.
  //iGauss_p: vector of internally generated gaussians index
  //massSize: size of reserved memory for each vector.
  //Returns the number of masses.
  int getCentroidsIntoRange(double *mass_p, int *px_p, int *iGauss_p, int massSize);
  
  //Determines the centers of groups of values whose distance does not exceed a given tolerance.
  //mass_p: vector masses.
  //px_p: vector of pixels.
  //iGauss_p: vector gaussians index
  //tolerance: maximum bin size.
  //Returns the number of centers detected (length of the centers_p array).
  int centers(double *mass_p, int *px_p, int *iGauss_p, int size, double tolerance);
  
  //generate de peak matrix with the centroids, their tolerance, and the number of support pixels.  
  int getCentroids();
  
  //Save the peak matrix to the file ".../tmpPeakMatrix.bin"
  int infoToFile(char *baseDir);
  
  int getSamplesPixelNumber();
    
  double *m_centers_p;
  int    *m_centersSize_p;
  int     m_nCentroids;
  
  int   m_totalPixels,
  m_nIons,
  m_nSamples;
  double m_massResolution, 
  m_mzHigh, 
  m_mzLow, 
  m_pxSupport,
  m_linkedPeaks;
  GAUSS_SP      *m_gaussians_p;
  MASS_RANGE    *m_massRange_p;
  ION_ENTRY     m_ionEntry;
  int           m_pixelsSample[MAX_SAMPLES];
  PIXEL_XY      *m_pixelsCoordinates_p;
  char *        m_baseDir;
};
#endif
