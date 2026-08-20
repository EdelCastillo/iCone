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
#include "kmeansR.h"
#include "segments.h"
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
  //massResolution: desired mass resolution.
  //         mzLow: lower  mass to consider
  //        mzHigh: higher mass to consider
  //     pxSupport: minimum percentage of pixels that must support an ion for it to be considered.
  //   massRange_p: array of structures with information on isolated mass segments.
  //   gaussians_p: array of structures with information about Gaussians.
  PeakMatrix(int totalPixels, double massResolution, double mzHigh, double mzLow, double pxSupport, MASS_RANGE *massRange_p, GAUSS_SP *gaussians_p, int nThreads);
  
  //destructor
  ~PeakMatrix();
  
  //Loads the file with the coordinates of each pixel
  //If there is more than one sample, all its pixels are integrated
  //Returns false if the load failed
  int loadPixelsCoordinates(char *fileName, int totalPixels);
  
  //mtSegmentation()
  //Kmeans segmentation in parallel processing.
  //For each mass segment, centroids are generated.
  //The information is stored in the array pointed to by m_ionEntry_p.
  void mtSegmentation(int thrIndex);
  
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
  List massRangeToCentroids(char *baseDir, MASS_RANGE massRange);
  
  //getCentroidsIntoRange()
  //Extracts the existing Gaussians within a mass range from the information in m_gaussians_p.
  //massRange: Mass range from which to extract the Gaussians.
  //gaussians_p: Requested Gaussians.
  //Returns the number of Gaussians.
  int getCentroidsIntoRange(MASS_RANGE massRange, double **gaussians_p, int size);
  
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
  
  // Determines the centers of groups of values whose distance does not exceed segmentSize.
  // mass_p: array of data to consider.
  // size: size of the mass_p array.
  // Centers_p: array of final centers.
  // centerSize_p: number of elements that make up each center.
  // Returns the number of centers detected (length of the centers_p array).
  int centers(double *mass_p, int size, double segmentSize, double* centers_p, int *centerSize_p);
  
  ION_ENTRY     **m_ionEntry_p;
  MASS_RANGE    *m_massRange_p;
  MASS_SEGMENT  m_massSegment;
  GAUSS_SP      *m_gaussians_p;
  bool          m_enable;
  int           m_nThreads,
                m_massRangeSize,
                m_totalPixels,
                m_totalIons[MAX_THREADS],
                m_pxSupport;
  std::thread  *m_thread_p[MAX_THREADS];                //thread that processes the spectrum.

  double         m_massResolution, 
                m_mzHigh, 
                m_mzLow; 
  PIXEL_XY      *m_pixelsCoordinates_p;
  int           m_pixelsSample[MAX_SAMPLES],
                m_nSamples;
};
#endif
