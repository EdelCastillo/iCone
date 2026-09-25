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
#include "noiseestimation.h"
#include "peakInfo.h"
#include "gmmPeak.h"
#include "kmeansR.h"

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
  PeakMatrix(int totalPixels, double massResolution, double mzLow, double mzHigh, double pxSupport, int *pxSamples_p, int nSamples, char *baseDir, int intMethod=MEAN, int nThreads=1);
    
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
  int getCentroidsIntoRange(double *mass_p, int *px_p, double *intensity_p, int massSize);
  
  //generate de peak matrix with the centroids, their tolerance, and the number of support pixels.  
  int getCentroids();
  
  //Save the peak matrix to the file ".../tmpPeakMatrix.bin"
  int infoToFile(char *baseDir);
  
  //sets the number of pixels per sample based on info in _pixelsCoord.bin
  //result stored in vector m_pixelsSample[]
  //returns the number of samples, or -1 on failure
  int getSamplesPixelNumber();
  
  //It establishes the task distribution lineal for each thread.
  //Each thread receives indices for the start and end segments.
  void setWorkThr(int vectorSize);
  
  
  //Parallel execution. Each mass segment is treated separately.
  //Based on the Gaussians, the centroids within a mass segment are determined.
  // a)An ordered vector containing the Gaussian centers of all pixels is created.
  // b)Iterative averaging, constrained by a tolerance threshold, is used to generate a proposal for centroids.
  // c)This estimate is used as the initialization centroid for the k-means algorithm.
  //   The dataset consists of the centers of Gaussians centered on the initialization value and delimited by the tolerance.
  // d)An intensity value is associated with each centroid for every pixel. 
  //   A pixel contributes the maximum value from the Gaussian nearest to the centroid, provided it falls within the tolerance. 
  //   If a pixel contributes to more than one Gaussian, either the maximum value or the mean value is extracted, depending on the `intMethod` argument.
  //   A representative intensity value is also generated for the centroid. It is the average or maximum value, depending on the intMethod parameter.
  //The centroid information is saved in a file with a name ending in _peakMatrix_n.bin, where n is the thread number.
  void mtSegments(int thrIdx);
  
  
  //Parallel operation. 
  //1) A vector is generated with spacing equal to one-tenth of the tolerance set for the lowest mass. 
  //2) One unit is accumulated in each cell for every Gaussian distribution coinciding with that mass interval, 
  //assuming a unified standard deviation equal to half the tolerance for the mass in question.
  //It serves to delineate regions of interest for subsequent parallel processing. Noise matters.
  void mtMassAxis(int thrIdx);
    
  
  int   m_totalPixels,
        m_nSamples,
        m_nCentroids;
  double 
        m_tolerance, 
        m_mzHigh, 
        m_mzLow, 
        m_pxSupport,
        m_noise;
  
  private:
    
  GAUSS_SP      *m_gaussians_p, 
                *m_centroids_p;
  GAUSS_PARAMS  *m_centroidsGaussians_p;
  MASS_RANGE    *m_massRange_p;
  int           *m_pxSamples_p;
  PIXEL_XY      *m_pixelsCoordinates_p;
  char *        m_baseDir;
  int           m_intMethod;

  MASS_RANGE  *m_MR_p,
              *m_MRthr_p;
  double      *m_massAxis_p,
              m_deltaMass;
  int         m_massAxisSize,
              m_nThreads,
              m_MRsize;
  SPECTRO     m_spectro;
  std::thread  *m_thread_p[MAX_THREADS];
  float       m_status[MAX_THREADS];
  
};
#endif
