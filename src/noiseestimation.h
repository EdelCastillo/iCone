/*********************************************************************************
 *     noiseEstimation
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
#ifndef NOISE_ESTIMATION
#define NOISE_ESTIMATION

#include "types.h"
#include "common_methods.h"

class NoiseEstimation
{
public:
  //Constructor.
  //method: 1=diff(differences), 2=sd(standard deviation), 3=mad(median of absolute deviations).
  //smoothing: 0=no, 1=yes.
  //winSize=window size for smoothing.
  NoiseEstimation(int method=3, int smoothing=1, int winSize=9);
  
  //destructor
  //The memory reserved for the filter is freed.
  ~NoiseEstimation();
  
  //getNoise()
  //spectro_p: Pointer to data.
  //size: Size of the data array.
  //Returns the estimated noise value of a spectrum.
  float getNoise(float *spectro_p, int size);
  
  //getSNR()
  //Sets the signal-to-noise ratio (SNR) of a spectrum.
  //spectro_p: Pointer to data.
  //size: Size of the data array.
  //Sets the result to the SNR_p pointer.
  //Returns the noise estimation
  float getSNR(float *spectro_p, int size, float *SNR_p);

private:  
  //getNoise_mad()
  //spectro_p: Pointer to data.
  //size: Size of the data array.
  //Returns the noise estimate using the MAD (median absolute deviation) method.
  float getNoise_mad(float *spectro_p, int size);
  
  //getNoise_sd()
  //spectro_p: Pointer to data.
  //size: Size of the data array.
  //Returns the noise estimate using the standard deviation method.
  float getNoise_sd(float *spectro_p, int size);
  
  //getNoise_diff()
  //spectro_p: Pointer to data.
  //size: Size of the data array.
  //Returns the noise estimate using the difference method.
  float getNoise_diff(float *spectro_p, int size);
  
  //gaussSmoothing()
  //Convolution of a Gaussian with the signal.
  //spectro_p: Pointer to data.
  //size: Size of the data array.
  //smoothing_p: Pointer to smoothed data.
  //Uses the data generated in the constructor for filtering.
  //Returns 0
  int gaussSmoothing(float *spectro_p, int size, float *smoothing_p);
    
  //gaussian()
  //Generates a Gaussian over the data x with parameters mean and sd.
  //The Gaussian is stored in the gauss pointer.
  //Returns -1 if the parameters are unacceptable; 0 if OK.
  int gaussian(double *x, int size, double mean=0, double sd=1, double *gauss=0);
      
  int     m_method,
          m_winSize,
          m_smoothing;
  double  *m_filter_p;
  bool    m_hit;
};
#endif
