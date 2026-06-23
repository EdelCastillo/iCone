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
#include <Rcpp.h>
#include "noiseestimation.h"

//Constructor.
//method: 1=diff(differences), 2=sd(standard deviation), 3=mad(median of absolute deviations).
//smoothing: 0=no, 1=yes.
//winSize=window size for smoothing.
NoiseEstimation::NoiseEstimation(int method, int smoothing, int winSize)
{
  m_hit=true;
  m_method=method;
  m_smoothing=smoothing;
  m_filter_p=0;
  m_winSize=winSize;

  if(method==2 || method==3) //the Gauss filter is required
  {
    if(!(winSize%2)) winSize++; //odd
    int radius=winSize/2;
    
    //memory for the filter
    m_filter_p=new double[winSize];
    double *radius_p=0;
    radius_p  =new double[winSize];
    
    for(int i=-radius; i<=radius; i++) //window for the Gaussian
      radius_p[i+radius]=(double)i;
    double sd=winSize/4; //(width %/% 2) / 2)  standard deviation
    float add=0;
    
    if(smoothing==1) //Gaussian filter
    {
      if(gaussian(radius_p, winSize, 0, sd, m_filter_p)!=0) //the Gaussian is formed
        {m_hit=false; return;}
      for(int i=0; i<winSize; i++) add+=m_filter_p[i];
      if(add==0) add=1e-30;
      for(int i=0; i<winSize; i++) {
        m_filter_p[i]/=add; //normalization
        }
    }
    else {printf("only gaussian filter is implemented\n"); m_hit=false;}
    
    m_winSize=winSize;
    if(radius_p) delete [] radius_p;
  }
}

//destructor
//The memory reserved for the filter is freed.
NoiseEstimation::~NoiseEstimation()
{
  if(m_filter_p) delete [] m_filter_p;  
}

//getNoise()
//spectro_p: Pointer to data.
//size: Size of the data array.
//Returns the estimated noise value of a spectrum.
float NoiseEstimation::getNoise(float *spectro_p, int size)
{
  float estNoise;
  switch(m_method) 
  {
    case 1: estNoise=getNoise_diff (spectro_p, size); break;
    case 2: estNoise=getNoise_sd   (spectro_p, size); break;
    case 3: estNoise=getNoise_mad  (spectro_p, size); break;
  }
  return estNoise;
}

//getSNR()
//Sets the signal-to-noise ratio (SNR) of a spectrum.
//spectro_p: Pointer to data.
//size: Size of the data array.
//Sets the result to the SNR_p pointer.
//Returns 0
float NoiseEstimation::getSNR(float *spectro_p, int size, float *SNR_p)
{
  float estNoise;
  
  switch(m_method) //function pointer initialization.
  {
    case 1: estNoise=getNoise_diff (spectro_p, size); break;
    case 2: estNoise=getNoise_sd   (spectro_p, size); break;
    case 3: estNoise=getNoise_mad  (spectro_p, size); break;
  }
  
  if(estNoise<=0)   for(int i=0; i<size; i++) SNR_p[i]=1e6;
  else
  {
    for(int i=0; i<size; i++) SNR_p[i]=spectro_p[i]/estNoise;
  }
  return estNoise;
}

//getNoise_mad()
//spectro_p: Pointer to data.
//size: Size of the data array.
//Returns the noise estimate using the MAD (median absolute deviation) method.
float NoiseEstimation::getNoise_mad(float *spectro_p, int size)
{
  Common tools;
  double add=0;
  float *smoothing_p=0, *ad_p=0;
  smoothing_p=new float[size];
  ad_p=new float[size];

  //Gaussian filter smoothing
  gaussSmoothing(spectro_p, size, smoothing_p);
 
  for(int i=0; i<size; i++)
    ad_p[i]=abs(spectro_p[i]-smoothing_p[i]); //differences
  
  float central=tools.medianF(ad_p, size); //median of the differences

  for(int i=0; i<size; i++) ad_p[i]=abs(ad_p[i]-central); //abs deviations from the median
  
  //median of the abs deviations from the central value.
  //adjust by a factor for asymptotically normal consistency.
  float noise=1.4826*tools.medianF(ad_p, size); 

  if(smoothing_p) delete [] smoothing_p;
  if(ad_p) delete [] ad_p;
  return noise;
}

//getNoise_sd()
//spectro_p: Pointer to data.
//size: Size of the data array.
//Returns the noise estimate using the standard deviation method.
float NoiseEstimation::getNoise_sd(float *spectro_p, int size)
{
  Common tools;
  double add=0;
  float *smoothing_p=0, *ad_p=0;
  smoothing_p=new float[size];
  ad_p=new float[size];

  //smoothed by Gaussian filter.
  gaussSmoothing(spectro_p, size, smoothing_p);
  
  for(int i=0; i<size; i++)
    ad_p[i]=abs(spectro_p[i]-smoothing_p[i]); //differences
  
  //estimated noise (standard deviation of the differences).
  float noise=(float)sqrt((double)tools.varF(ad_p, size));
  
  if(smoothing_p) delete [] smoothing_p;
  if(ad_p) delete [] ad_p;
  return noise;
}

//getNoise_diff()
//spectro_p: Pointer to data.
//size: Size of the data array.
//Returns the noise estimate using the difference method.
float NoiseEstimation::getNoise_diff(float *spectro_p, int size)
{
  Common tools;
  float *ad_p=0;
  ad_p=new float[size];
  if(size<2) return spectro_p[0];
  
  //differences between neighbors
  for(int i=0; i<size-1; i++)
    ad_p[i]=spectro_p[i+1]-spectro_p[i];
  
  //mean value of the differences
  float A=tools.meanF(ad_p, size-1);
  
  //absolute differences with the mean value
  for(int i=0; i<size; i++)
    ad_p[i]=abs(spectro_p[i]-A);
  
  //estimated noise (mean of the absolute differences)
  float noise=tools.meanF(ad_p, size);
  
  if(ad_p) delete [] ad_p;
  return noise;
}

//gaussSmoothing()
//Convolution of a Gaussian with the signal.
//spectro_p: Pointer to data.
//size: Size of the data array.
//smoothing_p: Pointer to smoothed data.
//Uses the data generated in the constructor for filtering.
//Returns 0
int NoiseEstimation::gaussSmoothing(float *spectro_p, int size, float *smoothing_p)
{
  double add;
  int index;
  int winSize2=(m_winSize-1)/2; //winSize/2 
  for(int i=0; i<size; i++)
  {
    add=0;
    for (int j=0; j<m_winSize; j++)
    {
      index=i-winSize2+j;
      if(index>=0 && index<size) //the extremes are removed
        add+=m_filter_p[j]*spectro_p[index];
    }
    smoothing_p[i]=add;
  }
  return 0;
}

//gaussian()
//Generates a Gaussian over the data x with parameters mean and sd.
//The Gaussian is stored in the gauss pointer.
//Returns -1 if the parameters are unacceptable; 0 if OK.
int NoiseEstimation::gaussian(double *x, int size, double mean, double sd, double *gauss)
{
  if(sd<=0 || gauss==0 || size==0) return -1;
  double sd2=2*sd*sd;
  double factor=1.0/(sqrt(M_PI*sd2)); 
  for(int i=0; i<size; i++)
  {
    gauss[i]=factor*exp(-(x[i]-mean)*(x[i]-mean)/sd2);
  }
  return 0;
}


