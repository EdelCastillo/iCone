/*************************************************************************
 *     types
 *     iCone - R package for MSI data processing
 *     Copyright (C) 2025 Esteban del Castillo Pérez
 *     esteban.delcastillo@urv.cat
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
 **************************************************************************/
#ifndef ICONE_TYPES
#define ICONE_TYPES
  
#include <thread>
#include <mutex>
#include "GMM.h"

#define MAX_THREADS 50
#define MAX_SAMPLES 50

    typedef struct{
        int x, y;
    }PIXEL_XY;                  //2D pixel

  typedef struct{
        PIXEL_XY *set;
        int      size;
    }GROUP_XY;                  //pixel array


     typedef struct xGROUP{
        int 	size;
        int 	*set;
        struct xGROUP 	*group;  //group of chained integer elements
    }GROUP;

    typedef struct GROUP_F{
        int 	size;
        float 	*set;
        struct GROUP_F 	*group; //group of chained real elements
    }GROUP_F;

    typedef struct {
        int 	*set;          //items
        int 	size;          //group size.
    }LGROUP;                 //group

    typedef struct
    {
      float   low,          //low mass
              high;         //high mass
      double  resolution;   //resolution with which the range is resolved.
      int     nGaussians;   //Gaussians in range
    }MASS_RANGE;
    
    
    typedef struct
    {
      int low, high;        //low and high indices to simple peaks that make up a compound peak.
    }PEAK_UNITED;           //information about united peaks
    
    typedef struct 
    {
      int low, high, max;   //indexes to elements left valley, right valley and maximum.
      bool confidence;      //true if the magnitude variation between valleys exceeds m_peaksGap.
    }ION_INDEX;             //indices to the magnitude_p[] array that delimit an ion.
    
    typedef struct
    {
      ION_INDEX *peakF_p;  //simple peak indices
      PEAK_UNITED *peakU_p;//compound peak indices
      int peakFsize, peakUsize; //sizes 
    }PEAK_F_GROUP;          //complete info on peaks.
    
    typedef struct
    {
      GAUSS_PARAMS *gauss_p;  //pointer to Gaussian parameters.
      int size;    //array size
    }GAUSS_SP;
    
    typedef struct SEGMENT
    {
      float   lowMass,        //lower mass of the segment
              highMass;       //upper mass of the segment
      int     nGaussians;     //number of Gaissians in the segment
      SEGMENT *next;          //pointer to the next link
    }SEGMENT;                 //chain link for segment info.
    
    typedef struct
    {
      MASS_RANGE *massRange_p;  //extreme masses
      SEGMENT    **segment_p;   //pointer to segment string
      int        *nGaussians_p; //pointer to the number of Gaussians in each link.
    }MASS_SEGMENT;
    
    typedef struct ION_ENTRY
    {
      float 	*set,          //items
              mass,
              massResolution;
      int 	  size;          //group size.
      struct  ION_ENTRY *group;
    }ION_ENTRY;             //chain link for ion info.
    
    typedef struct
    {
      float *int_p,         //intensity array
            *mass_p,        //mass array
            *SNR_p,         //signal-to-noise ratio array
            *tmpInt_p,      //array of temporal intensities
            *tmpMass_p,     //temporary mass array
            *tmpSNR_p,      //temporal signal-to-noise ratio array
            noise;
          
      int   *sort_p,        //indexes to ordered information
            size,           //spectrum size
            pixel;          //pixel associated with the spectrum.
      bool  flag;           
      std::thread *thread_p;                //thread that processes the spectrum.
      std::mutex  *mutexIn_p, *mutexOut_p;  //mutexes associated with the thread.
    }SPECTRO;               //info for a spectrum
    
    typedef struct          
    {                       //If it is a peak composed of a single simple peak, both entries coincide.
      int   peakLow,        //index to the lower peak of the compound peak.
            peakHigh;       //index to the upper peak of the composite peak.
    }UNITED_PEAK;          //info about a compound peak.
    
    typedef struct                  
    {
      ION_INDEX       *intPeak_p;//simple magnitude peak.
      int             nIntPeak;  //amount
      UNITED_PEAK    *uPeak_p;  //joined (compound) peaks
      int             nUPeak;    //amount
    }PEAK_LIST;                  //complete peak info.
    
    
    
#endif

