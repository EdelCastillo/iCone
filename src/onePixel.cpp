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
#include "onePixel.h"

/// R METHOD ////////////////////////////////////////////////////////////////////////

 //'
 //'  @name rGetPixelGaussians
 //'  @title It obtains spectrum information from an imzML file.
 //'  
 //'  @param ibdFname  absolute reference to the file with the ibd extension.
 //'  @param imzML  list with information extracted from the imzML file with import_imzML()
 //'  @param params specific parameters
 //'              "SNR": signal-to-noise ratio
 //'   "massResolution": mass resolution with which the spectra were acquired.
 //'      "noiseMethod": method for estimating noise.
 //' "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered.
 //'      "linkedPeaks": two peaks are considered linked if they are closer than the given standard deviation (by defect=3).   
 //'  @param mzLow  lower mass to consider
 //'  @param mzHigh higher mass to consider
 //'  @return List: gaussians, mass intensity, SNR, noise
 //'       gaussians: matrix with parameters for each Gaussian
 //'            mass: vector with the masses of the raw spectrum
 //'       intensity: intensity associated with each mass of the raw spectrum
 //'             SNR: signal-to-noise ratio associated with each mass of the raw spectrum.
 //'           noise: noise estimation
 //'     
 // [[Rcpp::export]]
List rGetPixelGaussians(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, int pixel)
{
  OnePixel onePixel(ibdFname, imzML, params, mzLow, mzHigh, pixel);  
  
  if(onePixel.m_exit==false && pixel+1>=onePixel.m_NPixels) 
  {
    printf("ERROR: the given pixel (%d) exceeds the maximum (%d).\n", pixel+1, onePixel.m_NPixels);
    return 0;
  }
  //loads data from a file and converts its peak into Gaussians.
  onePixel.rawToGaussians();

  return onePixel.getPixelGaussians(pixel, mzLow, mzHigh);
//  printf("getPixelsGaussians..........\n");
//  return ret;
}

//Returns a list with information about a spectrum
// gaussians: matrix with parameters for each Gaussian
// mass: vector with the masses of the raw spectrum
// intensity: intensity associated with each mass of the raw spectrum
// SNR: signal-to-noise ratio associated with each mass of the raw spectrum.
// noise: noise estimation
// require of rGetPixelGaussians()
List OnePixel::getPixelGaussians(int px, float mzLow, float mzHigh)
{
  if(px<0 || px>=m_NPixels) {return 0;}
  if(!m_gaussians.gauss_p)  {return 0;}
  
  NumericMatrix gaussians(m_gaussians.size, 3);
  for(int i=0 ; i<m_gaussians.size; i++) 
  {
    if(m_gaussians.gauss_p[i].mean>= mzLow && m_gaussians.gauss_p[i].mean<= mzHigh)
    {
      gaussians(i,0)=m_gaussians.gauss_p[i].mean;
      gaussians(i,1)=fabs(m_gaussians.gauss_p[i].sigma);
      gaussians(i,2)=m_gaussians.gauss_p[i].weight;
    }
    else
    {
      gaussians(i,0)=0;
      gaussians(i,1)=0;
      gaussians(i,2)=0;
    }
  }
 
  Common utils;
  int iLow=utils.nearestIndex(mzLow, m_spectro.mass_p, m_spectro.size);
  int iHigh=utils.nearestIndex(mzHigh, m_spectro.mass_p, m_spectro.size);
  int size=iHigh-iLow+1;

  NumericVector mass(size);
  NumericVector intensity(size);
  NumericVector SNR(size);
  
  for(int i=0, j=iLow; j<=iHigh; j++, i++) 
  {
    mass(i)=m_spectro.mass_p[j];
    intensity(i)=m_spectro.tmpInt_p [j];//original intensity data
//    intensity(i)=m_spectro.int_p [j]; //intensity data without lower SNR values.
    SNR(i)=m_spectro.SNR_p [j];
  }
  List ret=List::create(Named("gaussians")=gaussians, Named("mass")=mass, Named("intensity")=intensity, 
                        Named("SNR")=SNR, Named("noise")=m_noise);
  return ret;
}


//Constructor
//captures input information, allocates memory and initializes.
////////////////////////////////////////////////////////////////////////////////
OnePixel::OnePixel(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, int pixel)
{
  m_mzLow=mzLow;
  m_mzHigh=mzHigh;
  m_exit=true;
  m_getImzMLData_p=0;
  m_noiseEst_p=0;
  
  //information capture
  Rcpp::DataFrame df;
  df=imzML["run"];
  m_continuous=imzML["continuous_mode"];
  m_NPixels=df.nrows();
  
  //dimensioned
  NumericVector mzLength=df["mzLength"];
  m_spectro.mass_p=0;
  m_spectro.int_p=0;
  m_spectro.SNR_p=0;
  m_spectro.tmpMass_p=0;
  m_spectro.tmpInt_p=0;
  m_spectro.tmpSNR_p=0;
  m_spectro.sort_p=0;
  m_spectro.size=0;
  m_spectro.thread_p=0;
  
  //structure initialization for peak.
  m_peakFG.peakF_p=0;
  m_peakFG.peakU_p=0;
  m_peakFG.peakFsize=0; 
  m_peakFG.peakUsize=0;
  
  m_massSegment.massRange_p=0;
  m_massSegment.segment_p=0;
  m_massSegment.nGaussians_p=0;
  
  m_massRange_p=0;
  
  m_gaussians.gauss_p=0;
  m_gaussians.initGauss=0;
  
  
  printf("\tpixel=%d/%d; mzLow=%.4f; mzHigh=%.4f\n",pixel+1, m_NPixels, m_mzLow, m_mzHigh);
  if(pixel<0 || pixel+1>m_NPixels) {m_exit=false; return;}

  m_mzLength=mzLength[pixel];
  
  m_pixel=pixel;
  
  NumericVector nv;
  CharacterVector cv;
  
  nv=params["SNR"];
  m_SNR=nv[0];
  if(m_SNR<=0) m_SNR=1;
  
  nv=params["massResolution"];
  m_massResolution=nv[0];
  
  cv=params["noiseMethod"];
  String tmpStr=cv[0];
  const char* SNRmethod=tmpStr.get_cstring(); //conversion C
  
  if     (strcmp(SNRmethod, "estnoise_diff")==0) {m_SNRmethod=1; }
  else if(strcmp(SNRmethod, "estnoise_sd")  ==0) {m_SNRmethod=2; }
  else if(strcmp(SNRmethod, "estnoise_mad") ==0) {m_SNRmethod=3; }
  else {m_SNRmethod=0; printf("unknow noise method: %s so, estnoise_mad is used\n", SNRmethod);}
  m_noiseEst_p=new NoiseEstimation(m_SNRmethod, 1, 9);
  
  //memory and its initialization
  //class for accessing imzML files.
  m_getImzMLData_p= new GetImzMLData(ibdFname, imzML);
  
    m_enable=true; //terminates threads if false.

    printf("\ttotal data points of the spectrum=%d\n",m_mzLength);

    //keeps the info of a spectrum, along with the thread that processes it.
    //memory reservation and initialization.
      m_spectro.int_p      =new float[m_mzLength];
      m_spectro.mass_p     =new float[m_mzLength];
      m_spectro.SNR_p      =new float[m_mzLength];
      m_spectro.tmpMass_p  =new float[m_mzLength];
      m_spectro.tmpInt_p   =new float[m_mzLength];
      m_spectro.tmpSNR_p   =new float[m_mzLength];
      m_spectro.sort_p     =new int  [m_mzLength];
      m_spectro.size=0;
 //printf("...end constructor %d\n", m_mzLength);   
}

//destructor
//free reserved memory
OnePixel::~OnePixel()
{
  //  printf("init PeakMatrix destructor\n");
  if(m_gaussians.gauss_p) delete[] m_gaussians.gauss_p;
  
  //class for accessing imzML files
  if(m_getImzMLData_p) {delete m_getImzMLData_p; m_getImzMLData_p=0;}
  
  if(m_peakFG.peakF_p) 
  {delete [] m_peakFG.peakF_p; m_peakFG.peakF_p=0;}
  if(m_peakFG.peakU_p) 
  {delete [] m_peakFG.peakU_p; m_peakFG.peakU_p=0;}
  
  //spectrum info
  if(m_spectro.int_p)      {delete []m_spectro.int_p;    m_spectro.int_p=0;}
  if(m_spectro.mass_p)     {delete []m_spectro.mass_p;   m_spectro.mass_p=0;}
  if(m_spectro.SNR_p)      {delete []m_spectro.SNR_p;    m_spectro.SNR_p=0;}
  if(m_spectro.tmpMass_p)  {delete []m_spectro.tmpMass_p;m_spectro.tmpMass_p=0;}
  if(m_spectro.tmpInt_p)   {delete []m_spectro.tmpInt_p; m_spectro.tmpInt_p=0;}
  if(m_spectro.tmpSNR_p)   {delete []m_spectro.tmpSNR_p; m_spectro.tmpSNR_p=0;}
  if(m_spectro.sort_p)     {delete []m_spectro.sort_p;   m_spectro.sort_p=0;}
  //  printf("end PeakMatrix destructor\n");
  
}

//Loads the raw data from a file, separates the peaks, and establishes the Gaussian.
//results in internal structure. 
//return the number of Gaussians or -1 is KO
int OnePixel::rawToGaussians()
{
  IntensityPeak intPeak(m_SNR);
  Common common;
  GAUSS_PARAMS *gaussians_p=0;
  double *centroids_p=0;
  int *centroidsIndex_p=0;
 
    //Spectrum conditioning.
    //Full spectra are received and only part of it may be of interest.
    int iMzLow=-1, iMzHigh=-1;
    int massSize=m_getImzMLData_p->getPixelMassF(m_pixel, m_spectro.tmpMass_p);//mass vector
    m_spectro.size=massSize;
    if(massSize<=0) {return 0;}
    
    //raw info of intensities  
    m_getImzMLData_p->getPixelIntensityF(m_pixel, m_spectro.tmpInt_p);

    m_spectro.pixel=m_pixel;
    
    int spSize=m_spectro.size; //spectrum size
    //SNR
    m_noise=m_noiseEst_p->getSNR(m_spectro.tmpInt_p, m_spectro.size,  m_spectro.tmpSNR_p);
    m_spectro.noise=m_noise;
    
    //the spectrum is limited to the range of interest.
    iMzLow =common.nearestIndex(m_mzLow,  m_spectro.tmpMass_p, spSize); //low index
    iMzHigh=common.nearestIndex(m_mzHigh, m_spectro.tmpMass_p, spSize); //high index
    spSize=iMzHigh- iMzLow +1;
//fault control
    if(iMzLow<0 || iMzHigh>m_spectro.size-1 || iMzHigh-iMzLow>m_spectro.size || iMzHigh<iMzLow) 
    {
      return -1;
    }
    if(spSize<=1) 
    {
      printf("There is no information in the requested mass range (%.4f/%.4f) for the pixel %d.\nits mass range is %.4f/%.4f\n", 
             m_mzLow, m_mzHigh, m_pixel+1, m_spectro.tmpMass_p[0], m_spectro.tmpMass_p[massSize-1]);
      return 0;
    }
    //mass ordination.
    common.sortUpF(m_spectro.tmpMass_p+iMzLow, m_spectro.sort_p, spSize);

    //The part of interest is extracted from the mass, intensity and SNR vectors.
    for(int i=0; i<spSize; i++) 
    {
      int index=iMzLow+m_spectro.sort_p[i];
      m_spectro.mass_p[i]=m_spectro.tmpMass_p[index];
      m_spectro.int_p [i]=m_spectro.tmpInt_p [index];
      m_spectro.SNR_p [i]=m_spectro.tmpSNR_p [index];
    }
    for(int i=0; i<spSize; i++) 
      m_spectro.tmpInt_p [i]=m_spectro.int_p [i];//getPeakList() destroys intensity data.
    
    m_spectro.size=spSize; //useful range
    //-----------------------------------------------------------------    

    //conversion to Gaussians.
    
    int nPeak;
    
    //the peak are extracted from the spectrum (they are delimited by their indices).
    nPeak=intPeak.getPeakList(&m_spectro);
    
    m_gaussians.size=0;
    if(nPeak>0) //if there are peak to treat
    {
      //single peak info is saved.
      m_peakFG.peakF_p=new ION_INDEX[nPeak];
      for(int i=0; i<nPeak; i++) //copy peak
      {
        m_peakFG.peakF_p[i].low =intPeak.getSinglePeak(i).low;
        m_peakFG.peakF_p[i].max =intPeak.getSinglePeak(i).max;
        m_peakFG.peakF_p[i].high=intPeak.getSinglePeak(i).high;
      }
      m_peakFG.peakFsize=nPeak;//number of simple peak
      
      //the information of compound peak (simple joined-overlapping peak) is obtained.
      int nUPeak=intPeak.getCompoundPeakNumber(); 
     
      //Compound peak information is saved.
      //Each entry is a reference to the initial and final single peak of the merged peak.
      m_peakFG.peakU_p=new PEAK_UNITED[nUPeak];
      for(int i=0; i<nUPeak; i++)
      {
        m_peakFG.peakU_p[i].low =intPeak.getCompoundPeak(i).peakLow;
        m_peakFG.peakU_p[i].high=intPeak.getCompoundPeak(i).peakHigh;
      }
      m_peakFG.peakUsize=nUPeak;//number of compound peak
      
      //the peak are converted to Gaussians
      gaussians_p=new GAUSS_PARAMS[nPeak]; //temporal copy of Gaussians
      centroids_p=new double[nPeak];       //temporary mass copy
      centroidsIndex_p=new int[nPeak];     //indices to ordered masses
      
      int nGaussians=getGaussians(&m_spectro, gaussians_p);
     
      if(nGaussians>0 && nGaussians<=nPeak)//if Gaussians exist
      {
        //sorted from lowest to highest (some joined peak may be altered).
        for(int i=0; i<nGaussians; i++)
        {
          centroids_p[i]=gaussians_p[i].mean;
        }
        
        if(common.sortUp(centroids_p, centroidsIndex_p, nGaussians)==0) //the centers are ordered
        {
          //ordered copy to main structure
          for(int i=0; i<nGaussians; i++)
          {
            int index=centroidsIndex_p[i];
            m_gaussians.gauss_p[i].mean  =gaussians_p[index].mean;
            m_gaussians.gauss_p[i].sigma =fabs(gaussians_p[index].sigma);
            m_gaussians.gauss_p[i].weight=gaussians_p[index].weight;
          }
          m_gaussians.size=nGaussians;
        }
        
      }
      else m_gaussians.size=0;
      
    } //end of peak processing
    printf("\t#peaks=%d; #gaussians=%d; noise estimation=%.3f\n", nPeak, m_gaussians.size, m_spectro.noise);

    if(gaussians_p)     {delete [] gaussians_p;       gaussians_p=0;}
    if(centroids_p)     {delete [] centroids_p;       centroids_p=0;}
    if(centroidsIndex_p){delete [] centroidsIndex_p;  centroidsIndex_p=0;}
    return m_gaussians.size;
}
  

//getGaussians()
//Sets the Gaussians on the peak.
//Uses the peak separation information (m_peakFG_p).
//spectro: pointer to the spectrum.
//gaussians_p: pointer to the structure containing the Gaussians' parameters.
//Returns the number of Gaussians or a value < 0 on failure.
int OnePixel::getGaussians(SPECTRO *spectro_p, GAUSS_PARAMS *gaussians_p)
{
  float *intSpectrum_p=spectro_p->int_p, *massSpectrum_p=spectro_p->mass_p;
  int intSize=spectro_p->size;
  
  float minMeanPxMag=spectro_p->noise*m_SNR; //minimum value to consider a peak as valid.
  //class for conversion to Gaussians.
  GmmPeak gmmPeak(minMeanPxMag);
  
  int gaussIndex=0; //indices for each Gaussian.
  int nUPeak=m_peakFG.peakUsize; //#united peak
  
  //memory for predictable Gaussians.
  m_gaussians.gauss_p=new GAUSS_PARAMS[m_peakFG.peakFsize];
  
  //for each set of joined peak.
  for(int uPeak=0; uPeak<nUPeak; uPeak++)
  {
    int mPeakLow   =m_peakFG.peakU_p[uPeak].low;      //lower simple peak.
    int mPeakHigh  =m_peakFG.peakU_p[uPeak].high;     //upper simple peak.
    int lowMzIndex =m_peakFG.peakF_p[mPeakLow].low;   //lower  mass index.
    int highMzIndex=m_peakFG.peakF_p[mPeakHigh].high; //higher mass index.
    int mzSize=highMzIndex-lowMzIndex+1;
    
    int nPeak=mPeakHigh-mPeakLow+1; //number of simple peak into united peak.
    
    //The information of simple peak of magnitude that make up the segment is established.
    gmmPeak.setPeak(&m_peakFG.peakF_p[mPeakLow], nPeak); 
    
    //The magnitude information that must be adjusted is established.
    GROUP_F pxMag;
    pxMag.set=intSpectrum_p+lowMzIndex; 
    pxMag.size=mzSize;
    gmmPeak.setMagnitudes(&pxMag); //magnitude peak.
    
    //The simple Gaussians are formed whose sum reproduces the magnitude peak.
    if(gmmPeak.gmmDeconvolution()<0) //deconvolution
    {
      printf("ERROR: limits exceeded. The aim is to deconvolve a peak composed of more than %d Gaussians.\n", DECONV_MAX_GAUSSIAN);
      return -4;
    }
    if(intSize<mzSize) //if the dimensions are not correct.
    {
      //warning: the dimensions of mzAxis and data[magnitudes] must match. 
      printf("ERROR: \n");
      return -5;
    }
    
    int nGauss=gmmPeak.getDeconvNumber(); //# gaussians
    
    bool hit=true;
    GAUSSIAN gaussIn, gaussOut; //if it is required to adapt the mass axis.
    for(int g=0; g<nGauss; g++)
    {
      gaussIn=gmmPeak.getDeconv(g); //get a gaussian.
      
      //unit conversion (scans to Daltos).
      gmmPeak.gaussConversion(&gaussIn, &gaussOut, massSpectrum_p+lowMzIndex, mzSize); 
      
      if(gaussIndex>=m_peakFG.peakFsize) //control de error
      {
        //warning: memory overflow for gaussians 
        printf("warning: memory overflow for gaussians %d/%d\n", gaussIndex, m_peakFG.peakFsize); 
        hit=false; break;
      }
      
      //the Gaussians are saved.
      gaussians_p[gaussIndex].mean  =(double)gaussOut.mean;
      gaussians_p[gaussIndex].sigma =(double)gaussOut.sigma;
      gaussians_p[gaussIndex].weight=(double)gaussOut.yFactor*gaussOut.weight;
      gaussIndex++;
    }//end of gaussians
    
    if(hit==false) break;//fallo
  }//end uPeak loop
  return gaussIndex; //# gaussians
}
