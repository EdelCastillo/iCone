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
#include "peakMatrix.h"

int  gPeakCount=0, gSpectra=0, gProcSegments, gVez, gError=0; //observers.
float gMinMass, gMaxMass;
std::mutex gMutex;
//FILE *fp;

/// R METHOD ////////////////////////////////////////////////////////////////////////

//'
//'  @name peakMatrix
//'  @title converts the info in the imzML file into a peak matrix.
//'  
//'  @param ibdFname:  absolute reference to the file with the ibd extension.
//'  @param imzML:     list with information extracted from the imzML file with import_imzML()
//'  @param params:    specific parameters
//'              "SNR": signal-to-noise ratio
//'   "massResolution": mass resolution with which the spectra were acquired.
//'      "noiseMethod": method for estimating noise.
//' "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered.
//'      "linkedPeaks": two peaks are considered linked if they are closer than the given standard deviation (by defect=3).   
//'  @param mzLow:    lower mass to consider
//'  @param mzHigh:   higher mass to consider
//'  @param pxList:   list of pixels. First pixel=1. By default everyone.
//'  @param nThreads: number of threads suggested for parallel processing.
//'  @return lista: peakMatrix, massVector, massResolution, pixelsSupport
//'     peakMatrix: matrix of centroids and the intensity associated with each pixel.
//'     massResolution: 
//'     massVector: the mz associated with each column of peakmatrix.
//'     pixelsSupport: number of pixels with intensity >= minPixelsSupport
//'     
// [[Rcpp::export]]
List peakMatrix(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, Rcpp::NumericVector pxList, int nThreads)
{
  gPeakCount=0, gSpectra=0;
//printf("...%.0f\n", pxList[0]); 

  PeakMatrix pMatrix(ibdFname, imzML, params, pxList, mzLow, mzHigh, nThreads);
  if(pMatrix.m_hit==false) return 0;
 
  //Phase 1:
  //loads data from a file and converts its peak into Gaussians.
  int ret1=pMatrix.rawToGaussians(); //parallel processing
  if(ret1<0) 
    {return 0;}
  if(gPeakCount<=0) 
   {printf("No peaks are detected in the sample.\n"); return 0;}

  //Phase 2:
  float linkedPeaks=pMatrix.m_linkedPeaks;
  pMatrix.m_massRangeSize=pMatrix.getMassRanges(&linkedPeaks);
  if(pMatrix.m_massRangeSize==0) 
    {
    printf("There are no valid mass segments.\n"); 
    return 0;
    }
  if(pMatrix.m_massRangeSize<0) 
  {
    printf("The segmentation stage could not be resolved. There are segments larger than 3 Da. It is suggested to increase the SNR parameter.\n"); 
    return 0;
  }
  if(linkedPeaks!=pMatrix.m_linkedPeaks)
  {
    printf("\t\t\tsegmentation stage resolved with sigma factor=%.3f\n", linkedPeaks);
  }
  
  pMatrix.freeMemoryPeak(); //frees up unnecessary memory.

  MASS_RANGE massRange;
  massRange.low=mzLow;
  massRange.high=mzHigh;
  
  //Phase 3:
  //From the Gaussians, the centroids are obtained by kmeans segmentation.
  List ret=pMatrix.massRangeToCentroids(massRange); //parallel processing
  
  printf("\t#valid spectra=%d; #max Gaussians in a spectra=%d; #total spectrum peak in mass range=%d\n", 
         gSpectra, pMatrix.m_maxPxGaussians, gPeakCount);  
  printf("\tminimum/maximun sample mass=%.4f/%.4f (Da)\n", gMinMass, gMaxMass);
  
  return ret;
}


   //'  @name rGetAverageGaussianSpectrum
   //'  @title converts the info in the imzML file into two arrays: Gaussians average values versus masses
   //'  
   //'  @param ibdFname:  absolute reference to the file with the ibd extension.
   //'  @param imzML:     list with information extracted from the imzML file with import_imzML()
   //'  @param params:    specific parameters
   //'              "SNR": signal-to-noise ratio
   //'   "massResolution": mass resolution with which the spectra were acquired.
   //'      "noiseMethod": method for estimating noise.
   //'  @param mzLow:        lower  mass  to consider
   //'  @param mzHigh:       higher mass  to consider
   //'  @param pxList:       list of pixels. First pixel=1. By default everyone.
   //'  @param overSampling: interval between points on the mass axis = massResolution/overSampling.
   //'  @param nThreads:     number of threads suggested for parallel processing.
   //'  @return lista: averageMz and averageIntensity
   //'     averageMz: array of masses at intervals of 1/4 of the resolution
   //'     averageIntensity: array of average values with all Gaussians 
   //'     
// [[Rcpp::export]]
List rGetAverageGaussianSpectrum(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, Rcpp::NumericVector pxList, float overSampling, int nThreads)
  {
  gPeakCount=0, gSpectra=0;
  PeakMatrix pMatrix(ibdFname, imzML, params, pxList, mzLow, mzHigh, nThreads);  
  if(pMatrix.m_hit==false) return 0;
  
  //loads data from a file and converts its peak into Gaussians.
  int ret1=pMatrix.rawToGaussians(); //parallel processing
  if(ret1<0) 
  {return 0;}
  if(gPeakCount<=0) 
  {printf("No peaks are detected in the sample.\n"); return 0;}

  List ret=pMatrix.getMeanGaussianSpectrum(pMatrix.m_massResolution,(int)overSampling);
  return ret;
}
 
 
   //'  @name rGetAverageSpectrum
   //'  @title converts the info in the imzML file into two arrays: average values versus masses
   //'  
   //'  @param ibdFname:  absolute reference to the file with the ibd extension.
   //'  @param imzML:     list with information extracted from the imzML file with import_imzML()
   //'  @param params:    specific parameters not considered but must exist.
   //'  @param mzLow:        lower  mass  to consider
   //'  @param mzHigh:       higher mass  to consider
   //'  @param pxList:       list of pixels. First pixel=1. By default everyone.
   //'  @param overSampling: interval between points on the mass axis = massResolution/overSampling.
   //'  @return lista: averageMz and averageIntensity
   //'     averageMz: array of masses at intervals of mass resolution/overSampling
   //'     averageIntensity: array of average values with all Gaussians 
   //'     
   // [[Rcpp::export]]
   List rGetAverageSpectrum(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, Rcpp::NumericVector pxList, float overSampling)
   {
     PeakMatrix pMatrix(ibdFname, imzML, params, pxList, mzLow, mzHigh, 1);  
     if(pMatrix.m_hit==false) return 0;

     List ret=pMatrix.getMeanSpectrum(pMatrix.m_massResolution, (int)overSampling);
     return ret;
   }


//Constructor
//captures input information, allocates memory and initializes.
////////////////////////////////////////////////////////////////////////////////
PeakMatrix::PeakMatrix(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, Rcpp::NumericVector pxList, float mzLow, float mzHigh, int nThreads)
{
  m_hit=true;
  m_mzLow=mzLow;
  m_mzHigh=mzHigh;
  gMinMass=1e30;
  gMaxMass=0;
  m_getImzMLData_p=0;
  m_noiseEst_p=0;
  m_gaussians_p=0;
  m_maxMzLength=0;
  m_peakFG_p=0;
  m_massSegment.massRange_p=0;
  m_massSegment.segment_p=0;
  m_massSegment.nGaussians_p=0;
  m_ionEntry_p=0;
  m_maxPxGaussians=0;
  m_massRange_p=0;
  
  //information capture
  Rcpp::DataFrame df;
  df=imzML["run"];
  m_continuous=imzML["continuous_mode"];
  m_pxList=0;
  NumericVector nv;
  CharacterVector cv;
  m_pxMax=df.nrows()-1;
  m_pxMin=0;

  if(pxList.size()==1 && pxList[0]==-1) //by defect, all pixels
    {
    m_NPixels=df.nrows();
    m_pxList= new int[m_NPixels];
    for(int i=0; i<m_NPixels; i++) m_pxList[i]=i;
    }
  else //only pixels passed
    {
    m_NPixels=pxList.size();
    m_pxList= new int[m_NPixels];
    int pxMin=0x7FFFFFFF, pxMax=0;
    
    for(int i=0; i<m_NPixels; i++)
    {
      if(pxList[i]<pxMin) pxMin=pxList[i];
      if(pxList[i]>pxMax) pxMax=pxList[i];
      m_pxList[i]=pxList[i];
    }
    if(pxMin<0 || pxMax>m_pxMax)
    {
      printf("ERROR: some pixels are out of range(%d/%d).", 1, m_pxMax+1);
      m_hit=false;
      return;
    }
    m_pxMin=pxMin;
    m_pxMax=pxMax;
    }
  printf("pixel Minimun:%d; pixel Maximun:%d; total:%d\n", m_pxMin+1, m_pxMax+1, m_NPixels);
  
  //  m_NPixels=1000;
  printf("to process: #spectra=%d in mass range(Da)=%.4f to %.4f\n",m_NPixels, m_mzLow, m_mzHigh);

  //copy pixel coordinates
  m_x=0, m_y=0;
  m_x=new int[m_NPixels];
  m_y=new int[m_NPixels];
  NumericVector X=df["x"];
  NumericVector Y=df["y"];
  
  for(int i=0, px; i< m_NPixels; i++)
  {
    px=m_pxList[i];
    m_x[i]=X[px];
    m_y[i]=Y[px];
  }

  nv=params["SNR"];
  m_SNR=nv[0];
  if(m_SNR<=0) m_SNR=1;

  nv=params["massResolution"];
  m_massResolution=nv[0];

  nv=params["minPixelsSupport"];
  m_pxSupport=nv[0]*m_NPixels/100.0;
  
  cv=params["noiseMethod"];
  String tmpStr=cv[0];
  const char* SNRmethod=tmpStr.get_cstring(); //conversion C
  
  //Two peaks are considered linked if they are closer than the given standard deviation.
  nv=params["linkedPeaks"]; 
  m_linkedPeaks=nv[0];
  
  if     (strcmp(SNRmethod, "estnoise_diff")==0) {m_SNRmethod=1; }
  else if(strcmp(SNRmethod, "estnoise_sd")  ==0) {m_SNRmethod=2; }
  else if(strcmp(SNRmethod, "estnoise_mad") ==0) {m_SNRmethod=3; }
  else {m_SNRmethod=0; printf("unknow noise method: %s so, estnoise_mad is used\n", SNRmethod);}
  m_noiseEst_p=new NoiseEstimation(m_SNRmethod, 1, 9);

  //class for accessing imzML files.
  m_getImzMLData_p= new GetImzMLData(ibdFname, imzML);
  
  //memory and its initialization
  m_gaussians_p=new GAUSS_SP[m_NPixels];
  for(int px=0; px<m_NPixels; px++)
  {
    m_gaussians_p[px].gauss_p=0;
    m_gaussians_p[px].initGauss=0;
  }
  
  //dimensioned
  NumericVector mzLength=df["mzLength"];
  for(int i=0, px; i<m_NPixels; i++)
  {
    px=m_pxList[i];
    if(mzLength[px]>m_maxMzLength)m_maxMzLength=mzLength[px]; //maximum spectrum
  }

  //vectors to accommodate a spectrum and its masses.
  for(int i=0; i< MAX_THREADS; i++)
    {
      m_spectro[i].mass_p=0;
      m_spectro[i].int_p=0;
      m_spectro[i].SNR_p=0;
      m_spectro[i].tmpMass_p=0;
      m_spectro[i].tmpInt_p=0;
      m_spectro[i].tmpSNR_p=0;
      m_spectro[i].sort_p=0;
      m_spectro[i].size=0;
      m_spectro[i].thread_p=0;
      m_spectro[i].mutexIn_p=0;
      m_spectro[i].mutexOut_p=0;;
    }
  
  m_enable=true; //terminates threads if false.
  
  //parameter control.
  m_nThreads =thread::hardware_concurrency()-1; //a core is released
  if(nThreads<m_nThreads && nThreads>0)
    m_nThreads =nThreads;
  if(m_nThreads<=0) m_nThreads=1;
  if(m_nThreads>MAX_THREADS) m_nThreads=MAX_THREADS;
  if(m_nThreads> m_NPixels) m_nThreads=m_NPixels;
  
  printf("maximum data points of a spectrum: %d\n",m_maxMzLength);
  printf("Threads to use: %d\n", m_nThreads);

  //keeps the info of a spectrum, along with the thread that processes it.
  //memory reservation and initialization.
  for(int i=0; i<m_nThreads; i++)
  {
    m_spectro[i].int_p      =new float[m_maxMzLength];
    m_spectro[i].mass_p     =new float[m_maxMzLength];
    m_spectro[i].SNR_p      =new float[m_maxMzLength];
    m_spectro[i].tmpMass_p  =new float[m_maxMzLength];
    m_spectro[i].tmpInt_p   =new float[m_maxMzLength];
    m_spectro[i].tmpSNR_p   =new float[m_maxMzLength];
    m_spectro[i].sort_p     =new int  [m_maxMzLength];
    m_spectro[i].mutexIn_p  =new std::mutex;
    m_spectro[i].mutexOut_p =new std::mutex;
    m_spectro[i].size=0;
    m_spectro[i].mutexIn_p ->lock();
    m_spectro[i].mutexOut_p->lock();
    m_spectro[i].thread_p=new std::thread(&PeakMatrix::mtGetGaussians, this, i);
  }
  gMutex.unlock();
  gProcSegments=0;
  gVez=1;
  gError=0;
  
  for(int i=0; i<m_nThreads; i++)
    m_totalIons[i]=0;
    
  //structure initialization for peak.
  m_peakFG_p= new PEAK_F_GROUP[m_NPixels];
  for(int px=0; px<m_NPixels; px++)
  {
    m_peakFG_p[px].peakF_p=0;
    m_peakFG_p[px].peakU_p=0;
    m_peakFG_p[px].peakFsize=0; 
    m_peakFG_p[px].peakUsize=0;
  }

//  fp=fopen("kk.txt", "w");
}


//destructor
//free reserved memory
PeakMatrix::~PeakMatrix()
{
   //printf("init PeakMatrix destructor\n");
 if(m_gaussians_p)
  {
    for(int px=0; px<m_NPixels; px++)
    {
      if(m_gaussians_p[px].gauss_p) delete[] m_gaussians_p[px].gauss_p;
    }
    delete []m_gaussians_p;
  }

  freeMemoryPeak();

 for(int i=0; i<m_nThreads; i++)
 {
    //Error if the comment is removed.???
    //if(m_spectro[i].mutexIn_p)  {delete m_spectro[i].mutexIn_p; m_spectro[i].mutexIn_p=0;}
    if(m_spectro[i].mutexOut_p) {delete m_spectro[i].mutexOut_p; m_spectro[i].mutexOut_p=0;}
 }

  if(m_massSegment.massRange_p)     delete []m_massSegment.massRange_p;
  if(m_massSegment.nGaussians_p)    delete []m_massSegment.nGaussians_p;

  if(m_massRange_p) delete [] m_massRange_p;
  if(m_noiseEst_p)  delete m_noiseEst_p;
  if(m_x)           delete [] m_x;
  if(m_y)           delete [] m_y;
  if(m_pxList)      delete [] m_pxList;
  
  //if(fp) fclose(fp);
  
  //printf("end PeakMatrix destructor\n");
}

//freeing buffer.
void PeakMatrix::freeMemoryPeak()
{
  //class for accessing imzML files
  if(m_getImzMLData_p) {delete m_getImzMLData_p; m_getImzMLData_p=0;}

  if(m_peakFG_p) //if the peak structure exists
  {
    for(int px=0; px<m_NPixels; px++)
    {
      if(m_peakFG_p[px].peakF_p) 
        {delete [] m_peakFG_p[px].peakF_p; m_peakFG_p[px].peakF_p=0;}
      if(m_peakFG_p[px].peakU_p) 
        {delete [] m_peakFG_p[px].peakU_p; m_peakFG_p[px].peakU_p=0;}
    }
    delete []m_peakFG_p; m_peakFG_p=0;
  }
 
  //spectrum info
  for(int i=0; i<m_nThreads; i++)
   {
    if(m_spectro[i].int_p)      {delete []m_spectro[i].int_p;    m_spectro[i].int_p=0;}
    if(m_spectro[i].mass_p)     {delete []m_spectro[i].mass_p;   m_spectro[i].mass_p=0;}
    if(m_spectro[i].SNR_p)      {delete []m_spectro[i].SNR_p;    m_spectro[i].SNR_p=0;}
    if(m_spectro[i].tmpMass_p)  {delete []m_spectro[i].tmpMass_p;m_spectro[i].tmpMass_p=0;}
    if(m_spectro[i].tmpInt_p)   {delete []m_spectro[i].tmpInt_p; m_spectro[i].tmpInt_p=0;}
    if(m_spectro[i].tmpSNR_p)   {delete []m_spectro[i].tmpSNR_p; m_spectro[i].tmpSNR_p=0;}
    if(m_spectro[i].sort_p)     {delete []m_spectro[i].sort_p;   m_spectro[i].sort_p=0;}
   }

  //the memory reserved for the ion chain is released.
  if(m_ionEntry_p)
    {
      ION_ENTRY *ionEntry2_p, *ionEntry_p=m_ionEntry_p[0];
      while(ionEntry_p)
      {
        ionEntry2_p=ionEntry_p->group;
        if(ionEntry_p->set) {delete[] ionEntry_p->set;  ionEntry_p->set=0;}
        if(ionEntry_p)      {delete ionEntry_p;         ionEntry_p=0;}
        ionEntry_p=ionEntry2_p;
      }
      delete [] m_ionEntry_p; m_ionEntry_p=0;
    }
  
}

//rawToGaussians
//gets the intensity peak and converts them into Gaussians.
//The intensity and mass data adjusted to the range of interest are loaded from the imzML file.
//the SNR info is established for each point of the spectra.
//The generated information is stored in the m_spectro_p structure.
/////////////////////////////////////////////////////////////////////////////////////////////
int PeakMatrix::rawToGaussians()
{
  printf("Processing \n\tphase 1:\traw to gaussians(%%): 00 ");
  int vez=1;
  Common common;
  int spSize=0;
  bool hit=false;
  NumericMatrix SNR_Mx;
  NumericVector noiseSize;
  int px;
  
  //for all pixels (spectra)
  //For each spectrum, determine the intensity peak, the joined peak, and their Gaussians.
  for(int iPx=0; iPx<m_NPixels; ) 
    {
    //indication of the progress of the process (10% resolution)
    if((float)iPx/(float)m_NPixels>vez*0.1) {if(vez<10) printf("%d ", vez*10); vez++;}
    
    //We load as many spectra as threads are used.
    //Each spectrum is processed by a thread.
    hit=true;
    int nThr=0;
    for(int thr=0; thr<m_nThreads; thr++) //for each thread
      {
      
      while(iPx<m_NPixels) //iterates while the spectrum has length <=2 (eliminates empty spectra).
        {
        spSize=getRawInfo(iPx, thr);//capturing spectra from imzML file
        iPx++;
        //If the spectrum size <=2 its peak are not processed.
        if(spSize>2) break;//spectrum for analysis. There is no information of interest in the pixel spectrum <=3 peak.
        }
      if(spSize>2 && iPx<=m_NPixels) 
        {
        nThr++;
        {m_spectro[thr].mutexIn_p->unlock();} //This thread is allowed to run.
        }
      if(iPx==m_NPixels) break;
      }
    //wait for the conclusion of all threads.
    for(int thr=0; thr<nThr; thr++) //for each thread
    {
      m_spectro[thr].mutexOut_p->lock();
    }

    int tmp=0;
    gMutex.lock();
    tmp=gError; //error count
    gMutex.unlock();
    if(tmp>10) 
    {
      printf("Application aborted because there are too many warnings. It is suggested to increase the SNR parameter.\n"); 
      //thread removal
      m_enable=false; //threads end.
      for(int thr=0; thr<m_nThreads; thr++) //are unlocked for the conclusion.
        m_spectro[thr].mutexIn_p->unlock();
      
      for(int i=0; i<m_nThreads; i++)
      {
        m_spectro[i].thread_p->join();
        delete m_spectro[i].thread_p; m_spectro[i].thread_p=0;
      }
      return -1;
    }
  }

  //maximum amount of Gaussians in the spectra.
  for(int px=0; px<m_NPixels; px++)//for all pixels
  {
    if(m_gaussians_p[px].gauss_p && m_gaussians_p[px].size>m_maxPxGaussians) 
      m_maxPxGaussians=m_gaussians_p[px].size;
  }

  //mass range of interest
  MASS_RANGE mRange;
  mRange.low=m_mzLow;
  mRange.high=m_mzHigh;
    
  printf("100\n");
  
  //thread removal
  m_enable=false; //threads end.
  for(int thr=0; thr<m_nThreads; thr++) //are unlocked for the conclusion.
    m_spectro[thr].mutexIn_p->unlock();

  for(int i=0; i<m_nThreads; i++)
  {
    m_spectro[i].thread_p->join();
    delete m_spectro[i].thread_p; m_spectro[i].thread_p=0;
  }
  
  return 0;
}

//getRawInfo()
//Loads the full spectrum information associated with a pixel from an imzML file.
//The information is stored in the m_spectro structure, set to the range [m_mzLow, m_mzHigh].
//px: Pixel whose spectrum should be loaded.
//spIndex: Threads that manage it
//Returns the size of the spectrum.
int PeakMatrix::getRawInfo(int iPx, int spIndex)
{
  //raw info of mass
  int px=m_pxList[iPx];
  
  int massSize=m_getImzMLData_p->getPixelMassF(px, m_spectro[spIndex].tmpMass_p);//mass vector
  m_spectro[spIndex].size=massSize;
  if(massSize<=0) return 0;
  
  //raw info of intensities  
  int intSize=m_getImzMLData_p->getPixelIntensityF(px, m_spectro[spIndex].tmpInt_p);

  m_spectro[spIndex].pixel=iPx;
  return massSize;
} 

//mtGetGaussians()
//Parallel processing.
//Peak are delimited and their Gaussians are formed.
//This thread remains active, processing spectra until none remain.
//Each spectrum is converted into Gaussians that can overlap (join).
//spIndex: thread
//Returns -1 on failure, 0 = OK.
int PeakMatrix::mtGetGaussians(int spIndex)
  { 
  Common common;
  GAUSS_PARAMS *gaussians_p=0;
  double *centroids_p=0;
  int *centroidsIndex_p=0;
//  m_spectro[spIndex].mutexOut_p->unlock(); //end of spectrum processing.
//  return 0;
  
  //while execution is enabled.
  while(m_enable)
  {
    m_spectro[spIndex].mutexIn_p->lock(); //permission to continue (synchro)
    if(!m_enable) break; //The signal can be activated while waiting.

    IntensityPeak intPeak(m_SNR); //peak class
    
    //Spectrum conditioning.
    //Full spectra are received and only part of it may be of interest.
    
    int iMzLow=-1, iMzHigh=-1;
    int spSize, spSize_tmp;
    spSize=m_spectro[spIndex].size; //spectrum size
    spSize_tmp=spSize;
    //SNR
    m_spectro[spIndex].noise=m_noiseEst_p->getSNR(m_spectro[spIndex].tmpInt_p, m_spectro[spIndex].size,  m_spectro[spIndex].tmpSNR_p);

    //the spectrum is limited to the range of interest.
    iMzLow =common.nearestIndex(m_mzLow,  m_spectro[spIndex].tmpMass_p, spSize); //low index
    iMzHigh=common.nearestIndex(m_mzHigh, m_spectro[spIndex].tmpMass_p, spSize); //high index
    spSize=iMzHigh- iMzLow +1;

    //fault control
    if(iMzLow<0 || iMzHigh>spSize_tmp-1 || iMzHigh-iMzLow>spSize_tmp || iMzHigh<iMzLow) 
      {
        m_spectro[spIndex].mutexOut_p->unlock(); //end of spectrum processing.
        return -1;
      }
   //mass ordination.
    common.sortUpF(m_spectro[spIndex].tmpMass_p+iMzLow, m_spectro[spIndex].sort_p, spSize);

    //The part of interest is extracted from the mass, intensity and SNR vectors.
    for(int i=0; i<spSize; i++) 
      {
      int index=iMzLow+m_spectro[spIndex].sort_p[i];
        m_spectro[spIndex].mass_p[i]=m_spectro[spIndex].tmpMass_p[index];
        m_spectro[spIndex].int_p [i]=m_spectro[spIndex].tmpInt_p [index];
        m_spectro[spIndex].SNR_p [i]=m_spectro[spIndex].tmpSNR_p [index];
      }
    m_spectro[spIndex].size=spSize; //useful range
    //-----------------------------------------------------------------    
  
    //conversion to Gaussians.
    
    if(!m_enable) //end of thread?
      {m_spectro[spIndex].mutexOut_p->unlock(); return 0;}
    int nPeak;
    //the peak are extracted from the spectrum (they are delimited by their indices).
    nPeak=intPeak.getPeakList(&m_spectro[spIndex]);
    
    if(nPeak>0)
    {
    gMutex.lock();
    gPeakCount+=nPeak; //peak accumulation unsynchronized, therefore, it is not accurate.
    gSpectra++;    //processed spectra
    if(m_spectro[spIndex].tmpMass_p[0]<gMinMass) gMinMass=m_spectro[spIndex].tmpMass_p[0];
    if(m_spectro[spIndex].tmpMass_p[spSize_tmp-1]>gMaxMass) gMaxMass=m_spectro[spIndex].tmpMass_p[spSize_tmp-1];
    gMutex.unlock();
    }
   
    if(nPeak>0) //if there are peak to treat
      {
      int px=m_spectro[spIndex].pixel; //pixel
      //single peak info is saved.
      m_peakFG_p[px].peakF_p=new ION_INDEX[nPeak];
      for(int i=0; i<nPeak; i++) //copy peak
      {
        m_peakFG_p[px].peakF_p[i].low =intPeak.getSinglePeak(i).low;
        m_peakFG_p[px].peakF_p[i].max =intPeak.getSinglePeak(i).max;
        m_peakFG_p[px].peakF_p[i].high=intPeak.getSinglePeak(i).high;
      }

      m_peakFG_p[px].peakFsize=nPeak;//number of simple peak
      //the information of compound peak (simple joined-overlapping peak) is obtained.
      int nUPeak=intPeak.getCompoundPeakNumber(); 

      //Compound peak information is saved.
      //Each entry is a reference to the initial and final single peak of the merged peak.
      m_peakFG_p[px].peakU_p=new PEAK_UNITED[nUPeak];
      for(int i=0; i<nUPeak; i++)
      {
        m_peakFG_p[px].peakU_p[i].low =intPeak.getCompoundPeak(i).peakLow;
        m_peakFG_p[px].peakU_p[i].high=intPeak.getCompoundPeak(i).peakHigh;
      }
      m_peakFG_p[px].peakUsize=nUPeak;//number of compound peak
      
      //the peak are converted to Gaussians
      gaussians_p=new GAUSS_PARAMS[nPeak]; //temporal copy of Gaussians
      centroids_p=new double[nPeak];       //temporary mass copy
      centroidsIndex_p=new int[nPeak];     //indices to ordered masses
      
      int nGaussians=getGaussians(px, &m_spectro[spIndex], gaussians_p);
      if(nGaussians==-2) //too many gaussians in the spectra -> abort the iteration
      {
        if(gaussians_p)     {delete [] gaussians_p;       gaussians_p=0;}
        if(centroids_p)     {delete [] centroids_p;       centroids_p=0;}
        if(centroidsIndex_p){delete [] centroidsIndex_p;  centroidsIndex_p=0;}
        m_gaussians_p[px].gauss_p=0;
        m_gaussians_p[px].size=0;
        m_spectro[spIndex].mutexOut_p->unlock();
        continue;
      }

      m_gaussians_p[px].size=0;
 
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
            m_gaussians_p[px].gauss_p[i].mean  =gaussians_p[index].mean;
            m_gaussians_p[px].gauss_p[i].sigma =fabs(gaussians_p[index].sigma);
            m_gaussians_p[px].gauss_p[i].weight=gaussians_p[index].weight;
          }
          m_gaussians_p[px].size=nGaussians;
        }
       
      }
      else m_gaussians_p[px].size=0;

    } //end of peak processing

    if(gaussians_p)     {delete [] gaussians_p;       gaussians_p=0;}
    if(centroids_p)     {delete [] centroids_p;       centroids_p=0;}
    if(centroidsIndex_p){delete [] centroidsIndex_p;  centroidsIndex_p=0;}
    
    m_spectro[spIndex].mutexOut_p->unlock(); //end of spectrum processing.
  }

  return 0;
}


//getGaussians()
//Called from a thread.
//Sets the Gaussians on the peak.
//Uses the peak separation information (m_peakFG_p).
//px: pixel.
//spectro: pointer to the spectrum.
//gaussians_p: pointer to the structure containing the Gaussians' parameters.
//Returns the number of Gaussians or a value < 0 on failure.
int PeakMatrix::getGaussians(int px, SPECTRO *spectro_p, GAUSS_PARAMS *gaussians_p)
{
  float *intSpectrum_p=spectro_p->int_p, *massSpectrum_p=spectro_p->mass_p;
  int intSize=spectro_p->size;
  
  float minMeanPxMag=spectro_p->noise*m_SNR; //minimum value to consider a peak as valid.
  //class for conversion to Gaussians.
  GmmPeak gmmPeak(minMeanPxMag);
  
  int gaussIndex=0; //indices for each Gaussian.
  int nUPeak=m_peakFG_p[px].peakUsize; //#united peak
  
  //memory for predictable Gaussians.
  m_gaussians_p[px].gauss_p=new GAUSS_PARAMS[m_peakFG_p[px].peakFsize];
  
  //for each set of joined peak.
  for(int uPeak=0; uPeak<nUPeak; uPeak++)
  {
    int mPeakLow   =m_peakFG_p[px].peakU_p[uPeak].low;      //lower simple peak.
    int mPeakHigh  =m_peakFG_p[px].peakU_p[uPeak].high;     //upper simple peak.
    int lowMzIndex =m_peakFG_p[px].peakF_p[mPeakLow].low;   //lower  mass index.
    int highMzIndex=m_peakFG_p[px].peakF_p[mPeakHigh].high; //higher mass index.
    int mzSize=highMzIndex-lowMzIndex+1;
    
    int nPeak=mPeakHigh-mPeakLow+1; //number of simple peak into united peak.
    
    //The information of simple peak of magnitude that make up the segment is established.
    gmmPeak.setPeak(&m_peakFG_p[px].peakF_p[mPeakLow], nPeak); 

    //The magnitude information that must be adjusted is established.
    GROUP_F pxMag;
    pxMag.set=intSpectrum_p+lowMzIndex; 
    pxMag.size=mzSize;
    gmmPeak.setMagnitudes(&pxMag); //magnitude peak.
    
  //The simple Gaussians are formed whose sum reproduces the magnitude peak.
    int ret=gmmPeak.gmmDeconvolution(); //deconvolution
    if(ret==-1) //posible intensidad del espectro demasiado baja
      return -1;
    if(ret==-2)
    {
      gMutex.lock();
      printf("\nWarning: limits exceeded in pixel %4d. The aim is to deconvolve a mass segment composed of more than %d Gaussians.\n", spectro_p->pixel, DECONV_MAX_GAUSSIAN);
      gError++; //error count
      gMutex.unlock();
      return -2;
    }
    if(intSize<mzSize) //if the dimensions are not correct.
    {
      //warning: the dimensions of mzAxis and data[magnitudes] must match. 
      printf("Warning, match error detected in pixel %d. \n", spectro_p->pixel);
      return -1;
    }
    
    int nGauss=gmmPeak.getDeconvNumber(); //# gaussians
    
    bool hit=true;
    GAUSSIAN gaussIn, gaussOut; //if it is required to adapt the mass axis.
    for(int g=0; g<nGauss; g++)
    {
      gaussIn=gmmPeak.getDeconv(g); //get a gaussian.
      
      //unit conversion (scans to Daltos).
      gmmPeak.gaussConversion(&gaussIn, &gaussOut, massSpectrum_p+lowMzIndex, mzSize); 
      
      if(gaussIndex>=m_peakFG_p[px].peakFsize) //control de error
      {
        //warning: memory overflow for gaussians 
        printf("warning: memory overflow for gaussians %d/%d\n", gaussIndex, m_peakFG_p[px].peakFsize); 
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

//massRangeToCentroids()
//Parallel processing.
//1) Sets the number of independent mass ranges (joined peak).
//2) Establishes clusters within those ranges (kmeans segmentation).
//Requires preprocessing by getGaussians()
//Receives the total mass range to consider.
//Returns a list: peakMatrix, massVector, massResolution, pixelsSupport.
//peakMatrix: Matrix of centroids and the intensity associated with each pixel.
//massVector: The mz associated with each column of the peakMatrix.
//pixelsSupport: Number of pixels with intensity > minPixelsSupport.
List PeakMatrix::massRangeToCentroids(MASS_RANGE massRangeIn)
{
  Common common;
  int indexLow =common.nearestIndexMassRangeLow (massRangeIn.low,  m_massRange_p, m_massRangeSize);
  int indexHigh=common.nearestIndexMassRangeHigh(massRangeIn.high, m_massRange_p, m_massRangeSize);
  int Q =(indexHigh-indexLow+1)/m_nThreads; //partes iguales para cada thread
  int R =(indexHigh-indexLow+1)%m_nThreads; //partes iguales para cada thread
  int px;
  
  if(indexLow==-1 || indexHigh==-1) {return 0;}
  
  //Memory for necessary structures.
  //Holds mass range information.
  m_massSegment.massRange_p=  new MASS_RANGE[m_nThreads];
  m_massSegment.nGaussians_p= new int[m_nThreads];

  //Ion information with the intensities for each pixel.
  //Linked sequence.
  m_ionEntry_p=new ION_ENTRY*[m_nThreads];
  for(int i=0; i<m_nThreads; i++)
  {
    m_ionEntry_p[i]=new ION_ENTRY;
    m_ionEntry_p[i]->set=0;
    m_ionEntry_p[i]->size=0;
    m_ionEntry_p[i]->group=0;
  }
  
  m_enable=true; //allows the thread loop to operate.
  bool hit=false;
  int high=-1;
  //Mass ranges that each thread must process (linear distribution).
    for(int thr=0; thr<m_nThreads; thr++)
    {
      m_massSegment.massRange_p[thr].low =high+1; //Note: these are integer values over floats.
      if(thr<R)
      {
        m_massSegment.massRange_p[thr].high=high+Q+1;
      }
      else
      {
        m_massSegment.massRange_p[thr].high=high+Q;
      }
      high=m_massSegment.massRange_p[thr].high;
      if(high>=indexHigh) break;
    }
  //Thread
  //1) Mass ranges containing joined Gaussians are established.
  //2) Segmentation: Centroids and nearest pixels are established.
  //Threads are activated.
  m_enable=true;
  for(int i=0; i<m_nThreads; i++)
    m_spectro[i].thread_p=new std::thread(&PeakMatrix::mtSegmentation, this, i);
  
  printf("\tphase 3:\tsegments to centroids(%%): 00 ");

  //wait until the conclusion-
  for(int i=0; i<m_nThreads; i++)
  {
    m_spectro[i].thread_p->join();
    delete m_spectro[i].thread_p; m_spectro[i].thread_p=0;
  }
  printf("100\n");
    
  //adequacy of results for R.
  int totalIons=0;
  for(int i=0; i<m_nThreads; i++)
  {
    totalIons+=m_totalIons[i];
  }
  NumericMatrix peakMatrix(m_NPixels, totalIons);  //peak matrix.
  NumericVector massVector(totalIons);              //mass vector.
  NumericVector massResolution(totalIons);          //mass resolution.
  IntegerVector pixelsSupport(totalIons);           //#Support pixels.
 
  //Each thread contains the chained information of a consecutive part of the ions.
  ION_ENTRY *localIonEntry_p; //input
  int col=0;
  float lastMass=0;
  for(int thr=0; thr<m_nThreads; thr++)
   {
   localIonEntry_p=m_ionEntry_p[thr];  //entry to the info.
   while(localIonEntry_p->group) //as long as the next link is not zero.
   {
     //It may happen that the initial masses of a thread are lower than the last masses of 
     //the previous thread, and they should be discarded. 
     //This is a consequence of dealing with joined Gaussians.     
     if(localIonEntry_p->mass>lastMass)                 //the masses are ordered.
     {
       massVector(col)    =localIonEntry_p->mass;           //centroid
       massResolution(col)=localIonEntry_p->massResolution; //massResolution
       pixelsSupport(col) =localIonEntry_p->size;        //#pixels that support the ion
       for(int row=0; row<m_NPixels; row++)             //for each pixel
         peakMatrix(row, col)=localIonEntry_p->set[row]; //intensities
       col++;
     }
     if(localIonEntry_p->group->group==0) //the last link contains no info.
       lastMass=localIonEntry_p->mass;    //last mass of the thread
     localIonEntry_p=localIonEntry_p->group; //next item
   }
 }
  
  NumericMatrix coord(m_NPixels, 2);  //coordinates matrix.
  m_coordSize=0;
  for(int px=0; px<m_NPixels; px++)//for all pixels
  {
    if(m_gaussians_p[px].gauss_p==0) continue; //pixel without Gaussians
    else {
      coord(m_coordSize,0)=m_x[px];
      coord(m_coordSize,1)=m_y[px];
      m_coordSize++;
      }
  }
  
 List ret=List::create(Named("peakMatrix")=peakMatrix, Named("mass")=massVector, Named("massResolution")=massResolution, 
                       Named("pixelsSupport")=pixelsSupport, Named("coordinates")=coord);
 return ret;
 }

//getCentroidsIntoRange()
//Extracts the existing Gaussians within a mass range from the information in m_gaussians_p.
//massRange: Mass range from which to extract the Gaussians.
//gaussians_p: Requested Gaussians.
//size: size of reserved memory.
//Returns the number of Gaussians.
int PeakMatrix::getCentroidsIntoRange(MASS_RANGE massRange, float **gaussians_p, int size)
{
  int count=0;
  Common tools;
  bool hit=false;
  int iLow, iHigh, px;
  
  for(int px=0; px<m_NPixels; px++)//for all pixels
    {
    if(m_gaussians_p[px].gauss_p==0) continue; //pixel without Gaussians
    
    if(m_gaussians_p[px].size>33) //if long spectra
    {
      iLow =tools.nearestIndexGaussians(massRange.low,  m_gaussians_p[px].gauss_p, m_gaussians_p[px].size, 4|1); //ajuste a masa superior   
      iHigh=tools.nearestIndexGaussians(massRange.high, m_gaussians_p[px].gauss_p, m_gaussians_p[px].size, 4|2); //ajuste a masa inferior
      if(iLow==-1 || iHigh==-1 || iLow>iHigh) continue;

      for(int i=iLow; i<=iHigh; i++)
      {
        if(count>=size) //limits control
          {
          hit=true;
          break;
          }
        //copia
        gaussians_p[count][0]=m_gaussians_p[px].gauss_p[i].mean;
        gaussians_p[count][1]=fabs(m_gaussians_p[px].gauss_p[i].sigma);
        gaussians_p[count][2]=m_gaussians_p[px].gauss_p[i].weight;
        gaussians_p[count][3]=(float)px;
        count++;
      }
    if(hit) break;
    }
  else //short spectra
  {
    for(int i=0 ; i<m_gaussians_p[px].size; i++) 
      {
      
      if(m_gaussians_p[px].gauss_p[i].mean>= massRange.low && 
         m_gaussians_p[px].gauss_p[i].mean<= massRange.high)
        {
        if(count>=size) //limits control 
          {
          hit=true;
          break;
          }
        gaussians_p[count][0]=m_gaussians_p[px].gauss_p[i].mean;
        gaussians_p[count][1]=fabs(m_gaussians_p[px].gauss_p[i].sigma);
        gaussians_p[count][2]=m_gaussians_p[px].gauss_p[i].weight;
        gaussians_p[count][3]=(float)px;
        count++;
      }
    }
  }  
  }
  return count;
}

//getCentroidsNumberIntoRange()
//Returns the number of Gaussians in a mass range from the information in m_gaussians_p.
//massRange: Mass range from which to extract Gaussians.
//Returns the number of Gaussians.
int PeakMatrix::getCentroidsNumberIntoRange(MASS_RANGE massRange)
{
  int count=0, px;
  for(int px=0; px<m_NPixels; px++)//for all pixels
  {
    if(m_gaussians_p[px].gauss_p==0) continue; //pixel without Gaussians
    for(int i=0 ; i<m_gaussians_p[px].size; i++) //improve!!!!
    {
      if(m_gaussians_p[px].gauss_p[i].mean>= massRange.low && 
         m_gaussians_p[px].gauss_p[i].mean<= massRange.high)
        count++;
    }
  }
  return count;
}

//getMassRanges()
//establishes the mass segments where overlapping Gaussians exist.
//the information is stored in the array pointed to by m_massRange_p.
//returns the number of segments.
//procedure:
//A vector is generated containing information on the contributions of each Gaussian distribution across 
//all spectra. Overlapping portions of the Gaussian distributions are counted only once. Next, contiguous 
//and isolated peaks are extracted, and segments separated by noise levels below a certain threshold are 
//identified. The noise level is determined by the highest of three values: iterative averaging with saturation,
//the highest noise level, and the percentage of total pixels. If any segment exceeds 3 Da, the noise level is 
//increased until it is reached.
int PeakMatrix::getMassRanges(float* linkedPeak_p)
{
  printf("\tphase 2:\t#isolate mass segments: ");
  double highMass, lowMass, highMassOld, lowMassOld;
  int iLow, iHigh;
  double deltaMass=m_mzLow/(4*m_massResolution); //delta=1/4 of the minimum mass increment of the spectrometer
  int massAxisSize=1+(m_mzHigh-m_mzLow)/deltaMass;
  
  float *massAxis=0;
  massAxis= new float[massAxisSize];
  float linkedPeaks=m_linkedPeaks;
  int nUPeak=0, px;
  SPECTRO  spectro;
  spectro.SNR_p=0;
  spectro.SNR_p=new float[massAxisSize];
  bool hit;
  int iSegment=0;
  
    for(int i=0; i<massAxisSize; i++) massAxis[i]=0;
    
    double deltaMass2;
    
    //for all elements of the spectrum
    //vector con la cantidad de contribuciones de las gaussianas de cada espectro
    for(int px=0; px<m_NPixels; px++)
    {
      if(m_gaussians_p[px].gauss_p==0) continue; //if this entry does not contain info.
      lowMassOld=0; highMassOld=0;
      for(int i=0; i<m_gaussians_p[px].size; i++) 
      {
        //Very wide Gaussians (sigma>5*deltaMass) are discarded.
        deltaMass2=(m_gaussians_p[px].gauss_p[i].mean)/m_massResolution;
        if(fabs(m_gaussians_p[px].gauss_p[i].sigma)>5*deltaMass2) //sigma >>
        {
          continue;
        }
        //masses occupied by the Gaussian.
        //only add the overlapping part of the bell once.
        lowMass =m_gaussians_p[px].gauss_p[i].mean-linkedPeaks*m_gaussians_p[px].gauss_p[i].sigma;
        highMass=m_gaussians_p[px].gauss_p[i].mean+linkedPeaks*m_gaussians_p[px].gauss_p[i].sigma;
        if(lowMass<m_mzLow)   lowMass =m_mzLow;
        if(highMass>m_mzHigh) highMass=m_mzHigh;
        if(lowMass<highMassOld) //unidas
          lowMass=highMassOld;
        iLow =(lowMass -m_mzLow)/deltaMass;
        iHigh=(highMass-m_mzLow)/deltaMass;
        for(int j=iLow; j<=iHigh; j++) massAxis[j]+=1.0;
        if(highMass>highMassOld)
          highMassOld=highMass;
        }
      }
    
    //noise estimation
    NoiseEstimation noiseEst(3, 1, 9); //MAD, smoothig & windows=9
    float noise=noiseEst.getNoise(massAxis, massAxisSize);

    //The baseband is estimated (recursive averaging).
    //the iteration ends when the slope reaches a value less than 1% of the initial slope.
    //max iterations=10
    float acu=0, vMean;
    float value[2], minDiff=0;
    for(int j=0; j<massAxisSize; j++) acu+=massAxis[j]; vMean=acu/massAxisSize;
    for(int i=0; i<10; i++)
    {
      acu=0;
      for(int j=0; j<massAxisSize; j++)
        if(massAxis[j]<=vMean){acu+=massAxis[j];}
        else acu+=vMean;
      vMean=acu/massAxisSize;
      if(i==0) {value[0]=vMean;}
      else if(i==1) {value[1]=vMean; minDiff=abs(value[0]-value[1])*0.01; }
      else {if(abs(vMean-value[1])<minDiff) break; else value[1]=vMean;}
   }
    
    //decision of the value assigned to noise.
    vMean=noise>vMean?noise:vMean; // the greatest.
    float pxSupport=m_NPixels*0.0001; //minimum spectra that support.
    if(pxSupport<0.999) pxSupport=0.999; //at least one support (avoid unity).
    vMean=vMean>pxSupport?vMean: pxSupport; //the greatest.
    if(vMean<=0 && noise<1e-6) vMean=value[0]/100.0; //minimal noise
    if(vMean==0) vMean=1e-6;
    
    //A spectrum is created with massAxis info to obtain its peaks.
    spectro.int_p=massAxis;
    spectro.noise=vMean;
    spectro.size=massAxisSize;

    //the SNR is established.
    for(int i=0; i<massAxisSize; i++) 
    {
      if(spectro.int_p[i]<spectro.noise) 
      {spectro.int_p[i]=0; spectro.SNR_p[i]=0;} //Warning: The entire spectrum may be canceled.
      else
        {spectro.SNR_p[i]=spectro.int_p[i]/vMean;}
    }

    int maxSegments=1+(m_mzHigh-m_mzLow)/(m_mzLow/m_massResolution); //estimation
    maxSegments/=4;

    IntensityPeak intPeak(1); //simple and compound peaks are obtained with SNR=1.
    if(intPeak.getPeakList(&spectro)==-1) return 0; //no peak
    nUPeak=intPeak.getCompoundPeakNumber(); //#compound peaks.
    m_massRange_p=new MASS_RANGE[maxSegments];   //memory to store the joined mz ranges.

    //for each set of joined peak.
    int   pLow,  pHigh, iLowMass, iHighMass;
    float lowMz, highMz;
    hit=true;
    float localSNR=1;
    iSegment=0;
    for(int i=0; i<nUPeak; i++)
    {
      pLow =intPeak.getCompoundPeak(i).peakLow; //index to lower peak of the composite peak.
      pHigh=intPeak.getCompoundPeak(i).peakHigh;//index to upper peak of the composite peak.
      iLowMass =intPeak.getSinglePeak(pLow). low;
      iHighMass=intPeak.getSinglePeak(pHigh).high;
      lowMz =m_mzLow+iLowMass *deltaMass; //lower  mass.
      highMz=m_mzLow+iHighMass*deltaMass; //higher mass.
      
      //It is estimated that segments larger than 3 Da are not suitable for segmentation.
      if(highMz-lowMz <=3.0)
      {
        m_massRange_p[iSegment].low =lowMz;//package ends
        m_massRange_p[iSegment++].high=highMz;
      }
      else
      {
        localSNR=1;
        spectro.int_p=massAxis+iLowMass;
        spectro.size=iHighMass-iLowMass+1;
        spectro.noise=vMean;
        //the SNR is established.
        for(int i=0; i<=iHighMass-iLowMass+1; i++) 
        {
          if(spectro.int_p[i]<spectro.noise) 
            {spectro.SNR_p[i]=0;} //Warning: The entire spectrum may be canceled.
          else
            {spectro.SNR_p[i]=spectro.int_p[i]/vMean;}
        }
        
        int nUPeak2=0, iLowMass2, iHighMass2;
        //The base (noise) is raised one by one until the segment width < 3.0 Da.
        while(true)
        {
          localSNR+=1;
          IntensityPeak intPeak(localSNR); //simple and compound peaks are obtained with SNR=1.
          intPeak.getPeakList(&spectro); 
          nUPeak2=intPeak.getCompoundPeakNumber(); //#compound peaks.
          hit=true;
          for(int j=0; j<nUPeak2; j++)
          {
            pLow =intPeak.getCompoundPeak(j).peakLow; //index to lower peak of the composite peak.
            pHigh=intPeak.getCompoundPeak(j).peakHigh;//index to upper peak of the composite peak.
            iLowMass2 =intPeak.getSinglePeak(pLow). low +iLowMass;
            iHighMass2=intPeak.getSinglePeak(pHigh).high+iLowMass;
            lowMz =m_mzLow+iLowMass2 *deltaMass; //lower  mass.
            highMz=m_mzLow+iHighMass2*deltaMass; //higher mass.
          if(highMz-lowMz > 3.0) {hit=false;break;} //iterar para mejorar
          }
          if(hit) //solution found.
          {
            for(int j=0; j<nUPeak2; j++) 
            {
              pLow =intPeak.getCompoundPeak(j).peakLow; //index to lower peak of the composite peak.
              pHigh=intPeak.getCompoundPeak(j).peakHigh;//index to upper peak of the composite peak.
              iLowMass2 =intPeak.getSinglePeak(pLow). low +iLowMass;
              iHighMass2=intPeak.getSinglePeak(pHigh).high+iLowMass;
              lowMz =m_mzLow+iLowMass2 *deltaMass; //lower  mass.
              highMz=m_mzLow+iHighMass2*deltaMass; //higher mass.
              m_massRange_p[iSegment].low =lowMz;//package ends
              m_massRange_p[iSegment++].high=highMz;
            }
            break;
          }
        }
      }
    }
  *linkedPeak_p=linkedPeaks; //return info
  printf("%d\n", iSegment); //#mass segments
  
  if(massAxis)      delete [] massAxis;
  if(spectro.SNR_p) delete []spectro.SNR_p;
  return iSegment; 
 
}

  


//Obtains the average values of the Gaussians on an artificial mass axis.
//The mass axis is formed from the extreme masses to be considered and the desired resolution and overSampling.
//Returns a list with two arrays: averageMz and averageIntensity.
List PeakMatrix::getMeanGaussianSpectrum(float resolution, int overSampling)
  {
    double highMass, lowMass;
    int iLow, iHigh, px;
    Common tools;
    int minPx=0x7FFFFFFF, maxPx=0;
    for(int i=0; i<m_NPixels; i++)
    {
      if(m_pxList[i]<minPx) minPx=m_pxList[i];
      if(m_pxList[i]>maxPx) maxPx=m_pxList[i];
    }
    if(resolution>m_massResolution) resolution=m_massResolution;
    
    double deltaMass=(double)m_mzLow/(overSampling*resolution); //delta=1/4 of the minimum mass increment of the spectrometer
    int massAxisSize=1+(m_mzHigh-m_mzLow)/deltaMass;
    
    double *meanMassAxis_p=0;
    meanMassAxis_p= new double[massAxisSize];
    bool hit;
    int pxCount=0;
    for(int i=0; i<massAxisSize; i++) meanMassAxis_p[i]=0.0;
    
    double deltaMass2, gMean, gSigma, gIntensity, mz, den;
    double middle=deltaMass/2.0;
    
    //for all elements of the spectrum
    for(int px=0; px<m_NPixels; px++)
    {
      if(m_gaussians_p[px].gauss_p==0) continue; //if this entry does not contain info.
      for(int i=0; i<m_gaussians_p[px].size; i++) 
      {
        gSigma    =m_gaussians_p[px].gauss_p[i].sigma;
        gMean     =m_gaussians_p[px].gauss_p[i].mean;
        gIntensity=m_gaussians_p[px].gauss_p[i].weight;
        
        //masses occupied by the Gaussian.
        highMass=gMean+3*gSigma;
        lowMass =gMean-3*gSigma;
        if(lowMass <m_mzLow)   continue; 
        if(highMass>m_mzHigh)  continue; 
        
        iLow =(lowMass -m_mzLow)/deltaMass;
        iHigh=(highMass-m_mzLow)/deltaMass;
        den=2*gSigma*gSigma;
        for(int j=iLow; j<iHigh; j++) 
        {
          mz=j*deltaMass+m_mzLow+middle;
          meanMassAxis_p[j]+=gIntensity*exp(-(mz-gMean)*(mz-gMean)/den);
        }
      }
      pxCount++;
    }
    NumericVector meanMz(massAxisSize), meanInt(massAxisSize);
    for(int i=0; i<massAxisSize; i++)
    {
      meanMz[i]=i*deltaMass+m_mzLow;
      meanInt[i]=meanMassAxis_p[i]/pxCount;
    }
    if(meanMassAxis_p) delete [] meanMassAxis_p;
    
    List ret=List::create(Named("averageMz")=meanMz, Named("averageIntensity")=meanInt);
    return ret;
  }

//Obtains the average values of the intensities on an artificial mass axis.
//Noise is not taken into account.
//The mass axis is formed from the extreme masses to be considered and the desired resolution and overSampling.
//delta_mass=mzLow/(overSamplig*massResolution) where massResolution=mz/deltaMz
//Returns a list with two arrays: averageMz and averageIntensity.
List PeakMatrix::getMeanSpectrum(float resolution, int overSampling)
{
  double highMass, lowMass;
  int iLow, iHigh, px, massSize;
  Common tools;
  float *meanMassAxis_p=0;
  
  if(overSampling==1 && m_NPixels==1)
  {
    int massSize=m_getImzMLData_p->getPixelMassF(m_pxList[px], m_spectro[0].tmpMass_p);//mass vector
    if(massSize<=0) return 0;
    
    //raw info of intensities  
    m_getImzMLData_p->getPixelIntensityF(m_pxList[px], m_spectro[0].tmpInt_p);
    
    NumericVector meanMz(massSize), meanInt(massSize);
    for(int i=0; i<massSize; i++)
    {
      meanMz[i] =m_spectro[0].tmpMass_p[i];  //mass into bin
      meanInt[i]=m_spectro[0].tmpInt_p[i];//intensity into bin
    }
    
    List ret=List::create(Named("averageMz")=meanMz, Named("averageIntensity")=meanInt);
    return ret;
  }
  
  if(resolution>m_massResolution) resolution=m_massResolution;
  
  double deltaMass=(double)m_mzLow/(overSampling*resolution); 
  int massAxisSize=1+(m_mzHigh-m_mzLow)/deltaMass;
  
  meanMassAxis_p= new float[massAxisSize];

  int pxCount=0;
  for(int i=0; i<massAxisSize; i++) 
    {meanMassAxis_p[i]=0.0;}
  
  int idx0, idx1;
  double slope, newValue, I1, I0, M1, M0;
  double delta2, middle=deltaMass/2.0;
  
  //for all elements of the spectrum
  for(int px=0; px<m_NPixels; px++)
    {
    //raw info of masses 
    massSize=m_getImzMLData_p->getPixelMassF(m_pxList[px], m_spectro[0].tmpMass_p);//mass vector
    if(massSize<=0) continue;
   
    //raw info of intensities  
    m_getImzMLData_p->getPixelIntensityF(m_pxList[px], m_spectro[0].tmpInt_p);
    pxCount++;
    for(int i=1; i<massSize; i++) 
    {
      if(m_spectro[0].tmpMass_p[i-1]<m_mzLow) continue;
      idx0 =(m_spectro[0].tmpMass_p[i-1]-m_mzLow)/deltaMass; //low  mass point index
      idx1 =(m_spectro[0].tmpMass_p[i]  -m_mzLow)/deltaMass; //high mass point index
      if(idx0<0 || idx1>=massAxisSize) { break;}
      
      I0=m_spectro[0].tmpInt_p [i-1]; I1=m_spectro[0].tmpInt_p [i]; //intensity
      M0=m_spectro[0].tmpMass_p[i-1]; M1=m_spectro[0].tmpMass_p[i]; //mass
      if(M1<=M0) continue;
      slope=(I1-I0)/(M1-M0); //slope between two adjacent points
      for(int j=idx0; j<idx1; j++) //for all bins between points
      {
        delta2=j*deltaMass+middle+m_mzLow; //central mass
        if(delta2<M0) newValue=I0;
        else 
        {
          delta2-=M0; //incremento respecto a M0
          newValue=slope*delta2+I0; //new intensity
        }
        meanMassAxis_p[j]+=newValue; //accumulated into bin
      }
    }
  }
  NumericVector meanMz(massAxisSize), meanInt(massAxisSize);
  for(int i=0; i<massAxisSize; i++)
  {
    meanMz[i]=i*deltaMass+m_mzLow;  //mass into bin
    meanInt[i]=meanMassAxis_p[i]/pxCount;//intensity into bin
  }
  
  if(meanMassAxis_p) delete [] meanMassAxis_p;

  List ret=List::create(Named("averageMz")=meanMz, Named("averageIntensity")=meanInt);
  return ret;
}


//getMassRanges()
//establishes the mass segments where overlapping Gaussians exist.
//generates a Boolean mass axis with a resolution of 1/4 of the spectrometer's mass resolution.
//if any Gaussian invades its space, sets the corresponding +/-m_linkedPeaks*sigma checkboxes to true.
//Very wide Gaussians (sigma > 5*deltaMass) are discarded.
//the information is stored in the array pointed to by m_massRange_p.
//returns the number of segments.
int PeakMatrix::getMassRanges2()
{
  printf("\tphase 2:    #isolate mass segments: ");
  double highMass, lowMass;
  int iLow, iHigh, px;
  double deltaMass=m_mzLow/(4*m_massResolution); //delta=1/4 of the minimum mass increment of the spectrometer
  int massAxisSize=1+(m_mzHigh-m_mzLow)/deltaMass;
  
  bool *massAxis=0;
  massAxis= new bool[massAxisSize];
  for(int i=0; i<massAxisSize; i++) massAxis[i]=false;
  
  double deltaMass2;
  
  //for all elements of the spectrum
  for(int px=0; px<m_NPixels; px++)
  {
    if(m_gaussians_p[px].gauss_p==0) continue; //if this entry does not contain info.
    for(int i=0; i<m_gaussians_p[px].size; i++) 
    {
      //Very wide Gaussians (sigma>5*deltaMass) are discarded.
      deltaMass2=(m_gaussians_p[px].gauss_p[i].mean)/m_massResolution;
      if(fabs(m_gaussians_p[px].gauss_p[i].sigma)>5*deltaMass2) //sigma >>
      {
        continue;
      }
    //masses occupied by the Gaussian.
    highMass=m_gaussians_p[px].gauss_p[i].mean+m_linkedPeaks*m_gaussians_p[px].gauss_p[i].sigma;
    lowMass =m_gaussians_p[px].gauss_p[i].mean-m_linkedPeaks*m_gaussians_p[px].gauss_p[i].sigma;
    if(lowMass<m_mzLow) lowMass =m_mzLow;
    if(highMass>m_mzHigh) highMass=m_mzHigh;

    iLow =(lowMass -m_mzLow)/deltaMass;
    iHigh=(highMass-m_mzLow)/deltaMass;
    for(int j=iLow; j<=iHigh; j++) massAxis[j]=true;
    }
  }
  
  //separation into segments, observing unused spaces.
  bool into=false;
  int init=-1, end=-1;
  int count=0, index=0;
  for(int i=0; i<massAxisSize; i++)   //#packages (contiguous masses)
  {
    if(massAxis[i] && !into)          //home package
      {into=true; init=i;}            //inside the package
    if(!massAxis[i] && into)          //end of package
      {into=false; end=i-1; count++;} //out of the package
  }
  if(into && !count) //if all massAxis vector is true (only one segment)
    { 
    count=1;         
    m_massRange_p=new MASS_RANGE[1];
    m_massRange_p[index].low =m_mzLow;//package ends
    m_massRange_p[index].high=m_mzHigh;
    index=1;
    }
  else //there are more than one segment
  {
    m_massRange_p=new MASS_RANGE[count];
    
    into=false;
    for(int i=0; i<massAxisSize; i++)
    {
      if(massAxis[i] && !into) 
        {into=true; init=i;}
      if(!massAxis[i] && into) //stands out from the pack.
      {
        into=false; end=i-1; 
        m_massRange_p[index].low =m_mzLow+init*deltaMass;//package ends
        m_massRange_p[index].high=m_mzLow+end *deltaMass;
        if(m_massRange_p[index].high-m_massRange_p[index].low>1.0) 
          printf("%9.4f %9.4f\n", m_massRange_p[index].low, m_massRange_p[index].high);
        index++;
      }
    }
  }
  printf("%d\n", index); //#mass segments
  if(massAxis) delete [] massAxis;
  return index; 
}

//setGaussiansNumberIntoSegments()
//Determines the maximum number of Gaussians over the given mass intervals and all pixels.
//brute force.
//massRange: Mass range to consider.
//Returns the maximum value.
int PeakMatrix::setGaussiansNumberIntoSegments2(MASS_RANGE massRange)
{
  int maxGaussians=0, nGaussians, px;
  float massLow, massHigh;
  for(int mr=massRange.low; mr<=massRange.high; mr++) //for the entire mass range
  {
    massLow =m_massRange_p[mr].low;  //mass range of this segment
    massHigh=m_massRange_p[mr].high;
    nGaussians=0;

    for(int px=0; px<m_NPixels; px++)//for all pixels
    {
      if(m_gaussians_p[px].gauss_p==0) continue;  //no info
      for(int g=0; g<m_gaussians_p[px].size; g++)
      {
        if(m_gaussians_p[px].gauss_p[g].mean>=massLow && m_gaussians_p[px].gauss_p[g].mean<=massHigh)
          nGaussians++;
      }
      if(nGaussians>maxGaussians) maxGaussians=nGaussians;
    }
    m_massRange_p[mr].nGaussians=nGaussians; 
  }
  return maxGaussians;
}

//setGaussiansNumberIntoSegments()
//Determines the maximum number of Gaussians over the given mass intervals and all pixels.
//aproximaciones sucesivas
//massRange: Mass range to consider.
//Returns the maximum value.
int PeakMatrix::setGaussiansNumberIntoSegments(MASS_RANGE massRange)
{
  int count=0;
  Common tools;
  bool hit=false;
  int iLow, iHigh, px;
  int maxGaussians=0, nGaussians;
  
  float massLow, massHigh;
  for(int mr=massRange.low; mr<=massRange.high; mr++) //for the entire mass range
  {
    massLow =m_massRange_p[mr].low;  //mass range of this segment
    massHigh=m_massRange_p[mr].high;
    nGaussians=0;
    for(int px=0; px<m_NPixels; px++)//for all pixels
    {
      if(m_gaussians_p[px].gauss_p==0) continue; //pixel without Gaussians
      
      if(m_gaussians_p[px].size>33) //if long spectra
      {
        //index adjusted to the nearest; if out of range return -1.
        iLow =tools.nearestIndexGaussians(massLow,  m_gaussians_p[px].gauss_p, m_gaussians_p[px].size, 4);   
        iHigh=tools.nearestIndexGaussians(massHigh, m_gaussians_p[px].gauss_p, m_gaussians_p[px].size, 4);
        if(iLow==-1 && iHigh==-1) continue;
        if(iLow==-1 && iHigh!=-1) //iLow no content and iHigh content.
        {
          //index adjusted to the nearest above, even if it is out of range.
          iLow =tools.nearestIndexGaussians(massLow,  m_gaussians_p[px].gauss_p, m_gaussians_p[px].size, 1); //ajuste a masa superior   
         }
        if(iLow!=-1 && iHigh==-1) //iLow content and iHigh no content.
        {
          //index adjusted to the nearest below, even if it is out of range.
          iHigh =tools.nearestIndexGaussians(massHigh,  m_gaussians_p[px].gauss_p, m_gaussians_p[px].size, 2); //ajuste a masa superior   
          if(iHigh==-1) continue;
        }
        if(iLow>iHigh) continue;

        nGaussians+=iHigh-iLow+1; //accumulation of Gaussians in the segment.
      }
      else //short spectra
      {
        //brute force.
        for(int i=0 ; i<m_gaussians_p[px].size; i++) 
        {
          if(m_gaussians_p[px].gauss_p[i].mean>= massLow && m_gaussians_p[px].gauss_p[i].mean<= massHigh)
              nGaussians++;
        }
      } 
    }
    m_massRange_p[mr].nGaussians=nGaussians; 
    if(nGaussians>maxGaussians) maxGaussians=nGaussians;
  }
  return maxGaussians;
}

//mtSegmentation()
//Kmeans segmentation in parallel processing.
//For each mass segment, centroids are generated.
//The information is stored in the array pointed to by m_ionEntry_p.
void PeakMatrix::mtSegmentation(int thrIndex)
{
  float mass=0;
  float **gaussians_p=0; //pointers to info to return.
  int totalIons=0;
  Common common;
  MASS_RANGE massRange, massRange2;
  massRange.low=m_mzLow;
  massRange.high=m_mzHigh;
  int count=0;
  ION_ENTRY *localIonEntry_p=m_ionEntry_p[thrIndex];
  
  int mRange=0; 
  float maxRange=0;
  int maxRangeIndex=-1;
  int iLow =m_massSegment.massRange_p[thrIndex].low;
  int iHigh=m_massSegment.massRange_p[thrIndex].high;
  if(iHigh<iLow) return;
  
//Determines the maximum cluster number for the mass range to be evaluated.
//Determines the largest segment in the range.
//A maximum value for clusters is understood to be one in which all Gaussians are
//as close as the mass resolution allows.  
for(int i=iLow; i<=iHigh; i++) 
  {
    float range=m_massRange_p[i].high- m_massRange_p[i].low;
    if(range>maxRange) {maxRange=range; maxRangeIndex=i;} //greater range to evaluate
  }
  //minimum resolution value in this mass range (Da).
  double massRes=m_massRange_p[iLow].low/m_massResolution;//minimum resolution for the thread
  int maxClustersThr=ceil(maxRange/massRes); //maximum clusters

  int *tmpIndex_p=0, *sortedIndex_p=0;;
  float *prob_p=0, *mass_p=0, *mass2_p=0;
  double *tmpMass_p=0;
  
  tmpMass_p=new double[maxClustersThr];
  tmpIndex_p=new int[maxClustersThr];
  prob_p=new float[maxClustersThr];
  mass_p=new float[maxClustersThr];
  sortedIndex_p=new int [maxClustersThr];
  
  //maximum Gaussians in the given mass range.
  int maxGNumber, maxGNumber2;
  maxGNumber=setGaussiansNumberIntoSegments(m_massSegment.massRange_p[thrIndex]);
  gaussians_p=new float*[maxGNumber];
  for(int i=0; i<maxGNumber; i++)
  {
    gaussians_p[i]=0;
    gaussians_p[i]=new float[4];
  }
  float *massVect_p=0;
  massVect_p=new float[maxGNumber];
    
  double newMassRes;

  //main loop
  //Each isolated mass packet is segmented.
  /////////////////////////////////////////////////////////////////  
  for(int mrIndex=iLow; mrIndex<=iHigh; mrIndex++)
  {
    if(!m_enable) break;
    gMutex.lock();
    gProcSegments++;
    if(thrIndex==0) 
      if((float)gProcSegments/(float)m_massRangeSize>gVez*0.1) {if(gVez<10) printf("%d ", gVez*10); gVez++;}
    gMutex.unlock();

    massRange2.low =m_massRange_p[mrIndex].low; //mass segment to be evaluated
    massRange2.high=m_massRange_p[mrIndex].high;
    m_massRange_p[mrIndex].resolution=0;
    int rangeNCenters;
    //Capture of the centroids in the range of masses to be considered.
    rangeNCenters=getCentroidsIntoRange(massRange2, gaussians_p, maxGNumber);
    if(rangeNCenters<=m_pxSupport) continue; //minimum size control.

    //mass vector for the current range.
    for(int i=0; i<rangeNCenters; i++)
      massVect_p[i]=gaussians_p[i][0]; //centroids

    bool hit, hit2;
    int maxIter=40;
    float convergencia=1e-5;
    
    double massRes=massRange2.low/m_massResolution; //resolution (Da)
    newMassRes=massRes;
    double massResSqr=massRes*massRes; //Da*Da
    int nClusters, iter=0;
    hit=true;
    int maxMasterIter=10; //maximum coarse iterations.
    
    //It is assumed that the maximum number of clusters coincides with the case where there are 
    //Gaussians in all possible spaces.
    int maxClustersRange=ceil((massRange2.high-massRange2.low)/massRes);
    if(maxClustersRange<0) {continue;}
    else if(maxClustersRange==0) {maxClustersRange=1;} //single mass range case.

    float*centers_p=0;
    int *centersSize_p=0, nCenters;
    centers_p=new float[maxClustersRange+1];
    centersSize_p=new int[maxClustersRange+1];
//    centers_p=new float[rangeNCenters];
//    centersSize_p=new int[rangeNCenters];
    
   
    nCenters=centers(massVect_p, rangeNCenters, massRes, centers_p, centersSize_p);
//for(int i=0; i<nCenters; i++)
//  fprintf(fp, "[%3d]%9.4f %d\n",i, centers_p[i], centersSize_p[i]);
    if(nCenters<=0) {continue;}

      //indexes ordered from highest to lowest probability.
    common.sortDownI(centersSize_p, sortedIndex_p, nCenters);
  
    for(int i=0; i<nCenters; i++) //initialization masses for clusters.
      mass_p[i]=centers_p[sortedIndex_p[i]]; 
    
    //class for segmentation
    KmeansR kmeans(massVect_p, rangeNCenters, maxIter, convergencia);
    if(kmeans.m_error==1) {printf("###2\n"); continue;}
    mass2_p=mass_p;

    while(true) //Iterate until goals are reached, relaxing constraints if necessary.
    {
      //Iterates until all clusters reach the desired resolution or are small groups
      //with each iteration, the number of clusters increases by 1.
    for(nClusters=1; nClusters<=nCenters; nClusters++) //clusters search
      {
        if(nClusters>nCenters) mass_p=0;// aleatoriedad
        if(kmeans.getClusters(nClusters, mass_p)!=0) //clusters are obtained.
          {hit=false; printf("###3\n"); break;} //if algorithm fails
        hit=true;
        //The validity of the obtained clusters is evaluated:
        //All large clusters must have low dispersion.
        //Low dispersion if it is < mass resolution.
        for(int c=0; c<nClusters; c++) 
        {
          int cPixels=kmeans.m_kStruct.clusters_p[c].data.size; //px that support
          //dispersion measure.
          double res=kmeans.m_kStruct.clusters_p[c].withinss/cPixels; //mean square distance.
          if(res>massResSqr)
            {hit=false; break;}         //requires iteration to improve.
          }
        if(hit) {break;}                //all groups of size>minimum have low dispersion.
      } //end for() for nClusters
    
    //If the desired conditions are not met, the constraints are relaxed and iterates completely.
    if(!hit) 
    {
        iter++;
      if(iter>=maxMasterIter) //limit control. Ensures that it ends.
        {
          printf("Warning: clustering failed in mass range %.5f/%.5f with %d clusters (maxClusters=%d)\n", 
                 m_massRange_p[mrIndex].low, m_massRange_p[mrIndex].high, kmeans.m_kStruct.nClusters, maxClustersRange);
          break;
        }
        //The resolution is relaxed.
        if(iter<=4) newMassRes=massRes*(1.0+iter/2.0);  //massRes/2 increments
        else newMassRes+=massRes*(iter-4);              //exponential increases.
        massResSqr=newMassRes*newMassRes;
        continue;
    }
    else break;
  }//while end
    
    //The resolution with which the segmentation has been resolved is noted.
    m_massRange_p[mrIndex].resolution=newMassRes*1e6/m_massRange_p[mrIndex].low;
    
    if(iter>=maxMasterIter) {kmeans.freeClusters();continue;}//voided segment.
    
    //increasing mass ordering
    for(int c=0; c<kmeans.m_kStruct.nClusters; c++)
    {
      tmpMass_p[c]=kmeans.m_kStruct.clusters_p[c].center;
    }
    common.sortUp(tmpMass_p, tmpIndex_p, kmeans.m_kStruct.nClusters);

    //The results are saved.
    //Memory is created to store the cluster information about the peak array.
    //These are structures linked by pointers so they can grow.
    for(int c=0; c<kmeans.m_kStruct.nClusters; c++)
    {
      int sortIndex=tmpIndex_p[c]; //index to cluster
      
      //filter: discard poorly supported clusters -> "minPixelsSupport" parameter
      int size=kmeans.m_kStruct.clusters_p[sortIndex].data.size; //cluster size
      if(size<m_pxSupport) continue;
      
      localIonEntry_p->set=new float[m_NPixels]; //memory for intensities
      for(int i=0; i<m_NPixels; i++) //reset all intensities. 
        localIonEntry_p->set[i]=0.0;
      
      localIonEntry_p->mass=kmeans.m_kStruct.clusters_p[sortIndex].center; //centroide
      localIonEntry_p->size=size; //intensities
      localIonEntry_p->massResolution=massRange2.low/newMassRes; //resolution (Da)

      totalIons++; //input counter (centroides)
      
      for(int i=0; i<size; i++) //copy of intensities
      {
        int index=kmeans.m_kStruct.clusters_p[sortIndex].data.set[i]; //Gaussian index
        int px=round(gaussians_p[index][3]); //pixel of the Gaussian.
        localIonEntry_p->set[px]=gaussians_p[index][2]; //intensities
      }
      
      localIonEntry_p->group=new ION_ENTRY; //new entry for the ion.
      localIonEntry_p=localIonEntry_p->group; //the pointer is updated
      localIonEntry_p->group=0; //initializes the next element.
      localIonEntry_p->set=0;
      localIonEntry_p->size=0;
    }

    if(centers_p)     delete []centers_p;
    if(centersSize_p) delete []centersSize_p;
  
  } //end of analysis of a mass segment.
  
  m_totalIons[thrIndex]=totalIons;
  
  //reserved memory is freed.
  if(tmpMass_p) delete [] tmpMass_p;
  if(tmpIndex_p) delete [] tmpIndex_p;
  if(prob_p) delete [] prob_p;
  if(mass2_p) delete [] mass2_p;
  if(sortedIndex_p) delete []sortedIndex_p;
  
  if(massVect_p) delete [] massVect_p;
  if(gaussians_p)
  {
    for(int i=0; i<maxGNumber; i++)
      if(gaussians_p[i]) delete [] gaussians_p[i];
    delete []gaussians_p;
  }
}

// Determines the centers of groups of values whose distance does not exceed segmentSize.
// mass_p: array of data to consider.
// size: size of the mass_p array.
// Centers_p: array of final centers.
// centerSize_p: number of elements that make up each center.
// Returns the number of centers detected (length of the centers_p array).
int PeakMatrix::centers(float *mass_p, int size, float segmentSize, float* centers_p, int *centerSize_p)
{
  if(size<=0) return 0;
  
  int *sortedIndex_p=0;
  sortedIndex_p=new int[size];
  Common common;
  common.sortUpF(mass_p, sortedIndex_p, size); //increasing ordering of masses.
  
  centers_p[0]=mass_p[sortedIndex_p[0]];
  centerSize_p[0]=1;
  
  int k=0;
  for(int i=1; i<size; i++) //iterative averaging.
  {
    if(fabs(mass_p[sortedIndex_p[i]]-centers_p[k])<segmentSize)
    {
      centerSize_p[k]++;
      centers_p[k]+=(mass_p[sortedIndex_p[i]]-centers_p[k])/centerSize_p[k];
    }
    else //new bin
    {
      centers_p[++k]=mass_p[sortedIndex_p[i]];
      centerSize_p[k]=1;
    }
  }
  if(sortedIndex_p) delete [] sortedIndex_p;
  return k+1;
}
