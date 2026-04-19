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
#include "rawToGaussians.h"

 extern int gProcSegments, gVez; //observers.
 extern std::mutex gMutex;
 int  gPeakCount, gSpectra, gError;
 float gMinMass, gMaxMass;

/// R METHOD ////////////////////////////////////////////////////////////////////////

//'
 //'  @name rawToGaussiansR
 //'  @title converts the peaks of each spectrum in the imzML file into Gaussian waves.
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
 //'  @return number of pixels
 //'     
 // [[Rcpp::export]]
 int rawToGaussiansR(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, Rcpp::NumericVector pxList, int nThreads)
 {
   gPeakCount=0, gSpectra=0;

   RawToGaussians pMatrix(ibdFname, imzML, params, pxList, mzLow, mzHigh, nThreads);
   if(pMatrix.m_hit==false) return 0;
   
   //Phase 1:
   //loads data from a file and converts its peak into Gaussians.
   int ret1=pMatrix.rawToGaussians(); //parallel processing
   if(ret1<0) 
      {return 0;}
   if(gPeakCount<=0) 
      {printf("No peaks are detected in the sample.\n"); return 0;}
   bool hit=pMatrix.saveGaussians((char*)"tmpGaussians.bin", pMatrix.m_gaussians_p);
   return pMatrix.m_NPixels;
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
   RawToGaussians gauss(ibdFname, imzML, params, pxList, mzLow, mzHigh, nThreads);  
   if(gauss.m_hit==false) return 0;
   
   //loads data from a file and converts its peak into Gaussians.
   int ret1=gauss.rawToGaussians(); //parallel processing
   if(ret1<0) 
   {return 0;}
   if(gPeakCount<=0) 
   {printf("No peaks are detected in the sample.\n"); return 0;}
   
   List ret=gauss.getMeanGaussianSpectrum(gauss.m_massResolution,(int)overSampling);
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
   RawToGaussians gauss(ibdFname, imzML, params, pxList, mzLow, mzHigh, 1);  
   if(gauss.m_hit==false) return 0;
   
   List ret=gauss.getMeanSpectrum(gauss.m_massResolution, (int)overSampling);
   return ret;
 }
 
 
 //Constructor
 //captures input information, allocates memory and initializes.
 ////////////////////////////////////////////////////////////////////////////////
 RawToGaussians::RawToGaussians(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, Rcpp::NumericVector pxList, float mzLow, float mzHigh, int nThreads)
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
   m_maxPxGaussians=0;

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
//   m_NPixels=20;
   
   printf("pixel Minimun:%d; pixel Maximun:%d; total:%d\n", m_pxMin+1, m_pxMax+1, m_NPixels);
   
   printf("to process: #spectra=%d in mass range(Da)=%.4f to %.4f\n",m_NPixels, m_mzLow, m_mzHigh);
   
   //copy pixel coordinates
   NumericVector X=df["x"];
   NumericVector Y=df["y"];
   
   savePixelsCoordinates((char*)"tmpPixelsCoordinates.bin", X, Y);
   
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
//     m_gaussians_p[px].initGauss=0;
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
     m_spectro[i].thread_p=new std::thread(&RawToGaussians::mtGetGaussians, this, i);
   }
   gMutex.unlock();
   gProcSegments=0;
   gVez=1;
   gError=0;
   
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
 RawToGaussians::~RawToGaussians()
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
   

   if(m_noiseEst_p)  delete m_noiseEst_p;
   //printf("end PeakMatrix destructor\n");
 }
 
 //freeing buffer.
 void RawToGaussians::freeMemoryPeak()
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
   

 }
 
 //rawToGaussians
 //gets the intensity peak and converts them into Gaussians.
 //The intensity and mass data adjusted to the range of interest are loaded from the imzML file.
 //the SNR info is established for each point of the spectra.
 //The generated information is stored in the m_spectro_p structure.
 /////////////////////////////////////////////////////////////////////////////////////////////
 int RawToGaussians::rawToGaussians()
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
 int RawToGaussians::getRawInfo(int iPx, int spIndex)
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
 int RawToGaussians::mtGetGaussians(int spIndex)
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
//    nPeak=0; 
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
 int RawToGaussians::getGaussians(int px, SPECTRO *spectro_p, GAUSS_PARAMS *gaussians_p)
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
 

 
 //Obtains the average values of the Gaussians on an artificial mass axis.
 //The mass axis is formed from the extreme masses to be considered and the desired resolution and overSampling.
 //Returns a list with two arrays: averageMz and averageIntensity.
 List RawToGaussians::getMeanGaussianSpectrum(float resolution, int overSampling)
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
 List RawToGaussians::getMeanSpectrum(float resolution, int overSampling)
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
 
 // Saves the Gaussian data to the given file
 // Adds it to any existing data
 // Returns false if failed 
 bool RawToGaussians::saveGaussians(char *fileName, GAUSS_SP *gauss_p)
 {
   std::fstream fp;
   bool hit=true;
   Common tools;
   char txt[200];
   int size=0;
   
   fp.open(fileName, std::fstream::out | std::ios::binary | std::ios::app);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: %s file could not be opened\n.", fileName);
     throw std::runtime_error(txt);
   }

   int gaussSize=3*sizeof(double);
   fp.write((char*)&m_NPixels, sizeof(int));

   for(int px=0; px<m_NPixels; px++)
     {
     if(gauss_p[px].gauss_p)
      {
       fp.write((char*)&gauss_p[px].size, sizeof(int));
       fp.write((char*) gauss_p[px].gauss_p, gauss_p[px].size*gaussSize );
       if(fp.fail() || fp.bad()) //control de errores
        {
         char txt[200];
         sprintf(txt, "Error: %s file could not be written completely\n.", fileName);
         throw std::runtime_error(txt);
         hit =false; break;
        }
       }
     else
     {
       fp.write((char*)&size, sizeof(int));
     }
   }
   fp.close();
   return hit;  
 }
 
 //Save the coordinate (XY) information of each pixel of m_pxList[] to the fileName file.
 //It is added to any existing data.
 //returns false if it failed.
 bool RawToGaussians::savePixelsCoordinates(char *fileName, NumericVector X, NumericVector Y)
 {
   std::fstream fp;
   bool hit=true;
   int coordX, coordY, px;
   
   fp.open(fileName, std::fstream::out | std::ios::binary | std::ios::app);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: %s file could not be opened\n.", fileName);
     throw std::runtime_error(txt);
   }
   
   fp.write((char*)&m_NPixels, sizeof(int)); //# de pixeles
   
   for(int i=0; i<m_NPixels; i++)
   {
     px=m_pxList[i];
     coordX=X[px];
     coordY=Y[px];
     fp.write((char*)&coordX, sizeof(int));
     fp.write((char*)&coordY, sizeof(int));
     if(fp.fail() || fp.bad()) //control de errores
       {
         char txt[200];
         sprintf(txt, "Error: %s file could not be written completely\n.", fileName);
         throw std::runtime_error(txt);
         hit =false; break;
       }
   }
   fp.close();
   return hit;  
 }
 