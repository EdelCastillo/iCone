#include "toleranceEstimation.h"

extern int  gProcSegments, gVez; //observers.
extern std::mutex gMutex, gMutex2;
extern int  gPeakCount, gSpectra, gError;
extern double gMinMass, gMaxMass;
#define MAX_PX  1000

//Constructor
//captures input information, allocates memory and initializes.
////////////////////////////////////////////////////////////////////////////////
ToleranceEstimation::ToleranceEstimation(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, int nThreads)
{
  m_hit=true;
  gMinMass=1e30;
  gMaxMass=0;
  m_getImzMLData_p=0;
  m_noiseEst_p=0;
  m_maxMzLength=0;
  m_pxLen_p=0;
  m_pxLenIdx_p=0;
  m_peaksFG_p=0;
  
  //information capture
  Rcpp::DataFrame df;
  df=imzML["run"];
  m_continuous=imzML["continuous_mode"];
  NumericVector nv;
  CharacterVector cv;
  m_pxMax=df.nrows()-1;
  m_pxMin=0;
  try
    {  
    m_NPixels=df.nrows();
    m_pxLen_p   =new int[m_NPixels];
    m_pxLenIdx_p=new int[m_NPixels];
    //dimensioned
    NumericVector mzLength=df["mzLength"];
    for(int i=0; i<m_NPixels; i++)
    {
      m_pxLen_p[i]=mzLength[i];
      if(mzLength[i]>m_maxMzLength)m_maxMzLength=mzLength[i]; //maximum spectrum
    }
    Common common;
    common.sortDownI(m_pxLen_p, m_pxLenIdx_p, m_NPixels);
    
    if(m_NPixels>MAX_PX) m_NPixels=MAX_PX;
    m_peaksFG_p=new PEAKS_FG[m_NPixels];
      for(int i=0; i<m_NPixels; i++)
      {
        m_peaksFG_p[i].size=0;
        m_peaksFG_p[i].peaks_p=0;
      }
         //m_NPixels=20;
      
      printf("to process: #spectra=%d\n",m_NPixels);
      
      nv=params["SNR"];
      m_SNR=nv[0];
      if(m_SNR<=0) m_SNR=1;
      
      cv=params["noiseMethod"];
      String tmpStr=cv[0];
      const char* SNRmethod=tmpStr.get_cstring(); //conversion to C
      
      if     (strcmp(SNRmethod, "estnoise_diff")==0) {m_SNRmethod=1; }
      else if(strcmp(SNRmethod, "estnoise_sd")  ==0) {m_SNRmethod=2; }
      else if(strcmp(SNRmethod, "estnoise_mad") ==0) {m_SNRmethod=3; }
      else {m_SNRmethod=0; printf("unknow noise method: %s so, estnoise_mad is used\n", SNRmethod);}
      m_noiseEst_p=new NoiseEstimation(m_SNRmethod, 1, 9);
      
      //class for accessing imzML files.
      m_getImzMLData_p= new GetImzMLData(ibdFname, imzML);
      
      //memory and its initialization
      
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
        m_spectro[i].mass_p  =new double[m_maxMzLength];
        m_spectro[i].int_p   =new double[m_maxMzLength];
        m_spectro[i].SNR_p   =new double[m_maxMzLength];
        m_spectro[i].sort_p     =new int  [m_maxMzLength];
        m_spectro[i].mutexIn_p  =new std::mutex;
        m_spectro[i].mutexOut_p =new std::mutex;
        m_spectro[i].size=0;
        m_spectro[i].mutexIn_p ->lock();
        m_spectro[i].mutexOut_p->unlock();
        m_spectro[i].thread_p=new std::thread(&ToleranceEstimation::mtToleranceEstimation, this, i);
        
      }
      gMutex.unlock();
      gMutex2.unlock();
      gProcSegments=0;
      gVez=1;
      gError=0;
      
    }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return;
  }
}

ToleranceEstimation::~ToleranceEstimation()
{
  if(m_pxLenIdx_p)  delete [] m_pxLenIdx_p;
  if(m_pxLen_p)     delete [] m_pxLen_p;
  if(m_peaksFG_p)   delete [] m_peaksFG_p;
  if(m_noiseEst_p)  delete m_noiseEst_p;
  if(m_getImzMLData_p) delete m_getImzMLData_p;
  
  for(int i=0; i<m_nThreads; i++)
  {
    if(m_spectro[i].int_p)      delete []m_spectro[i].int_p; 
    if(m_spectro[i].mass_p)     delete []m_spectro[i].mass_p;
    if(m_spectro[i].SNR_p)      delete []m_spectro[i].SNR_p; 
    if(m_spectro[i].sort_p)     delete []m_spectro[i].sort_p; 
  }
}

//parallel processing
//Estimates tolerance by counting the mass points that make up each peak
//Analyzes the 1,000 spectra with the highest number of mass points
//Returns a matrix containing the estimate for 10 consecutive mass ranges.
NumericMatrix ToleranceEstimation::getTolerance()
{
  printf("Processing \n\tphase 1:   from raw data to peaks (%%): 00 ");
  int vez=1;
  Common common;
  int spSize=0;
  bool hit=false;
  NumericMatrix SNR_Mx;
  NumericVector noiseSize;
  int px;

  //For each spectrum, determine the intensity peak, the joined peak, and their Gaussians.
  for(int iPx=0; iPx<m_NPixels; ) 
  {
    //indication of the progress of the process (10% tolerance)
    if((double)iPx/(double)m_NPixels>vez*0.1) {if(vez<10) printf("%d ", vez*10); vez++;}
    if(iPx==m_NPixels) break;
    
    //We load as many spectra as threads are used.
    //Each spectrum is processed by a thread.
    for(int thr=0; thr<m_nThreads && iPx<m_NPixels; thr++) //for each thread
    {
      if(m_spectro[thr].mutexOut_p->try_lock()) //if this thread is free
      {
        while(iPx<m_NPixels) //iterates while the spectrum has length <=2 (eliminates empty spectra).
        {
          spSize=getRawInfo(m_pxLenIdx_p[iPx], thr);//capturing spectra from imzML file
          if(spSize>2) break;//spectrum for analysis. There is no information of interest in the pixel spectrum <=3 peak.
          else
            iPx++;
        }
        if(spSize>2 && iPx<m_NPixels) 
          m_spectro[thr].mutexIn_p->unlock(); //This thread is allowed to run.
        iPx++;
      }
    }
    //wait the conclusion
    while(true)  
    {
      hit=false;
      for(int thr=0; thr<m_nThreads; thr++) //for each thread
      {
        if(m_spectro[thr].mutexOut_p->try_lock()) //if it was free
          m_spectro[thr].mutexOut_p->unlock(); 
        else 
          hit=true; 
      }
      if(!hit) break;
    }
  }
  printf("100\n");
  
  //mass range
  int minIdx, maxIdx;
  double minMz=0x7FFFFFFF, maxMz=0;
  for(int px=0; px<m_NPixels; px++)
  {
    int lastIdx=m_peaksFG_p[px].size-1;
    if(m_peaksFG_p[px].size==0) continue;
    if(m_peaksFG_p[px].peaks_p[0].lowMz<minMz) 
      minMz=m_peaksFG_p[px].peaks_p[0].lowMz;
    if(m_peaksFG_p[px].peaks_p[lastIdx].highMz>maxMz) 
      maxMz=m_peaksFG_p[px].peaks_p[lastIdx].highMz;
  }
  printf("\t\t\tmass range: %.4f %.4f\n", minMz, maxMz);
  printf("\tphase 2:   from peaks to tolerance(%%): 00 ");
  
  //oversampled mass axis * 20
  //pixels where peaks exist are accumulated along the axis
  double deltaMz=(maxMz-minMz)/(20.0*m_maxMzLength);
  int idxSize=(maxMz-minMz)/deltaMz;
  double *idxAxis_p=0;
  idxAxis_p=new double[idxSize];

  for(int i=0; i<idxSize; i++) idxAxis_p[i]=0;
  int idxLow, idxHigh;

double acuMz=0, acuM=0, mzCenter;
for(int px=0; px<m_NPixels; px++) //for all pixels
    {
    for(int i=0; i<m_peaksFG_p[px].size; i++) //for all peaks
      {
      if(m_peaksFG_p[px].size==0) continue;
      idxLow =(m_peaksFG_p[px].peaks_p[i].lowMz -minMz)/deltaMz;
      idxHigh=(m_peaksFG_p[px].peaks_p[i].highMz-minMz)/deltaMz;
      if(idxLow<0 || idxHigh>=idxSize)  continue;
      for(int k=idxLow; k<=idxHigh; k++) //pixel accumulation
        idxAxis_p[k]+=1.0;
      }
    }
  //the peaks are obtained on the oversampled mass axis
  SPECTRO spectro;
  spectro.int_p=idxAxis_p;
  spectro.size=idxSize;
  try{
    spectro.SNR_p=new double[idxSize];
    }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }
  
  //NoiseEstimation noiseEst_p(1, 1, 9); //sd
  //spectro.noise=noiseEst_p.getSNR(spectro.int_p, spectro.size, spectro.SNR_p);
  
  spectro.noise=0.1*m_NPixels;
  for(int i=0; i<idxSize; i++)
    spectro.SNR_p[i]=spectro.int_p[i]/spectro.noise;
  
  //IntensityPeak intPeak(m_SNR, false); //peak class
  IntensityPeak intPeak(1, false); //SNR=1, united peaks=false
  int nPeak=intPeak.getPeakList(&spectro, 0.15);  //intensity to zero if intensity<SNR*0.15

  //tolerance estimate for each peak
  double *tmpTolerance_p=0;
  double *tmpMass_p=0;
  double *pointsAcu_p=0; 
  int *pointsIdx_p=0; 
  double massAcu=0;
  try{  
    tmpTolerance_p=new double[nPeak];
    tmpMass_p     =new double[nPeak];
    pointsAcu_p   =new double[nPeak];
    pointsIdx_p   =new int[nPeak];
   }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }

  int pkCount=0;
  int lowIdx, highIdx;
  double lowMz, highMz;
  for(int i=0; i<nPeak; i++) {pointsAcu_p[i]=0;}
  
  for(int pk=0; pk<nPeak; pk++)
  {
    lowIdx =intPeak.getSinglePeak(pk).low;
    maxIdx =intPeak.getSinglePeak(pk).max;
    highIdx=intPeak.getSinglePeak(pk).high;
    lowMz =lowIdx *deltaMz+minMz;
    maxMz =maxIdx *deltaMz+minMz;
    highMz=highIdx*deltaMz+minMz;
    pkCount=0;
    massAcu=0;
    for(int px=0; px<m_NPixels; px++)
    {
      PEAKS_F *peaks_p=m_peaksFG_p[px].peaks_p;
      if(m_peaksFG_p[px].size==0) continue;
      int idx=nearestIndexPeak(maxMz, peaks_p, m_peaksFG_p[px].size, 0);
      int size=peaks_p[idx].highIdx-peaks_p[idx].lowIdx+1;
      if(peaks_p[idx].lowMz<=maxMz && peaks_p[idx].highMz>=maxMz && size<30) //into
        {
        pointsAcu_p[pk]+=size; //scans count
        massAcu+=peaks_p[idx].highMz-peaks_p[idx].lowMz; //mz cumulated
        pkCount++;//peak count
        }
    }
    double deltaMass=massAcu/pointsAcu_p[pk]; //average mass of the scan interval
    double pointInPeak=(double)pointsAcu_p[pk]/pkCount; //average scans at a peak
    tmpTolerance_p[pk]=1e6*pointInPeak*deltaMass/(2.5*maxMz);
    tmpMass_p[pk]=maxMz;
  }
  
  //Tolerance estimation across 10 mass segments.
  int segmentSize=nPeak/10;
  double *tolerance_p=0;
  try{
    tolerance_p=new double[segmentSize];
    }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }
  double mzInit, mzEnd;
  int initIdx, endIdx=-1, maxPointIdx=-1;
  NumericMatrix tol(10, 4);
  //For each segment, the tolerance associated with the peak exhibiting the 
  //greatest number of averaged points is extracted.
  for(int seg=0; seg<10; seg++) //10 segments
  {
    initIdx=endIdx+1;
    endIdx+=segmentSize;
    mzInit=tmpMass_p[initIdx];
    mzEnd=tmpMass_p[endIdx];
    double maxPoint=0;
    for(int i=initIdx; i<endIdx; i++)
    if(pointsAcu_p[i]>maxPoint) {maxPoint=pointsAcu_p[i]; maxPointIdx=i;}
    tolerance_p[seg]=tmpTolerance_p[maxPointIdx];
    //info a matrix
    tol(seg, 0)=mzInit;
    tol(seg, 1)=mzEnd;
    tol(seg, 2)=tmpTolerance_p[maxPointIdx];
    tol(seg, 3)=1e6/tmpTolerance_p[maxPointIdx]; 
  }
  printf("100\n");
  printf("\t\t\tpeaks under consideration:%d\n", nPeak);
  
  if(idxAxis_p)   delete [] idxAxis_p;
  if(pointsAcu_p) delete [] pointsAcu_p;
  if(pointsIdx_p) delete [] pointsIdx_p;
  if(tmpTolerance_p) delete []tmpTolerance_p;
  if(tmpMass_p)   delete []tmpMass_p;
  if(tolerance_p) delete []tolerance_p;
  
  return tol;
}
  
int ToleranceEstimation::mtToleranceEstimation(int spIndex)
  { 
    while(m_enable)
    {
      m_spectro[spIndex].mutexIn_p->lock(); //permission to continue (synchro)
      if(!m_enable) break; //The signal can be activated while waiting.
      IntensityPeak intPeak(m_SNR, false); //peak class
      
      //Spectrum conditioning.
      //Full spectra are received and only part of it may be of interest.
      
      int iMzLow=-1, iMzHigh=-1;
      int spSize, spSize_tmp, nPeak;
      int lowIdx, highIdx, maxIdx, px;
      spSize=m_spectro[spIndex].size; //spectrum size
      spSize_tmp=spSize;
      
      //SNR
      m_spectro[spIndex].noise=m_noiseEst_p->getSNR(m_spectro[spIndex].int_p, m_spectro[spIndex].size,  m_spectro[spIndex].SNR_p);
      
      if(!m_enable) //end of thread?
        {m_spectro[spIndex].mutexOut_p->unlock();  return 0;}
      
      //the peak are extracted from the spectrum (they are delimited by their indices).
      nPeak=intPeak.getPeakList(&m_spectro[spIndex], 0.75);  //intensity to zero if intensity<SNR*0.5

      if(nPeak>0)
      {
        gMutex.lock();
        gPeakCount+=nPeak; //peak accumulation synchronized.
        gSpectra++;    //processed spectra
        if(m_spectro[spIndex].mass_p[0]<gMinMass) gMinMass=m_spectro[spIndex].mass_p[0];
        if(m_spectro[spIndex].mass_p[spSize_tmp-1]>gMaxMass) gMaxMass=m_spectro[spIndex].mass_p[spSize_tmp-1];
        
        for(int i=0; i<m_NPixels; i++) if(m_pxLenIdx_p[i]==m_spectro[spIndex].pixel) {px=i; break;}
        gMutex.unlock();
        try{
          m_peaksFG_p[px].peaks_p=new PEAKS_F[nPeak];
          }
        catch(const std::bad_alloc& e)
        {
          printf("Error reserving memory: %s\n",e.what());
          return 0;
        }
        m_peaksFG_p[px].size=nPeak;
          
        for(int i=0; i<nPeak; i++) //copy peak
        {
          lowIdx =intPeak.getSinglePeak(i).low;
          maxIdx =intPeak.getSinglePeak(i).max;
          highIdx=intPeak.getSinglePeak(i).high;
          m_peaksFG_p[px].peaks_p[i].lowIdx =lowIdx;
          m_peaksFG_p[px].peaks_p[i].maxIdx =maxIdx;
          m_peaksFG_p[px].peaks_p[i].highIdx=highIdx;
          m_peaksFG_p[px].peaks_p[i].lowMz = m_spectro[spIndex].mass_p[lowIdx]; 
          m_peaksFG_p[px].peaks_p[i].maxMz = m_spectro[spIndex].mass_p[maxIdx]; 
          m_peaksFG_p[px].peaks_p[i].highMz= m_spectro[spIndex].mass_p[highIdx]; 
          m_peaksFG_p[px].peaks_p[i].maxIntensity= m_spectro[spIndex].int_p[maxIdx];
        }
      }
    m_spectro[spIndex].mutexOut_p->unlock(); //end of spectrum processing.
    }
    return 0;
  }

//getRawInfo()
//Loads the full spectrum information associated with a pixel from an imzML file.
//The information is stored in the m_spectro structure, set to the range [m_mzLow, m_mzHigh].
//px: Pixel whose spectrum should be loaded.
//spIndex: Threads that manage it
//Returns the size of the spectrum.
int ToleranceEstimation::getRawInfo(int iPx, int spIndex)
{
  int px=iPx;
  int massSize=m_getImzMLData_p->getPixelMassF(px, m_spectro[spIndex].mass_p);//mass vector
  m_spectro[spIndex].size=massSize;
  if(massSize<=0) return 0;
  
  //raw info of intensities  
  int intSize=m_getImzMLData_p->getPixelIntensityF(px, m_spectro[spIndex].int_p);
  
  m_spectro[spIndex].pixel=iPx;
  return massSize;
} 

//Returns the index of data.maxMz closest to value
//If nearest bits 1:0 == 00, returns the nearest
//If nearest bits 1:0 == 01, returns the nearest above
//If nearest bits 1:0 == 10, returns the nearest below
//If nearest bit    2 == 0 & value is out of range returns the nearest
//If nearest bit    2 == 1 & value is out of range returns -1
int ToleranceEstimation::nearestIndexPeak(float value, PEAKS_F *data, int size, int nearest)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  if(value<data[indexLow].maxMz)
  {
    if(nearest&4) return -1;
    else return 0;
  }
  if(value>data[indexHigh].maxMz)
  {
    if(nearest&4) return -1;
    else return indexHigh;
  }
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1)
  {
    if((nearest&3)==1)      return indexHigh;
    else if((nearest&3)==2) return indexLow;
    else
    {
      if(value-data[indexLow].maxMz <= data[indexHigh].maxMz-value) return indexLow;
      else return indexHigh;
    }
  }
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter].maxMz) return indexCenter;
    if(value<data[indexCenter].maxMz) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow) return indexLow;
    else if(indexHigh==indexLow+1)
    {
      if((nearest&3)==1)      return indexHigh;
      else if((nearest&3)==2) return indexLow;
      else
      {
        if(value-data[indexLow].maxMz <= data[indexHigh].maxMz-value) return indexLow;
        else return indexHigh;
      }
    }
  }
}
