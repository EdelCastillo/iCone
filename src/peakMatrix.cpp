/*********************************************************************************
 *     Segments
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



//Constructor
 //'  @param totalPixels: total pixels for all samples
 //'  @param    baseDir: directory for report files 
 //'  @param  massResolution: mass resolution for binning  stage
 //'  @param      mzLow: lower  mass to consider
 //'  @param     mzHigh: higher mass to consider
 //'  @param  pxSupport: minimum percentage of pixels that must support an ion for it to be considered.
 //'  @param  pxSamples: vector containing the number of pixels for each sample
 //'  @param   nSamples: Number of samples
 //'  @param    baseDir: Directory for data files
 //'  @param  intMethod: intensity values for the binning stage: mean, max
 //'  @param   nThreads: number of threads for parallel processing
  PeakMatrix::PeakMatrix(int totalPixels, double tolerance, double mzLow, double mzHigh, double pxSupport, int *pxSamples_p, int nSamples, char *baseDir, int intMethod, int nThreads)
{
  m_totalPixels=totalPixels; //It is also set in the loadGaussians() function.
  m_nSamples=nSamples;
  m_pxSamples_p=pxSamples_p;
  m_tolerance=tolerance;
  m_mzHigh=mzHigh;
  m_mzLow=mzLow;
  m_baseDir=baseDir;
  m_gaussians_p=0;
  m_massRange_p=0;
  m_pxSupport=pxSupport;
  m_pixelsCoordinates_p=0;
  m_intMethod=intMethod;
  m_centroidsGaussians_p=0;
  m_MRthr_p=0;
  m_MRthr_p=new MASS_RANGE[nThreads];
  m_MR_p=0;
  m_massAxis_p=0;
  m_spectro.SNR_p=0;
  m_nThreads=nThreads;
  m_centroids_p=new GAUSS_SP[nThreads];
  for(int i=0; i<nThreads; i++) {m_centroids_p[i].gauss_p=0; m_centroids_p[i].size=0;}
  for(int i=0; i<MAX_THREADS; i++)
  {
    m_thread_p[i]=0;
    m_status[i]=0; 
  }
}

//Destructor
PeakMatrix::~PeakMatrix()
{
  //printf("PeakMatrix destructor init\n");

  if(m_gaussians_p)
  {
    for(int i=0; i<m_totalPixels; i++)
    {
      if(m_gaussians_p[i].gauss_p) delete [] m_gaussians_p[i].gauss_p;
    }
    delete []m_gaussians_p;
  }
  if(m_centroidsGaussians_p) delete []m_centroidsGaussians_p;
  if(m_pixelsCoordinates_p) delete[] m_pixelsCoordinates_p;
 
  if(m_centroids_p)
  {
    for(int i=0; i<m_nThreads; i++) if(m_centroids_p[i].gauss_p) delete [] m_centroids_p[i].gauss_p;
    delete [] m_centroids_p;
  }
  if(m_spectro.SNR_p) delete [] m_spectro.SNR_p;
  if(m_MRthr_p)     delete [] m_MRthr_p;
  if(m_MR_p)        delete [] m_MR_p;
  if(m_massAxis_p)  delete [] m_massAxis_p;
  
  //printf("PeakMatrix destructor finish\n");
}

//Load the file with information about the Gaussian curves associated with each pixel
//If there is more than one sample, all its pixels are integrated
//fileName: temporary file name (".../_gaussians.bin")
//Returns false if the load failed
int PeakMatrix::loadGaussians(char *fileName)
{
  std::fstream fp;
  fp.open(fileName, std::fstream::in | std::ios::binary);
  if(!fp.is_open())
  {
    char txt[200];
    sprintf(txt, "Error: %s file could not be opened.\n", fileName);
    throw std::runtime_error(txt);
  }
  int nSamplePixels, nPxGauss;
  int gaussSize=3*sizeof(double);
  bool hit=true;
  
  int pxTotal=0, pxIndex=0;
//  int pixelsSample[MAX_SAMPLES];
  while(true) //first reading to obtain information.
  {
    fp.read((char*)&nSamplePixels, sizeof(int)); //#pixels into the sample
    if(fp.eof()) break;
//    pixelsSample[nSamples++]=nSamplePixels;
    pxTotal+=nSamplePixels;
    for(int pxSample=0; pxSample<nSamplePixels; pxSample++) //for each pixel of the sample 
    {
      fp.read((char*)&nPxGauss, sizeof(int)); //#gaussianas into the pixel
      if(fp.eof()) break;
      fp.seekg(gaussSize*nPxGauss, std::ios_base::cur); //to next pixel
    }
  }
  
  fp.close(); //If you exit with an error, you need to close the file.
  
  //second reading to obtain the coordinates.
  fp.open(fileName, std::fstream::in | std::ios::binary); //reopen
  try{
      m_gaussians_p=new GAUSS_SP[pxTotal]; //array of structs
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }
  
  for(int i=0; i<pxTotal; i++) {m_gaussians_p[i].gauss_p=0; m_gaussians_p[i].size=0;}  
  //  int totalGaussians=0;
  
  while(!fp.eof())
  {
    fp.read((char*)&nSamplePixels, sizeof(int)); //#pixels of sample
    
    for(int pxSample=0; pxSample<nSamplePixels; pxSample++) //for each pixel of the sample 
    {
      fp.read((char*)&nPxGauss, sizeof(int)); //#gaussians into the pixel
      if(fp.eof()) break;
      m_gaussians_p[pxIndex].size=nPxGauss; 
      if(nPxGauss==0) {pxIndex++; continue;}//px without conten
      
      try{
      m_gaussians_p[pxIndex].gauss_p=new GAUSS_PARAMS[nPxGauss]; //memory
      }
      catch(const std::bad_alloc& e)
      {
        printf("Error reserving memory: %s\n",e.what());
        return 0;
      }
      
      fp.read((char*)m_gaussians_p[pxIndex].gauss_p, nPxGauss*gaussSize ); //load
      
      if(fp.eof() || fp.fail() || fp.bad()) //boundary control
      {
        char txt[200];
        sprintf(txt, "Error: %s file could not be read completely.", fileName);
        throw std::runtime_error(txt);
        hit =false; break;
      }
      
      if(fp.eof() || !hit) break;
      pxIndex++;
    }
  }
  fp.close();
  m_totalPixels=pxIndex;
  return pxIndex;
}

//getCentroidsIntoRange()
//Extracts information from all existing Gaussians in m_gaussians_p.
//     mass_p: vector of internally generated masses.
//       px_p: vector of internally generated pixels.
//intensity_p: vector of intensities
//   massSize: size of reserved memory for each vector.
//Returns the number of masses.
int PeakMatrix::getCentroidsIntoRange(double *mass_p, int *px_p, double *intensity_p, int massSize)
{
  int count=0;

  for(int px=0; px<m_totalPixels; px++)//for all pixels
  {
    if(m_gaussians_p[px].gauss_p==0) continue; //pixel without Gaussians
    for(int i=0; i<m_gaussians_p[px].size; i++) 
    {
      if(count>=massSize) break;
      else 
      {
        mass_p[count]=m_gaussians_p[px].gauss_p[i].mean;
        if(mass_p[count]==0) {continue;}
        px_p[count]=px;
        intensity_p[count]=m_gaussians_p[px].gauss_p[i].weight;
        count++;
      }
    }
  }
  return count;
}


//Loads the file with the coordinates of each pixel
//If there is more than one sample, all its pixels are integrated
//update m_nSamples variable
//Returns false if the load failed
int PeakMatrix::loadPixelsCoordinates(char *fileName, int totalPixels)
{
  std::fstream fp;
  fp.open(fileName, std::fstream::in | std::ios::binary);
  if(!fp.is_open())
  {
    char txt[200];
    sprintf(txt, "Error: %s file could not be opened\n", fileName);
    throw std::runtime_error(txt);
  }
  int nSamplePixels, coordSize=2*sizeof(int);
  bool hit=true;
  
//  int nSamples=0;
  try{
      m_pixelsCoordinates_p=new PIXEL_XY[totalPixels]; //coordinates
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }

  int pxTotal=0;
  while(true)
  {
    fp.read((char*)&nSamplePixels, sizeof(int)); //#pixels into the sample
    if(fp.eof()) break;
//    m_pixelsSample[nSamples++]=nSamplePixels;
    for(int pxSample=0; pxSample<nSamplePixels; pxSample++) //for each pixel into the sample
    {
      fp.read((char*)&m_pixelsCoordinates_p[pxSample], coordSize); //#gaussianas into pixel
      if(fp.fail() || fp.bad()) //fault control
      {
        char txt[200];
        sprintf(txt, "Error: %s file could not be read completely.", fileName);
        throw std::runtime_error(txt);
        hit =false; break;
      }
      
      if(fp.eof() || !hit) break;
      pxTotal++;
    }
  }
  fp.close();
  if(0)//if(nSamples != m_nSamples)
  {
    printf("ERROR: inconsistency in the number of samples.\n");
    return -1;
  }
  return pxTotal;
}

//sets the number of pixels per sample based on info in _pixelsCoord.bin
//result stored in vector m_pixelsSample[]
//returns the number of samples, or -1 on failure
int PeakMatrix::getSamplesPixelNumber()
{
  std::fstream fp;
  char fileName[200];
  strcpy(fileName, m_baseDir);
  strcat(fileName, (char*)"_pixelsCoord.bin");

  fp.open(fileName, std::fstream::in | std::ios::binary);
  if(!fp.is_open())
  {
    char txt[200];
    sprintf(txt, "Error: %s file could not be opened\n", fileName);
    throw std::runtime_error(txt);
  }
  int nSamplePixels;
  bool hit=true;
  
  int nSamples=0;
  while(true)
  {
    fp.read((char*)&nSamplePixels, sizeof(int)); //#pixels into the sample
    if(fp.eof()) break;
    m_pxSamples_p[nSamples++]=nSamplePixels;
    fp.seekg((std::streampos)(nSamplePixels*sizeof(PIXEL_XY)), std::ios_base::cur);
    if(fp.fail() || fp.bad()) //fault control
      {
        char txt[250];
        sprintf(txt, "Error: %s file could not be read completely.", fileName);
        throw std::runtime_error(txt);
        hit =false; break;
      }
      
    if(fp.eof() || !hit) break;
  }
  fp.close();
  if(!hit) return -1;
  return nSamples;
}


//getCentroids()

//Based on the Gaussians, the centroids within a mass segment are determined.
//1) A histogram is generated based on Gaussian centers and added uncertainty. 
// It serves to delineate regions of interest for subsequent parallel processing. Noise matters.
//2)Each mass segment is treated separately.
// a)An ordered vector containing the Gaussian centers of all pixels is created.
// b)Iterative averaging, constrained by a tolerance threshold, is used to generate a proposal for centroids.
// c)This estimate is used as the initialization centroid for the k-means algorithm.
//   The dataset consists of the centers of Gaussians centered on the initialization value and delimited by the tolerance.
// d)An intensity value is associated with each centroid for every pixel. 
//   A pixel contributes the maximum value from the Gaussian nearest to the centroid, provided it falls within the tolerance. 
//   If a pixel contributes to more than one Gaussian, either the maximum value or the mean value is extracted, depending on the `intMethod` argument.
//   A representative intensity value is also generated for the centroid. It is the average or maximum value, depending on the intMethod parameter.
//The centroid information is saved in a file with a name ending in _peakMatrix_n.bin, where n is the thread number.
//These files are destroyed upon completion, after all the information has been consolidated into a new file (_peakMatrix.bin).
int PeakMatrix::getCentroids()
{
  Common tools;

  m_deltaMass=m_tolerance*m_mzLow/(10*1e6); //delta=1/4 of the minimum mass increment of the spectrometer
  m_massAxisSize=1+(m_mzHigh-m_mzLow)/m_deltaMass;
  try{
      m_massAxis_p= new double[m_massAxisSize];
      m_spectro.SNR_p=new double[m_massAxisSize];
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }
  
  //First parallel stage: the mass axis is generated via addition based on Gaussian centers. Oversampling = 10.
  //Each center of mass has an associated uncertainty equal to half the tolerance.
  ///////////////////////////////////////////////////////////////////////////////
  
  setWorkThr(m_massAxisSize); //Task distribution for the threads (equal number of elements).
  try{  
      for(int i=0; i<m_nThreads; i++)
        m_thread_p[i]=new std::thread(&PeakMatrix::mtMassAxis, this, i);
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }
  
  //wait until the conclusion-
  for(int i=0; i<m_nThreads; i++)
  {
    m_thread_p[i]->join();
    delete m_thread_p[i]; m_thread_p[i]=0;
  }

  //Noise is estimated and isolated segments are detected for subsequent parallel processing.
  //////////////////////////////////////////////////////////////////////////////////////////

  m_spectro.int_p=m_massAxis_p;
  m_spectro.size=m_massAxisSize;
  
  NoiseEstimation noiseEst_p(2, 1, 9); //sd
  m_noise=noiseEst_p.getSNR(m_spectro.int_p, m_spectro.size, m_spectro.SNR_p);
  m_spectro.noise=m_noise;

  //the occupied mass segments (delimited by noise) are determined.
  bool state=0;
  int count=0;
  //required memory
  for(int i=0; i<m_massAxisSize; i++)
  {
    if(state==0 && m_massAxis_p[i]>=m_noise) {state=1;}
    else if(state==1 && m_massAxis_p[i]<m_noise) {state=0; count++;}
  }
  if(state==1) count++; //el último 
  try{
      m_MR_p=new MASS_RANGE[count+1];
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }
  
  state=0;
  count=0;

  for(int i=0; i<m_massAxisSize; i++)
  {
    if     (state==0 && m_massAxis_p[i]>=m_noise) {state=1; m_MR_p[count].low=i;}
    else if(state==1 && m_massAxis_p[i]<m_noise)  {state=0; m_MR_p[count++].high=i;}
  }
  
  if(state==1) m_MR_p[count++].high=m_massAxisSize-1; //the last 
  m_MRsize=count;
//  printf("segments:%d noise:%.2f\n", count, m_noise);  
 
  //Second parallel stage: peaks are detected and centroids are estimated.
  //Each thread generates a file containing its own results (fileName_peakMatrix_n, where n=thread).
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  setWorkThr(m_MRsize); //Task distribution for the threads (same number of segments).
  try{  
      for(int i=0; i<m_nThreads; i++)
        m_thread_p[i]=new std::thread(&PeakMatrix::mtSegments, this, i);
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return 0;
  }
  //Progress information for the user.
  float status=0, vez=1;
  while(status<0.9)
  {
    status=0;
    for(int i=0; i<m_nThreads; i++)
    {
      status+=m_status[i];
    }
    status/=m_nThreads;
    if(status>vez*0.1) {if(vez<10) printf("%.0f ", vez*10); vez+=1.0;}
  }
 
    //wait until the conclusion-
  for(int i=0; i<m_nThreads; i++)
  {
    m_thread_p[i]->join();
    delete m_thread_p[i]; m_thread_p[i]=0;
  }


 //The information from each threads file is consolidated into a single file.
 ///////////////////////////////////////////////////////////////////////////

char file[250];
sprintf(file,"%s_peakMatrix.bin", m_baseDir);
std::fstream fpOut, fpIn;
fpOut.open(file, std::fstream::out | std::ios::binary | std::ios::trunc);
if(!fpOut.is_open())
{
  char txt[200];
  sprintf(txt, "Error: The internal %s file could not be created.\n The peak matrix cannot be saved.\n", file);
  throw std::runtime_error(txt);
}

std::streampos initPos, finalPos;
std::streampos ionPos, colSize, offset;
int totalSamples, totalPx, totalIons=0, nIonsThr;
Common common;
int pxSample[MAX_SAMPLES];

//Each thread generates a file containing partial centroids. Here, they are aggregated into a single file.
for(int thrIdx=0; thrIdx<m_nThreads; thrIdx++)
  {
  sprintf(file,"%s_peakMatrix_%d.bin", m_baseDir, thrIdx);
  fpIn.open(file, std::fstream::in | std::ios::binary);
    if(!fpIn.is_open())
    {
      char txt[300];
      sprintf(txt, "Error: The internal %s file could not be created.\n The peak matrix cannot be saved.\n", file);
      throw std::runtime_error(txt);
    }
    //info común (read)
      //metadata
      fpIn.read((char*)&totalSamples, sizeof(int)); //samples number
      fpIn.read((char*)&totalPx,      sizeof(int)); //total pixels in all samples
      fpIn.read((char*)&nIonsThr,    sizeof(int)); //ions number
      if(fpIn.fail()) 
        {
        printf("Error while reading in %s", file);
        fpIn.close(); return 0;
        }
      fpIn.read((char*)pxSample,  totalSamples*sizeof(int)); //px in each sample
      
      //info común (write)
    if(thrIdx==0)
      {
      fpOut.write((char*)&totalSamples, sizeof(int)); //samples number
      fpOut.write((char*)&totalPx,      sizeof(int)); //total pixels in all samples
      fpOut.write((char*)&nIonsThr,    sizeof(int)); //ions number
      if(fpOut.fail()) 
        {
        printf("Error while writing in %s", file);
        fpOut.close(); return 0;
        }
      fpOut.write((char*)pxSample,  totalSamples*sizeof(int)); //px in each sample
      }
    totalIons+=nIonsThr;
    
    //The data on disk is organized by centroids.
    //Each centroid is accompanied by a list of pixels from all samples in pixel-intensity format.
    //To separate them, the pixel range associated with each sample must be known.    
    //The minimum and maximum pixels in each sample are delimited.
    PIXEL_XY samplesPxLimit[totalSamples]; //x -> low; y -> high
    samplesPxLimit[0].x=0;
    samplesPxLimit[0].y=pxSample[0]-1;
    for(int i=1; i<totalSamples; i++)
    {
      samplesPxLimit[i].x=samplesPxLimit[i-1].y+1;
      samplesPxLimit[i].y=samplesPxLimit[i].x+pxSample[i]-1;
    }
    //The information is copied to the peak matrix file.
    double tmpMass, tmpIntensity, tmpTolerance;
    int tmpPx, nPx;
    
    for(int ion=0; ion<nIonsThr; ion++) //para cada centroide
    {
      //Info on this centroid.
      fpIn.read((char*)&tmpMass, sizeof(double));
      fpIn.read((char*)&tmpIntensity, sizeof(double));
      fpIn.read((char*)&tmpTolerance, sizeof(double));
      fpIn.read((char*)&nPx, sizeof(int));
      if(fpIn.fail()) 
      {
        printf("Error while reading in %s", file);
        fpIn.close(); return 0;
      }

    //Copy the centroid to the destination file. 
    fpOut.write((char*)&tmpMass, sizeof(double));
    fpOut.write((char*)&tmpIntensity, sizeof(double));
    fpOut.write((char*)&tmpTolerance, sizeof(double));
    fpOut.write((char*)&nPx, sizeof(int));
    if(fpOut.fail()) 
    {
      printf("Error while writing in %s", file);
      fpOut.close(); return 0;
    }
    
    //pixels and their intensities.
    for(int i=0; i< nPx; i++)
    {
        fpIn.read((char*)&tmpPx, sizeof(int));
        fpIn.read((char*)&tmpIntensity, sizeof(double));
        if(fpIn.fail()) 
        {
          printf("Error while reading in %s", file);
          fpIn.close(); return 0;
        }
      fpOut.write((char*)&tmpPx, sizeof(int));
      fpOut.write((char*)&tmpIntensity, sizeof(double));
      if(fpOut.fail()) 
      {
        printf("Error while reading in %s", file);
        fpOut.close(); return 0;
      }
    }
  }

    fpIn.close();
    sprintf(file,"rm %s_peakMatrix_%d.bin", m_baseDir, thrIdx);
    
    system(file); //The thread's temporary files are deleted.
  }
fpOut.seekg((std::streampos)(2*sizeof(int)), std::ios_base::beg);
fpOut.write((char*)&totalIons, sizeof(int));

fpOut.close();


printf("100\n");
printf("\t\t\ttotal centroids:%d\n", totalIons);
return totalIons; 
}


//Parallel operation. 
//1) A vector is generated with spacing equal to one-tenth of the tolerance set for the lowest mass. 
//2) One unit is accumulated in each cell for every Gaussian distribution coinciding with that mass interval, 
//assuming a unified standard deviation equal to half the tolerance for the mass in question.
//It serves to delineate regions of interest for subsequent parallel processing. Noise matters.
void PeakMatrix::mtMassAxis(int thrIdx)
{
  double tolerance=m_tolerance, deltaMass=m_deltaMass, mzLow=m_mzLow, mzHigh=m_mzHigh;
  double maxMass, sigma, deltaMass2, lowMass, highMass;
  int iLow, iHigh, massAxisSize=m_massAxisSize;
  
  int lowSegIdx =round(m_MRthr_p[thrIdx].low); 
  int highSegIdx=round(m_MRthr_p[thrIdx].high); 
  if(highSegIdx<lowSegIdx) return;
  double lowSegMass =mzLow+lowSegIdx *deltaMass; //extreme masses for the thread
  double highSegMass=mzLow+highSegIdx*deltaMass;
  double *massAxis_p=m_massAxis_p;

  for(int i=lowSegIdx; i<=highSegIdx; i++) {massAxis_p[i]=0;} //init
  
  //for all elements of the spectrum
  for(int px=0; px<m_totalPixels; px++)
  {
    if(m_gaussians_p[px].gauss_p==0 || m_gaussians_p[px].size==0) {continue;} //if this entry does not contain info.
    
    for(int i=0; i<m_gaussians_p[px].size; i++) 
    {
      maxMass=m_gaussians_p[px].gauss_p[i].mean; //gaussian central mass
      
      //Some ions may be zero in mean/sigma/weight. They are discarded.
      if(maxMass<mzLow || maxMass>mzHigh) 
        continue;
      deltaMass2=tolerance*maxMass/1e6;
      sigma=deltaMass2/2.0;
      //sigma=m_gaussians_p[px].gauss_p[i].sigma;

      //Very wide Gaussians (sigma>5*deltaMass) are discarded.
      if(fabs(m_gaussians_p[px].gauss_p[i].sigma)>5*deltaMass2) //sigma >>
        continue;
      
      //vector points occupied by the peak.
      lowMass =maxMass-sigma;
      highMass=maxMass+sigma;
      if(lowMass<lowSegMass || highMass>highSegMass) continue;
      iLow =(lowMass -mzLow)/deltaMass;
      iHigh=(highMass-mzLow)/deltaMass;
      for(int j=iLow; j<=iHigh && j<massAxisSize; j++) 
        massAxis_p[j]+=1.0;
    }
  }
}


//It establishes the task distribution lineal for each thread.
//Each thread receives indices for the start and end segments.
/////////////////////////////////////////////////////////////
void PeakMatrix::setWorkThr(int vectorSize)
{
  int nSegThr=vectorSize/m_nThreads;
  if(nSegThr==0) //vectorSize<m_nThreads => A single thread bears the entire load.
    {
    m_MRthr_p[0].low=0;  //integer over double
    m_MRthr_p[0].high=vectorSize-1;
    for( int i=1; i<m_nThreads; i++) //The remaining threads are cancelled.
      {
      m_MRthr_p[i].low=1;  
      m_MRthr_p[i].high=0; //high<low
      }
    }
  else //one or more segments per thread
    {
    for(int i=0; i<m_nThreads-1; i++)
      {
      m_MRthr_p[i].low=i*nSegThr;  //integer over double
      m_MRthr_p[i].high=(i+1)*nSegThr-1;
      }
    m_MRthr_p[m_nThreads-1].low=(m_nThreads-1)*nSegThr;
    m_MRthr_p[m_nThreads-1].high=vectorSize-1;
    }
}


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
void PeakMatrix::mtSegments(int thrIdx)
{
  Common tools;
  double tolerance=m_tolerance;
  int lowSegIdx =round(m_MRthr_p[thrIdx].low); 
  int highSegIdx=round(m_MRthr_p[thrIdx].high);
  if(highSegIdx<lowSegIdx) return;
  
  int lowMassAxisIdx =round(m_MR_p[lowSegIdx].low);
  int highMassAxisIdx=round(m_MR_p[highSegIdx].high);
  int segSize=highMassAxisIdx-lowMassAxisIdx+1;
  
  double localDeltaMass=m_deltaMass;
  double mzLow =m_mzLow+lowMassAxisIdx *localDeltaMass;
  double mzHigh=m_mzLow+highMassAxisIdx*localDeltaMass;

  //determines the number of Gaussians in the interval for subsequent memory allocation.
  int gCount=0;
  for(int i=0; i<m_totalPixels; i++)
  {
    if(m_gaussians_p[i].gauss_p && m_gaussians_p[i].size>0)
    {
      int idxLow =tools.nearestIndexGaussians(mzLow,  m_gaussians_p[i].gauss_p, m_gaussians_p[i].size);
      int idxHigh=tools.nearestIndexGaussians(mzHigh, m_gaussians_p[i].gauss_p, m_gaussians_p[i].size);
      gCount+=idxHigh-idxLow+1;
    }
  }
  //memory allocation for the centers of the Gaussians.
  double *tmpMassAxis_p=0, *massAxis_p=0;
  int    *tmpMassAxisIdx_p=0;
  try{
      tmpMassAxis_p   =new double[gCount];
      massAxis_p      =new double[gCount];
      tmpMassAxisIdx_p=new int[gCount];
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return;
  }
  
  gCount=0;
  //storage of the Gaussian centers
  for(int i=0; i<m_totalPixels; i++)
  {
    if(m_gaussians_p[i].gauss_p && m_gaussians_p[i].size>0)
    {
      int idxLow =tools.nearestIndexGaussians(mzLow,  m_gaussians_p[i].gauss_p, m_gaussians_p[i].size);
      int idxHigh=tools.nearestIndexGaussians(mzHigh, m_gaussians_p[i].gauss_p, m_gaussians_p[i].size);
      for(int j=idxLow; j<=idxHigh; j++)
      {
        tmpMassAxis_p[gCount]=m_gaussians_p[i].gauss_p[j].mean;
          gCount++;
      }
    }
  }
  //ascending order
  //massAxis_p stores the sorted central masses of the Gaussians for all pixels.
  tools.sortUp(tmpMassAxis_p, tmpMassAxisIdx_p, gCount); 

  for(int i=0; i<gCount; i++)
    massAxis_p[i]=tmpMassAxis_p[tmpMassAxisIdx_p[i]];
  int massAxisSize=gCount;

  if(tmpMassAxisIdx_p)delete [] tmpMassAxisIdx_p;

  //centroid estimation by iterative averaging.
  /////////////////////////////////////////////
  double *centers_p=0;
  int *centersSize_p=0;
  try{
      centers_p=new double[massAxisSize];
      centersSize_p=new int[massAxisSize];
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return;
  }
  
  for(int i=0; i<massAxisSize; i++) centersSize_p[i]=0;
  
  int idxLow =tools.nearestIndex(mzLow,  massAxis_p, massAxisSize); //lower mass index
  int idxHigh=tools.nearestIndex(mzHigh, massAxis_p, massAxisSize);
  int nIons=0;
  centers_p[nIons]=massAxis_p[idxLow];
  double segmentSize=centers_p[nIons]*tolerance/1e6;
  double massDiff;

  centersSize_p[0]=1;
  for(int i=idxLow+1; i<=idxHigh; i++)//iterative averaging over the interval, taking tolerance into account.
  {
    massDiff=fabs(massAxis_p[i]-centers_p[nIons]);
    if(massDiff<segmentSize) //within tolerance
      {
      centers_p[nIons]=(1.0/((double)centersSize_p[nIons]+1.0))*((double)centersSize_p[nIons]*centers_p[nIons]+massAxis_p[i]);
      centersSize_p[nIons]++; //Number of elements in the bin.
      }
    // new ion if the pixels provide sufficient support.
    else //out of tolerance
      {
      if(centersSize_p[nIons]>=m_pxSupport) nIons++; //new ion (Initial approach to offloading work.)
      
      centers_p[nIons]=massAxis_p[i]; //ion fail
      centersSize_p[nIons]=1;
      segmentSize=massAxis_p[i]*tolerance/1e6;
      }
  }
  //The possible final peak is coming to an end.
  if(massDiff<segmentSize && centersSize_p[nIons]>=m_pxSupport) 
    nIons++;  

  
  //k-means Clustering
  int clusterSize;
  int maxIter=40; //maximum iterations
  double convergenceValue=1e-6;
  double mass, deltaMass, lowMass, highMass;
  int lowMassIdx, highMassIdx;
  try{
      m_centroids_p[thrIdx].gauss_p=new GAUSS_PARAMS[nIons];
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return;
  }
  
  for(int i=0; i<nIons; i++)
  {
    //The set of participating central masses is determined.
    mass=centers_p[i];
    deltaMass=tolerance*mass/1e6;
    lowMass =mass-deltaMass/2;
    highMass=mass+deltaMass/2;
    lowMassIdx =tools.nearestIndex(lowMass,  massAxis_p, massAxisSize);
    highMassIdx=tools.nearestIndex(highMass, massAxis_p, massAxisSize);
    clusterSize=0;
    for(int j=lowMassIdx; j<=highMassIdx; j++)
      tmpMassAxis_p[clusterSize++]=massAxis_p[j]; 
    
    KmeansR kmeans(tmpMassAxis_p, clusterSize, 1, maxIter, convergenceValue); //k-means constructor
    kmeans.getClusters(1, &mass); 
    
    //results
    m_centroids_p[thrIdx].gauss_p[i].mean=kmeans.m_kStruct.clusters_p[0].center;
    m_centroids_p[thrIdx].gauss_p[i].sigma=sqrt(kmeans.m_kStruct.clusters_p[0].withinss);
    m_centroids_p[thrIdx].gauss_p[i].weight=1;

  }
  m_centroids_p[thrIdx].size=nIons;
  
  if(tmpMassAxis_p)   delete [] tmpMassAxis_p;
  if(centers_p)       delete [] centers_p;
  /////////////////////////////////////////7

  //info to file
  std::fstream fp;
  char fileName[200];
  sprintf(fileName,"%s_peakMatrix_%d.bin", m_baseDir, thrIdx);
  fp.open(fileName, std::fstream::out | std::ios::binary | std::ios::trunc);
  if(!fp.is_open())
  {
    char txt[300];
    sprintf(txt, "Error: The internal file %s could not be created.\n The peak matrix cannot be saved.\n", fileName);
    throw std::runtime_error(txt);
  }
  std::streampos initPos, finalPos;
  
  fp.write((char*)&m_nSamples, sizeof(int)); //samples number
  fp.write((char*)&m_totalPixels, sizeof(int)); //space for matrix rows
  initPos=fp.tellg();
  fp.write((char*)&m_nCentroids, sizeof(int)); //space for matrix cols
  for(int i=0; i<m_nSamples; i++)
    fp.write((char*)&m_pxSamples_p[i], sizeof(int)); //space for pixel number into each sample
  
  double dispersion, max;
  int pxCount, peakCount=0;
  int *nPx_p=0;
  try{
      nPx_p=new int[m_totalPixels];
      }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return;
  }
  
  int vez=1;
  double maxMass, sigma;
  m_status[thrIdx]=1; //avoids problems in case nIons=0

  //management of intensities and dispersion
  for(int i=0; i<nIons; i++)
  {
    pxCount=0;
    dispersion=0;
    max=m_centroids_p[thrIdx].gauss_p[i].mean;    
    double intensity, maxInt=0, meanInt=0;
    double tmpTolerance=max*m_tolerance/1e6;
    
    for(int px=0; px<m_totalPixels; px++)
    {
      if(m_gaussians_p[px].gauss_p==0 || m_gaussians_p[px].size==0) {continue;} //if this entry does not contain info.
      //Some ions may be zero in mean/sigma/weight. They are discarded.
      
      int gIdx=tools.nearestIndexGaussians(max, m_gaussians_p[px].gauss_p, m_gaussians_p[px].size);
      maxMass=m_gaussians_p[px].gauss_p[gIdx].mean;
      double averageValue=0, maxValue=0;
      
      //mass range for this pixel within tolerance.
      int idxLow=gIdx, idxHigh=gIdx;
      if(gIdx>0)  //limits control
        for(idxLow=gIdx-1; idxLow>=0; idxLow--)
        {
          if(fabs(m_gaussians_p[px].gauss_p[idxLow].mean-maxMass)>tmpTolerance) 
              {idxLow ++; break;}
        }
      if(gIdx<m_gaussians_p[px].size-1) //limits control
        for(idxHigh=gIdx+1; idxHigh<m_gaussians_p[px].size; idxHigh++)
        {
          if(fabs(m_gaussians_p[px].gauss_p[idxHigh].mean-maxMass)>tmpTolerance) 
              {idxHigh--; break;}
       }
      //valid indices
      for(int idx=idxLow; idx<=idxHigh; idx++)
      {
        double value=m_gaussians_p[px].gauss_p[idx].weight;
        averageValue+=value;
        if(value>maxValue) maxValue=value;
      }
      if(idxHigh>=idxLow)
        averageValue/=(idxHigh-idxLow+1.0);
      else   averageValue=0;
      
      if(maxMass<m_mzLow || maxMass>m_mzHigh) {continue;}
      //sigma=maxMass*tolerance/2e6; //SD from mass
      sigma=m_gaussians_p[px].gauss_p[gIdx].sigma; //SD from gaussians
      if(fabs(max-maxMass)<=3*sigma)
      {
        if(m_intMethod==MAX) //max value
          intensity=maxValue;
        else{
          intensity=averageValue; //mean value
          meanInt+=intensity;
        }
        massAxis_p[pxCount]=intensity; //Note: memory is reused.
        if(intensity>maxInt) maxInt=intensity;
        dispersion+=(max-maxMass)*(max-maxMass);
        nPx_p[pxCount]=px;
        pxCount++;
      }
    }
    if(pxCount>=m_pxSupport) //ions that are sufficiently supported.
    {
      peakCount++;
      //centroid header in the file
      dispersion=sqrt(dispersion/pxCount);
      meanInt/=pxCount;
      fp.write((char*)&max, sizeof(double)); //centroid mass
      if(m_intMethod==MAX)
        fp.write((char*)&maxInt, sizeof(double));      //centroid intensity
      else
        fp.write((char*)&meanInt, sizeof(double));      //centroid intensity
      double centroidTolerance=1e6*dispersion/max;  
      fp.write((char*)&centroidTolerance, sizeof(double)); //centroid tolerance (ppm)
      fp.write((char*)&pxCount, sizeof(int));         //px number support
      
      //pixel-intensity pairs
      for(int i=0; i<pxCount; i++)
      {
        fp.write((char*)&nPx_p[i], sizeof(int)); //px 
        fp.write((char*)&massAxis_p[i], sizeof(double)); //px intensity
      }
    }
    m_status[thrIdx]=(i+1.0)/(float)nIons; //para el observador: evolución de tareas
  }
  fp.seekg(initPos, std::ios_base::beg);
  fp.write((char*)&peakCount, sizeof(int)); //update nIons
  
  fp.close();
  if(nPx_p)           delete [] nPx_p;
  if(massAxis_p)      delete [] massAxis_p;
}

