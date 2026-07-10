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

int  gProcSegments, gVez; //observers.
std::mutex gMutex;

/// R METHODS ////////////////////////////////////////////////////////////////////////

//' @name rGetMetaDataFromFile
//' @title returns generic information about the peak matrix located in a file.
//' @param file -> file name with peak matrix (tmpPeakMatrix.bin)
//' @return a list:
//'   nSamples     -> number of samples analyzed.
//'   totalPx      -> total number of pixels (cumulative of each sample).
//'   nIons        -> number of columns in the matrix.
//'   pixelsSample -> vector with the pixels in each sample.

// [[Rcpp::export]]
List rGetMetaDataFromFile(const char* file)
{
  std::fstream fp;
  std::streampos pos;
  int nSamples, totalPx, nIons, colSize;
  
  fp.open(file, std::fstream::in | std::ios::binary);
  if(!fp.is_open())
  {
    char txt[200];
    sprintf(txt, "Error: The internal file %s could not be created.\n The peak matrix cannot be saved.\n", file);
    throw std::runtime_error(txt);
  }
  int px;
  fp.read((char*)&nSamples, sizeof(int)); //samples number
  fp.read((char*)&totalPx,  sizeof(int)); //samples number
  fp.read((char*)&nIons,    sizeof(int)); //samples number
  
  NumericVector pxSample(nSamples);
  for(int i=0; i<nSamples; i++)
  {
    fp.read((char*)&px, sizeof(int)); //samples number
    pxSample[i]=px;
  }
  fp.close();
  List ret=List::create(Named("nSamples")=nSamples, Named("totalPx")=totalPx, Named("nIons")=nIons, 
                        Named("pixelsSample")=pxSample);
  return ret;
}

//' @name rGetCoordinatesFromFile()
//' @title returns a matrix with the coordinates of all pixels (X/Y).
//' If there are multiple samples, they appear sequentially; that is, the matrix has as many rows 
//' as the cumulative number of pixels in each sample and two columns.
//' @param file   -> file name with pixels coordinates (tmpPixelsCoordinates.bin)
//' @param sample -> just download the pixels from this sample.
//'                 if sample < 0, all sample coordinates are returned
//' @return a matrix with the coordinates (X/Y) of pixels.
//' 
// [[Rcpp::export]]
NumericMatrix rGetCoordinatesFromFile(const char* fileName, int sample)
{
  std::fstream fp;
  std::streampos initPos;
  fp.open(fileName, std::fstream::in | std::ios::binary);
  if(!fp.is_open())
  {
    char txt[200];
    sprintf(txt, "Error: %s file could not be opened\n", fileName);
    throw std::runtime_error(txt);
  }
  int totalSamples=0, nSamples, nSamplePixels, nPixels, coordSize=2*sizeof(int);
  bool hit=true;
  int pxTotal=0, pxIndex=0;
  int pixelsSample[MAX_SAMPLES];
  while(true) //first reading to obtain information.
  {
    fp.read((char*)&nSamplePixels, sizeof(int)); //#pixels into the sample
    if(fp.eof()) break;
    fp.seekg((std::streampos)nSamplePixels*2*sizeof(int), std::ios_base::cur);
    
    pixelsSample[totalSamples++]=nSamplePixels;
    pxTotal+=nSamplePixels;
  }

  sample--; //to C++
  if(sample>totalSamples)
  {
    printf("warning: max samples in file are %d\n", nSamples);
    fp.close();
    return 0;
  }
  
  //get pixels range from sample parameter
  int initPx=0;
  if(sample<0) //all samples
  {
    nSamples=totalSamples; 
    initPos=0;
    printf("The %d samples are considered.\n", nSamples);
  }
  else //only one sample
  {
    nSamples=1;
    for(int i=0; i<sample; i++)
    {
      initPos+=(1+pixelsSample[i]*2)*(std::streampos)sizeof(int); //#px, X0/Y0, X1/Y1, ...
    }
    pxTotal=pixelsSample[sample];
  }

  NumericMatrix pxCoord(pxTotal, 2);
  fp.close(); //If you exit with an error, you need to close the file.
  
  //second reading to obtain the coordinates.
  fp.open(fileName, std::fstream::in | std::ios::binary); //reopen
  fp.seekg(initPos, std::ios_base::beg); //file position
  
  int xy[2];
  
  for(int sample=0; sample<nSamples; sample++) //for all samples
  {
    fp.read((char*)&nPixels, sizeof(int)); //#pixels into the sample
    for(int pxSample=0; pxSample<nPixels; pxSample++) //for each pixel into the sample
      {
      fp.read((char*)xy, 2*sizeof(int)); //XY coordinate
      pxCoord(pxIndex,0)  =xy[0];
      pxCoord(pxIndex++,1)=xy[1];
      if(fp.fail() || fp.bad()) //fault control
        {
          char txt[200];
          sprintf(txt, "Error: %s file could not be read completely.", fileName);
          throw std::runtime_error(txt);
          hit =false; break;
        }
        
        if(fp.eof() || !hit) break;
      }
  }
  
  if(pxIndex!=pxTotal) printf("warning: not all coordinates were loaded.\n");
  fp.close();
  return pxCoord;
}



 //' @name rGetMassVectorFromFile()
 //' @title returns a vector with all the masses in peak matrix 
 //' @param file -> file name with peak matrix (tmpPeakMatrix.bin)
 //' @return mass vector
 //' 
 // [[Rcpp::export]]
NumericVector rGetMassVectorFromFile(const char* file)
 {
   std::fstream fp;
   std::streampos ionPos, colSize, offset;
   int nSamples, totalPx, nIons, massIndex;
   int pxSample[nSamples];
   Common common;
   
   fp.open(file, std::fstream::in | std::ios::binary);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: The internal %s file could not be created.\n The peak matrix cannot be saved.\n", file);
     throw std::runtime_error(txt);
   }
   fp.read((char*)&nSamples, sizeof(int)); //samples number
   fp.read((char*)&totalPx,  sizeof(int)); //total pixels in all samples
   fp.read((char*)&nIons,    sizeof(int)); //ions number
   fp.read((char*)pxSample,  nSamples*sizeof(int)); //px in each sample
   
   colSize=((totalPx+2)*sizeof(float)+sizeof(int));
   offset=3*sizeof(int)+nSamples*sizeof(int);

   //get matrix column from mass parameters
   NumericVector massVect(nIons);
   float tmpMass=0;
   for(int i=0; i<nIons; i++)
     {
       ionPos=offset+(i*colSize);
       fp.seekg(ionPos, std::ios_base::beg);  
       fp.read((char*)&tmpMass, sizeof(float));
       massVect(i)=tmpMass;
     }
   fp.close();
   return massVect;
 }

//' @name rGetIntensityFromFile()
 //' @title returns a column whit pixels intensities from the peak matrix: 
 //' 
 //' @param file -> file name with peak matrix (tmpPeakMatrix.bin)
 //' @param column -> desired column (first = 1)
 //' @return a vector with intensity of each pixel.
 
 // [[Rcpp::export]]
 NumericVector rGetIntensityFromFile(const char* file, int column)
 {
   std::fstream fp;
   std::streampos pos, colSize;
   int nSamples, totalPx, nIons;
   
   fp.open(file, std::fstream::in | std::ios::binary);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: The internal %s file could not be created.\n The peak matrix cannot be saved.\n", file);
     throw std::runtime_error(txt);
   }
   fp.read((char*)&nSamples, sizeof(int)); //samples number
   fp.read((char*)&totalPx,  sizeof(int)); //total pixels in all samples
   fp.read((char*)&nIons,    sizeof(int)); //ions number
   
   colSize=((totalPx+2)*sizeof(float)+sizeof(int));
   pos=3*sizeof(int)+nSamples*sizeof(int)+(column-1)*colSize;
   
   NumericVector ret(totalPx+3);
   float intensity;
   int pxSupport;
   fp.seekg(pos, std::ios_base::beg);
   for(int i=0; i<totalPx+3; i++)
   {
     fp.read((char*)&intensity, sizeof(float));
     ret[i]=intensity;
   }
   fp.close();
   return ret;  
 }

//' @name rGetMassColumFromFile()
 //' @title returns a column information of the peak matrix: 
 //' 
 //' @param file     -> file name with peak matrix (tmpPeakMatrix.bin)
 //' @param mass     -> reference to the desired initial column of the peak matrix (Da).
 //' @param sample   -> just download the pixels from this sample.
 //'                 if sample < 0, all sample coordinates are returned
 //' @return a list:
 //'     intensity: vector of intesities 
 //'          mass: mass associated with the column of peakMatrix.
 //'massResolution: final mass resolution at centroid.
 //' pixelsSupport: number of pixels in column with non-zero intensity. 

 // [[Rcpp::export]]
 List rGetMassColumnFromFile(const char* file, float mass, int sample)
 {
   std::fstream fp;
   std::streampos ionPos, colSize, offset;
   int nSamples, totalPx, nIons, massIndex;
   Common common;
   
   fp.open(file, std::fstream::in | std::ios::binary);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: The internal %s file could not be created.\n The peak matrix cannot be saved.\n", file);
     throw std::runtime_error(txt);
   }

   fp.read((char*)&nSamples, sizeof(int)); //samples number
   fp.read((char*)&totalPx,  sizeof(int)); //total pixels in all samples
   fp.read((char*)&nIons,    sizeof(int)); //ions number
   int pxSample[nSamples];
   fp.read((char*)pxSample,  nSamples*sizeof(int)); //px in each sample
   
   if(sample>nSamples)
   {
     printf("warnimg: max samples in file are %d\n", nSamples);
     fp.close();
     return 0;
   }
   
   colSize=(totalPx+3)*sizeof(float); //size of column
   offset=3*sizeof(int)+nSamples*sizeof(int); //file header
   
   //get the mass vector.
   float* tmpMassVect=0;
   tmpMassVect=new float[nIons];
   int massVectSize;

   for(int i=0; i<nIons; i++)
     {
       ionPos=offset+(i*colSize);
       fp.seekg(ionPos, std::ios_base::beg);  
       fp.read((char*)(tmpMassVect+i), sizeof(float));
     }
   
   if(mass<tmpMassVect[0] || mass>tmpMassVect[nIons-1])
     printf("warning: mass is out of range (%.4f/%.4f)\n", tmpMassVect[0], tmpMassVect[nIons-1]);
   
   //mass index
   massIndex =common.nearestIndex(mass, tmpMassVect, nIons);
   if(tmpMassVect) delete []tmpMassVect;
   
   //get pixels range from sample parameter
   sample--; //to C++
   int initPx=0, finalPx;
   if(sample<0) //all samples
    {
     initPx=0; finalPx=totalPx;
     printf("The %d samples are considered.\n", nSamples);
     }
   else //only one sample
   {
   for(int i=0; i<sample; i++)
     initPx+=pxSample[i];
   finalPx=initPx+pxSample[sample];
   }
   int nPixels=finalPx-initPx;

   //object for R link
   NumericVector intVect(nPixels);
   float matrixMass, massRes, intensity, pxSupport;
   
   ionPos=offset+(massIndex*colSize);
   fp.seekg(ionPos, std::ios_base::beg); //file position
   
   //read column data
   fp.read((char*)&matrixMass, sizeof(float));
   fp.read((char*)&massRes,  sizeof(float));
   fp.read((char*)&pxSupport,  sizeof(float));
   ionPos=offset+(massIndex*colSize)+initPx*(std::streampos)sizeof(float);

   for(int i=0; i<nPixels; i++)
   {
     fp.read((char*)&intensity, sizeof(float));
     intVect[i]=intensity;
   }
   fp.close();
   List ret=List::create(Named("intensity")=intVect, Named("mass")=matrixMass, Named("massResolution")=massRes, 
                               Named("pixelsSupport")=pxSupport);
   return ret;  
 }



/// R METHOD ////////////////////////////////////////////////////////////////////////

//'
 //'  @name peakMatrixR
 //'  @title construct the peak matrix. It requires the prior contribution of Class RawToGaussians.
 //'  
 //'  @param totalPixels: Number of pixels to be evaluated. Sum of contributions from each previously analyzed sample.
 //'  @param params:    specific parameters
 //'   "massResolution": desired mass resolution.
 //' "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered.
 //'      "linkedPeaks": two peaks are considered linked if they are closer than the given standard deviation (by defect=3).   
 //'  @param mzLow:    lower  mass to consider
 //'  @param mzHigh:   higher mass to consider
 //'  @param nThreads: number of threads suggested for parallel processing.
 //'  @return a list:
 //'    peakMatrix: Matrix of peak (centroids) rows = pixels, columns = intensity of each pixel.
 //'          mass: Vector with the masses associated with each column of peakMatrix.
 //'massResolution: Vector with final mass resolution at each centroid.
 //' pixelsSupport: Vector with the number of pixels in each column with non-zero intensity. 
 //'   coordinates: Two columns matrix with pixel coordinates.
 //'  pixelsSample: Vector with number of pixels in each sample analyzed.
 //'     
 // [[Rcpp::export]]
List peakMatrixR(Rcpp::String baseDir, Rcpp::List params, float mzLow, float mzHigh, int nThreads)
{
  NumericVector nv;
//  double *Z=new double[1000000000];
  nv=params["massResolution"];
  float massResolution=nv[0];
  
  //Two peaks are considered linked if they are closer than the given standard deviation.
  nv=params["linkedPeaks"]; 
  float linkedPeaks=nv[0];

  //class for segmentation of the mass axis.
  Segments segments(massResolution, mzLow, mzHigh, linkedPeaks);
  
  //conversion Rcpp::String to char*
  char *fileName=new char[200];
  strcpy(fileName, (char*)baseDir.get_cstring());
  strcat(fileName, (char*)"tmpGaussians.bin");

  //Loading the Gaussians generated by the RawToGaussians class from a temporary file.
  int totalPixels=segments.loadGaussians(fileName);

  nv=params["minPixelsSupport"];
  float pxSupport=nv[0]*totalPixels/100.0;
  
  float linked=linkedPeaks;
  //Segments from mass axis
  int size=segments.getMassRanges(&linked);

  MASS_RANGE massRange; //mass range to evaluate
  massRange.low=mzLow;
  massRange.high=mzHigh;
  
  //class for peak matrix
  PeakMatrix peakMatrix(totalPixels, massResolution, mzLow, mzHigh, pxSupport, segments.m_massRange_p, segments.m_gaussians_p, nThreads);
  peakMatrix.m_massRangeSize=size;
  
  strcpy(fileName, (char*)baseDir.get_cstring());
  //build the peak matrix
  List ret=peakMatrix.massRangeToCentroids(fileName, massRange);
  
  delete [] fileName;
  return ret;
}

//Constructor
//Generates the peak matrix; final step. It requires the prior contribution of classes RawToGaussians and Segments.
//   totalPixels: cumulative value of the pixels of each sample analyzed.
//massResolution: desired mass resolution.
//         mzLow: lower  mass to consider
//        mzHigh: higher mass to consider
//     pxSupport: minimum percentage of pixels that must support an ion for it to be considered.
//   massRange_p: array of structures with information on isolated mass segments.
//   gaussians_p: array of structures with information about Gaussians.
PeakMatrix::PeakMatrix(int totalPixels, float massResolution, float mzLow, float mzHigh, float pxSupport, MASS_RANGE *massRange_p, GAUSS_SP *gaussians_p, int nThreads)
{
  m_totalPixels=totalPixels;
  m_massResolution=massResolution;
  m_mzHigh=mzHigh;
  m_mzLow=mzLow;
  m_pxSupport=pxSupport;
  m_massRange_p=massRange_p;
  
  m_nThreads =thread::hardware_concurrency()-1; //a core is released
  if(nThreads<m_nThreads && nThreads>0)
    m_nThreads =nThreads;
  if(m_nThreads<=0) m_nThreads=1;
  if(m_nThreads>MAX_THREADS) m_nThreads=MAX_THREADS;
  if(m_nThreads> m_totalPixels) m_nThreads=m_totalPixels;
  
  m_gaussians_p=gaussians_p;
  for(int i=0; i<MAX_THREADS; i++)
    m_thread_p[i]=0;
  m_pixelsCoordinates_p=0;
}

//Destructor
PeakMatrix::~PeakMatrix()
{
  return;
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
  if(m_pixelsCoordinates_p) delete[] m_pixelsCoordinates_p;
  
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
List PeakMatrix::massRangeToCentroids(char *baseDir, MASS_RANGE massRangeIn)
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
    m_thread_p[i]=new std::thread(&PeakMatrix::mtSegmentation, this, i);
  
  printf("\tphase 3:\tsegments to centroids(%%): 00 ");
  
  //wait until the conclusion-
  for(int i=0; i<m_nThreads; i++)
  {
    m_thread_p[i]->join();
    delete m_thread_p[i]; m_thread_p[i]=0;
  }
  printf("100\n");
  

  //adequacy of results for R.
  int totalIons=0;
  for(int i=0; i<m_nThreads; i++)
  {
    totalIons+=m_totalIons[i];
  }
  std::fstream fp;
  
  char fileName[200];
  strcpy(fileName, baseDir);
  strcat(fileName, (char*)"tmpPixelsCoordinates.bin");
  int nPx=loadPixelsCoordinates(fileName, m_totalPixels);
  if(nPx!=m_totalPixels){printf("ERROR: pixel consistency error (%d/%d)\n", nPx, m_totalPixels); return 0;}
  
//The peak matrix is saved to tmpPeakMatrix.bin file and no information is returned.
  
  strcpy(fileName, baseDir);
  strcat(fileName, (char*)"tmpPeakMatrix.bin");
  fp.open(fileName, std::fstream::out | std::ios::binary | std::ios::trunc);
    if(!fp.is_open())
    {
      char txt[300];
      sprintf(txt, "Error: The internal file %s could not be created.\n The peak matrix cannot be saved.\n", fileName);
      throw std::runtime_error(txt);
    }
    int col=0;
    fp.write((char*)&m_nSamples, sizeof(int)); //samples number
    fp.write((char*)&col, sizeof(int)); //space for matrix rows
    fp.write((char*)&col, sizeof(int)); //space for matrix cols
    
    for(int i=0; i<m_nSamples; i++) 
    {
      fp.write((char*)&m_pixelsSample[i], sizeof(int)); //size of each sample
    }
    
    ION_ENTRY *localIonEntry_p; //input
    float lastMass=0;
    float tmp;
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
          fp.write((char*)&localIonEntry_p->mass, sizeof(float));
          fp.write((char*)&localIonEntry_p->massResolution, sizeof(float));
          tmp=(float)localIonEntry_p->size;
          fp.write((char*)&tmp, sizeof(float));
          for(int row=0; row<m_totalPixels; row++)             //for each pixel
            fp.write((char*)&localIonEntry_p->set[row], sizeof(float));
          col++;
        }
        if(localIonEntry_p->group->group==0) //the last link contains no info.
          lastMass=localIonEntry_p->mass;    //last mass of the thread
        localIonEntry_p=localIonEntry_p->group; //next item
      }
    }
    fp.seekg(sizeof(int), std::ios_base::beg);
    fp.write((char*)&m_totalPixels, sizeof(int));
    fp.write((char*)&col, sizeof(int));

    fp.close();

  return 0;
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
  
  for(int px=0; px<m_totalPixels; px++)//for all pixels
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
  for(int px=0; px<m_totalPixels; px++)//for all pixels
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
        
        localIonEntry_p->set=new float[m_totalPixels]; //memory for intensities
        for(int i=0; i<m_totalPixels; i++) //reset all intensities. 
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
    
    for(int px=0; px<m_totalPixels; px++)//for all pixels
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
    for(int px=0; px<m_totalPixels; px++)//for all pixels
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

//Loads the file with the coordinates of each pixel
//If there is more than one sample, all its pixels are integrated
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
  int nSamplePixels, nPixels, nPxGauss, coordSize=2*sizeof(int);
  bool hit=true;
  
  m_nSamples=0;
  m_pixelsCoordinates_p=new PIXEL_XY[totalPixels]; //coordinates
  int pxTotal=0;
  while(true)
  {
    fp.read((char*)&nSamplePixels, sizeof(int)); //#pixels into the sample
    if(fp.eof()) break;
    m_pixelsSample[m_nSamples++]=nSamplePixels;
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
  return pxTotal;
}

 
