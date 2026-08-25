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
//'                 if sample < 1, all sample coordinates are returned
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
   
   colSize=((totalPx+2)*sizeof(double)+sizeof(int));
   offset=3*sizeof(int)+nSamples*sizeof(int);

   //get matrix column from mass parameters
   NumericVector massVect(nIons);
   double tmpMass=0;
   for(int i=0; i<nIons; i++)
     {
       ionPos=offset+(i*colSize);
       fp.seekg(ionPos, std::ios_base::beg);  
       fp.read((char*)&tmpMass, sizeof(double));
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
   if(column<1)
   {printf("warning: 'column' must be an integer value greater than zero.\n"); return 0;}
   
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
   
   colSize=((totalPx+2)*sizeof(double)+sizeof(int));
   pos=3*sizeof(int)+nSamples*sizeof(int)+(column-1)*colSize;
   
   NumericVector ret(totalPx+3);
   double intensity;
   int pxSupport;
   fp.seekg(pos, std::ios_base::beg);
   fp.read((char*)&intensity, sizeof(double)); //mass
   ret[0]=intensity;
   fp.read((char*)&intensity, sizeof(double)); //massRes
   ret[1]=intensity;
   fp.read((char*)&pxSupport, sizeof(int)); //nPx
   ret[2]=pxSupport;
   for(int i=3; i<totalPx+3; i++)
   {
     fp.read((char*)&intensity, sizeof(double));
     ret[i]=intensity;
   }
   fp.close();
   return ret;  
 }

//' @name rGetColumFromFile()
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
 List rGetColumnFromFile(const char* file, double mass, int sample)
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
   
   colSize=(totalPx+2)*sizeof(double)+sizeof(int); //size of column
   offset=3*sizeof(int)+nSamples*sizeof(int); //file header
   
   //get the mass vector.
   double* tmpMassVect=0;
   tmpMassVect=new double[nIons];
   int massVectSize;

   for(int i=0; i<nIons; i++)
     {
       ionPos=offset+(i*colSize);
       fp.seekg(ionPos, std::ios_base::beg);  
       fp.read((char*)(tmpMassVect+i), sizeof(double));
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
   double matrixMass, massRes, intensity;
   int pxSupport;
   
   ionPos=offset+(massIndex*colSize);
   fp.seekg(ionPos, std::ios_base::beg); //file position
   
   //read column data
   fp.read((char*)&matrixMass, sizeof(double));
   fp.read((char*)&massRes,  sizeof(double));
   fp.read((char*)&pxSupport,  sizeof(int));
   ionPos=offset+(massIndex*colSize)+initPx*(std::streampos)sizeof(double);

   for(int i=0; i<nPixels; i++)
   {
     fp.read((char*)&intensity, sizeof(double));
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
List peakMatrixR(Rcpp::String baseDir, Rcpp::List params, double mzLow, double mzHigh, int nPixels)
{
  NumericVector nv;
  nv=params["tolerance"];
  double tolerance=nv[0];

  //minimum percentage of pixels that must support a centroid.
  nv=params["minPixelsSupport"]; 
  double pxSupport=nv[0]*nPixels/100.0;
  
  printf("\tphase 2:\tbinning (with max bin tolerance=%.3f ppm)(%%): 00 ", tolerance);
  
  PeakMatrix peakMatrix(nPixels, 1e6/tolerance, mzLow, mzHigh, pxSupport);
  
  //conversion Rcpp::String to char*
  char *fileName=new char[200];
  strcpy(fileName, (char*)baseDir.get_cstring());
  strcat(fileName, (char*)"tmpGaussians.bin");

  //Loading the Gaussians generated by the RawToGaussians class from a temporary file.
  int totalPixels=peakMatrix.loadGaussians(fileName);
  
  int nCentroids=peakMatrix.getCentroids();
  
  peakMatrix.infoToFile((char*)baseDir.get_cstring());
  printf("100\n");
  
  printf("total centroids:%d\n", nCentroids);
 return 0;
}

//Constructor
//Generates the peak matrix; final step. It requires the prior contribution of classes RawToGaussians and Segments.
//   totalPixels: cumulative value of the pixels of each sample analyzed.
//massResolution: desired mass resolution (mz/delta_mz).
//         mzLow: lower  mass to consider
//        mzHigh: higher mass to consider
//     pxSupport: minimum percentage of pixels that must support an ion for it to be considered.
PeakMatrix::PeakMatrix(int totalPixels, double massResolution, double mzLow, double mzHigh, double pxSupport)
{
  m_totalPixels=0; //It is set in the loadGaussians() function.
  m_massResolution=massResolution;
  m_mzHigh=mzHigh;
  m_mzLow=mzLow;
  m_gaussians_p=0;
  m_massRange_p=0;
  m_centers_p=0;
  m_centersSize_p=0;
  m_pxSupport=pxSupport;
  m_ionEntry.set=0;
  m_ionEntry.size=0;
  m_ionEntry.group=0;
  m_pixelsCoordinates_p=0;
}

//Destructor
PeakMatrix::~PeakMatrix()
{
  //printf("PeakMatrix destructor init\n");
  if(m_centers_p) delete []m_centers_p;
  if(m_centersSize_p) delete []m_centersSize_p;
 
  if(m_gaussians_p)
  {
    for(int i=0; i<m_totalPixels; i++)
    {
      if(m_gaussians_p[i].gauss_p) delete [] m_gaussians_p[i].gauss_p;
    }
    delete []m_gaussians_p;
  }

  if(m_pixelsCoordinates_p) delete[] m_pixelsCoordinates_p;
 
  //the memory reserved for the ion chain is released.
  if(m_ionEntry.group)
  {
    ION_ENTRY *ionEntry2_p, *ionEntry_p=m_ionEntry.group;
    while(ionEntry_p)
    {
      ionEntry2_p=ionEntry_p->group;
      if(ionEntry_p->set) {delete[] ionEntry_p->set;}
      if(ionEntry_p)      {delete ionEntry_p;}
      
      ionEntry_p=ionEntry2_p;
    }
    if(m_ionEntry.set) delete [] m_ionEntry.set;
  }
  //printf("PeakMatrix destructor finish\n");
}

//generate de peak matrix with the centroids, their tolerance, and the number of support pixels.  
int PeakMatrix::getCentroids()
{
  int nCenters=0;
  
  for(int px=0; px<m_totalPixels; px++)//for all pixels
  {
    if(m_gaussians_p[px].gauss_p!=0) 
      nCenters+=m_gaussians_p[px].size;
  }
  
  double *mass_p=new double[nCenters];
  int *px_p=new int[nCenters];
  int *iGauss_p=new int[nCenters];
  getCentroidsIntoRange(mass_p, px_p, iGauss_p, nCenters);
  
  double deltaMass=m_mzLow/(m_massResolution); //delta=1/2 of the minimum mass increment of the spectrometer
  int massAxisSize; 
  massAxisSize=nCenters;
  
  double tolerance=1e6/m_massResolution;
  
  m_centers_p=new double[massAxisSize];
  m_centersSize_p=new int[massAxisSize];
  
  m_nCentroids=centers(mass_p, px_p, iGauss_p, nCenters, tolerance);
  
  if(px_p) delete [] px_p;
  if(iGauss_p) delete [] iGauss_p;
  delete mass_p;
  return m_nCentroids;
}


//Load the file with information about the Gaussian curves associated with each pixel
//If there is more than one sample, all its pixels are integrated
//fileName: temporary file name (".../tmpGaussians.bin")
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
  int nSamplePixels, nPixels, nPxGauss;
  int gaussSize=3*sizeof(double);
  bool hit=true;
  
  int pxTotal=0, pxIndex=0, nSamples=0;
  int pixelsSample[MAX_SAMPLES];
  while(true) //first reading to obtain information.
  {
    fp.read((char*)&nSamplePixels, sizeof(int)); //#pixels into the sample
    if(fp.eof()) break;
    pixelsSample[nSamples++]=nSamplePixels;
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
  
  m_gaussians_p=new GAUSS_SP[pxTotal]; //array of structs
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
      if(nPxGauss==0) {pxIndex++; continue;}

      m_gaussians_p[pxIndex].gauss_p=new GAUSS_PARAMS[nPxGauss]; //memory
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

//Determines the centers of groups of values whose distance does not exceed a given tolerance.
//mass_p: vector masses.
//px_p: vector of pixels.
//iGauss_p: vector gaussians index
//tolerance: maximum bin size.
// Returns the number of centers detected (length of the centers_p array).
int PeakMatrix::centers(double *mass_p, int *px_p, int *iGauss_p, int size, double tolerance)
{
  if(size<=0) return 0;
  //tolerance*=2.0;
  ION_ENTRY *localIonEntry_p=&m_ionEntry;
  localIonEntry_p->set=0;
  localIonEntry_p->group=0;
  
  int *sortedIndex_p=0;
  sortedIndex_p=new int[size];
  int *massIndex_p=new int[size];
  Common common;
  common.sortUp(mass_p, sortedIndex_p, size); //increasing ordering of masses.
  
  m_centers_p[0]=mass_p[sortedIndex_p[0]];
  m_centersSize_p[0]=1;
  massIndex_p[0]=sortedIndex_p[0];
  double segmentSize=m_centers_p[0]*tolerance/1e6;

  int k=0, mi=0;
  double massDiff;
  int count=0;
  for(int i=1; i<size; i++) //iterative averaging.
  {
    massDiff=fabs(mass_p[sortedIndex_p[i]]-m_centers_p[k]);
    if(massDiff<segmentSize)
    {
      m_centers_p[k]=(1.0/((double)m_centersSize_p[k]+1.0))*((double)m_centersSize_p[k]*m_centers_p[k]+mass_p[sortedIndex_p[i]]);
      massIndex_p[mi]=sortedIndex_p[i];
      m_centersSize_p[k]++;
      mi++;
    }
    else //new bin
    {
      if(m_centersSize_p[k]>=m_pxSupport)     
      {
        //The results for the ion are noted:
        //Analysis of centroid average intensity and dispersion.
        double intensity=0, dispersion=0;
        int iMass;

        for(int j=0; j<m_centersSize_p[k]; j++)
        {
          iMass=massIndex_p[j];
          intensity+=m_gaussians_p[px_p[iMass]].gauss_p[iGauss_p[iMass]].weight; //intensity.
          dispersion+=(mass_p[iMass] - m_centers_p[k])*(mass_p[iMass] - m_centers_p[k]);
        }
        dispersion=sqrt(dispersion/m_centersSize_p[k]);
        intensity/=m_centersSize_p[k];
        
        localIonEntry_p->set=new double[m_totalPixels]; //memory for px intensities
        for(int i=0; i<m_totalPixels; i++) //reset all intensities. 
          localIonEntry_p->set[i]=0.0;
        
        //centroids & tolerance (ppm)
        localIonEntry_p->mass=m_centers_p[k]; //centroide
        localIonEntry_p->size=m_centersSize_p[k]; //masses into centroid
        localIonEntry_p->massResolution=1e6*dispersion/m_centers_p[k]; //resolution (ppm)
        
        //px intensities
        for(int i=0; i<m_centersSize_p[k]; i++) //copy of intensities
        {
          iMass=massIndex_p[i];
          int px=px_p[iMass]; //pixel of the Gaussian.
          localIonEntry_p->set[px]=m_gaussians_p[px_p[iMass]].gauss_p[iGauss_p[iMass]].weight; //px intensities
        }
        
        localIonEntry_p->group=new ION_ENTRY; //new entry for the ion.
        localIonEntry_p=localIonEntry_p->group; //the pointer is updated
        localIonEntry_p->group=0; //initializes the next element.
        localIonEntry_p->set=0;
        localIonEntry_p->size=0;
      
        k++;
      }
      //new centroid
      mi=1;
      segmentSize=mass_p[sortedIndex_p[i]]*tolerance/1e6;
      m_centers_p[k]=mass_p[sortedIndex_p[i]];
      m_centersSize_p[k]=1;
      massIndex_p[0]=sortedIndex_p[i];
    }
  }
  if(sortedIndex_p) delete [] sortedIndex_p;
  if(massIndex_p) delete [] massIndex_p;
  m_nIons=k;
  
  return k;
}


//getCentroidsIntoRange()
//Extracts the existing Gaussians within a mass range from the information in m_gaussians_p.
//massRange: Mass range from which to extract the Gaussians.
//mass_p: vector of internally generated masses.
//px_p: vector of internally generated pixels.
//iGauss_p: vector of internally generated gaussians index
//massSize: size of reserved memory for each vector.
//Returns the number of masses.
int PeakMatrix::getCentroidsIntoRange(double *mass_p, int *px_p, int *iGauss_p, int massSize)
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
        px_p[count]=px;
        iGauss_p[count++]=i;
      }
    }
  }
  return count;
}

//generate de peak matrix with the centroids, their tolerance, and the number of support pixels.  
int PeakMatrix::infoToFile(char *baseDir)
{
  //adequacy of results for R.
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
  double lastMass=0;
  int tmp;
  localIonEntry_p=&m_ionEntry;  //entry to the info.
  while(localIonEntry_p->group) //as long as the next link is not zero.
  {
    //It may happen that the initial masses of a thread are lower than the last masses of 
    //the previous thread, and they should be discarded. 
    //This is a consequence of dealing with joined Gaussians.     
    if(localIonEntry_p->mass>lastMass)                 //the masses are ordered.
    {
      fp.write((char*)&localIonEntry_p->mass, sizeof(double));
      fp.write((char*)&localIonEntry_p->massResolution, sizeof(double));
      tmp=localIonEntry_p->size;
      fp.write((char*)&tmp, sizeof(int));
      for(int row=0; row<m_totalPixels; row++)             //for each pixel
        fp.write((char*)&localIonEntry_p->set[row], sizeof(double));
      col++;
    }
    if(localIonEntry_p->group->group==0) //the last link contains no info.
      lastMass=localIonEntry_p->mass;    //last mass of the thread
    localIonEntry_p=localIonEntry_p->group; //next item
  }
  
  fp.seekg(sizeof(int), std::ios_base::beg);
  fp.write((char*)&m_totalPixels, sizeof(int));
  fp.write((char*)&col, sizeof(int));
  
  fp.close();
  return 0;
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


