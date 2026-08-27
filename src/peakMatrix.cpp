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


//' @name rGetCoordinatesFromFile()
//' @title returns a matrix with the coordinates of all pixels (X/Y).
//' If there are multiple samples, they appear sequentially; that is, the matrix has as many rows 
//' as the cumulative number of pixels in each sample and two columns.
//' @param file   -> file name with pixels coordinates (_pixelsCoord.bin)
//' @param sample -> just download the pixels from this sample.
//'                 if sample < 1, all sample coordinates are returned
//' @return a matrix with the coordinates (X/Y) of pixels.
//' 
// [[Rcpp::export]]
NumericMatrix rGetPixelsCoordinates(const char* fileName, int sample)
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


 //' @name rGetCentroidsFromFile()
 //' @title returns info of centroids from the peak matrix file: 
 //' 
 //' @param file -> file name with peak matrix (_peakMatrix.bin)
 //' @return a list:
 //'      mass: vector of centroids
 //'      intensity: vector of average intensities of all pixels associated to centroid
 //'      tolerance: vector of tolerances od aech centroid
 //'      pixelsSupport: vector of pixels that support the centroid.
 //'      pixelsSample: vector of number of pixels into samples
 
 // [[Rcpp::export]]
 List rGetBasic(const char* file)
 {
   std::fstream fp;
   std::streampos pos, colSize;
   int nSamples, totalPx, nIons;
   bool hit=true;
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
   
   NumericVector mass(nIons);
   NumericVector intensity(nIons);
   NumericVector tolerance(nIons);
   IntegerVector pxSupport(nIons);
   IntegerVector pixelsSample(nSamples);
   
   int pxSize;
   for(int i=0; i<nSamples; i++)
   {
     fp.read((char*)&pxSize, sizeof(int));
     pixelsSample[i]=pxSize;
   }
   
   int index=0, tmpPxSupport;
   double tmpMass, tmpIntensity, tmpTolerance;
   std::streampos ionPos=(3+nSamples)*sizeof(int);
   
   while(index<nIons)
   {
     if(fp.eof()) break;
     fp.seekg(ionPos, std::ios_base::beg);
     if(fp.fail() || fp.bad()) //fault control
     {
       char txt[200];
       sprintf(txt, "Error: %s file could not be read completely.", file);
       throw std::runtime_error(txt);
       hit =false; break;
     }
     fp.read((char*)&tmpMass, sizeof(double));
     fp.read((char*)&tmpIntensity, sizeof(double));
     fp.read((char*)&tmpTolerance, sizeof(double));
     fp.read((char*)&tmpPxSupport, sizeof(int));
     mass[index]=tmpMass;
     intensity[index]=tmpIntensity;
     tolerance[index]=tmpTolerance;
     pxSupport[index]=tmpPxSupport;
     index++;
     ionPos=fp.tellg()+(std::streampos)(tmpPxSupport*(sizeof(int)+sizeof(double)));
   }
   fp.close();
   List ret=List::create(Named("mass")=mass, Named("intensity")=intensity,  Named("tolerance")=tolerance, 
                         Named("pixelsSupport")=pxSupport, Named("pixelsSample")=pixelsSample);
   return ret;
 }

 //' @name rGetColumFromFile()
 //' @title returns a column information of the peak matrix: 
 //' 
 //' @param file     -> file name with peak matrix (_peakMatrix.bin)
 //' @param mass     -> reference to the desired initial column of the peak matrix (Da).
 //' @param sample   -> just download the pixels intensity from this sample.
 //'                    if sample is out of limits, pixels from all samples are returned
 //' @return a list:
 //'     intensity: vector of intesities 
 //'         pixel: pixel associated with intensity.
 //'          mass: mass of ion
 //'     tolerance: centroid tolerance (ppm).

 // [[Rcpp::export]]
 List rGetCentroid(const char* file, double mass, int sample)
 {
   std::fstream fp;
   std::streampos ionPos, colSize, offset;
   int nSamples, totalPx, nIons, massIndex;
   Common common;
   double *massAxis_p=0;
   int *colSize_p=0, ion;
   bool allSamples=false;
   
   fp.open(file, std::fstream::in | std::ios::binary);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: The internal %s file could not be created.\n The peak matrix cannot be saved.\n", file);
     throw std::runtime_error(txt);
   }
    //metadata
   fp.read((char*)&nSamples, sizeof(int)); //samples number
   fp.read((char*)&totalPx,  sizeof(int)); //total pixels in all samples
   fp.read((char*)&nIons,    sizeof(int)); //ions number
   int pxSample[nSamples];
   fp.read((char*)pxSample,  nSamples*sizeof(int)); //px in each sample
   
   //px range to load
   int pxLow=0, pxHigh=-1;
   if(sample<0 || sample>=nSamples) //all column pixel
   {
    for(int i=0; i<nSamples; i++)
       pxHigh+=pxSample[i];
   }   
   else
    {
    for(int i=0; i<nSamples; i++)
    {
     pxLow=pxHigh+1;
     pxHigh=pxLow+pxSample[i]-1;
     if(i==sample) break;
    }
    }
   offset=fp.tellg(); //file position

   //mass axis and size of columns
    massAxis_p=new double[nIons];
    colSize_p=new int[nIons];
     double tmpIntensity, tmpTolerance;
     int nPx;
     for(int ion=0; ion<nIons; ion++)
     {
        fp.read((char*)&massAxis_p[ion], sizeof(double));
        fp.read((char*)&tmpIntensity, sizeof(double));
        fp.read((char*)&tmpTolerance, sizeof(double));
        fp.read((char*)&nPx, sizeof(int));
        colSize_p[ion]=3*sizeof(double)+sizeof(int)+nPx*(sizeof(int)+sizeof(double)); //columns size

        fp.seekg((std::streampos)(nPx*(sizeof(int)+sizeof(double))), std::ios_base::cur);
     }
     //nearest index to mass
      massIndex =common.nearestIndex(mass, massAxis_p, nIons); //sort up

     //column position into file
     ionPos=offset; //first column input
     for(int i=0; i<massIndex; i++)
       ionPos+=(std::streampos)(colSize_p[i]); //offset to column
     ionPos+=(std::streampos)(2*sizeof(double));
     fp.seekg(ionPos, std::ios_base::beg); //positioning
     fp.read((char*)&tmpTolerance, sizeof(double)); //size of column
     fp.read((char*)&nPx, sizeof(int)); //size of column
  
    double *intensity_p=0;
    int *pixel_p=0;
    intensity_p=new double[nPx];
    pixel_p=new int[nPx];
    
    int tmpPixel, pxSize=0;
    //data load
    for(int i=0; i<nPx; i++) //data load
     {
      fp.read((char*)&tmpPixel, sizeof(int));
      fp.read((char*)&tmpIntensity, sizeof(double));
      if(tmpPixel>=pxLow && tmpPixel<=pxHigh) //pixel into range?
        {
        pixel_p[pxSize]=tmpPixel-pxLow+1; //+1 to R
        intensity_p[pxSize++]=tmpIntensity;
        }
     }
    //link to R
    NumericVector intensity(pxSize);
    IntegerVector pixel(pxSize);
    for(int i=0; i<pxSize; i++)
    {
      intensity[i]=intensity_p[i];
      pixel[i]=pixel_p[i];
    }
    
    fp.close();
    List ret=List::create(Named("mass")=massAxis_p[massIndex], Named("tolerance")=tmpTolerance, 
                          Named("intensity")=intensity, Named("pixel")=pixel);
    if(massAxis_p)    delete []massAxis_p;
    if(colSize_p)     delete []colSize_p;
    if(pixel_p)       delete []pixel_p;
    if(intensity_p)   delete []intensity_p;
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
List peakMatrixR(Rcpp::String baseDir, Rcpp::List params, double mzLow, double mzHigh, int nPixels, int nSamples)
{
  NumericVector nv;
  nv=params["tolerance"];
  double tolerance=nv[0];

  //minimum percentage of pixels that must support a centroid.
  nv=params["minPixelsSupport"]; 
  double pxSupport=nv[0]*nPixels/100.0;
  
  //conversion Rcpp::String to char*
  char *baseDir2=new char[200];
  char *fileName=new char[200];
  strcpy(baseDir2, (char*)baseDir.get_cstring());
  
  printf("\tphase 2:\tbinning (with max bin tolerance=%6.3f ppm)(%%): 00 ", tolerance);
  
  PeakMatrix peakMatrix(nPixels, 1e6/tolerance, mzLow, mzHigh, pxSupport, nSamples, baseDir2);
  
  //Loading the Gaussians generated by the RawToGaussians class from a temporary file.
  strcpy(fileName, baseDir2);
  strcat(fileName, (char*)"_gaussians.bin");
  int totalPx=peakMatrix.loadGaussians(fileName);
  if(totalPx<0 || totalPx!=nPixels) //Some pixels might be empty
  {
    printf("ERROR: fail of consistence in pixels number\n");
    return -1;
  }
  int nCentroids=peakMatrix.getCentroids();//get centroids a save it to file
  
  //peakMatrix.infoToFile(baseDir2);
  printf("100\n");
  
  printf("total centroids:%d\n", nCentroids);
  if(baseDir2) delete []baseDir2;
  if(fileName) delete []fileName;
 return nCentroids;
}

//Constructor
//Generates the peak matrix; final step. It requires the prior contribution of classes RawToGaussians and Segments.
//   totalPixels: cumulative value of the pixels of each sample analyzed.
//massResolution: desired mass resolution (mz/delta_mz).
//         mzLow: lower  mass to consider
//        mzHigh: higher mass to consider
//     pxSupport: minimum percentage of pixels that must support an ion for it to be considered.
PeakMatrix::PeakMatrix(int totalPixels, double massResolution, double mzLow, double mzHigh, double pxSupport, int nSamples, char *baseDir)
{
  m_totalPixels=totalPixels; //It is also set in the loadGaussians() function.
  m_nSamples=nSamples;
  m_massResolution=massResolution;
  m_mzHigh=mzHigh;
  m_mzLow=mzLow;
  m_baseDir=baseDir;
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
      if(nPxGauss==0) {pxIndex++; continue;}//px without conten

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
//Formato: 
//cabecera: totalSamples|totalPx|totalMass|pxSample 0, 1,...|
//para cada ion: mass|intensity|tolerance|pxSupport|px1-int1, px2-int2,...|
int PeakMatrix::centers(double *mass_p, int *px_p, int *iGauss_p, int size, double tolerance)
{
  if(size<=0) return 0;
  int nSamples=getSamplesPixelNumber(); //info de px en cada muestra
  if(nSamples==-1) return -1;
  
  std::fstream fp;
  char fileName[200];
  strcpy(fileName, m_baseDir);
  strcat(fileName, (char*)"_peakMatrix.bin");
  fp.open(fileName, std::fstream::out | std::ios::binary | std::ios::trunc);
  int tmp=0;
  double kk;
  
  if(!fp.is_open())
  {
    char txt[250];
    sprintf(txt, "Error: The internal file %s could not be created.\n The peak matrix cannot be saved.\n", fileName);
    throw std::runtime_error(txt);
  }
  fp.write((char*)&m_nSamples, sizeof(int)); //samples number
  fp.write((char*)&tmp, sizeof(int)); //space for matrix rows
  fp.write((char*)&tmp, sizeof(int)); //space for matrix cols
  for(int i=0; i<m_nSamples; i++)
    fp.write((char*)&tmp, sizeof(int)); //space for pixel number into each sample

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
  double segmentSize=m_centers_p[0]*tolerance/1e6, centroidTolerance;

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
        double intensity=0, dispersion=0, pxIntensity;
        int iMass;

        for(int j=0; j<m_centersSize_p[k]; j++)
        {
          iMass=massIndex_p[j];
          intensity+=m_gaussians_p[px_p[iMass]].gauss_p[iGauss_p[iMass]].weight; //intensity.
          dispersion+=(mass_p[iMass] - m_centers_p[k])*(mass_p[iMass] - m_centers_p[k]);
        }
        dispersion=sqrt(dispersion/m_centersSize_p[k]);
        intensity/=m_centersSize_p[k];
        
        fp.write((char*)&m_centers_p[k], sizeof(double)); //centroid mass
        fp.write((char*)&intensity, sizeof(double));      //centroid intensity
        centroidTolerance=1e6*dispersion/m_centers_p[k];  
        fp.write((char*)&centroidTolerance, sizeof(double)); //centroid tolerance (ppm)
        fp.write((char*)&m_centersSize_p[k], sizeof(int)); //px number support
        //px intensities
        for(int i=0; i<m_centersSize_p[k]; i++) //copy of intensities
        {
          iMass=massIndex_p[i];
          int px=px_p[iMass]; //pixel of the Gaussian.
          pxIntensity=m_gaussians_p[px_p[iMass]].gauss_p[iGauss_p[iMass]].weight; //px intensities
          fp.write((char*)&px, sizeof(int)); //px 
          fp.write((char*)&pxIntensity, sizeof(double)); //px intensity
        }
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
  fp.seekg(sizeof(int), std::ios_base::beg);
  fp.write((char*)&m_totalPixels, sizeof(int));
  fp.write((char*)&k, sizeof(int));
  fp.write((char*)m_pixelsSample, nSamples*sizeof(int));
  fp.close();
  
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
  strcat(fileName, (char*)"_pixelsCoord.bin");
  int nPx=loadPixelsCoordinates(fileName, m_totalPixels);
  if(nPx!=m_totalPixels){printf("ERROR: pixel consistency error (%d/%d)\n", nPx, m_totalPixels); return 0;}
  
  //The peak matrix is saved to _peakMatrix.bin file and no information is returned.
  
  strcpy(fileName, baseDir);
  strcat(fileName, (char*)"_peakMatrix.bin");
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
  int nSamplePixels, nPixels, nPxGauss, coordSize=2*sizeof(int);
  bool hit=true;
  
  int nSamples=0;
  m_pixelsCoordinates_p=new PIXEL_XY[totalPixels]; //coordinates
  int pxTotal=0;
  while(true)
  {
    fp.read((char*)&nSamplePixels, sizeof(int)); //#pixels into the sample
    if(fp.eof()) break;
    m_pixelsSample[nSamples++]=nSamplePixels;
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
  if(nSamples != m_nSamples)
  {
    printf("ERROR: inconsistency in the number of samples.\n");
    return -1;
  }
  return pxTotal;
}

//establece la cantidad de píxeles en cada muestra a partir de info en _pixelsCoord.bin
//resultado en el vector m_pixelsSample[]
//retorna la cantidad de muestras o -1 si fallo
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
    m_pixelsSample[nSamples++]=nSamplePixels;
    fp.seekg((std::streampos)(nSamplePixels*sizeof(PIXEL_XY)), std::ios_base::cur);
    if(fp.fail() || fp.bad()) //fault control
      {
        char txt[200];
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

