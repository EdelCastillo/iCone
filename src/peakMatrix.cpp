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
    if(fp.fail()) 
    {
      printf("Error while reading in %s", fileName);
      fp.close(); return 0;
    }
    
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
  if(fp.fail()) 
  {
    printf("Error while reading in %s", fileName);
    fp.close(); return 0;
  }
  
  int xy[2];
  
  for(int sample=0; sample<nSamples; sample++) //for all samples
  {
    fp.read((char*)&nPixels, sizeof(int)); //#pixels into the sample
    for(int pxSample=0; pxSample<nPixels; pxSample++) //for each pixel into the sample
      {
      fp.read((char*)xy, 2*sizeof(int)); //XY coordinate
      if(fp.fail()) 
      {
        printf("Error while reading in %s", fileName);
        fp.close(); return 0;
      }
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
     if(fp.fail()) 
     {
       printf("Error while reading in %s", file);
       fp.close(); return 0;
     }
     mass[index]=tmpMass;
     intensity[index]=tmpIntensity;
     tolerance[index]=tmpTolerance;
     pxSupport[index]=tmpPxSupport;
//printf("[%3d]%.4f\n", index, tmpMass);
     index++;
     ionPos=fp.tellg()+(std::streampos)(tmpPxSupport*(sizeof(int)+sizeof(double)));
   }
   fp.close();
   List ret=List::create(Named("mass")=mass, Named("intensity")=intensity,  Named("tolerance")=tolerance, 
                         Named("pixelsSupport")=pxSupport, Named("pixelsSample")=pixelsSample);
   return ret;
 }

 //' @name rGetCentroid()
 //' @title returns a column information of the peak matrix. only non-zero pixels. 
 //' 
 //' @param file     -> file name with peak matrix (_peakMatrix.bin)
 //' @param mass     -> reference to the desired initial column of the peak matrix (Da).
 //' @param sample   -> just download the pixels intensity from this sample.
 //' @param expand   -> if true,  it returns the intensity values for all the pixels in the sample.
 //'                    if false, it returns the intensity values only for the sample pixels with non-zero values.
 //'                       in this case, return two vectors: intensity and pixel
 //'                    if sample is out of limits, pixels from all samples are returned
 //' @return a list:
 //'     intensity: vector of intesities 
 //'         pixel: pixel associated with intensity.
 //'          mass: mass of ion
 //'     tolerance: centroid tolerance (ppm).

 // [[Rcpp::export]]
 List rGetCentroid(const char* file, double mass, int sample, bool expand)
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
   if(fp.fail()) 
   {
     printf("Error while reading in %s", file);
     fp.close(); return 0;
   }
   
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
        if(fp.fail()) 
        {
          printf("Error while reading in %s", file);
          fp.close(); return 0;
        }
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
     if(fp.fail()) 
     {
       printf("Error while reading in %s", file);
       fp.close(); return 0;
     }
     
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
      if(fp.fail()) 
      {
        printf("Error while reading in %s", file);
        fp.close(); return 0;
      }
      if(tmpPixel>=pxLow && tmpPixel<=pxHigh) //pixel into range?
        {
        pixel_p[pxSize]=tmpPixel-pxLow;
        intensity_p[pxSize++]=tmpIntensity;
        }
     }
    List ret;
    //link to R
    if(expand)
    {
      int pxIntoSample=pxSample[sample];

      NumericVector intensity(pxIntoSample);
      for(int i=0; i<pxIntoSample; i++) intensity[i]=0;
      for(int i=0; i<pxSize; i++)
      {
        intensity[pixel_p[i]]=intensity_p[i];
      }
      ret=List::create(Named("mass")=massAxis_p[massIndex], Named("tolerance")=tmpTolerance, 
                            Named("intensity")=intensity);
    }
    else
    {
      NumericVector intensity(pxSize);
      IntegerVector pixel(pxSize);
      for(int i=0; i<pxSize; i++)
      {
        intensity[i]=intensity_p[i];
        pixel[i]=pixel_p[i];
      }
      ret=List::create(Named("mass")=massAxis_p[massIndex], Named("tolerance")=tmpTolerance, 
                            Named("intensity")=intensity, Named("pixel")=pixel+1);  //+1 to R
    }
    fp.close();
    if(massAxis_p)    delete []massAxis_p;
    if(colSize_p)     delete []colSize_p;
    if(pixel_p)       delete []pixel_p;
    if(intensity_p)   delete []intensity_p;
    return ret;
 }


 //' @name rGetMatrix()
 //' @title returns the intensity matrix associated with a given sample. 
 //' 
 //' @param file      -> file name with peak matrix (_peakMatrix.bin)
 //' @param sample    -> just download the matrix intensity from this sample.
 //' @return a matrix -> row = pixels; column=centroids

 // [[Rcpp::export]]
 NumericMatrix rGetMatrix(const char* file, int sample)
 {
   std::fstream fp;
   std::streampos ionPos, colSize, offset;
   int totalSamples, totalPx, totalIons, massIndex;
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
   fp.read((char*)&totalSamples, sizeof(int)); //samples number
   fp.read((char*)&totalPx,  sizeof(int)); //total pixels in all samples
   fp.read((char*)&totalIons,    sizeof(int)); //ions number
   if(fp.fail()) 
   {
     printf("Error while reading in %s", file);
     fp.close(); return 0;
   }
   if(sample<0 || sample>=totalSamples) 
   {
     printf("sample out of limits [%d/%d]\n", 0, totalSamples-1);
     fp.close(); return 0;
   }
   int pxSample[totalSamples];
   fp.read((char*)pxSample,  totalSamples*sizeof(int)); //px in each sample

   //la info en disco se organiza por centriodes.
   //cada centroide está acompañado de una lista de píxeles de todas las muestras en formato pixel-intensidad
   //para separarlos, se deben conocer en rango de píxeles asociado a cada muestra
   
   //se delimitan los pixeles mínimo y máximo en cada muestra
   PIXEL_XY samplesPxLimit[totalSamples]; //x -> low; y -> high
   int pxLow, pxHigh;
   samplesPxLimit[0].x=0;
   samplesPxLimit[0].y=pxSample[0]-1;
   pxLow=pxSample[0];
   for(int i=1; i<totalSamples; i++)
   {
     samplesPxLimit[i].x=samplesPxLimit[i-1].y+1;
     samplesPxLimit[i].y=samplesPxLimit[i].x+pxSample[i]-1;
   }
   //se copia la info a la matriz
   NumericMatrix pkMat(pxSample[sample], totalIons); //filas, columnas
   double tmpMass, tmpIntensity, tmpTolerance;
   int tmpPx, nPx;
   
   for(int ion=0; ion<totalIons; ion++) //para cada centroide
   {
     for(int i=0; i<pxSample[sample]; i++) //puesta cero de la columna
       pkMat(i, ion)=0;
     
     //info de este centroide
     fp.read((char*)&tmpMass, sizeof(double));
     fp.read((char*)&tmpIntensity, sizeof(double));
     fp.read((char*)&tmpTolerance, sizeof(double));
     fp.read((char*)&nPx, sizeof(int));
     if(fp.fail()) 
     {
       printf("Error while reading in %s", file);
       fp.close(); return 0;
     }
     pxLow=samplesPxLimit[sample].x; //rango de píxeles de interés
     pxHigh=samplesPxLimit[sample].y;
     for(int i=0; i< nPx; i++)
     {
       fp.read((char*)&tmpPx, sizeof(int));
       fp.read((char*)&tmpIntensity, sizeof(double));
       if(fp.fail()) 
       {
         printf("Error while reading in %s", file);
         fp.close(); return 0;
       }
       if(tmpPx>=pxLow && tmpPx<=pxHigh)
         pkMat(tmpPx-pxLow, ion)=tmpIntensity;
     }
   }
   fp.close();
   return (pkMat);
 }


/// R METHOD ////////////////////////////////////////////////////////////////////////

//'
 //'  @name peakMatrixR
 //'  @title construct the peak matrix. It requires the prior contribution of Class RawToGaussians.
 //'  
 //'  @param  "baseDir": directory for report files 
 //'  @param     params: specific parameters
 //'        "tolerance": desired mass tolerance for binning.
 //' "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered.
 //'  @param      mzLow: lower  mass to consider
 //'  @param     mzHigh: higher mass to consider
 //'  @param    nPixels: total pixels
 //'  @param  pxSamples: vector containing the number of pixels for each sample
 //'  @param   nSamples: Number of samples
 //'  @return   Number of centroids
 //'     
 // [[Rcpp::export]]
List peakMatrixR(Rcpp::String baseDir, Rcpp::List params, double mzLow, double mzHigh, int nPixels, IntegerVector pxSamples, int nSamples)
{
  NumericVector nv;
  nv=params["tolerance"];
  double tolerance=nv[0];

  //minimum percentage of pixels that must support a centroid.
  nv=params["minPixelsSupport"]; 
  double pxSupport=nv[0]*nPixels/100.0;

  CharacterVector cv;
  cv=params["intMethod"];
  int intMethod=0;
  if(cv[0]=="max")
    intMethod=1;
  else
    intMethod=0; //mean
  
  
  //conversion Rcpp::String to char*
  char *baseDir2=new char[200];
  char *fileName=new char[200];
  strcpy(baseDir2, (char*)baseDir.get_cstring());
 
  int *pxSamples_p=0;
  pxSamples_p=new int[nSamples];
  for(int i=0; i<nSamples; i++)
  {
    pxSamples_p[i]=pxSamples[i];
  }
  
  printf("\tphase 2:   binning (with bin tolerance=%6.2f ppm)(%%): 00 ", tolerance);
  
  PeakMatrix peakMatrix(nPixels, 1e6/tolerance, mzLow, mzHigh, pxSupport, pxSamples_p, nSamples, baseDir2, intMethod);

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
  
  if(baseDir2) delete []baseDir2;
  if(fileName) delete []fileName;
  if(pxSamples_p) delete []pxSamples_p;
  return nCentroids;
}

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
  PeakMatrix::PeakMatrix(int totalPixels, double massResolution, double mzLow, double mzHigh, double pxSupport, int *pxSamples_p, int nSamples, char *baseDir, int intMethod)
{
  m_totalPixels=totalPixels; //It is also set in the loadGaussians() function.
  m_nSamples=nSamples;
  m_pxSamples_p=pxSamples_p;
  m_massResolution=massResolution;
  m_mzHigh=mzHigh;
  m_mzLow=mzLow;
  m_baseDir=baseDir;
  m_gaussians_p=0;
  m_massRange_p=0;
  m_centers_p=0;
  m_centersSize_p=0;
  m_pxSupport=pxSupport;
  m_pixelsCoordinates_p=0;
  m_intMethod=intMethod;
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
 
  //printf("PeakMatrix destructor finish\n");
}

//generate de peak matrix with the centroids, their tolerance, and the number of support pixels.  
int PeakMatrix::getCentroids()
{
  int nCenters=0; //número de mzMax(gaussianas) sobre toda la muestra
  
  for(int px=0; px<m_totalPixels; px++)//for all pixels
  {
    if(m_gaussians_p[px].gauss_p!=0) 
      nCenters+=m_gaussians_p[px].size;
  }
  
  double *mass_p=new double[nCenters];
  int *px_p=new int[nCenters];
  double *intensity_p=new double[nCenters];
//  for(int i=0; i<nCenters; i++) intensity_p[i]=0;
  
  //gaussians to vectors
  nCenters=getCentroidsIntoRange(mass_p, px_p, intensity_p, nCenters);
  
  double deltaMass=m_mzLow/(m_massResolution); //delta=1/2 of the minimum mass increment of the spectrometer

  int massAxisSize; 
  massAxisSize=nCenters;
  
  double tolerance=1e6/m_massResolution;
  
  m_centers_p=new double[massAxisSize];
  m_centersSize_p=new int[massAxisSize];
  
  //get all centroids (binning)
  m_nCentroids=centers(mass_p, px_p, intensity_p, nCenters, tolerance);
  
  if(mass_p)      delete [] mass_p;
  if(px_p)        delete [] px_p;
  if(intensity_p) delete [] intensity_p;
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
//     mass_p: vector masses.
//       px_p: vector of pixels.
//intensity_p: vector of intensities
//       size: vectors length
// tolerance: maximum bin size (ppm).
// Returns the number of centers detected (length of the centers_p array).
//Formato: 
//cabecera: totalSamples|totalPx|totalMass|pxSample 0, 1,...|
//para cada ion: mass|intensity|tolerance|pxSupport|px1-int1, px2-int2,...|
int PeakMatrix::centers(double *mass_p, int *px_p, double *intensity_p, int size, double tolerance)
{ 
  if(size<=0) return 0;
  if(m_nSamples<=0) return -1;
  
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
  fp.write((char*)&m_totalPixels, sizeof(int)); //space for matrix rows
  fp.write((char*)&tmp, sizeof(int)); //space for matrix cols
  for(int i=0; i<m_nSamples; i++)
  {
    fp.write((char*)&m_pxSamples_p[i], sizeof(int)); //space for pixel number into each sample
  }

  int *sortedIndex_p=0;
  sortedIndex_p=new int[size];
  int *massIndex_p=new int[size];
  
  int *tmpPx_p=new int[size];
  int *tmpPxIndex_p=new int[size];

  Common common;
  common.sortUp(mass_p, sortedIndex_p, size); //increasing ordering of masses.

  m_centers_p[0]=mass_p[sortedIndex_p[0]];

  m_centersSize_p[0]=1;
  massIndex_p[0]=sortedIndex_p[0];
  double segmentSize=m_centers_p[0]*tolerance/1e6;
  double centroidTolerance;

  int nIons=0, ionSize=0, mi=0;
  double massDiff;
  int totalPixels=0;
  int repesCount=0, tmpCount=0;
  bool repes=false;
  std::streampos initPos, finalPos;
  int vez=1;
  
  for(int i=1; i<size; i++) //iterative averaging.
  {
    //avance
    if((double)i/(double)size>vez*0.1) {if(vez<10) printf("%d ", vez*10); vez++;}
    
    massDiff=fabs(mass_p[sortedIndex_p[i]]-m_centers_p[nIons]);
    if(massDiff<segmentSize)
    {
      m_centers_p[nIons]=(1.0/((double)m_centersSize_p[nIons]+1.0))*((double)m_centersSize_p[nIons]*m_centers_p[nIons]+mass_p[sortedIndex_p[i]]);
      massIndex_p[mi]=sortedIndex_p[i];
      m_centersSize_p[nIons]++; //masas en el bin
      mi++;
    }
    else //new bin
    {
      if(m_centersSize_p[nIons]>=m_pxSupport)     
      {
        //The results for the ion are noted:
        //Analysis of centroid average intensity and dispersion.
        double intensity=0, dispersion=0, pxIntensity;
        int iMass;
        for(int j=0; j<m_centersSize_p[nIons]; j++)
        {
          iMass=massIndex_p[j];
          dispersion+=(mass_p[iMass] - m_centers_p[nIons])*(mass_p[iMass] - m_centers_p[nIons]);
        }
        dispersion=sqrt(dispersion/m_centersSize_p[nIons]);
        
        fp.write((char*)&m_centers_p[nIons], sizeof(double)); //centroid mass
        initPos=fp.tellg();
        fp.write((char*)&intensity, sizeof(double));      //centroid intensity
        centroidTolerance=1e6*dispersion/m_centers_p[nIons];  
        fp.write((char*)&centroidTolerance, sizeof(double)); //centroid tolerance (ppm)
        fp.write((char*)&m_centersSize_p[nIons], sizeof(int)); //px number support
        
        //px intensities
        for(int k=0; k<m_centersSize_p[nIons]; k++) //copy of intensities
        {
          iMass=massIndex_p[k];
          int px=px_p[iMass]; //pixel of the Gaussian.
          tmpPx_p[k]=px;
        }
        //Puede que dentro de un mismo bin exista más de una gaussiana para el mismo pixel
        //se genera una intensidad representativa
        common.sortUpI(tmpPx_p, tmpPxIndex_p, m_centersSize_p[nIons]); //increasing ordering of pixeles.
        bool into=false;
        int first=-1, last=-1;
        for(int k=0; k<m_centersSize_p[nIons]-1; k++) //para todas las masas del bin
        {
          int px=tmpPx_p[tmpPxIndex_p[k]]; //px asociado a 
          if(px==tmpPx_p[tmpPxIndex_p[k+1]]&& into==false) 
            {into=true; first=k; last=-1;}
          if(into==true && px!=tmpPx_p[tmpPxIndex_p[k+1]])
            {into=false; last=k;}
          if(into==true && k+2==m_centersSize_p[nIons])
            {into=false; last=k+1;}

         if(into==false)
          {
           double pxIntensity=0;
           if(first!=-1 && last!=-1) //px coinciden en la secuencia
              {
             repesCount+=last-first+1;
              int px=tmpPx_p[tmpPxIndex_p[k]];
              if(m_intMethod==MEAN) //mean  
              {
                for(int j=first; j<=last; j++)
                {
                  int pxIdx_a=tmpPxIndex_p[j];
                  int pxIdx_b=massIndex_p[pxIdx_a];
                  pxIntensity+=intensity_p[pxIdx_b];
                }
                             
                pxIntensity/=(last-first+1);
              }
              else if(m_intMethod==MAX) //max value
              {
                double intValue;
                pxIntensity=0;
                for(int j=first; j<=last; j++)
                {
                  int pxIdx_a=tmpPxIndex_p[j];
                  int pxIdx_b=massIndex_p[pxIdx_a];
                  intValue=intensity_p[pxIdx_b];
                  if(intValue>pxIntensity) pxIntensity=intValue; 
                }
              }
              //intensidad representativa del grupo
              fp.write((char*)&px, sizeof(int)); //px 
              fp.write((char*)&pxIntensity, sizeof(double)); //px intensity
              first=-1; last=-1;
              ionSize++;
            }
           
            else //px difieren
            {
              px=tmpPx_p[k];
              pxIntensity=intensity_p[massIndex_p[k]];
              fp.write((char*)&px, sizeof(int)); //px 
              fp.write((char*)&pxIntensity, sizeof(double)); //px intensity
              ionSize++;
            }
          }
          
        }
          //save info to disk
          finalPos=fp.tellg();
          if(m_intMethod==MAX)
          {
            intensity=0;
            for(int i=0; i<ionSize; i++)
              if(intensity_p[i]>intensity) intensity=intensity_p[i];
            fp.seekg(initPos, std::ios_base::beg);
            fp.write((char*)&intensity, sizeof(double)); //px in centroid
            fp.seekg(sizeof(double), std::ios_base::cur); //saltamos la tolerancia
          }
          else
          //if(intMethos==MEAN)
          {
            intensity=0;
            for(int i=0; i<ionSize; i++)
                intensity+=intensity_p[i];
            if(ionSize>0)
              intensity/=(double)ionSize;
            else intensity=0;
            fp.seekg(initPos, std::ios_base::beg);
            fp.write((char*)&intensity, sizeof(double)); //px in centroid
            fp.seekg(sizeof(double), std::ios_base::cur); //saltamos la tolerancia
          }
          totalPixels+=ionSize;
          fp.write((char*)&ionSize, sizeof(int)); //px in centroid
          fp.seekg(finalPos, std::ios_base::beg);
          ionSize=0;
          nIons++;
      }
      //new centroid
      mi=1;
      segmentSize=mass_p[sortedIndex_p[i]]*tolerance/1e6;
      m_centers_p[nIons]=mass_p[sortedIndex_p[i]];
      m_centersSize_p[nIons]=1;
      massIndex_p[0]=sortedIndex_p[i];
    }
           
  }
  printf("100\n");
  printf("\t\t\ttotal centroids:%d\n", nIons);
  printf("\t\t\ttotal number of Gaussians in the same pixel within a bin=%d; \n", repesCount);
  fp.seekg(sizeof(int), std::ios_base::beg);
  fp.write((char*)&m_totalPixels, sizeof(int));
  fp.write((char*)&nIons, sizeof(int));
  fp.close();
  
  if(sortedIndex_p) delete [] sortedIndex_p;
  if(massIndex_p)   delete [] massIndex_p;
  if(tmpPx_p)       delete [] tmpPx_p;
  if(tmpPxIndex_p)  delete [] tmpPxIndex_p;
  m_nIons=nIons;
  
  return nIons;
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
  int nSamplePixels, nPixels, nPxGauss, coordSize=2*sizeof(int);
  bool hit=true;
  
//  int nSamples=0;
  m_pixelsCoordinates_p=new PIXEL_XY[totalPixels]; //coordinates
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
    m_pxSamples_p[nSamples++]=nSamplePixels;
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

