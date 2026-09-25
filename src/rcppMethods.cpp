/*********************************************************************************
 *     rcppMethods.cpp
 *     
 *     iCone - R package for MSI data processing
 *     Copyright (C) 2026 Esteban del Castillo Pérez (esteban.delcastillo@urv.cat)
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
#include "peakMatrix.h"

extern int  gPeakCount, gSpectra;

/// R METHODS ////////////////////////////////////////////////////////////////////////

 //'  @name rawToGaussiansR
 //'  @title converts the peaks of each spectrum in the imzML file into Gaussian waves.
 //'  
 //'  @param ibdFname:  absolute reference to the file with the ibd extension.
 //'  @param imzML:     list with information extracted from the imzML file with import_imzML()
 //'  @param params:    specific parameters
 //'              "SNR": signal-to-noise ratio
 //'          tolerance: desired tolerance for the centroids (ppm)
 //'      "noiseMethod": method for estimating noise.
 //' "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered.
 //'  @param mzLow:    lower mass to consider
 //'  @param mzHigh:   higher mass to consider
 //'  @param pxList:   list of pixels. First pixel=1. By default everyone.
 //'  @param nThreads: number of threads suggested for parallel processing.
 //'  @return number of pixels
 //'     
 // [[Rcpp::export]]
 int rawToGaussiansR(Rcpp::String baseDir, const char* ibdFname, Rcpp::List imzML, Rcpp::List params, double mzLow, double mzHigh, Rcpp::NumericVector pxList, int nThreads)
 {
   gPeakCount=0, gSpectra=0;
   char *directory=new char[200];
   strcpy(directory, (char*)baseDir.get_cstring());
   
   RawToGaussians pMatrix(directory, ibdFname, imzML, params, pxList, mzLow, mzHigh, nThreads);
   if(pMatrix.m_hit==false) return 0;
   
   //Phase 1:
   //loads data from a file and converts its peak into Gaussians.
   int ret1=pMatrix.rawToGaussians(); //parallel processing
   if(ret1<0) 
   {return 0;}
   if(gPeakCount<=0) 
   {printf("No peaks are detected in the sample.\n"); return 0;}
   
   //conversion Rcpp::String to char*
   char *fileName=new char[200];
   strcpy(fileName, (char*)baseDir.get_cstring());
   strcat(fileName, (char*)"_gaussians.bin");
   
   int nPeaks=pMatrix.saveGaussians(fileName, pMatrix.m_gaussians_p);
   printf("\t\t\ttotal peaks in the samples: %d\n", nPeaks);
   delete []fileName;
   return pMatrix.m_NPixels;
 }
 
 //'
 //'  @name peakMatrixR
 //'  @title construct the peak matrix into a file. It requires the prior contribution of Class RawToGaussians.
 //'  
 //'  @param  "baseDir": directory for report files 
 //'  @param     params: specific parameters
 //'        "tolerance": desired mass tolerance for binning.
 //'         "nThreads": number of threads for parallel processing
 //'       "intMethods": intensity values for the binning stage: mean (defect), max
 //' "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered.
 //'  @param      mzLow: lower  mass to consider
 //'  @param     mzHigh: higher mass to consider
 //'  @param    nPixels: total pixels
 //'  @param  pxSamples: vector containing the number of pixels for each sample
 //'  @param   nSamples: Number of samples
 //'  @return   Number of centroids
 //'     
 // [[Rcpp::export]]
 int peakMatrixR(Rcpp::String baseDir, Rcpp::List params, double mzLow, double mzHigh, int nPixels, IntegerVector pxSamples, int nSamples)
 {
   NumericVector nv;
   nv=params["tolerance"];
   double tolerance=nv[0];
   
   //minimum percentage of pixels that must support a centroid.
   nv=params["minPixelsSupport"]; 
   double pxSupport=nv[0]*nPixels/100.0;
   
   nv=params["nThreads"];
   int tmpThreads=(int)nv[0];
   //parameter control.
   
   int nThreads =thread::hardware_concurrency()-1; //a core is released
   if(tmpThreads<nThreads && tmpThreads>0)
     nThreads =tmpThreads;
   if(nThreads<=0) nThreads=1;
   if(nThreads>MAX_THREADS) nThreads=MAX_THREADS;
   if(nThreads> nPixels) nThreads=nPixels;
   
   CharacterVector cv;
   cv=params["intMethod"];
   int intMethod=0;
   if(cv[0]=="max")
     intMethod=1;
   else
     intMethod=0; //mean
   
   
   //conversion Rcpp::String to char*
   int *pxSamples_p=0;
   char *baseDir2=0;
   char *fileName=0;
   try{
     baseDir2=new char[200];
     fileName=new char[200];
     strcpy(baseDir2, (char*)baseDir.get_cstring());
     
     pxSamples_p=new int[nSamples];
   }
   catch(const std::bad_alloc& e)
   {
     printf("Error reserving memory: %s\n",e.what());
     return 0;
   }
   
   for(int i=0; i<nSamples; i++)
   {
     pxSamples_p[i]=pxSamples[i];
   }
   
   printf("\tphase 2:   from peaks to centroids(%%): 00 ");
   
   //constructor
   PeakMatrix peakMatrix(nPixels, tolerance, mzLow, mzHigh, pxSupport, pxSamples_p, nSamples, baseDir2, intMethod, nThreads);
   
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
   
   if(baseDir2)    delete []baseDir2;
   if(fileName)    delete []fileName;
   if(pxSamples_p) delete []pxSamples_p;
   return nCentroids;
   //  return ret;
 }
 
 
 //information capture
 //'  @name rGetRawBasic
 //'  @title returns basic information about dimensions from the imzML file.
 //'  
 //'  @param ibdFname:  absolute reference to the file with the ibd extension.
 //'  @param imzML:     list with information extracted from the imzML file with import_imzML()
 //'  @param pxList:    list of pixels. First pixel=1. By default everyone.
 //'  @return list:     minPixel, maxPixel, minMz, maxMz
 //'     
 // [[Rcpp::export]]
 List rGetRawBasic(const char* ibdFname, Rcpp::List imzML, Rcpp::NumericVector pxList)
 {
   //class for accessing imzML files.
   GetImzMLData *getImzMLData_p=0;
   double *rawMzData=0;
   Rcpp::DataFrame df;
   df=imzML["run"];
   bool continuous=imzML["continuous_mode"], hit;
   int nPixels, minPx, maxPx, maxMzLength=0, minMzLength=0x7FFFFFFF, massSize, px;
   double minMz=1e32, maxMz=0;
   
   NumericVector mzLength=df["mzLength"];
   
   if(pxList.size()==1 && pxList[0]==-1) //by defect, all pixels
   {
     hit=false;
     nPixels=df.nrows();
     minPx=1; maxPx=nPixels;
     for(int px=0; px<nPixels; px++)
     {
       if(mzLength[px]>maxMzLength) maxMzLength=mzLength[px]; //maximum spectrum length
       if(mzLength[px]<minMzLength) minMzLength=mzLength[px]; //minimum spectrum length
     }
     try{
       rawMzData=new double[maxMzLength]; 
       getImzMLData_p= new GetImzMLData(ibdFname, imzML);   
     }
     catch(const std::bad_alloc& e)
     {
       printf("Error reserving memory: %s\n",e.what());
       return 0;
     }
     
     for(int px=0; px<nPixels; px++)
     {
       massSize=getImzMLData_p->getPixelMassF(px, rawMzData);//mass vector
       if(massSize<=0)
       {
         printf("ERROR reading from the %s file\n", ibdFname);
         hit=true;
         break;
       }
       massSize--;
       if(rawMzData[0]<minMz)        minMz=rawMzData[0];
       if(rawMzData[massSize]>maxMz) maxMz=rawMzData[massSize];
     }
     
   }
   else //only pixels passed
   {
     hit=false; //all OK
     nPixels=pxList.size();
     minPx=0x7FFFFFFF, maxPx=0;
     
     for(int i=0; i<nPixels; i++)
     {
       if(pxList[i]<minPx) minPx=pxList[i];
       if(pxList[i]>maxPx) maxPx=pxList[i];
     }
     if(minPx<0 || minPx>maxPx)
     {
       printf("ERROR: some pixels are out of range(%d/%d)\n", minPx, maxPx);
       hit=true;
       return true; //warning
     }
     
     for(int i=0; i<nPixels; i++)
     {
       px=pxList[i];
       if(mzLength[px]>maxMzLength) maxMzLength=mzLength[px]; //maximum spectrum length
       if(mzLength[px]<minMzLength) minMzLength=mzLength[px]; //minimum spectrum length
     }
     try{
       rawMzData=new double[maxMzLength]; 
       getImzMLData_p= new GetImzMLData(ibdFname, imzML);   
     }
     catch(const std::bad_alloc& e)
     {
       printf("Error reserving memory: %s\n",e.what());
       return 0;
     }
     
     for(int i=0; i<nPixels; i++)
     {
       px=pxList[i];
       massSize=getImzMLData_p->getPixelMassF(px, rawMzData);//mass vector
       if(massSize<=0)
       {
         printf("ERROR reading from the %s file\n", ibdFname);
         hit=true;
         break;
       }
       massSize--;
       if(rawMzData[0]<minMz) minMz=rawMzData[0];
       if(rawMzData[massSize]>maxMz) maxMz=rawMzData[massSize];
     }
     minPx++; maxPx++; //for R
   }
   
   if(getImzMLData_p) delete getImzMLData_p;
   if(rawMzData)      delete [] rawMzData;
   if(!hit){
     List ret=List::create(Named("minPixel")=minPx, Named("maxPixel")=maxPx, Named("minMz")=minMz, Named("maxMz")=maxMz);
     return ret;
   }
   return 0;
 }
 
 //'  @name rGetAverageGaussianSpectrum
 //'  @title converts the info in the imzML file into two arrays: Gaussians average values versus masses
 //'  
 //'  @param ibdFname:  absolute reference to the file with the ibd extension.
 //'  @param imzML:     list with information extracted from the imzML file with import_imzML()
 //'  @param params:    specific parameters
 //'              "SNR": signal-to-noise ratio
 //'          tolerance: desired tolerance for the centroids (ppm)
 //'      "noiseMethod": method for estimating noise.
 //'  @param mzLow:        lower  mass  to consider
 //'  @param mzHigh:       higher mass  to consider
 //'  @param pxList:       list of pixels. First pixel=1. By default everyone.
 //'  @param oversampling: interval between points on the mass axis = tolerance/oversampling.
 //'  @param nThreads:     number of threads suggested for parallel processing.
 //'  @return lista: averageMz and averageIntensity
 //'     averageMz: array of masses at intervals of 1/4 of the tolerance
 //'     averageIntensity: array of average values with all Gaussians 
 //'     
 // [[Rcpp::export]]
 List rGetAverageGaussianSpectrum(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, double mzLow, double mzHigh, Rcpp::NumericVector pxList, double oversampling, int nThreads)
 {
   gPeakCount=0, gSpectra=0;
   char dir[200];
   Common common;
   common.getFileNameDirectory(ibdFname, dir);
   
   NumericVector nv;
   nv=params["tolerance"];
   double tolerance=nv[0];
   
   RawToGaussians gauss(dir, ibdFname, imzML, params, pxList, mzLow, mzHigh, nThreads);  
   if(gauss.m_hit==false) return 0;
   
   //loads data from a file and converts its peak into Gaussians.
   int ret1=gauss.rawToGaussians(); //parallel processing
   if(ret1<0) 
   {return 0;}
   if(gPeakCount<=0) 
      {printf("No peaks are detected in the sample.\n"); return 0;}
   List ret=gauss.getMeanGaussianSpectrum(tolerance,(int)oversampling);
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
 //'  @param oversampling: interval between points on the mass axis = tolerance/oversampling.
 //'  @return lista: averageMz and averageIntensity
 //'     averageMz: array of masses at intervals of mass tolerance/oversampling
 //'     averageIntensity: array of average values with all Gaussians 
 //'     
 // [[Rcpp::export]]
 List rGetAverageSpectrum(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, double mzLow, double mzHigh, Rcpp::NumericVector pxList, double oversampling, int nThreads)
 {
   char dir[200];
   Common common;
   common.getFileNameDirectory(ibdFname, dir);
   
   NumericVector nv;
   nv=params["tolerance"];
   double tolerance=nv[0];
   
   RawToGaussians gauss(dir, ibdFname, imzML, params, pxList, mzLow, mzHigh, nThreads);  
   if(gauss.m_hit==false) return 0;
   
   List ret=gauss.getMeanSpectrum(tolerance, (int)oversampling);
   return ret;
 }
 
 //'  @name rSaveMassRange
 //'  @title save the mz range information to file _massRange.bin
 //'  
 //'  @param mzLow    -> low  m/z
 //'  @param mzHigh   -> high m/z
 //'  @param fileName ->absolute path to file
 //'  @return TRUE if the file could not be opened.
 //'  
 // [[Rcpp::export]]
 bool rSaveMassRange(const char* fileName, double mzLow, double mzHigh, NumericVector pixelSize)
 {
   std::fstream fp;
   fp.open(fileName, std::fstream::out | std::ios::binary | std::ios::trunc);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: %s file could not be opened\n", fileName);
     throw std::runtime_error(txt);
     return TRUE;
   }
   
   fp.write((char*)&mzLow,  sizeof(double)); 
   fp.write((char*)&mzHigh, sizeof(double)); 
   int pxSize=pixelSize.size();
   double tmp;
   fp.write((char*)&pxSize, sizeof(int));
   for(int i=0; i<pxSize; i++)
   {
     tmp=pixelSize[i];
     fp.write((char*)&tmp, sizeof(double)); 
   }
   
   fp.close();
   return FALSE;
 }
 
 //'  @name rLoadMassRange
 //'  @title load the mz range information from file _massRange.bin
 //'  @param fileName ->absolute path to file
 //'  @return a vector whit mzLow, mzHigh. Zero if error
 //'  
 // [[Rcpp::export]]
 NumericVector rLoadMassRange(const char* fileName)
 {
   std::fstream fp;
   fp.open(fileName, std::fstream::in | std::ios::binary);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: %s file could not be opened\n", fileName);
     throw std::runtime_error(txt);
     return 0;
   }
   double mzLow, mzHigh, pixelSize;
   fp.read((char*)&mzLow,  sizeof(double)); 
   fp.read((char*)&mzHigh, sizeof(double)); 
   int pxSize;
   double tmp;
   fp.read((char*)&pxSize, sizeof(int)); 
   NumericVector ret(2+pxSize);
   ret[0]=mzLow;
   ret[1]=mzHigh;
   for(int i=0;  i<pxSize; i++)
   {
     fp.read((char*)&tmp, sizeof(double)); 
     ret[2+i]=tmp;
   }
   fp.close();
   return ret;
 }
 
 //'  @name rGetDirectory
 //'  @title 
 //'  @param fileName ->absolute path to file
 //'  @return a vector 
 //'  
 // [[Rcpp::export]]
 CharacterVector rGetDirectory(const char* path)
 {
   Common common;
   char fileName[200];
   common.getFileNameDirectory(path, fileName);
   CharacterVector vect=fileName;//(strlen(fileName));
   return vect;
 }
 

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
   int totalSamples=0, nSamples, nSamplePixels, nPixels;
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
     printf("warning: max samples in file are %d\n", totalSamples);
     fp.close();
     return 0;
   }
   
   //get pixels range from sample parameter
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
 
 
 //' @name rGetBasic()
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
   //   bool hit=true;
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
       break;
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
   int *colSize_p=0;
   
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
   try{
     massAxis_p=new double[nIons];
     colSize_p=new int[nIons];
   }
   catch(const std::bad_alloc& e)
   {
     printf("Error reserving memory: %s\n",e.what());
     return 0;
   }
   
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
   try{
     intensity_p=new double[nPx];
     pixel_p=new int[nPx];
   }
   catch(const std::bad_alloc& e)
   {
     printf("Error reserving memory: %s\n",e.what());
     return 0;
   }
   
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
   int totalSamples, totalPx, totalIons;
   Common common;
   
   fp.open(file, std::fstream::in | std::ios::binary);
   if(!fp.is_open())
   {
     char txt[200];
     sprintf(txt, "Error: The internal %s file could not be created.\n The peak matrix cannot be saved.\n", file);
     throw std::runtime_error(txt);
   }
   //metadata
   fp.read((char*)&totalSamples, sizeof(int));  //samples number
   fp.read((char*)&totalPx,  sizeof(int));      //total pixels in all samples
   fp.read((char*)&totalIons,    sizeof(int));  //ions number
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
   
   //Data on disk is organized by centroids.
   //Each centroid is accompanied by a list of pixels—with non-zero intensity values—from all samples, in pixel-intensity format.
   //To separate them, the pixel range associated with each sample must be known.
   
   //the minimum and maximum pixels in each sample are delimited
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
     for(int i=0; i<pxSample[sample]; i++) //zeroing of the column
       pkMat(i, ion)=0;
     
     //Info on this centroid
     fp.read((char*)&tmpMass, sizeof(double));
     fp.read((char*)&tmpIntensity, sizeof(double));
     fp.read((char*)&tmpTolerance, sizeof(double));
     fp.read((char*)&nPx, sizeof(int));
     if(fp.fail()) 
     {
       printf("Error while reading in %s", file);
       fp.close(); return 0;
     }
     pxLow=samplesPxLimit[sample].x; //range of pixels of interest
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
 
 
