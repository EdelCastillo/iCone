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
#include "segments.h"

//constructor
//class: extracts mass segments where overlapping Gaussians exist. Over all pixels.
//The source information is located in a temporary file.
Segments::Segments(double massResolution, double mzLow, double mzHigh, double linkedPeaks)
{
  
  m_totalPixels=0; //It is set in the loadGaussians() function.
  m_massResolution=massResolution;
  m_mzHigh=mzHigh;
  m_mzLow=mzLow;
  m_linkedPeaks=linkedPeaks;
  m_gaussians_p=0;
  m_massRange_p=0;
}

//destructor
Segments::~Segments()
{
  if(m_gaussians_p)
  {
    for(int i=0; i<m_totalPixels; i++)
    {
      if(m_gaussians_p[i].gauss_p) delete [] m_gaussians_p[i].gauss_p;
    }
  delete []m_gaussians_p;
  }
  if(m_massRange_p) delete [] m_massRange_p;
  //printf("~Segments\n");
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
int Segments::getMassRanges(double* linkedPeak_p)
{
  printf("\tphase 2:\t#isolate mass segments: ");
  double highMass, lowMass, highMassOld, lowMassOld;
  int iLow, iHigh;
  double deltaMass=m_mzLow/(4*m_massResolution); //delta=1/4 of the minimum mass increment of the spectrometer
  int massAxisSize=1+(m_mzHigh-m_mzLow)/deltaMass;
  
  double *massAxis=0;
  massAxis= new double[massAxisSize];
  double linkedPeaks=m_linkedPeaks;
  int nUPeak=0, px;
  SPECTRO  spectro;
  spectro.SNR_p=0;
  spectro.SNR_p=new double[massAxisSize];
  bool hit;
  int iSegment=0;
  
  for(int i=0; i<massAxisSize; i++) massAxis[i]=0;
  double deltaMass2;
  

  //for all elements of the spectrum
  //vector con la cantidad de contribuciones de las gaussianas de cada espectro
  for(int px=0; px<m_totalPixels; px++)
  {
    if(m_gaussians_p[px].gauss_p==0 || m_gaussians_p[px].size==0) {/*printf("...%d\n", px);*/ continue;} //if this entry does not contain info.
    lowMassOld=0; highMassOld=0;
    for(int i=0; i<m_gaussians_p[px].size; i++) 
    {
      //Some ions may be zero in mean/sigma/weight. They are discarded.
      if(m_gaussians_p[px].gauss_p[i].mean<m_mzLow || m_gaussians_p[px].gauss_p[i].mean>m_mzHigh) 
          continue;
      //Very wide Gaussians (sigma>5*deltaMass) are discarded.
      deltaMass2=(m_gaussians_p[px].gauss_p[i].mean)/m_massResolution;
      if(fabs(m_gaussians_p[px].gauss_p[i].sigma)>5*deltaMass2) //sigma >>
        continue;
      
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
      
      for(int j=iLow; j<=iHigh; j++) 
        if(j<massAxisSize) massAxis[j]+=1.0;
      
      if(highMass>highMassOld)
        highMassOld=highMass;
    }
  }

  //noise estimation
  NoiseEstimation noiseEst(3, 1, 9); //MAD, smoothig & windows=9
  double noise=noiseEst.getNoise(massAxis, massAxisSize);
  
  //The baseband is estimated (recursive averaging).
  //the iteration ends when the slope reaches a value less than 1% of the initial slope.
  //max iterations=10
  double acu=0, vMean;
  double value[2], minDiff=0;

  for(int j=0; j<massAxisSize; j++) acu+=massAxis[j]; vMean=acu/massAxisSize;

  for(int i=0; i<10; i++)
  {
    acu=0;
    for(int j=0; j<massAxisSize; j++)
    {
      if(massAxis[j]<=vMean){acu+=massAxis[j];}
      else acu+=vMean;
    }
    vMean=acu/massAxisSize;
    if(i==0) {value[0]=vMean;}
    else if(i==1) {value[1]=vMean; minDiff=abs(value[0]-value[1])*0.01; }
    else {if(abs(vMean-value[1])<minDiff) break; else value[1]=vMean;}
  }
 
  //decision of the value assigned to noise.
  vMean=noise>vMean?noise:vMean; // the greatest.
  double pxSupport=m_totalPixels*0.0001; //minimum spectra that support.
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
  double lowMz, highMz;
  hit=true;
  double localSNR=1;
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
      if(iSegment>=maxSegments) 
      {
        printf("Warning: The limit of planned segments (%d) has been reached.\nIt is suggested to increase the value of the linkedPeaks argument.", maxSegments);
        break;   //no more space
      }
      
      unsigned int maxG=0;
      for(int i=iLowMass; i<iHighMass; i++)
      {
      //  if(massAxis[i]>maxG)maxG=massAxis[i];
          acu+=massAxis[i];
      }
      //acu*=4.0/(iHighMass-iLowMass);
      //if(maxG*4>=m_pxSupport)
      {
        m_massRange_p[iSegment].low =lowMz;//package ends
        m_massRange_p[iSegment].high=highMz;
        m_massRange_p[iSegment++].nGaussians=acu;
      }
      //else printf("...\n");
    }
    else
    {
      //It is analyzed in terms of the range of conflicting masses.
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
        IntensityPeak intPeak(localSNR); //simple and compound peaks are obtained with SNR update.
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
          if(highMz-lowMz > 3.0) {hit=false;break;} //iterate to improve
        }
        if(hit) //solution found.
        {
          if(iSegment>=maxSegments)
          {
            printf("Warning: The limit of planned segments (%d) has been reached.\nIt is suggested to increase the value of the linkedPeaks argument.", maxSegments);
            break;   //no more space
          }
          int acu=0;
          for(int j=0; j<nUPeak2; j++) 
            {
              pLow =intPeak.getCompoundPeak(j).peakLow; //index to lower peak of the composite peak.
              pHigh=intPeak.getCompoundPeak(j).peakHigh;//index to upper peak of the composite peak.
              iLowMass2 =intPeak.getSinglePeak(pLow). low +iLowMass;
              iHighMass2=intPeak.getSinglePeak(pHigh).high+iLowMass;
              lowMz =m_mzLow+iLowMass2 *deltaMass; //lower  mass.
              highMz=m_mzLow+iHighMass2*deltaMass; //higher mass.
              acu=0;
              unsigned int maxG=0;
              for(int i=iLowMass2; i<iHighMass2; i++)
              {
                //if(massAxis[i]>maxG)maxG=massAxis[i];
                acu+=massAxis[i];
              }
              //if(maxG*4>=m_pxSupport)
              {
                m_massRange_p[iSegment].low =lowMz;//package ends
                m_massRange_p[iSegment].high=highMz;
                m_massRange_p[iSegment++].nGaussians=acu;
              }
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

//Load the file with information about the Gaussian curves associated with each pixel
//If there is more than one sample, all its pixels are integrated
//fileName: temporary file name
//totalPixels: total pixels in the file
//Returns false if the load failed
int Segments::loadGaussians(char *fileName)
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
      fp.read((char*)&nPxGauss, sizeof(int)); //#gaussianas into the pixel
      if(fp.eof()) break;
      m_gaussians_p[pxIndex].size=nPxGauss;
      if(nPxGauss==0) {pxIndex++; continue;}
//  totalGaussians+=nPxGauss;
      
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
//  printf("...%d\n", totalGaussians);
  m_totalPixels=pxIndex;
  return pxIndex;
}


