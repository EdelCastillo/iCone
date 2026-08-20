/*************************************************************************
 *     Deconvolution of peak
 *     Copyright (C) 2023 Esteban del Castillo Pérez
 *     esteban.delcastillo@urv.cat
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
 **************************************************************************/

#include "gmmPeak.h"

//Constructor:
//Argument:
//minPeakMagnitude -> minimum magnitude for a peak to be deconvolved
GmmPeak::GmmPeak(double minPeakMagnitude)
{
    m_minPeakMagnitude=minPeakMagnitude;
    m_deconv_p=0;
    m_nDeconv=0;
    m_sGmm.x=0;
    m_sGmm.y=0;
    m_mag_p=0;
}

//Destructor: frees reserved memory.
GmmPeak::~GmmPeak()
{
//printf("...ini GmmPeak destructor\n");
    if(m_deconv_p)  delete [] m_deconv_p;
    if(m_sGmm.x) delete []m_sGmm.x;
    if(m_sGmm.y) delete []m_sGmm.y;
    //printf("...end GmmPeak destructor\n");

}

//magnitudes associated with all the peak to be treated
//Arguments:
//mag_p -> structure with magnitude information associated with each scan
void GmmPeak::setMagnitudes(GROUP_F *mag_p)
{
    m_mag_p=mag_p->set;
    m_nMag=mag_p->size;
}

//peak to consider
//Arguments:
//cnvPeak_p -> pointer to array of the set of simple peak that constitute the composite peak
//nUnitedPeak -> #simple peak
void GmmPeak::setPeak(ION_INDEX *intPeak_p, int nIntPeak)
{
  m_intPeak_p=intPeak_p;
  m_nIntPeak=nIntPeak;

    if(m_sGmm.x) {delete []m_sGmm.x; m_sGmm.x=0;}
    if(m_sGmm.y) {delete []m_sGmm.y; m_sGmm.y=0;}
    int nScans=0;
    nScans=m_intPeak_p[nIntPeak-1].high-m_intPeak_p[0].low+1; //#composite peak scans
    nScans+=2;   //to accommodate the extreme zeros
    m_sGmm.x=0;
    m_sGmm.y=0;
    try
        {
        m_sGmm.x=new double [nScans]; //mz
        m_sGmm.y=new double [nScans]; //magnitude
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return;
        }

    m_sGmm.size=nScans;
    m_sGmm.maxIter=100;     //maximum iterations of the GMM algorithm
    m_sGmm.relChange=0.001; //minimum change between iterations to stop the algorithm
    m_sGmm.quality=0.9;     //fit quality
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//prepare the necessary information before calling the GMM algorithm: initialize the m_sGmm structure
//Requires the reference to the composite peak to be deconvolved and the associated magnitudes
//Return:
//  -1 if the average magnitude within the peak does not reach a minimum (m_minPeakMagnitude)
//  -2 if gaussian number>DECONV_MAX_GAUSSIAN
int GmmPeak::iniGMM()
{
    int nGauss=0;
    m_etpHits=0;
    int etpLow, etpHigh, mzSize;
    double sumaMagEpeak, sumaMagMpeak, maxEpeak;
    ION_INDEX *magPeak_p;

    int lowPeakIndex, highPeakIndex, lowMzIndex, highMzIndex;
    lowPeakIndex =0;  //reference to the lower simple peak
    highPeakIndex=m_nIntPeak-1;//reference to the upper simple peak

    if(lowPeakIndex>highPeakIndex) return -1; //limit control

    lowMzIndex =m_intPeak_p[lowPeakIndex] .low;  //lower mass index
    highMzIndex=m_intPeak_p[highPeakIndex].high; //upper mass index
    mzSize=highMzIndex-lowMzIndex+1;

    if(lowMzIndex>=highMzIndex) return -1; //limit control

    //sum of the magnitudes within the magnitude peak
    sumaMagMpeak=0;
    double maxValue=0, tmpMag;
    for(int i=0; i< mzSize; i++)
        {
        tmpMag=m_mag_p[i];
        sumaMagMpeak+=tmpMag;
        if(tmpMag > maxValue) maxValue=tmpMag;
        }
    if(maxValue < m_minPeakMagnitude)
//    if(sumaMagMpeak < m_minPeakMagnitude)
      {m_sGmm.quality=-1; nGauss=0; return -1;}
    
    //The array of magnitudes and the array of masses are created
    //A zero is artificially added to the data on each side of the peak: it improves the fit.
     m_sGmm.x[0]=0.0;//To the left of the peak
    m_sGmm.y[0]=0.0;
    int j=1;
    for(int i=0; i< mzSize; i++, j++)
        {
        m_sGmm.x[j]=i; //magnitude peak
        m_sGmm.y[j]=m_mag_p[i];
        }
    //A la derecha del pico
    m_sGmm.x[j]=highMzIndex+1;
    m_sGmm.y[j]=0.0;
    m_sGmm.size=mzSize+2; //data size increased by two
    nGauss=0; //num of gausianas
    int iCount=-1;
    m_nDeconv=0;
    
    //in each single magnitude peak there may be several entropy peak, or none.
    //are made to fit within the limits of the peak magnitude
    for(int peak=lowPeakIndex; peak<=highPeakIndex; peak++)
      {
      magPeak_p=&m_intPeak_p[peak];
        //the parameters of the Gaussian are initialized. Same as before only with magnitude values
      if(magPeak_p->low > magPeak_p->high || magPeak_p->max<magPeak_p->low || magPeak_p->max>magPeak_p->high)
      {m_sGmm.params[nGauss].weight=0; nGauss=0; break;}
      
        m_sGmm.params[nGauss].mean=magPeak_p->max-lowMzIndex;
        m_sGmm.params[nGauss].sigma=((double)magPeak_p->high-magPeak_p->low)/2.0;
        //m_sGmm.params[nGauss].weight=1; //sum_of_each_ion/sum_of_peak
        sumaMagEpeak=0; maxEpeak=0;
        for(int j=magPeak_p->low; j<=magPeak_p->high; j++)//isolated peak
        {
          tmpMag=m_mag_p[j-lowMzIndex];
          sumaMagEpeak+=tmpMag;
          //if(sumaMagEpeak>maxEpeak)maxEpeak=sumaMagEpeak;
          if(tmpMag>maxEpeak) maxEpeak=tmpMag;
        }
        //if(sumaMagEpeak<m_minPeakMagnitude) {printf("...\n"); continue;}
        if(maxEpeak<m_minPeakMagnitude) {
          //printf("....%.4f/%.4f\n",maxEpeak, m_minPeakMagnitude); 
          continue;}
        else m_sGmm.params[nGauss].weight=sumaMagEpeak/sumaMagMpeak;
        
        //the limits of the parameters of the Gaussian are established
        m_sGmm.limits[nGauss].minMean=magPeak_p->low-lowMzIndex-2;
        if(m_sGmm.limits[nGauss].minMean<0) m_sGmm.limits[nGauss].minMean=0;
        m_sGmm.limits[nGauss].maxMean=magPeak_p->high-lowMzIndex+2;
        double minSigma=0.5;
        m_sGmm.limits[nGauss].minSigma=minSigma;
        m_sGmm.limits[nGauss].maxSigma=(m_sGmm.limits[nGauss].maxMean-m_sGmm.limits[nGauss].minMean)*1.0;
        if(m_sGmm.params[nGauss].sigma<=minSigma)  m_sGmm.params[nGauss].sigma=minSigma;
        m_sGmm.limits[nGauss].minWeight=0.0;
        m_sGmm.limits[nGauss].maxWeight=1.0;
        m_etpHits=1;
        nGauss++;
        m_nDeconv++;
        if(nGauss>=DECONV_MAX_GAUSSIAN) return -2;
        }
    m_sGmm.maxIter=nGauss*40; //maximum iterations for adjustment in GMM
    return nGauss; //number of Gaussians that will make up the peak magnitude
}

//Deconvolution
//From a compound magnitude peak and its entropy peak,
//obtains the Gaussians that compose it: EM Algorithm
//establishes the array m_deconv_p[] with m_nDeconv elements with info of each Gaussian ordered from lowest to highest mass
//if an element of m_deconv_p[] has its values null -> non-deconvolved entropy peak.
//If a peak could not be deconvolved, only one structure appears, its field nGauss=0 and its field quality=-1.
//The mean and sigma information of the Gaussians is in units of scans
//return:
//Number of gaussians
int GmmPeak::gmmDeconvolution()
{
    Cgmm cGmm; //Class for GMM

    int nGauss=0;

    nGauss=iniGMM(); //structure initialization for GMM

    if(nGauss>0) //if it can be decomposed
        {
        //the Gaussians that make up the magnitude peak are constructed
        m_sGmm.nGauss=nGauss;
        cGmm.gmm(&m_sGmm);
        }
    else if(nGauss<0)//if the average magnitude within the peak does not reach a minimum (m_minPeakMagnitude)
        {
        m_sGmm.quality=-1;
        m_nDeconv=nGauss;//a structure appears in the file without deconvolution info
        if(m_deconv_p) {delete [] m_deconv_p; m_deconv_p=0;} //old memory is released
        return nGauss;
        }
    
    if(isnan(m_sGmm.quality)) m_sGmm.quality=-1; //limit control

    if(m_deconv_p) {delete [] m_deconv_p; m_deconv_p=0;} //old memory is released
    try
        {
        m_deconv_p=new GAUSSIAN[m_nDeconv];   //new memory is reserved
        //m_deconv_p=new GAUSSIAN[nGauss];   //new memory is reserved
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }
    //m_nDeconv=nGauss;

    //saves deconvolution info.
    //m_sGmm returns the constructed Gaussians but this function returns all of them; 
    //that is to say not all entropy peak may generate Gaussians. Those who don't
    //they set all their values to zero
//    for(int g=0; g<nGauss; g++)
    for(int g=0, k=0; g<m_nDeconv; g++)
      {
      if(1) //if((m_etpHits|((unsigned long)1<<g))!=0)          //if the gaussian exist
        {
        if(isnan(m_sGmm.params[k].mean))      {m_sGmm.params[k].mean=0;   m_sGmm.quality=-1;}
        if(isnan(m_sGmm.params[k].sigma))     {m_sGmm.params[k].sigma=0;  m_sGmm.quality=-1;}
        if(isnan(m_sGmm.params[k].weight))    {m_sGmm.params[k].weight=0; m_sGmm.quality=-1;}
        m_deconv_p[g].mean=m_sGmm.params[k].mean;       //mean value
        m_deconv_p[g].sigma=m_sGmm.params[k].sigma;     //standard deviation
        m_deconv_p[g].weight=m_sGmm.params[k].weight;   //weighted weight
        m_deconv_p[g].yFactor=m_sGmm.yFactor;           //factor for reconstruction
        k++;
//        if(m_deconv_p[g].sigma<0) printf("..........%.5f\n",m_deconv_p[g].mean);
        }
      else
        {
        m_deconv_p[g].mean=0;       
        m_deconv_p[g].sigma=0;     
        m_deconv_p[g].weight=0;   
        m_deconv_p[g].yFactor=0;           
        }
      }
      
    m_quality=m_sGmm.quality;           //fit quality
    
//return nGauss; //all good
return m_nDeconv; //all good
}

//unit conversion.
//mean and sigma have units of scans and must be converted to Daltons.
//yFactor must adapt to the new units.
//proceeds to linear interpolation.
//Arguments:
// deconvIn_p -> pointer to structure with original units (scans).
// deconOut_p -> pointer to structure with units of Daltons.
// mzAxis_p   -> pointer to array with Dalton units.
//               Each new value is a new scan. Association between mz axes and magnitudes.
// mzAxisSize -> mz array size. It must coincide with the size of the magnitudes to be analyzed.
//NOTA: si los scans son muy irregulares, puede no hacer una buena conversión!!!!
//////////////////////////////////////////
void GmmPeak::gaussConversion(GAUSSIAN *deconvIn_p, GAUSSIAN *deconvOut_p, double *mzAxis_p, int mzAxisSize)
{
  double delta, delta1, delta2, offset, mean, sigma;
  int tmp;
  if(mzAxisSize<=1 || deconvIn_p->mean<=0)
  {
    deconvOut_p->mean=0;
    deconvOut_p->sigma=0;
    deconvOut_p->yFactor=0;
    return;
  }
  //mean adjustment
  mean=deconvIn_p->mean;
  tmp=(int)deconvIn_p->mean; //lower integer value.
  //If tmp is very large, it may become negative when truncated to integer.
  
  //interval in Daltons between two consecutive neighboring scans of mean.
  if(tmp < mzAxisSize-1 && tmp>0) //if it is within range.
  {
    delta1=mzAxis_p[tmp+1]-mzAxis_p[tmp];//posterior delta.
    delta2=mzAxis_p[tmp]-mzAxis_p[tmp-1];//anterior delta.
//    delta=delta1<delta2?delta1:delta2;   //mínima distancia
    delta=delta1;
    delta2=(delta1+delta2)/2.0;
  }
  else if(tmp==0) //origen eje de masas
  {
    delta=mzAxis_p[tmp+1]-mzAxis_p[tmp];//previous delta.
    delta2=delta;
  }
  else if(tmp==mzAxisSize-1) //fin del eje de masas
  {
    delta=mzAxis_p[tmp]-mzAxis_p[tmp-1];//previous delta.
    delta2=delta;
  }
  else //error
  {
    deconvOut_p->mean=0;
    deconvOut_p->sigma=0;
    deconvOut_p->yFactor=0;
    return;
  }

  offset=delta*(mean-(double)tmp); //displacement with respect to the origin of the compound peak.
  //  if(offset<0) return;
  deconvOut_p->mean=mzAxis_p[tmp]+offset; //converted value.
  //sigma adjustment, as in mean.
  sigma=deconvIn_p->sigma;
  tmp=(int)sigma;
  offset=delta*(sigma-(double)tmp);
//  deconvOut_p->sigma=tmp*delta+offset;
  deconvOut_p->sigma=tmp*delta2+offset;
  
  //Adjustment in Y. yFactor is increased by the quotient factor between the sigma in Daltons and in scans.
  if(deconvIn_p->sigma<1e-10) 
  {
    deconvOut_p->yFactor=0;
  }
  else
  {
    deconvOut_p->yFactor = deconvIn_p->yFactor*deconvOut_p->sigma/deconvIn_p->sigma;
  }
  deconvOut_p->yFactor/=(deconvOut_p->sigma*2.50662827); //normalizing factor (area=1.0).
  deconvOut_p->weight=deconvIn_p->weight; //The weight remains the same.
}


//Returns information about the sirem peak that are candidates for generating Gaussians
unsigned long GmmPeak::getEtpHits()
    {return m_etpHits;}

//Returns the number of Gaussians that integrate a compound magnitude peak
int GmmPeak::getDeconvNumber()
    {return m_nDeconv;}

//Returns the Gaussian indexed by 'index' housed in a compound magnitude peak
GAUSSIAN GmmPeak::getDeconv(int index)
    {return m_deconv_p[index];}

//Returns the direction of a Gaussian indexed by 'index' housed in a compound magnitude peak
GAUSSIAN *GmmPeak::getDeconv_p(int index)
    {return &m_deconv_p[index];}

//Returns the quality of the Gaussian fit with the composite magnitude info
double GmmPeak::getQuality()
    {return m_quality;}

//Returns a pointer to the magnitude information handled by the composite peak
double *GmmPeak::getMagnitude()
    {return m_mag_p;}

//Returns an element of the magnitude information handled by the composite peak
double GmmPeak::getMagnitude(int index)
    {return m_mag_p[index];}

//Returns the number of magnitude elements in the composite peak
int GmmPeak::getMagnitudeNumber()
    {return m_nMag;}
