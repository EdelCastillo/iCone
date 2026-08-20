/*************************************************************************
 *     simple and compound peaks
 *     iCone - R package for MSI data processing
 *     Copyright (C) octubre 2023 Esteban del Castillo Pérez
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

#include "peakInfo.h"

//Constructor: initialization
IntensityPeak::IntensityPeak(double SNR)
{
    m_SNR=SNR;
    m_peakList.intPeak_p=0;
    m_peakList.uPeak_p=0;
    m_intPeak_p=0;
}

//destructor: free reserved memory
IntensityPeak::~IntensityPeak()
{
//printf("...ini intPeak destructor\n");
    if(m_peakList.intPeak_p)    delete [] m_peakList.intPeak_p;
    if(m_peakList.uPeak_p)      delete [] m_peakList.uPeak_p;
    if(m_intPeak_p) delete m_intPeak_p;

//printf("...end intPeak destructor\n");
}

//Get a list of simple peak associated with each magnitude simple peak
//Get the list of compound peak: simple peak joined by their valleys
//makes use of the magnitude information
//Return -1 if unable to do so
int IntensityPeak::getPeakList(SPECTRO *spectro_p)
{
  int nIntPeak=0;
  for(int i=0; i<spectro_p->size; i++)
    if(spectro_p->SNR_p[i]< m_SNR) spectro_p->int_p[i]=0;
    
  try
        {
        //get the list of peak from the magnitude
        //each block in the list contains three indices: mzLow, mzMax and mzHigh
         m_intPeak_p=new Peak(spectro_p, m_SNR);
        nIntPeak=m_intPeak_p->get(0, spectro_p->size-1); //#magnitude peak

        if(nIntPeak<=0) return -1; //no peak in data

        m_peakList.intPeak_p=new ION_INDEX [nIntPeak];
        for(int i=0; i<nIntPeak; i++)
        {
          m_peakList.intPeak_p[i].low =m_intPeak_p->m_mzIndex_p[i].low;
          m_peakList.intPeak_p[i].max =m_intPeak_p->m_mzIndex_p[i].max;
          m_peakList.intPeak_p[i].high=m_intPeak_p->m_mzIndex_p[i].high;

        }
        m_peakList.nIntPeak=nIntPeak;
        m_peakList.uPeak_p=new UNITED_PEAK[nIntPeak];
        
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }
  m_nUnitedPeak=unitedPeak(&m_peakList, spectro_p);
  m_peakList.nUPeak=m_nUnitedPeak;
  m_nIntPeak=nIntPeak;
  return nIntPeak;
}

//The array of CONVOLVED_PEAK structures, of simple peak, is analyzed in search of possible compound peak.
//A composite peak is a grouping of consecutive simple peak joined by a valley of magnitude greater than m_noiseLevel.
//The info is in the array m_unitedPeak_p[] of size m_nUnitedPeak. The size is returned.
//Arguments:
// cnvPeak_p -> array of structures with simple peak info
// nMagPeak -> array size
//Returns the number of compound peak
int IntensityPeak::unitedPeak(PEAK_LIST *peakList_p, SPECTRO *spectro_p)
{
    int nUPeak=0, nIntPeak=peakList_p->nIntPeak;
    bool hit;
    int iPeak=0;
    while(true)
        {
        //if the index is less than the maximum and is linked to the next peak
      if(iPeak<(nIntPeak-1) && (peakList_p->intPeak_p[iPeak].high==peakList_p->intPeak_p[iPeak+1].low) &&
           spectro_p->SNR_p[peakList_p->intPeak_p[iPeak].high]>m_SNR)
            {
          //se toma nota
            peakList_p->uPeak_p[nUPeak].peakLow =iPeak;
            peakList_p->uPeak_p[nUPeak].peakHigh=iPeak+1;
            iPeak++;
            //perhaps there are more peak linked in a chain
            hit=false;
            while(true)
                {
              //if the index is less than the maximum and is linked to the next peak.
                if(iPeak<(nIntPeak-1) && (peakList_p->intPeak_p[iPeak].high==peakList_p->intPeak_p[iPeak+1].low) &&
                   spectro_p->SNR_p[peakList_p->intPeak_p[iPeak].high]>m_SNR)
                    {
                  //the end of the joint peak is updated.
                    peakList_p->uPeak_p[nUPeak].peakHigh=iPeak+1;
                    iPeak++;
                    }
                //If you have reached the last peak and it is linked to the previous one: end.
                else if(iPeak>=(nIntPeak-1))
                    {
                  iPeak++; nUPeak++;
                    hit=true;
                    break;
                    }
                //if you have not reached the end and are not linked to the next one.
                else {iPeak++; nUPeak++; break;} //two separate peak.
                }
            if(hit) break;
            }

        else  //single peak
            {
            peakList_p->uPeak_p[nUPeak].peakLow =iPeak;
            peakList_p->uPeak_p[nUPeak].peakHigh=iPeak;
            nUPeak++;
            if(iPeak >= nIntPeak-1) break;
            else iPeak++;
            }
        }
  return nUPeak;
}

//conversión de valores en un sistema no lineal
//Argumentos:
//  dataIn_p  -> viene en unidades en el rango [0:size-1]
//  dataOut_p -> son los datos convertidos
//  newUnit_p -> es un array con las unidades de destino
//  size      -> es el tamaño de todos los arrays
void IntensityPeak::unitConversion(double *dataIn_p, double *dataOut_p, double *newUnit_p, int size)
{
  double data, delta, offset;
  int tmp;

  //ajuste de mean
  for(int i=0; i<size; i++)
  {
    data=dataIn_p[i];
    tmp=(int)data;
    //delta se adecúa a la diferencial de los datos finales
    if(tmp+1 < size) //si dentro de rango
      delta=newUnit_p[tmp+1]-newUnit_p[tmp];//delta posterior
    else
      delta=newUnit_p[tmp]-newUnit_p[tmp-1];//delta anterior

    offset=delta*(data-(double)tmp); //desplazamiento respecto al origen del pico compuesto
    dataOut_p[i]=newUnit_p[tmp]+offset; //valor convertido

  }
}

//Returns the information of a simple magnitude peak.
ION_INDEX IntensityPeak::getSinglePeak(int iPeak)
{return m_peakList.intPeak_p[iPeak];}

//Returns the number of simple peak.
int IntensityPeak::getSinglePeakNumber()
    {return m_nIntPeak;}

//Returns the simple low and high peak that make up a composite peak.
UNITED_PEAK IntensityPeak::getCompoundPeak(int uPeak)
    {return m_peakList.uPeak_p[uPeak];}

//Returns the number of compound peak.
int IntensityPeak::getCompoundPeakNumber()
    {return m_nUnitedPeak;}


