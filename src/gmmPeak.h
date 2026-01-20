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
#ifndef GMM_PEAK_H
#define GMM_PEAK_H

#include "stdio.h"
#include <stdexcept>
#include "peak.h"
#include "GMM.h"
#include "types.h"


//Decomposes a single composite peak, possibly formed by several simple peak, into Gaussians.
class GmmPeak
{
public:
    //Constructor:
    //Argument:
    //minMeanPeakMagnitude -> minimum averaged magnitude for a peak to be deconvolved
    GmmPeak(float minMeanPeakMagnitude);

    //Destructor: frees reserved memory.
    ~GmmPeak();

    //magnitudes associated with all the peak to be treated
    //Arguments:
    //mag_p -> structure with magnitude information associated with each scan
    void setMagnitudes(GROUP_F *mag_p);

    //peak to consider
    //Arguments:
    //cnvPeak_p -> pointer to array of the set of simple peak that constitute the composite peak
    //nUnitedPeak -> #simple peak
    void setPeak(ION_INDEX *intPeak_p, int nIntPeak);

    //Deconvolution
    //From a compound magnitude peak and its entropy peak,
    //obtains the Gaussians that compose it: EM Algorithm
    //establishes the array m_deconv_p[] with m_nDeconv elements with info of each Gaussian ordered from lowest to highest mass
    //If a peak could not be deconvolved, only one structure appears, its field nGauss=0 and its field quality=-1.
    //The mean and sigma information of the Gaussians is in units of scans
    //return:
    //Number of gaussians
    int gmmDeconvolution();

    //conversión de unidades
    //mean y sigma tienen unidades de scans y deben convertirse a Daltons
    //yFactor debe adaptarse a las nuevas unidades
    //se procede a una interpolación lineal
    //Argumentos:
    //  deconvIn_p -> puntero a estructura con unidades originales (scans)
    //  deconOut_p -> puntero a estructura con unidades de Daltons
    //  mzAxis_p   -> puntero a array con unidades de Dalton.
    //                Cada nuevo valor es un nuevo scan. Asociación entre ejes de mz y magnitudes
    //  mzAxisSize -> tamaño del array de mz. Debe coincidir con el tamaño de las magnitudes a analizar
    //NOTA: si los scans son muy irregulares, puede no hacer una buena conversión!!!!
    void gaussConversion(GAUSSIAN *deconvIn_p, GAUSSIAN *deconvOut_p, float *mzAxis_p, int mzAxisSize);

    //Returns information about the sirem peak that are candidates for generating Gaussians
    unsigned long getEtpHits();

    //Returns the number of Gaussians that integrate a compound magnitude peak
    int getDeconvNumber();

    //Returns the Gaussian indexed by 'index' housed in a compound magnitude peak
    GAUSSIAN getDeconv(int index);

    //Retorna la dirección de una gausiana de un pico de magnitud compuesto
    GAUSSIAN *getDeconv_p(int index);

    //Returns the quality of the Gaussian fit with the composite magnitude info
    float getQuality();

    //Returns a pointer to the magnitude information handled by the composite peak
    float *getMagnitude();

    //Returns an element of the magnitude information handled by the composite peak
    float getMagnitude(int index);

    //Returns the number of magnitude elements in the composite peak
    int getMagnitudeNumber();

private:
    //prepare the necessary information before calling the GMM algorithm: initialize the m_sGmm structure
    //Requires the reference to the composite peak to be deconvolved and the associated magnitudes
    //Return -1 if the average magnitude within the peak does not reach a minimum (m_minMeanPeakMagnitude)
    int iniGMM();

    float           *m_mag_p,       //magnitude array
                    m_minPeakMagnitude; //peak of lower magnitude are discarded
    int             m_nMag;         //magnitude array size
    GMM_STRUCT      m_sGmm;         //structure for deconvolution of a peak
    int             m_nUnitedPeak; //number of entries in m_unitedPeak_p
    GAUSSIAN        *m_deconv_p,    //pointer to gausians
                    m_meanDeconv[32]; //array de gausianas con valores promediados
    int             m_nDeconv;      //deconvolution number of gaussians
    unsigned long   m_etpHits;      //sirem peak array
    float           m_quality;      //quality of the Gaussian fit
    ION_INDEX       *m_intPeak_p;
    int             m_nIntPeak;
};

#endif
