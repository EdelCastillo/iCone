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
#ifndef SIREM_PEAK_H
#define SIREM_PEAK_H

#include "stdio.h"
#include <stdexcept>
#include "peak.h"
#include "types.h"


class IntensityPeak
{
    public:
        //Constructor: initialization
        IntensityPeak(double SNR, bool unitedPeak=true);

        //destructor: frees reserved memory.
        ~IntensityPeak();

        //Get a list of simple peak associated with each magnitude simple peak.
        //Gets the list of compound peak: simple peak joined by their valleys.
        //makes use of the magnitude information.
        //Return -1 if unable to do so
        int getPeakList(SPECTRO *spectro_p);

        //conversión de valores en un sistema no lineal
        //Argumentos:
        //  dataIn_p  -> viene en unidades en el rango [0:size-1]
        //  dataOut_p -> son los datos convertidos
        //  newUnit_p -> es un array con las unidades de destino
        //  size      -> es el tamaño de todos los arrays
        void unitConversion(double *dataIn_p, double *dataOut_p, double *newUnit_p, int size);

        //Returns the information of a simple magnitude peak.
        ION_INDEX getSinglePeak(int mPeak);

        //Returns the number of simple peak.
        int getSinglePeakNumber();

        //Returns the simple low and high peak that make up a composite peak.
        UNITED_PEAK getCompoundPeak(int peak);

        //Returns the number of compound peak.
        int getCompoundPeakNumber();

    private:

        //The array of CONVOLVED_PEAK structures, of simple peak, is analyzed in search of possible compound peak.
        //A composite peak is a grouping of consecutive simple peak joined by a valley of magnitude greater than m_noiseLevel.
        //The info is in the array m_unitedPeak_p[] of size m_nUnitedPeak. The size is returned.
        //Arguments:
        // cnvPeak_p -> array of structures with simple peak info
        // nMagPeak -> array size
        //Returns the number of compound peak
        int unitedPeak(PEAK_LIST *peakList_p, SPECTRO *spectro_p);

        Peak           *m_intPeak_p;  //pointer to Peak class for magnitude.
        int             m_nUnitedPeak; //number of entries in m_unitedPeak_p.
        bool            m_unitedPeak; //true if connected peaks are considered.
        PEAK_LIST       m_peakList;
        int             m_nIntPeak;
        double           m_SNR;
};

#endif
