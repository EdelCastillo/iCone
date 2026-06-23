/*********************************************************************************
 *     GetImzMLData
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
#include <Rcpp.h>
#include "common_methods.h"

using namespace Rcpp;

#include "imzMLBin.h"

class GetImzMLData
{
public:
  //Constructor
  //ibdFname: Absolute path to the imzML file.
  //imzML: Information about the file's contents.
  GetImzMLData(const char* ibdFname, Rcpp::List imzML);
  
  //Destructor
  //The class that reads from the imzML file is released.
  ~GetImzMLData();
  
  //getPixelSpectrum()
  //pixel: Reference to the desired spectrum.
  //Returns the intensity vector along with the mass vector.
  DataFrame getPixelSpectrum(int pixel);
  
  //getPixelMass()
  //pixel: Reference to the desired spectrum.
  //Returns the mass vector.
  NumericVector getPixelMass(int pixel);
  
  //getPixelIntensity()
  //pixel: Reference to the desired spectrum.
  //Returns the intensity vector.
  NumericVector getPixelIntensity(int pixel);
  
  //getPixelMassF()
  //pixel: Reference to the desired spectrum.
  //data_p: Pointer to the mass data.
  //Returns the size of the returned data.
  int getPixelMassF(int pixel, float *data_p);
  
  //getPixelIntensityF()
  //pixel: Reference to the desired spectrum.
  //data_p: Pointer to the intensity data.
  //Returns the size of the returned data.
  int getPixelIntensityF(int pixel, float *data_p);
  

  ImzMLBinRead *m_myReader_p;  
  Rcpp::IntegerVector m_mzLength,
                      //m_mzOffset,
                      m_intLength;
                      //m_intOffset;
  Rcpp::NumericVector m_mzOffset,
                      m_intOffset;
  
  Rcpp::String        m_mz_dataType, 
                      m_int_dataType;
  int m_NPixels;
  bool m_continuous;
};