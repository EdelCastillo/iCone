/*************************************************************************
 *     common_methods 
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
 **************************************************************************/


#ifndef RMSI_COMMON_METHODS_H
#define RMSI_COMMON_METHODS_H

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <string>
#include <string.h>
#include "types.h"

using namespace std; 

  
//Parse the UUID from a XML file to get just the hex representation in a string (without dashes and {})
std::string parse_xml_uuid(std::string uuid);

class Common
{
public:
  Common();
  //TIC Normalization
  bool TICnormalization(float *vect, int size);
  
  //ascending ordering of doubles values.
  //If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
  int sortUp  (double *bufferIn, int *sort, int size);
  
  //descending ordering of doubles values.
  //If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
  int sortDown(double *bufferIn, int *sort, int size);
  
  //ascending ordering of floats values.
  //If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
  int sortUpF  (float *bufferIn, int *sort, int size);
  
  //descending ordering of floats values.
  //If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
  int sortDownF(float* bufferIn, int *sort, int size);
  
  //ascending ordering of integer values.
  //If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
  int sortUpI  (int *bufferIn, int *sort, int size);
  
  //descending ordering of integer values.
  //If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
  int sortDownI(int *bufferIn, int *sort, int size);
  
  //returns the index to data closest to value.
  int nearestIndex(int value, int *data, int size);
  
  //returns the index to data closest to value.
  int nearestIndex(float value, float *data, int size);
  
  //returns the index to data closest to value.
  int nearestIndex(double value, double *data, int size);
  
  //Returns the index of data closest to value
  //If nearest== 1, returns the nearest above
  //If nearest==-1, returns the nearest below
  //If nearest== 0, returns the nearest
  int nearestIndexGaussians(float value, GAUSS_PARAMS *data, int size, int nearest);
  
  //returns the index to data.low closest to value.
  int nearestIndexMassRangeLow(float value, MASS_RANGE *data, int size);
  
  //returns the index to data.high closest to value.
  int nearestIndexMassRangeHigh(float value, MASS_RANGE *data, int size);
  
  //Returns the median.
  //If the data has an even length, the average of the two middle values is returned.
  float median(double *data, int size);
  
  //Returns the median.
  //If the data has an even length, the average of the two middle values is returned.
  float medianF(float *data, int size);
  
  //average value of an array of floats.
  float meanF(float *A, int size);
  
  //average value of an array of doubles
  double mean(double *A, int size);
  
  //average value of an array of floats.
  float meanF(float *prob, float *data, int size);
  
  //average value of an array of doubles
  double mean(double *prob, double *data, int size);
  
  //variance of an array of floats.
  float varF(float *A, int size, float *mean=NULL);
  
  //variance of an array of double
  double var(double *A, int size, double *mean=NULL);
  
  //variance of an array of doubles
  double var(double *prob, double *data, int size, double *mean_p=NULL);
  
  //variance of an array of floats.
  float varF(float *prob, float *data, int size, float *mean_p=NULL);

  //Extracts the file extension from a filename.
  //Returns true if the extension is not found.
  bool getFileExtension(const char *fileName, char* extension);
    
  //Extracts the file name from a path, includes the possible extension.
  //Returns true if fileName is null.
  bool getFileNameWithExt(const char *fileName, char* name);

  //Extracts the file name from a path, excludes the possible extension.
  //Returns true if fileName is null.
  bool getFileNameWithoutExt(const char *fileName, char* name);

  //Extracts the directory from a path.
  //Returns true if fileName is null.
  bool getFileNameDirectory(const char *fileName, char* directory);
    
};
  
#endif
