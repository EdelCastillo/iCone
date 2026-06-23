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

#include "common_methods.h"

//for data ordenation
struct
{
  bool operator()(int a, int b) const { return a > b; }
}
cmp_Down;
struct
{
  bool operator()(int a, int b) const { return a < b; }
}
cmp_Up;

struct
{
  bool operator()(float a, float b) const { return a > b; }
}
cmp_DownF;
struct
{
  bool operator()(float a, float b) const { return a < b; }
}
cmp_UpF;

struct
{
  bool operator()(double a, double b) const { return a > b; }
}
cmp_DownD;
struct
{
  bool operator()(double a, double b) const { return a < b; }
}
cmp_UpD;


std::string parse_xml_uuid(std::string uuid)
{
  std::size_t ipos = uuid.find('{');
  if( ipos != std::string::npos )
  {
    uuid.erase(ipos, 1);
  }
  do
  {
    ipos = uuid.find('-'); 
    if( ipos != std::string::npos )
    {
      uuid.erase(ipos, 1);
    } 
  } while ( ipos != std::string::npos);
  ipos = uuid.find('}');
  if( ipos != std::string::npos )
  {
    uuid.erase(ipos, 1);
  }
  for( unsigned int i=0; i < uuid.length(); i++)
  {
    uuid[i] = toupper(uuid[i]);
  }
  return uuid;
}

Common::Common(){}

//TIC Normalization
bool Common::TICnormalization(float *vect, int size)
{
  double suma=0;
  for(int i=0; i<size; i++)
    suma+=vect[i];
  suma/=1e6;
  if(suma<1e-10) return(false); //error control
  for(int i=0; i<size; i++)
    vect[i]/=suma;
  return(true);
}

//ascending ordering of doubles values.
//If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
int Common::sortUp(double *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] < bufferIn[i2];});
  }
  else
    stable_sort(bufferIn, bufferIn+size, cmp_UpD /*[&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] < bufferIn[i2];}*/);
  return 0;
}

//descending ordering of doubles values.
//If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
int Common::sortDown(double *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] > bufferIn[i2];});
  }
  else
    stable_sort(bufferIn, bufferIn+size, cmp_DownD);
  return 0;
}

//ascending ordering of floats values.
//If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
int Common::sortUpF(float *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] < bufferIn[i2];});
  }
  else    stable_sort(bufferIn, bufferIn+size, cmp_UpF /*[&bufferIn](size_t i1, size_t i2) {return ((float*)bufferIn)[i1] < ((float*)bufferIn)[i2];}*/);
  return 0;
}

//descending ordering of floats values.
//If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
int Common::sortDownF(float *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] > bufferIn[i2];});
  }
  else
    stable_sort(bufferIn, bufferIn+size, cmp_DownF);
  return 0;
}

//ascending ordering of integer values.
//If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
int Common::sortUpI(int *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] < bufferIn[i2];});
  }
  else    stable_sort(bufferIn, bufferIn+size, cmp_Up /*[&bufferIn](size_t i1, size_t i2) {return ((float*)bufferIn)[i1] < ((float*)bufferIn)[i2];}*/);
  return 0;
}

//descending ordering of integer values.
//If index=null, the sorted information is returned in bufferIn; otherwise, bufferIn remains unchanged.
int Common::sortDownI(int *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] > bufferIn[i2];});
  }
  else
    stable_sort(bufferIn, bufferIn+size, cmp_Down);
  return 0;
}


//returns the index to data closest to value.
int Common::nearestIndex(int value, int *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1) //the nearest one returns.
  {
    if(value-data[indexLow] <= data[indexHigh]-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow])       {/*printf(".%d\n", indexLow); */return indexLow;}
  else if(value>=data[indexHigh]) {/*printf("..%d\n", indexHigh);*/return indexHigh;}
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter]) {/*printf("...%d\n", indexCenter);*/return indexCenter;}
    if(value<data[indexCenter]) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow] <= data[indexHigh]-value) {/*printf("....%d\n", indexLow);*/return indexLow;}
      else return indexHigh;
    }
  }
}

//returns the index to data closest to value.
int Common::nearestIndex(float value, float *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1) //the nearest one returns.
  {
    if(value-data[indexLow] <= data[indexHigh]-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow])       {/*printf(".%d\n", indexLow); */return indexLow;}
  else if(value>=data[indexHigh]) {/*printf("..%d\n", indexHigh);*/return indexHigh;}
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter]) {/*printf("...%d\n", indexCenter);*/return indexCenter;}
    if(value<data[indexCenter]) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow] <= data[indexHigh]-value) {/*printf("....%d\n", indexLow);*/return indexLow;}
      else return indexHigh;
    }
  }
}


//returns the index to data closest to value.
int Common::nearestIndex(double value, double *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1) //the nearest one returns.
  {
    if(value-data[indexLow] <= data[indexHigh]-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow])       return indexLow;
  else if(value>=data[indexHigh]) return indexHigh;
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter]) return indexCenter;
    if(value<data[indexCenter]) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow] <= data[indexHigh]-value) return indexLow;
      else return indexHigh;
    }
  }
}

//Returns the index of data closest to value
//If nearest bits 1:0 == 00, returns the nearest
//If nearest bits 1:0 == 01, returns the nearest above
//If nearest bits 1:0 == 10, returns the nearest below
//If nearest bit    2 == 0 & value is out of range returns the nearest
//If nearest bit    2 == 1 & value is out of range returns -1
int Common::nearestIndexGaussians(float value, GAUSS_PARAMS *data, int size, int nearest)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  if(value<data[indexLow].mean)
  {
    if(nearest&4) return -1;
    else return 0;
  }
  if(value>data[indexHigh].mean)
  {
    if(nearest&4) return -1;
    else return indexHigh;
  }
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1)
  {
    if((nearest&3)==1)      return indexHigh;
    else if((nearest&3)==2) return indexLow;
    else
    {
      if(value-data[indexLow].mean <= data[indexHigh].mean-value) return indexLow;
      else return indexHigh;
    }
  }
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter].mean) return indexCenter;
    if(value<data[indexCenter].mean) indexHigh=indexCenter; 
    else indexLow=indexCenter;
  if(indexHigh==indexLow) return indexLow;
  else if(indexHigh==indexLow+1)
    {
      if((nearest&3)==1)      return indexHigh;
      else if((nearest&3)==2) return indexLow;
      else
      {
        if(value-data[indexLow].mean <= data[indexHigh].mean-value) return indexLow;
        else return indexHigh;
      }
    }
  }
}


//returns the index to data.low closest to value.
int Common::nearestIndexMassRangeLow(float value, MASS_RANGE *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;

  if(size<=0) return -1;
  if(indexHigh==indexLow) //the nearest one returns.
  {
    if(value-data[indexLow].low <= data[indexHigh].low-value) return indexLow;
    else return indexHigh;
  }
  if(value<=data[indexLow].low)       return indexLow;
  else if(value>=data[indexHigh].low) return indexHigh;

  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter].low) return indexCenter;
    if(value<data[indexCenter].low) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow].low <= data[indexHigh].low-value) return indexLow;
      else return indexHigh;
    }
  }
}

//returns the index to data.high closest to value.
int Common::nearestIndexMassRangeHigh(float value, MASS_RANGE *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(size<=0) return -1;
  if(indexHigh==indexLow) //the nearest one returns
  {
    if(value-data[indexLow].high <= data[indexHigh].high-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow].high)       return indexLow;
  else if(value>=data[indexHigh].high) return indexHigh;
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter].high) return indexCenter;
    if(value<data[indexCenter].high) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow].high <= data[indexHigh].high-value) return indexLow;
      else return indexHigh;
    }
  }
}

//Returns the median.
//If the data has an even length, the average of the two middle values is returned.
float Common::median(double *data, int size)
{
  if (size<=0) return -1;
  double median;
  int half=size/2; //If even, it keeps the lower value.
  double *tmpData=0;
  tmpData=new double[size];
  for(int i=0; i<size; i++) //copy, since sort() destroys the origin.
    tmpData[i]=data[i];
  sortUp(tmpData, 0, size); //increasing order.
  
  if(!(size%2)) {//if even
    if(half>0)
    {median=(tmpData[half]+tmpData[half-1])/2.0;} //average value
    else {median=tmpData[0];}
  }
  else
  {median=tmpData[half];}
  if(tmpData) delete []tmpData;
  return median;
}

//Returns the median.
//If the data has an even length, the average of the two middle values is returned.
float Common::medianF(float *data, int size)
{
  if (size<=0) return -1;
  float median;
  int half=size/2; //If even, it keeps the lower value.
  float *tmpData=0;
  tmpData=new float[size];
  for(int i=0; i<size; i++) //copy, since sort() destroys the origin.
    tmpData[i]=data[i];
  
  sortUpF(tmpData, 0, size); //increasing order.

  if(!(size%2)) {//if even
    if(half>0)
    {median=(tmpData[half]+tmpData[half-1])/2.0;} //average value
    else {median=tmpData[0];}
  }
  else
    {median=tmpData[half];}
   if(tmpData) delete []tmpData;
  return median;
}

//average value of an array of floats.
float Common::meanF(float *A, int size)
{
  float S=0;
  if(size<=0) return -1;
  for(int i=0; i<size; i++)
    S+=A[i];
  return S/size;
}

//average value of an array of doubles
double Common::mean(double *A, int size)
{
  double S=0;
  if(size<=0) return -1;
  for(int i=0; i<size; i++)
    S+=A[i];
  return S/size;
}

//average value of an array of doubles
double Common::mean(double *prob, double *data, int size)
{
  double S=0;   
  if(size<=0) return -1;
  for(int i=0; i<size; i++)
    S+=prob[i]*data[i];
  return S;
}

//average value of an array of floats.
float Common::meanF(float *prob, float *data, int size)
{
  float S=0;   
  if(size<=0) return -1;
  for(int i=0; i<size; i++)
    S+=prob[i]*data[i];
  return S;
}

//variance of an array of floats.
float Common::varF(float *A, int size, float *mean_p)
{
  float m;
  float S=0;
  if(size<=0) return -1;
  if(mean_p) m=*mean_p;
  else m=meanF(A, size);
  
  for(int i=0; i<size; i++)
    S+=(A[i]-m)*(A[i]-m);
  return S/size;
}

//variance of an array of double
double Common::var(double *A, int size, double *mean_p)
{
  double m;
  double S=0;
  if(size<=0) return -1;
  
  if(mean_p) m=*mean_p;
  else m=mean(A, size);
  
  for(int i=0; i<size; i++)
    S+=(A[i]-m)*(A[i]-m);
  return S/size;
}

//variance of an array of doubles
double Common::var(double *prob, double *data, int size, double *mean_p)
{
  double m;
  double S=0;
  if(size<=0) return -1;
  
  if(mean_p) m=*mean_p;
  else m=mean(prob, data, size);
  
  for(int i=0; i<size; i++)
    S+=(data[i]-m)*(data[i]-m)*prob[i];
  return S;
}

//variance of an array of floats.
float Common::varF(float *prob, float *data, int size, float *mean_p)
{
  float m;
  float S=0;
  if(size<=0) return -1;
  if(mean_p) m=*mean_p;
  else m=meanF(prob, data, size);
  
  for(int i=0; i<size; i++)
    S+=(data[i]-m)*(data[i]-m)*prob[i];
  return S;
}

//Extracts the file extension from a filename.
//Returns true if the extension is not found.
bool Common::getFileExtension(const char *fileName, char* extension)
{
  int index=0;
  for(int i=strlen(fileName); i; i--)
  {
    if(fileName[i]=='.') {index=i; break;} 
  }
  if(index) 
  {
    strcpy(extension, fileName+index+1);
    return false;
  }
  
  return 0;
}

//Extracts the file name from a path, includes the possible extension.
//Returns true if fileName is null.
bool Common::getFileNameWithExt(const char *fileName, char* name)
{
  if(strlen(fileName)==0) return true;
  int index=0;
  for(int i=strlen(fileName); i; i--)
  {
    if(fileName[i]=='/' || fileName[i]=='\\') {index=i; break;} 
  }
  strcpy(name, fileName+index+1);
  
  return 0;
}

//Extracts the file name from a path, excludes the possible extension.
//Returns true if fileName is null.
bool Common::getFileNameWithoutExt(const char *fileName, char* name)
{
  if(strlen(fileName)==0) return true;
  int index1=0, index2=0;
  for(int i=strlen(fileName); i; i--)
  {
    if(fileName[i]=='.') index1= i;
    if(fileName[i]=='/' || fileName[i]=='\\') {index2=i+1; break;} 
  }
  for(int i=index2, k=0; i<index1; i++, k++)
    name[k]=fileName[i];
  name[index1]=0;
  return 0;
}

//Extracts the directory from a path.
//Returns true if fileName is null.
bool Common::getFileNameDirectory(const char *fileName, char* directory)
{
  if(strlen(fileName)==0) return true;
  int index=0;
  for(int i=strlen(fileName); i; i--)
  {
    if(fileName[i]=='/' || fileName[i]=='\\') {index=i; break;} 
  }
  if(index==0) strcpy(directory, ".");
  else
    for(int i=0; i<=index; i++)
      directory[i]=fileName[i];
  directory[index+1]=0;
  
  return 0;
}
