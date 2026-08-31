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
#include "rGetImzMLData.h"
using namespace Rcpp;

//Constructor
//ibdFname: Absolute path to the imzML file.
//imzML: Information about the file's contents.
GetImzMLData::GetImzMLData(const char* ibdFname, Rcpp::List imzML)
{
  Rcpp::DataFrame df;
  df=imzML["run"];
  m_continuous=imzML["continuous_mode"];
  m_mz_dataType=(const char*)imzML["mz_dataType"];
  m_int_dataType=(const char*)imzML["int_dataType"];
  m_NPixels=df.nrows();
  m_mzLength=df["mzLength"];
  m_mzOffset=df["mzOffset"];
  m_intLength=df["intLength"];
  m_intOffset=df["intOffset"];
  m_myReader_p=0;
  //clase para leer del fichero
//  printf("---%s %s\n", m_mz_dataType.get_cstring(), m_int_dataType.get_cstring());
  m_myReader_p=new ImzMLBinRead(ibdFname, m_NPixels, m_mz_dataType, m_int_dataType, m_continuous);
}

//Destructor
//The class that reads from the imzML file is released.
GetImzMLData::~GetImzMLData()
{
  if(m_myReader_p) delete m_myReader_p;
}

//getPixelSpectrum()
//pixel: Reference to the desired spectrum.
//Returns the intensity vector along with the mass vector.
DataFrame GetImzMLData::getPixelSpectrum(int pixel)
{
//  printf("...%d %d\n", m_mzLength[pixel], m_intLength[pixel]);
  NumericVector mass(m_mzLength[pixel]), intensity(m_intLength[pixel]);
  m_myReader_p->readMzData ((std::streampos)m_mzOffset [pixel], m_mzLength [pixel], mass.begin());
  m_myReader_p->readIntData((std::streampos)m_intOffset[pixel], m_intLength[pixel], intensity.begin());
  DataFrame df=DataFrame::create(Named("mass") =clone(mass), Named("int") =clone(intensity));
  return df;
}

//getPixelMass()
//pixel: Reference to the desired spectrum.
//Returns the mass vector.
NumericVector GetImzMLData::getPixelMass(int pixel)
{
  NumericVector mass(m_mzLength[pixel]);
  m_myReader_p->readMzData ((std::streampos)m_mzOffset[pixel], m_mzLength[pixel], mass.begin());
  return mass;
}

//getPixelIntensity()
//pixel: Reference to the desired spectrum.
//Returns the intensity vector.
NumericVector GetImzMLData::getPixelIntensity(int pixel)
{
  NumericVector intensity(m_intLength[pixel]);
  m_myReader_p->readIntData ((std::streampos)m_intOffset[pixel], m_intLength[pixel], intensity.begin());
  return intensity;
}

//getPixelMassF()
//pixel: Reference to the desired spectrum.
//data_p: Pointer to the mass data.
//Returns the size of the returned data.
int GetImzMLData::getPixelMassF(int pixel, double *data_p)
{
  int mzSize=m_mzLength[pixel];
  NumericVector mass(mzSize);
  m_myReader_p->readMzData ((std::streampos)m_mzOffset[pixel], mzSize, mass.begin());
  if(data_p==0) return -1;
  for(int i=0; i<mzSize; i++)
    data_p[i]=(double)mass[i];
  return mzSize;
}

//getPixelIntensityF()
//pixel: Reference to the desired spectrum.
//data_p: Pointer to the intensity data.
//Returns the size of the returned data.
int GetImzMLData::getPixelIntensityF(int pixel, double *data_p)
{
  int intSize=m_intLength[pixel];
  NumericVector intensity(intSize);
  m_myReader_p->readIntData ((std::streampos)m_intOffset[pixel], intSize, intensity.begin());
  if(data_p==0) return -1;
  for(int i=0; i<intSize; i++)
  {
    data_p[i]=(double)intensity[i];
  }
  return intSize;
}

//getPixelSpectrum() links to R.
//pixel: Reference to the desired spectrum.
//ibdFname: Absolute path to the imzML file.
//imzML: Information about the file's contents.
//Returns the intensity vector along with the mass vector for a given pixel.
//Applies sorting from lowest to highest, since the raw data may be jumbled.// [[Rcpp::export]]
DataFrame getPixelSpectrum(int pixel, const char* ibdFname, Rcpp::List imzML)
{
  Common tools;
  double *tmp_p=0;
  int *sort_p=0;
  GetImzMLData myReader(ibdFname, imzML);
  int size=myReader.m_mzLength [pixel];
  tmp_p=new double[size];
  sort_p=new int[size];
  
  NumericVector mass(myReader.m_mzLength[pixel]), intensity(myReader.m_intLength[pixel]);
//  myReader.m_myReader_p->readMzData (myReader.m_mzOffset [pixel], myReader.m_mzLength [pixel], mass.begin());
  myReader.m_myReader_p->readMzData ((std::streampos)myReader.m_mzOffset [pixel], size, tmp_p);
  tools.sortUp(tmp_p, sort_p, size); //ordenación
  for(int i=0; i<size; i++) mass[i]=tmp_p[sort_p[i]];
  
//  myReader.m_myReader_p->readIntData(myReader.m_intOffset[pixel], myReader.m_intLength[pixel], intensity.begin());
  myReader.m_myReader_p->readIntData((std::streampos)myReader.m_intOffset[pixel], myReader.m_intLength[pixel], tmp_p);
  for(int i=0; i<size; i++) intensity[i]=tmp_p[sort_p[i]];
  
  if(tmp_p)  delete [] tmp_p;
  if(sort_p) delete [] sort_p;
  
  DataFrame df=DataFrame::create(Named("mass") =clone(mass), Named("intensity") =clone(intensity));
  return df;
}

