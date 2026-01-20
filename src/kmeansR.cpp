/*************************************************************************
 *     kMeansR
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
#include "kmeansR.h"

//Constructor.
//Implements one-dimensional segmentation in the style of kmeans.
//data_p: Data to segment.
//size: Array size.
//maxIter: Maximum number of iterations allowed.
//convergenceValue: Average of squared deviations for convergence.
//Sets m_error = true if the array is invalid.
KmeansR::KmeansR(float *data_p, int size, int maxIter, double convergenceValue)
{
  m_error=0;
  m_convergenceValue=convergenceValue;
  m_size=size;
  m_maxIter=maxIter;
  m_data_p=data_p;
  
  m_kStruct.clusters_p=0;
  m_kStruct.totss=0;
  m_kStruct.totWithinss=0;
  m_kStruct.betweenss=0;
  m_kStruct.iter=0;
  m_kStruct.ifault=0;
  if(!data_p || size<=0) {m_error=true; return;}
}


//destructor
//frees reserved memory.  
KmeansR::~KmeansR()
{
  //printf("KmeansR destructor init\n");
  freeClusters();
  //printf("KmeansR destructor end\n");
}

//frees reserved memory.  
void KmeansR::freeClusters()
{
  if(m_kStruct.clusters_p)
  {
    for(int i=0; i<m_kStruct.nClusters; i++)
    {
      if(m_kStruct.clusters_p[i].data.set)
      {
        delete [] m_kStruct.clusters_p[i].data.set;
        m_kStruct.clusters_p[i].data.set=0;
      }
    }
    delete []m_kStruct.clusters_p; 
    m_kStruct.clusters_p=0;
  } 
}

//Segmentation.
//Initialization can operate in two ways: based on passed data or randomly.
//Arguments.
//nCluster: Number of desired segments.
//massInit_p: Pointer to initial data (if null, initialize randomly).
//Sets m_kStruct.ifault to 1 if the algorithm does not converge. =-1 if the algorithm fails.
//Returns 0 if OK.
int KmeansR::getClusters(int nClusters, float *massInit_p)
{
  if(nClusters<=0 || m_size<=0) return -1;
  int *binIndex_p=0;
  
  bool random=false;
  freeClusters(); //frees any previously reserved memory.
  if(!massInit_p) //if no initial masses are provided -> randomness.
  {
    binIndex_p=new int[nClusters]; //indexes to bins
    
    //random initialization mode
    random=true;
    for(int bin=0; bin<nClusters; bin++)
    {
        int ran;
        while(1) //avoid repeated values
            {
            bool hit=false;
            ran=rand()%m_size; 
            for(int i=0; i<bin; i++)
                if(ran==binIndex_p[i]) {hit=true; break;} //if it already exists
            if(!hit) break;
            }    
        binIndex_p[bin]=ran; //indexes  
    }
  }
  //initialization
  double distance;
  double *centralMass_p =0; 
  double *centralMass2_p=0; 
  centralMass_p = new double[nClusters]; //centroids
  centralMass2_p= new double[nClusters];
  
  //memory reservation and its initialization
  m_kStruct.clusters_p=new CLUSTER[nClusters];
  m_kStruct.nClusters=nClusters;  
  for(int i=0; i<nClusters; i++)
  {
    m_kStruct.clusters_p[i].data.set=0;
    m_kStruct.clusters_p[i].data.set=new int[m_size];
    m_kStruct.clusters_p[i].data.size=0;
    m_kStruct.clusters_p[i].withinss=0;
    m_kStruct.clusters_p[i].acu=0;
    }
  m_kStruct.totss=0;
  m_kStruct.totWithinss=0;
  
  //initial centroids
  for(int iCenter=0; iCenter<nClusters; iCenter++)
  {
    if(random)
      centralMass_p[iCenter]=binIndex_p[iCenter];
    else
      centralMass_p[iCenter]=massInit_p[iCenter]; //past values
  }

  int iter=0;
  if(binIndex_p) {delete []binIndex_p; binIndex_p=0;}  
  
  //iterations until convergence
  while(1)
  {
    //Assigning elements to centroids (segmentation).
    //The closest elements are assigned to a centroid (one-dimensional distance).
    for(int iMass=0; iMass<m_size; iMass++)//for all the elements to consider
    {
      int minIndex=-1;
      double minSqrDistance=1e10;
      for(int iCenter=0; iCenter<nClusters; iCenter++) //for all desired segments
      {
        double sqrDistance=(centralMass_p[iCenter]-m_data_p[iMass])*(centralMass_p[iCenter]-m_data_p[iMass]);
        if(sqrDistance<minSqrDistance) 
          {minSqrDistance=sqrDistance; minIndex=iCenter;} //minimum distance
      }
      if(minIndex==-1) return -1; //should never occur, unless nClusters==0 and is controlled.
      //data
      m_kStruct.clusters_p[minIndex].data.set[m_kStruct.clusters_p[minIndex].data.size++]=iMass;
      m_kStruct.clusters_p[minIndex].acu+=m_data_p[iMass]; //to average later
      m_kStruct.clusters_p[minIndex].withinss+=minSqrDistance; //measure of dispersion
    }
  
  //New centroids are established (averaged across elements).
  //The deviation of the centroids from the previous iteration is determined.
  //The deviation of each segment is weighted according to its number of elements.
  double varAcu=0;
  for(int iCenter=0; iCenter<nClusters; iCenter++)
  {
    //new center
    int size=m_kStruct.clusters_p[iCenter].data.size;
    if(size==0 ) centralMass2_p[iCenter]=0;
    else
      centralMass2_p[iCenter]=m_kStruct.clusters_p[iCenter].acu/size; //new centroid (average)
    
    //quadratic variation with respect to the previous centroid.
    double var=(centralMass_p[iCenter]-centralMass2_p[iCenter])*(centralMass_p[iCenter]-centralMass2_p[iCenter]);
    var*=(double)size/m_size; //is weighted based on its weight.
    varAcu+=var; //accumulated
  }
  varAcu/=nClusters; //mean square deviation.
  
  //If the average deviation is large and the iteration limit was not reached, iterate again.
  if(varAcu>m_convergenceValue && ++iter<m_maxIter)
  {
    for(int iCenter=0; iCenter<nClusters; iCenter++)
    {
      centralMass_p[iCenter]=centralMass2_p[iCenter]; //new centroids for the next iteration.
      m_kStruct.clusters_p[iCenter].data.size=0;
      m_kStruct.clusters_p[iCenter].acu=0;
      m_kStruct.clusters_p[iCenter].withinss=0;
    }
  }
  else //convergence or iterations have been exhausted.
    {
    //the data is noted
    for(int iCenter=0; iCenter<nClusters; iCenter++)
        {
      m_kStruct.clusters_p[iCenter].center=centralMass2_p[iCenter]; 
      m_kStruct.totWithinss+=m_kStruct.clusters_p[iCenter].withinss;
      m_kStruct.betweenss=0;
      m_kStruct.totss=0;
      m_kStruct.iter=iter;
      if(iter>=m_maxIter)
        m_kStruct.ifault=1; //non-convergence
        }
    break;
    }
  }
  //reserved memory is freed
  if(centralMass_p)  delete [] centralMass_p;
  if(centralMass2_p) delete [] centralMass2_p;
  
  return m_kStruct.ifault; //returns the error
}


