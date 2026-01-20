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
#ifndef KMEANS_R 
#define KMEANS_R

#include "types.h"
#include <cmath>
#include <stdlib.h>


typedef struct
{
  GROUP   data;     //indexes to elements associated with each segment.
  double  center,   //controide
          withinss, //cumulative square distances.
          acu;      //data summation (internal use)
}CLUSTER;

typedef struct
{
  CLUSTER *clusters_p; //particular information of each segment.
  double  
          totss,        //The total sum of squares.
          totWithinss,  //Total within-cluster sum of squares, i.e. sum(withinss).
          betweenss;    //The between-cluster sum of squares, i.e. totss-tot.withinss.
   int
          nClusters,    //number of clusters
          iter,         //The number of (outer) iterations.
          ifault;       //integer: indicator of a possible algorithm problem.
}K_MEANS;


class KmeansR
{
public:
  //Constructor.
  //Implements one-dimensional segmentation in the style of kmeans.
  //data_p: Data to segment.
  //size: Array size.
  //maxIter: Maximum number of iterations allowed.
  //convergenceValue: Average of squared deviations for convergence.
  //Sets m_error = true if the array is invalid.
  KmeansR(float *data_p, int size, int maxIter, double convergenceValue=1e-4);
  
  //Segmentation.
  //Initialization can operate in two ways: based on passed data or randomly.
  //Arguments.
  //nCluster: Number of desired segments.
  //massInit_p: Pointer to initial data (if null, initialize randomly).
  //Sets m_kStruct.ifault to 1 if the algorithm does not converge. =-1 if the algorithm fails.
  //Returns 0 if OK.
  int getClusters(int nClusters, float *massInit_p);
  
  //destructor
  //frees reserved memory.  
  ~KmeansR();
  
  //frees reserved memory.
  void freeClusters();
  
  K_MEANS  m_kStruct;
  bool      m_error;
  
private:  
  
  float     *m_data_p;
  int       m_size,
            m_nClusters,
            m_maxIter;
  double    m_convergenceValue;
};

#endif
