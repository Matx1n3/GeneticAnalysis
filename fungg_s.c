/*
   CA - OpenMP
   fungg_s.c
   Routines used in gengroups_s.c program

   TO BE COMPLETED
***************************************************************/


#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "definegg.h"		// definition of constants


/* 1 - Function to calculate the genetic distance; Euclidean distance between two elements.
       Input:   two elements of NFEAT characteristics (by reference)
       Output:  distance (double)
***************************************************************************************************/

double geneticdistance (float *elem1, float *elem2)
{
    double distance = 0.0;
    for (int i = 0; i < 40; i++){
        distance += pow(elem1[i] - elem2[i], 2);
    }
    return sqrt(distance);

}



/* 2 - Function to calculate the closest group (closest centroid) for each element.
   Input:   nelems   number of elements, int
            elems    matrix, with the information of the elements, of size MAXELE x NFEAT, by reference
            cent    matrix, with the centroids, of size NGROUPS x NFEAT, by reference
   Output:  grind   vector of size MAXELE, by reference, closest group for each element
***************************************************************************************************/

void closestgroup (int nelems, float elems[][NFEAT], float cent[][NFEAT], int *grind)
{
   // TO DO
   // grind: closest group/centroid for each element
   int minGroup;
   double minDistance;
   double actDistance;
    for (int i = 0; i < nelems; i++){
        minGroup = 0;
        minDistance = geneticdistance(elems[i], cent[0]);
        for (int j = 1; j < ngroups; j++){
            actDistance = geneticdistance(elems[i], cent[j]);
            if (actDistance < minDistance) {
                minDistance = actDistance;
                minGroup = j;
            }
        }
        grind[i] = minGroup;
    }
}

/* 3 - Function to validate the classification: group compactness and centroids compactness 
       Calculate the CVI index

       Input:  elems    elements (matrix of size MAXELE x NFEAT, by reference)
               iingrs   indices of the elements in each group (vector with information for each group)
               cent     matrix, with the centroids, by reference
       Output: cvi index
               compact  compactness of each group (vector of size NGROUPS, by reference) 
******************************************************************************************/

double validation (float elems[][NFEAT], struct ginfo *iingrs, float cent[][NFEAT], float *compact)
{

   // TO DO

   // Calculate group compactness: mean distance between elements of each group
   // ==========================================================================

          
   // Calculate cluster distances: mean distance between centroids 
   // ===========================================================================
   
  
 
   // Calculate CVI index
   // ===================

}




/* 4 - Function to analyse diseases 
   Input:  iingrs   indices of the elements in each group (matrix with MAXELE elements per group, by reference)
           dise     information about the diseases
   Output: disepro  analysis of the diseases: maximum, minimum of the medians and groups
***************************************************************************************************/

void diseases (struct ginfo *iingrs, float dise[][TDISEASE], struct analysis *disepro)
{

  // TO DO
   // process the information about diseases. 
   // for each group obtain the median value of each disease. To calculate the median order the
   // probabilities of the diseases of each group (using a simple algorithm for instance bubble
   // sort) and take the value in the middle position.
   // for each disease obtain the maximum and the group where the maximum is found 
   // the minimum and the group where the minimum is found (for each disease)
}




// TWO OTHER FUNCTIONS IN THE MAIN PROGRAM
// ========================================


/* 5 - Initial values for centroids
**************************************************************/

void firstcentroids (float cent[][NFEAT])
{
  int  i, j;


  srand (147);
  for (i=0; i<ngroups; i++)
  for (j=0; j<NFEAT/2; j++)
  {
     cent[i][j] = (rand() % 10000) / 100.0;
     cent[i][j+NFEAT/2] = cent[i][j];
  }
}



/* 6 - New centroids  
***************************************************************/

int newcentroids (float elems[][NFEAT], float cent[][NFEAT], int grind[], int nelems)
{
  int     i, j, finish;
  float   newcent[ngroups][NFEAT];
  double  discent;
  double  additions[ngroups][NFEAT+1];


  // calculate new centroids for each group:  average of each dimension or feature
  // additions: to accumulate the values for each feature and cluster. Last value: number of elements in the group

  for (i=0; i<ngroups; i++)
  for (j=0; j<NFEAT+1; j++)
    additions[i][j] = 0.0;

  for (i=0; i<nelems; i++)
  {
    for (j=0; j<NFEAT; j++)
      additions[grind[i]][j] += elems[i][j];
    additions[grind[i]][NFEAT] ++;
  }


  finish = 1;
  for (i=0; i<ngroups; i++)
  {
    if (additions[i][NFEAT] > 0) //the group is not empty
    {
      for (j=0; j<NFEAT; j++) newcent[i][j] = additions[i][j] / additions[i][NFEAT];

      // decide if the process needs to be finished

      discent = geneticdistance (&newcent[i][0], &cent[i][0]);
      if (discent > DELTA1) finish = 0;        // there is change at least in one of the dimensions; continue with the process

      // copy new centroids
      for (j=0; j<NFEAT; j++) cent[i][j] = newcent[i][j];
     }
   }
   
   return (finish);
}
