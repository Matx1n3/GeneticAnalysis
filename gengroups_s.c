/* 
    CA - practical work OpenMP
    gengroups_s.c SERIAL VERSION

    Processing genetic characteristics to discover information about diseases
    Classify in NGROUPS groups, elements of NFEAT features, according to "distances"    

    Input:  dbgen.dat 	   input file with genetic information
            dbdise.dat     input file with information about diseases
    Output: results_s.out  centroids, number of group members and compactness, and diseases

    Compile with module fungg_s.c and include option -lm
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "definegg.h"
#include "fungg.h"


float            elems[MAXELE][NFEAT];	// matrix to keep information about every element
struct ginfo     iingrs[NGROUPSMAX];	// vector to store information about each group: members and size

float            dise[MAXELE][TDISEASE];// probabilities of diseases (from dbdise.dat)
struct analysis  disepro[TDISEASE];	// vector to store information about each disease (max, min, group...)

int  ngroups = 35;			// initial number of groups



// Main programa 
// ================

void main (int argc, char *argv[])
{
  float   cent[NGROUPSMAX][NFEAT], newcent[NGROUPSMAX][NFEAT]; //centroid and new centroid
  float   compact[NGROUPSMAX];	 	// compactness of each group or cluster

  int     i, j, nelems, group, count;
  int     grind[MAXELE];   		//group assigned to each element 
  int     finish = 0, niter = 0, finish_classif; 
  double  cvi, cvi_old, diff;

  FILE    *f1, *f2;
  struct timespec  t1, t2, t3, t4, t5;
  double  t_read, t_clus, t_anal, t_write;


  if ((argc < 3) || (argc > 4)) {
    printf ("ATTENTION:  progr file1 (elems) file2 (dise) [num elems])\n");
    exit (-1);
  }


  printf ("\n >> Serial execution\n");
  clock_gettime (CLOCK_REALTIME, &t1);

  
  // read data from files: elems[i][j] and dise[i][j]
  // ================================================
  
  f1 = fopen (argv[1], "r");
  if (f1 == NULL) {
    printf ("Error opening file %s \n", argv[1]);
    exit (-1);
  }

  fscanf (f1, "%d", &nelems);
  if (argc == 4) nelems = atoi(argv[3]);

  for (i=0; i<nelems; i++) 
  for (j=0; j<NFEAT; j++) 
    fscanf (f1, "%f", &(elems[i][j]));
  
  fclose (f1);

  f1 = fopen (argv[2], "r");
  if (f1 == NULL) {
    printf ("Error opening file %s \n", argv[1]);
    exit (-1);
  }

  for (i=0; i<nelems; i++) 
  for (j=0; j<TDISEASE; j++) 
    fscanf (f1, "%f", &dise[i][j]);
  
  fclose (f1);

  clock_gettime (CLOCK_REALTIME, &t2);


  // PHASE 1: iterative process to classify elements in ngroups groups 
  // until a minimum difference (DELTA2) in the CVI index
  // =================================================================

  finish_classif = 0; 
  cvi_old = -1;

  while ((ngroups < NGROUPSMAX) && (finish_classif == 0))
  {
    // select randomly the first centroids
    firstcentroids (cent);


    // A. Classification process, ngroups
    // ===================================

    niter = 0; 
    finish = 0;
    while ((finish == 0) && (niter < MAXIT))
    {
      // Obtain the closest group or cluster for each element
      closestgroup (nelems, elems, cent, grind);

      // Calculate new centroids and decide to finish or not depending on DELTA
 
      finish = newcentroids (elems, cent, grind, nelems);

      niter ++;
    } 


    // B. Evaluation of the partition
    // ===================================
   
    for (i=0; i<ngroups; i++) iingrs[i].size = 0;

    // number of elements and elements of each group

    for (i=0; i<nelems; i++) 
    {
      group = grind[i];
      count = iingrs[group].size;
      iingrs[group].members[count] = i;
      iingrs[group].size ++; 
    }

    // validation process and convergence

     cvi = validation (elems, iingrs, cent, compact);

     diff = cvi - cvi_old;
     if (diff < DELTA2) finish_classif = 1;
     else {
       ngroups += 10;
       cvi_old = cvi;
     }
   }

   clock_gettime (CLOCK_REALTIME, &t3);


  // PHASE 2: analyse diseases
  // =========================

    diseases (iingrs, dise, disepro);

    clock_gettime (CLOCK_REALTIME, &t4);



  // write results in a file
  // =======================
  
  f2 = fopen ("results_s.out", "w");
  if (f2 == NULL) {
    printf ("Error when opening file results_s.outs \n");
    exit (-1);
  }
  
  fprintf (f2, " >> Centroids of groups \n\n");
  for (i=0; i<ngroups; i++) {
    for (j=0; j<NFEAT; j++) fprintf (f2, "%7.3f", cent[i][j]);
    fprintf (f2,"\n");
  }

  fprintf (f2, "\n >> Size of the groups: %d groups \n\n", ngroups);
  for (i=0; i<ngroups/10; i++) {
    for (j=0; j<10; j++) fprintf (f2, "%9d", iingrs[10*i+j].size);
    fprintf(f2, "\n");
  }
  for (i=i*10; i<ngroups; i++) fprintf (f2, "%9d", iingrs[i].size);
  fprintf (f2, "\n");

  fprintf (f2, "\n >> Group compactness \n\n");
  for (i=0; i<ngroups/10; i++) {
    for (j=0; j<10; j++) fprintf (f2, "%9.2f", compact[10*i+j]);
    fprintf(f2, "\n");
  }
  for (i=i*10; i<ngroups; i++) fprintf (f2, "%9.2f", compact[i]);
  fprintf (f2, "\n");

  fprintf (f2, "\n\n Analysis of diseases (medians)\n\n");
  fprintf (f2, "\n Dise.  M_max - Group   M_min - Group");
  fprintf (f2, "\n ==================================\n");
  for (i=0; i<TDISEASE; i++)
    fprintf (f2, "  %2d     %4.2f - %2d      %4.2f - %2d\n", i, disepro[i].mmax,
 			disepro[i].gmax, disepro[i].mmin, disepro[i].gmin);

  fclose (f2);

  clock_gettime (CLOCK_REALTIME, &t5);


  // print some results in the screen
  // =================================

  t_read    = (t2.tv_sec-t1.tv_sec) + (t2.tv_nsec-t1.tv_nsec) / (double)1e9;
  t_clus    = (t3.tv_sec-t2.tv_sec) + (t3.tv_nsec-t2.tv_nsec) / (double)1e9;
  t_anal    = (t4.tv_sec-t3.tv_sec) + (t4.tv_nsec-t3.tv_nsec) / (double)1e9;
  t_write   = (t5.tv_sec-t4.tv_sec) + (t5.tv_nsec-t4.tv_nsec) / (double)1e9;

  printf ("\n    Number of iterations: %d", niter);
  printf ("\n    T_read:    %6.3f s", t_read);
  printf ("\n    T_clus:    %6.3f s", t_clus);
  printf ("\n    T_anal:    %6.3f s", t_anal);
  printf ("\n    T_write:   %6.3f s", t_write);
  printf ("\n    ========================");
  printf ("\n    T_total:  %6.3f s\n\n", t_read + t_clus + t_anal + t_write);


  printf ("\n >> Centroids 0, 30 and 60 \n ");
  for (j=0; j<NFEAT; j++) printf ("%7.3f", cent[0][j]);
  printf("\n");
  for (j=0; j<NFEAT; j++) printf ("%7.3f", cent[20][j]);
  printf("\n");
  for (j=0; j<NFEAT; j++) printf ("%7.3f", cent[40][j]);
  printf("\n");

  printf ("\n >> Size of the groups: %d groups \n\n", ngroups);
  for (i=0; i<ngroups/10; i++) {
    for (j=0; j<10; j++) printf ("%9d", iingrs[10*i+j].size);
    printf("\n");
  }
  for (i=i*10; i<ngroups; i++) printf ("%9d", iingrs[i].size);
  printf ("\n");

  printf ("\n >> Group compactness \n\n");
  for (i=0; i<ngroups/10; i++) {
    for (j=0; j<10; j++) printf ("%9.2f", compact[10*i+j]);
    printf("\n");
  }
  for (i=i*10; i<ngroups; i++) printf ("%9.2f", compact[i]);
  printf ("\n");


  printf ("\n\n Analysis of diseases (medians)\n\n");
  printf ("\n Dise.  M_max - Group   M_min - Group");
  printf ("\n ==================================\n");
  for (i=0; i<TDISEASE; i++)
    printf ("  %2d     %4.2f - %2d      %4.2f - %2d\n", i, disepro[i].mmax,
 			disepro[i].gmax, disepro[i].mmin, disepro[i].gmin);

  printf("\n");
}

