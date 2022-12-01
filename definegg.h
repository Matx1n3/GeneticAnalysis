/*
   definegg.h
   Constants and structs used in files gengrops_s.c and fungg_s.c
*****************************************************************/

#define MAXELE      230000	//number of elements (samples)
#define NGROUPSMAX  100	  	//number of clusters
#define NFEAT       40		//features of each instance
#define TDISEASE    18		//types of disease

#define DELTA1      0.01	//convergence: minimum change in centroids
#define DELTA2      0.01	//convergence: classification (cvi)
#define MAXIT       1000 	//convergence: maximum number of iterations


extern int ngroups;


struct ginfo               // information about groups
{
 int  members[MAXELE];     // members
 int  size;                // number of elements
};


struct analysis            // analysis of diseases
{
 float  mmax, mmin;        // median maximum and minimum for each disease
 int    gmax, gmin;        // grous for maximums and minimums
};

