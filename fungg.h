/*
   fungg.h
   headers of the functions used in gengroups_s.c
***************************************************************/

extern void	firstcentroids
		(float cent[][NFEAT]);

extern int	newcentroids
                (float elems[][NFEAT],
                 float cent[][NFEAT],
                 int   *grind,
                 int   nelems
                );

extern double	geneticdistance
                (float *elem1,
                 float *elem2
                );

extern void	closestgroup
                (int   nelems,
                 float elems[][NFEAT],
                 float cent[][NFEAT],
                 int   *grind
                );

extern double	validation
                (float  elems[][NFEAT],
                 struct ginfo *iingrs,
                 float  cent[][NFEAT],
                 float  *compact
                );

extern void	diseases
                (struct ginfo *iingrs,
                 float  dise[][TDISEASE],
                 struct analysis *disepro
                );


