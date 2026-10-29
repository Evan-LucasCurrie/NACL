#ifndef BAC_ABD_H
#define BAC_ABD_H
#include "linear.h"
typedef struct BacAbd {
    size_t p, m, nblocks;
    double *top, *blocks, *bottom;
    size_t *pivots;
} BacAbd;
int bac_abd_factor(const BacBand *, BacAbd **);
int bac_abd_solve(const BacAbd *, double *);
int bac_abd_clone(BacAbd **, const BacAbd *);
void bac_abd_free(BacAbd *);
/* Returns source LAMDEC workspace length; optionally copies that workspace. */
size_t bac_abd_workspace(const BacLu *, double *);
#endif
