#ifndef BAC_LINEAR_H
#define BAC_LINEAR_H
#include <stddef.h>
typedef struct {
    size_t n, lower, upper, stride;
    double *values;
    size_t npde, kcol, nint;
} BacBand;
struct BacAbd;
typedef struct {
    size_t n, lower, upper, stride;
    /* Generic factors retain the pivot-expanded column workspace. Active
     * lengths trim trailing exact zeros without copying the factor matrix. */
    double *values;
    size_t *pivots;
    struct BacAbd *abd;
    size_t lower_count, upper_count;
    size_t *lower_lengths, *upper_lengths;
    int dense;
} BacLu;
int bac_band_new(BacBand *, size_t n, size_t lower, size_t upper);
void bac_band_free(BacBand *);
double bac_band_get(const BacBand *, size_t row, size_t col);
void bac_band_set(BacBand *, size_t row, size_t col, double);
void bac_band_add(BacBand *, size_t row, size_t col, double);
void bac_band_row_norms(const BacBand *, const double *weights, double *out);
int bac_band_factor(const BacBand *, BacLu *);
int bac_dense_factor(const double *, size_t n, BacLu *);
int bac_lu_solve(const BacLu *, double *rhs);
void bac_lu_free(BacLu *);
int bac_lu_clone(BacLu *, const BacLu *);
size_t bac_abd_workspace(const BacLu *, double *);
#endif
