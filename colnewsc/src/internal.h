#ifndef COLNEWSC_INTERNAL_H
#define COLNEWSC_INTERNAL_H
#include "colnewsc.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CN_MIN(a, b) ((a) < (b) ? (a) : (b))
#define CN_MAX(a, b) ((a) > (b) ? (a) : (b))
#define CN_TRY(x)                                                                                  \
    do {                                                                                           \
        colnewsc_status cn_e = (x);                                                                \
        if (cn_e)                                                                                  \
            return cn_e;                                                                           \
    } while (0)
typedef struct {
    size_t k;
    double nodes[7], a[49], b[7], integrated[7][8];
} cn_basis;
struct colnewsc_solution {
    size_t count, dimension;
    double *mesh, *values, *stages, *extra;
    cn_basis basis;
    colnewsc_statistics stats;
};
typedef struct {
    size_t n, lower, upper, stride;
    int band;
    double *data;
    size_t *pivots;
    void *storage;
    size_t data_capacity, pivot_capacity;
    int borrowed;
} cn_lu;
static inline void *cn_alloc(size_t n, size_t width) {
    return n && width > SIZE_MAX / n ? NULL : calloc(n ? n : 1, width);
}
static inline int cn_finite(const double *a, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (!isfinite(a[i]))
            return 0;
    return 1;
}
void cn_basis_new(cn_basis *, size_t);
void cn_basis_weights(const cn_basis *, double, double *);
colnewsc_status cn_checked_rhs(const colnewsc_problem *, double, const double *, double *);
colnewsc_status cn_prepare(const colnewsc_problem *, colnewsc_solution *);
void cn_evaluate_interval(const colnewsc_solution *, size_t, double, double *);
colnewsc_status cn_band_new(cn_lu *, size_t, size_t, size_t);
colnewsc_status cn_band_reset(cn_lu *, size_t, size_t, size_t);
void cn_lu_free(cn_lu *);
void cn_band_set(cn_lu *, size_t, size_t, double);
colnewsc_status cn_lu_factor(cn_lu *);
colnewsc_status cn_stage_factor(double *, size_t, cn_lu *);
colnewsc_status cn_stage_factor_reuse(const double *, size_t, cn_lu *);
colnewsc_status cn_lu_solve(const cn_lu *, double *, size_t);
#endif
