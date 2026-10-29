#ifndef NA_BVP_INTERNAL_H
#define NA_BVP_INTERNAL_H
#include "bvp.h"
#include <float.h>
#include <math.h>
#include <setjmp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define CLAMP(a, l, h) MAX((l), MIN((a), (h)))
#define SWAP(T, a, b)                                                                              \
    do {                                                                                           \
        T tmp_ = (a);                                                                              \
        (a) = (b);                                                                                 \
        (b) = tmp_;                                                                                \
    } while (0)
typedef struct allocation {
    void *ptr;
    size_t bytes;
    struct allocation *next;
} allocation;
typedef struct {
    const bvp_problem *p;
    size_t n, q, left;
    double a;
    double *endpoint_inverse, *projection;
    bvp_statistics stats;
    bvp_status error;
    jmp_buf *trap;
    allocation *memory;
    allocation *scratch_cache;
    size_t cached_bytes;
} context;
void *bvp_alloc(context *, size_t, size_t);
void bvp_release(context *, allocation *mark);
void bvp_detach(context *, void *);
void bvp_fail(context *, bvp_status);
void bvp_finite(context *, const double *, size_t);
double bvp_maximum(const double *, size_t);
typedef struct {
    size_t s, extra;
    double c[9], v[9], b[9], x[9][9], ce[3], ve[3], xe[3][9], defect[2];
} tableau;
void bvp_tableau(int method, tableau *);
typedef struct {
    size_t n, q, left, intervals, lower, upper, width, count;
    int colrow;
    double *a;
    size_t *pivot, *ends;
} matrix;
matrix bvp_matrix(context *, size_t q, size_t left, size_t points);
double *bvp_matrix_row(matrix *, size_t row, size_t first);
void bvp_factor(context *, matrix *);
void bvp_linear(context *, const matrix *, double *, int transpose);
double *bvp_pinv(context *, const double *, size_t);
#endif
