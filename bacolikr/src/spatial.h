#ifndef BAC_SPATIAL_H
#define BAC_SPATIAL_H
#include "linear.h"
typedef enum { BAC_LOI, BAC_SCI } BacEstimator;
typedef struct {
    size_t start, len;
    double b[13], bx[13], bxx[13];
} BacBasis;
typedef struct {
    double ratio, distribution_ratio;
    double *indicators, *component_ratios;
    size_t intervals, npde;
    BacEstimator estimator;
} BacSpatialError;
typedef struct BacEstimatorGeometry BacEstimatorGeometry;
typedef struct BacViewCache BacViewCache;
typedef struct {
    double x, weight;
    BacBasis basis;
} BacQuadrature;
typedef struct {
    size_t nint, kcol, n, nknots;
    double *mesh, *knots, *collocation;
    BacBasis *collocation_basis;
    BacLu interpolation;
    BacEstimatorGeometry *geometry[2];
    BacBasis *flux_basis;
    double *flux_derivatives;
    BacQuadrature *quadrature;
    BacViewCache *view_cache;
    double *boundary_lift;
    double lift_values[2][2], lift_derivatives[2][2];
} BacGrid;
int bac_gauss(size_t n, int weighted, double *points, double *weights);
int bac_grid_new(BacGrid *, const double *mesh, size_t intervals, size_t kcol);
void bac_grid_free(BacGrid *);
int bac_grid_basis(const BacGrid *, double x, BacBasis *);
int bac_grid_cached_basis(const BacGrid *, double x, BacBasis *);
const BacQuadrature *bac_grid_quadrature(const BacGrid *);
int bac_grid_flux_geometry(BacGrid *);
void bac_basis_sum(const BacBasis *, const double *coeff, size_t npde, int derivative, double *out);
int bac_grid_evaluate(const BacGrid *, const double *coeff, size_t npde, double x, int derivative,
                      double *out);
int bac_grid_interpolate(const BacGrid *, size_t npde, void (*function)(void *, double, double *),
                         void *ctx, double *coeff);
int bac_grid_project(const BacGrid *target, const BacGrid *source, const double *coeff, size_t npde,
                     double *out);
int bac_grid_estimate(BacGrid *, const double *coeff, size_t npde, const double *atol,
                      const double *rtol, BacEstimator, BacSpatialError *);
void bac_spatial_error_free(BacSpatialError *);
int bac_grid_refined_mesh(const BacGrid *, const BacSpatialError *, size_t max, size_t prescribed,
                          int capacity_trial, double **mesh, size_t *intervals);
int bac_grid_projection_error(const BacGrid *target, const BacGrid *source, const double *original,
                              const double *projected, size_t npde, const double *atol,
                              const double *rtol, double *ratio, double *maximum);
typedef struct {
    double ratio, maximum;
    size_t intervals;
    double *interval_maximum, *interval_squared;
} BacProjectionError;
int bac_grid_projection_measure(const BacGrid *, const BacGrid *, const double *, const double *,
                                size_t, const double *, const double *, BacProjectionError *);
void bac_projection_error_free(BacProjectionError *);
int bac_grid_refine_projection(const BacGrid *, const BacGrid *, const BacProjectionError *, size_t,
                               double **, size_t *);
typedef struct {
    const BacGrid *source, *target;
    size_t npde;
    BacBasis *basis;
    BacLu factor;
} BacSourceProjection;
int bac_source_projection_new(BacSourceProjection *, const BacGrid *target, const BacGrid *source,
                              size_t);
int bac_source_projection_apply(const BacSourceProjection *, const double *, double *);
void bac_source_projection_free(BacSourceProjection *);
#endif
