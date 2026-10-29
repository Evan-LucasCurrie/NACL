#ifndef BACOLIKR_H
#define BACOLIKR_H
#include "spatial.h"
#include "integrator.h"
typedef struct {
    size_t components;
    void *ctx;
    void (*initial)(void *, double, double *);
    void (*rhs)(void *, double, double, const double *, const double *, const double *, double *);
    void (*boundary)(void *, int, double, const double *, const double *, double *);
    int (*rhs_jacobian)(void *, double, double, const double *, const double *, const double *,
                        double *, double *, double *);
    int (*boundary_jacobian)(void *, int, double, const double *, const double *, double *,
                             double *, double *);
    int (*conservation_form)(void *, double, double, const double *, const double *, double *,
                             double *);
    int (*conservation_jacobian)(void *, double, double, const double *, const double *, double *,
                                 double *, double *, double *);
    int (*source_rhs_jacobian)(void *, double, double, const double *, const double *,
                               const double *, double *, double *, double *);
    int (*source_boundary_jacobian)(void *, int, double, const double *, const double *, double *,
                                    double *, double *);
} BacProblem;
typedef struct {
    const double *mesh;
    size_t intervals, kcol;
    const double *atol, *rtol;
    size_t atol_len, rtol_len;
    double absolute, relative;
    BacEstimator estimator;
    size_t max_intervals, max_remeshes;
    int spatial_adaptivity;
    TimeOptions time;
} BacOptions;
typedef struct {
    IntegrationStatistics time;
    size_t remeshes, warm_remeshes, cold_remeshes, spatial_rejections, intervals;
    double spatial_error_ratio;
} BacStatistics;
typedef struct {
    double time;
    const BacGrid *grid;
    size_t npde;
    const double *y, *yp;
} BacView;
typedef double (*BacEventFunction)(void *, const BacView *);
typedef struct {
    double time;
    int event;
    size_t count, *indices;
    int *directions;
} BacAdvanceResult;
typedef struct {
    BacProblem problem;
    BacOptions options;
    BacGrid grid;
    BdfIntegrator integrator;
    double *atol, *rtol;
    BacStatistics statistics;
    char error[256];
    size_t event_count;
    BacEventFunction *events;
    void **event_contexts;
    int *event_directions;
} BacSolver;
BacOptions bac_options(void);
int bac_solver_new(BacSolver *, const BacProblem *, const BacOptions *, double initial_time);
void bac_solver_free(BacSolver *);
int bac_solver_advance(BacSolver *, double tout, BacAdvanceResult *);
int bac_solver_reinitialize(BacSolver *);
BacView bac_solver_view(const BacSolver *);
BacStatistics bac_solver_statistics(const BacSolver *);
int bac_view_evaluate(const BacView *, double x, double *u, double *ux, double *uxx, double *ut);
double bac_view_integrate(const BacView *,
                          double (*)(void *, double, const double *, const double *, const double *,
                                     const double *),
                          void *);
int bac_solver_add_event(BacSolver *, int direction, BacEventFunction, void *);
void bac_solver_clear_events(BacSolver *);
void bac_advance_result_free(BacAdvanceResult *);
#endif
