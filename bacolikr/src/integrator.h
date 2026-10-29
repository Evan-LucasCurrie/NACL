#ifndef BAC_INTEGRATOR_H
#define BAC_INTEGRATOR_H
#include "linear.h"
#include <stddef.h>

/* Implicit DAE integration. All arrays are owned
 * by the integrator, except callback context and input option tolerances. */
typedef struct {
    size_t n;
    void *ctx;
    /* cj is NAN for an ordinary (unscaled) residual evaluation. */
    void (*residual)(void *, double, const double *, const double *, double, double *);
    /* Return 1 if supplied, 0 for finite differences, -1 on callback failure.
     * band is non-NULL for compact assembly; otherwise dense is row major. */
    int (*jacobian)(void *, double, const double *, const double *, double, BacBand *, double *);
    int banded;
    size_t lower, upper, npde, kcol, nint;
} DaeSystem;
typedef DaeSystem BacDae;
enum { BI_VERIFIED = 0, BI_DASKR = 1, BI_MAXIMUM = 0, BI_RMS = 1 };
typedef struct {
    const double *atol, *rtol;
    size_t atol_len, rtol_len;
    double absolute, relative, initial_step, max_step, min_step;
    size_t max_steps, max_order;
    int method, error_norm;
} TimeOptions;
typedef struct {
    size_t accepted_steps, rejected_steps, residual_evaluations, jacobian_evaluations;
    size_t newton_iterations, linear_solves, error_test_failures, convergence_failures,
        maximum_order;
} IntegrationStatistics;
typedef struct {
    double old_time, new_time, step_size;
    size_t order;
    double error_norm;
} StepOutcome;
struct BiNewtonCache;
typedef struct {
    size_t n;
    double t, *y, *yp, *phi, *source_states;
    size_t source_count;
    TimeOptions options;
    double *atol_owned, *rtol_owned;
    IntegrationStatistics statistics;
    double psi[7], alpha[7], beta[7], gamma[7], sigma[7];
    double previous_cj, newton_convergence_factor, h, hold, last_start;
    size_t order, last_order, stable_steps;
    int startup, initialized, has_last_start;
    double published_time, *published_y, *published_yp;
    int published;
    double root_cursor, last_event_time, initial_root_time, initial_root_h;
    int has_root_cursor, has_last_event, has_initial_root;
    double *initial_root_y, *initial_root_yp;
    size_t *last_event_indices, last_event_count;
    /* Immutable shared cache; matrix/factor below are borrowed views. */
    struct BiNewtonCache *newton_cache;
    BacBand matrix;
    BacLu factor;
    /* Private work arrays, retained across steps and never shared by clones. */
    double *step_work, *newton_work;
    size_t step_work_len, newton_work_len;
    int has_matrix, matrix_banded;
    size_t matrix_lower, matrix_upper, matrix_step;
    double matrix_cj;
    char error[256];
} BdfIntegrator;
enum { BI_ANY = 0, BI_RISING = 1, BI_FALLING = -1 };
typedef struct {
    double time, *state, *derivative;
    size_t count, *indices;
    int *directions;
} RootEvent;
typedef int (*BiRootFunction)(void *, double, const double *, const double *, double *, size_t);
typedef int (*BiProjection)(void *, const double *, size_t, double *, size_t);
TimeOptions bi_options(double atol, double rtol);
TimeOptions bi_daskr_options(double atol, double rtol);
int bi_validate_options(const TimeOptions *, size_t n, char *error, size_t error_size);
int bi_new(BdfIntegrator *, double t, const double *y, const double *yp, size_t n,
           const TimeOptions *);
void bi_free(BdfIntegrator *);
int bi_clone(BdfIntegrator *destination, const BdfIntegrator *source);
/* Replace an initialized (or zero-filled) receiver, reusing compatible buffers. */
int bi_clone_from(BdfIntegrator *destination, const BdfIntegrator *source);
int bi_reset(BdfIntegrator *, double, const double *, const double *, size_t);
void bi_invalidate_jacobian(BdfIntegrator *);
int bi_initialize_consistent(BdfIntegrator *, const DaeSystem *, const unsigned char *differential);
int bi_step(BdfIntegrator *, const DaeSystem *, double bound, int clip_at_bound, StepOutcome *);
int bi_interpolate(const BdfIntegrator *, double, double *y, double *yp);
int bi_publish_dense_output(BdfIntegrator *, double);
void bi_clear_published_output(BdfIntegrator *);
double bi_time(const BdfIntegrator *);
const double *bi_state(const BdfIntegrator *);
const double *bi_derivative(const BdfIntegrator *);
void bi_clear_event_history(BdfIntegrator *);
/* Returns 1 for an event, 0 for none, -1 for error. */
int bi_locate_roots(BdfIntegrator *, BiRootFunction, void *, size_t, const int *, double output,
                    RootEvent *);
void bi_root_free(RootEvent *);
/* Returns 1 on transfer, 0 without accepted history, -1 on error. */
int bi_project_history(BdfIntegrator *, const BdfIntegrator *, size_t, const TimeOptions *,
                       BiProjection, void *);
int bi_restore_source_spare(BdfIntegrator *, const double *, size_t);
#endif
