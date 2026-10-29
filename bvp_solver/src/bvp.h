/* Native C11 port of BVP_SOLVER; upstream notices are in LICENSE. */
#ifndef NA_BVP_H
#define NA_BVP_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    BVP_OK = 0,
    BVP_INVALID_INPUT,
    BVP_NONFINITE,
    BVP_SINGULAR,
    BVP_NEWTON_FAILED,
    BVP_MESH_LIMIT,
    BVP_REFINEMENT_LIMIT,
    BVP_OUT_OF_DOMAIN,
    BVP_IO_ERROR,
    BVP_NO_MEMORY
} bvp_status;
typedef enum { BVP_MIRK2 = 2, BVP_MIRK4 = 4, BVP_MIRK6 = 6 } bvp_method;
typedef enum {
    BVP_DEFECT = 0,
    BVP_GLOBAL,
    BVP_DEFECT_AND_GLOBAL,
    BVP_HYBRID,
    BVP_DEFECT_AND_DERIVATIVE
} bvp_error_control;
/* All arrays are row-major. Callbacks must fill every output. Optional Jacobian
 * callbacks return nonzero if supplied, zero to request forward differences. */
typedef struct {
    size_t dimension, parameter_count, left_conditions;
    int is_affine;
    void *user;
    void (*rhs)(void *, double, const double *, const double *, double *);
    void (*boundary)(void *, const double *, const double *, const double *, double *, double *);
    int (*jacobian)(void *, double, const double *, const double *, double *);
    int (*boundary_jacobian)(void *, const double *, const double *, const double *, double *,
                             double *);
    const double *singular_term; /* optional dimension by dimension matrix */
} bvp_problem;
typedef struct {
    bvp_method method;
    double tolerance;
    bvp_error_control error_control;
    size_t max_intervals, max_newton_iterations, max_refinements;
    int estimate_error, estimate_condition, trace;
} bvp_options;
typedef struct {
    size_t intervals, newton_iterations, mesh_refinements, rhs_evaluations, jacobian_evaluations;
    double defect_norm, derivative_error, boundary_residual, scaled_boundary_residual;
    double global_error, richardson_error, higher_order_error, deferred_correction_error;
    double condition_estimate, conditioning_error; /* unavailable diagnostics are NAN */
} bvp_statistics;
typedef struct {
    size_t points, dimension, parameter_count;
    double *mesh, *values, *parameters;
} bvp_initial_guess;
typedef struct {
    size_t points, dimension, parameter_count, stage_count;
    bvp_method method;
    double *mesh, *values, *parameters, *augmented, *stages;
    bvp_statistics statistics;
} bvp_solution;
void bvp_default_options(bvp_options *);
const char *bvp_status_string(bvp_status);
bvp_status bvp_linspace(double a, double b, size_t points, double *out);
bvp_status bvp_guess_from_values(size_t points, size_t dimension, size_t parameter_count,
                                 const double *mesh, const double *values, const double *parameters,
                                 bvp_initial_guess *);
bvp_status bvp_guess_from_function(size_t points, size_t dimension, size_t parameter_count,
                                   const double *mesh, const double *parameters,
                                   void (*guess)(void *, double, double *), void *user,
                                   bvp_initial_guess *);
void bvp_guess_free(bvp_initial_guess *);
bvp_status bvp_solve(const bvp_problem *, const bvp_initial_guess *, const bvp_options *,
                     bvp_solution *);
void bvp_solution_free(bvp_solution *);
bvp_status bvp_evaluate(const bvp_solution *, double x, double *out);
bvp_status bvp_derivative(const bvp_solution *, double x, double *out);
bvp_status bvp_solution_initial_guess(const bvp_solution *, bvp_initial_guess *);
bvp_status bvp_extend(const bvp_solution *, double left, double right, int order,
                      bvp_initial_guess *);
bvp_status bvp_extend_with_values(const bvp_solution *, double left, double right,
                                  const double *left_values, const double *right_values,
                                  bvp_initial_guess *);
bvp_status bvp_save(const bvp_solution *, const char *path);
bvp_status bvp_load(const char *path, bvp_solution *);
size_t bvp_weights(bvp_method, double tau, double weights[8], double derivatives[8]);
#ifdef __cplusplus
}
#endif
#endif
