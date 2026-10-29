#ifndef COLNEWSC_H
#define COLNEWSC_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Native C11 adaptive Gaussian collocation solver.
   All matrices are row major. Callback output arrays must be fully initialized.
   Optional jacobian/gradient return 0 to request centered finite differences. */
typedef enum {
    COLNEWSC_OK = 0,
    COLNEWSC_INVALID_INPUT,
    COLNEWSC_NON_FINITE,
    COLNEWSC_SINGULAR_JACOBIAN,
    COLNEWSC_NEWTON_NOT_CONVERGED,
    COLNEWSC_MESH_LIMIT,
    COLNEWSC_REFINEMENT_LIMIT,
    COLNEWSC_OUTSIDE_INTERVAL,
    COLNEWSC_ALLOCATION_FAILED
} colnewsc_status;
typedef struct {
    size_t dimension, left_boundary_conditions;
    double a, b;
    void *context;
    void (*rhs)(void *, double, const double *, double *);
    double (*boundary)(void *, size_t, const double *);
    void (*initial_guess)(void *, double, double *);
    int (*jacobian)(void *, double, const double *, double *);
    int (*boundary_gradient)(void *, size_t, const double *, double *);
} colnewsc_problem;
typedef struct {
    double tolerance;
    const double *component_tolerances;
    size_t component_tolerance_count;
    size_t collocation_points, initial_subintervals, max_subintervals;
    size_t max_newton_iterations, max_refinements;
    int trace;
    const double *initial_mesh;
    size_t initial_mesh_count;
    const double *fixed_points;
    size_t fixed_point_count;
} colnewsc_options;
typedef struct {
    size_t newton_iterations, refinements, subintervals;
    double *estimated_errors;
    double estimated_error_ratio, boundary_residual;
} colnewsc_statistics;
typedef struct colnewsc_solution colnewsc_solution;
colnewsc_options colnewsc_default_options(void);
const char *colnewsc_status_string(colnewsc_status);
colnewsc_status colnewsc_solve(const colnewsc_problem *, const colnewsc_options *,
                               colnewsc_solution **);
colnewsc_status colnewsc_solve_with_initial_solution(const colnewsc_problem *,
                                                     const colnewsc_options *,
                                                     const colnewsc_solution *,
                                                     colnewsc_solution **);
void colnewsc_free(colnewsc_solution *);
size_t colnewsc_dimension(const colnewsc_solution *);
size_t colnewsc_mesh_count(const colnewsc_solution *);
const double *colnewsc_mesh(const colnewsc_solution *);
const double *colnewsc_mesh_values(const colnewsc_solution *);
const colnewsc_statistics *colnewsc_stats(const colnewsc_solution *);
colnewsc_status colnewsc_evaluate(const colnewsc_solution *, double, double *, size_t);
#ifdef __cplusplus
}
#endif
#endif
