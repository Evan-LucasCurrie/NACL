#ifndef COLNEWSC_REFERENCE_H
#define COLNEWSC_REFERENCE_H
#include "colnewsc.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Upstream-compatible COLNEWSC policies. Analytic Jacobians are deliberately
 * unused; the original one-sided finite differences are retained.
 * The problem and callback context must outlive the returned solution, whose
 * interpolation can re-evaluate the RHS. Separate solutions own separate state.
 * A single solution is mutable during evaluation and requires external locking. */
typedef struct {
    colnewsc_problem base;
    void (*reference_guess)(void *, double, double *, double *);
    void (*reference_rhs)(void *, double, const double *, double *);
    double (*reference_boundary)(void *, size_t, const double *);
    int preserve_non_finite_interpolation;
} colnewsc_reference_problem;
typedef struct {
    size_t non_finite_interpolation_evaluations;
    int has_first_non_finite_interpolation_point;
    double first_non_finite_interpolation_point;
} colnewsc_reference_diagnostics;
typedef struct colnewsc_reference_solution colnewsc_reference_solution;
colnewsc_status colnewsc_reference_solve(const colnewsc_reference_problem *,
                                         const colnewsc_options *, colnewsc_reference_solution **);
void colnewsc_reference_free(colnewsc_reference_solution *);
size_t colnewsc_reference_mesh_count(const colnewsc_reference_solution *);
const double *colnewsc_reference_mesh(const colnewsc_reference_solution *);
const int *colnewsc_reference_mesh_history(const colnewsc_reference_solution *, size_t *count);
const double *colnewsc_reference_representation(const colnewsc_reference_solution *, size_t *count);
colnewsc_reference_diagnostics
colnewsc_reference_get_diagnostics(const colnewsc_reference_solution *);
colnewsc_status colnewsc_reference_evaluate(colnewsc_reference_solution *, double, double *,
                                            size_t);
#ifdef __cplusplus
}
#endif
#endif
