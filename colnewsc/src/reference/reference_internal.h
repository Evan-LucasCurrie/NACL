#ifndef COLNEWSC_REFERENCE_INTERNAL_H
#define COLNEWSC_REFERENCE_INTERNAL_H
#include "../colnewsc_reference.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef double *R;
typedef int *I;
#define GET(a, i) ((a)[(i) - 1])
#define SET(a, i, x) ((a)[(i) - 1] = (x))
#define GET2(a, i, j, rows) GET((a), (i) + ((j) - 1) * (rows))
#define SET2(a, i, j, rows, x) SET((a), (i) + ((j) - 1) * (rows), (x))
#define VIEW(a, i) ((a) + (i) - 1)
#define CN_MIN(a, b) ((a) < (b) ? (a) : (b))
#define CN_MAX(a, b) ((a) > (b) ? (a) : (b))
typedef struct cnref_allocation {
    void *data;
    struct cnref_allocation *next;
} cnref_allocation;
typedef struct cnref_context {
    const colnewsc_reference_problem *problem;
    double precis;
    int iprint;
    R rho, coef;
    int k, ncomp, mstar, kd, mmax;
    I m;
    int n, nold, nmax, nz, ndmz, mshflg, mshnum, mshlmt, mshalt;
    R zeta;
    double aleft, aright;
    int izeta, idum, nonlin, iter, limit, icare, iguess;
    R wgtmsh, wgterr, tolin, root;
    I jtol, ltol;
    int ntol;
    R b, acol, asave;
    int switch_;
    R y_meshnew, y_meshold;
    int ykstore, yknew;
    R zold, dmzold, meshold, y_super;
    R sk_discrete_super, sk_interp_super, c_star, v_star, x_star;
    int s_star, imaxe;
    int *mesh_history;
    size_t mesh_history_count, mesh_history_capacity;
    int invalid_evaluation;
    colnewsc_reference_diagnostics diagnostics;
    size_t maximum_corrections, maximum_refinements;
    cnref_allocation *allocations;
    int allocation_failed;
} cnref_context;
/* `switch` is a C keyword, so translated accessors use switch_. */
struct colnewsc_reference_solution {
    cnref_context *context;
    R fspace;
    I ispace;
    R mesh, output;
    size_t mesh_count;
};
void cnref_fsub(cnref_context *, double, R, R);
void cnref_fsub_interpolant(cnref_context *, double, R, R);
double cnref_gsub(cnref_context *, int, R);
void cnref_guess(cnref_context *, double, R, R);
void cnref_dfdy(cnref_context *, double, R, R);
void cnref_dbcdy(cnref_context *, int, R, R);
void cnref_history(cnref_context *, int);
int cnref_valid_mesh(R, int);
double cnref_powi(double, int);
int cnref_all_finite(cnref_context *, R, int);
int cnref_all_pair_finite(cnref_context *, R, R, int);
int cnref_check_mesh(cnref_context *, R, int);
double cnref_mesh_density(cnref_context *, double, double, int, double, double);
void cnref_dgefa(cnref_context *, R, int, int, I, int *);
void cnref_dgesl(cnref_context *, R, int, int, I, R, int);
void cnref_daxpy(cnref_context *, int, double, R, int, R, int);
double cnref_ddot(cnref_context *, int, R, int, R, int);
void cnref_dscal(cnref_context *, int, double, R, int);
int cnref_idamax(cnref_context *, int, R, int);
void cnref_dcopy(cnref_context *, int, R, int, R, int);
int cnref_plain_finite(R, int);
void cnref_interval(cnref_context *, int, R, double, int *);
#include "methods.h"
#endif
