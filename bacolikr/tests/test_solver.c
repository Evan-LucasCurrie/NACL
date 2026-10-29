#include "bacolikr.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI 3.14159265358979323846
static void initial(void *ctx, double x, double *v) {
    (void)ctx;
    v[0] = sin(PI * x);
}
static void rhs(void *ctx, double t, double x, const double *u, const double *ux, const double *uxx,
                double *f) {
    (void)ctx;
    (void)t;
    (void)x;
    (void)u;
    (void)ux;
    f[0] = uxx[0];
}
static void bc(void *ctx, int side, double t, const double *u, const double *ux, double *f) {
    (void)ctx;
    (void)side;
    (void)t;
    (void)ux;
    f[0] = u[0];
}
static int jac(void *ctx, double t, double x, const double *u, const double *ux, const double *uxx,
               double *a, double *b, double *c) {
    (void)ctx;
    (void)t;
    (void)x;
    (void)u;
    (void)ux;
    (void)uxx;
    a[0] = b[0] = 0;
    c[0] = 1;
    return 1;
}
static int bcjac(void *ctx, int side, double t, const double *u, const double *ux, double *a,
                 double *b, double *c) {
    (void)ctx;
    (void)side;
    (void)t;
    (void)u;
    (void)ux;
    a[0] = 1;
    b[0] = c[0] = 0;
    return 1;
}
static double root(void *ctx, const BacView *v) {
    (void)ctx;
    double u;
    bac_view_evaluate(v, .5, &u, NULL, NULL, NULL);
    return u - .5;
}
static void poly(void *ctx, double x, double *out) {
    (void)ctx;
    out[0] = 1 + x + x * x;
}
static void spatial_test(void) {
    double mesh[] = {-.2, 0, .3, 1};
    for (size_t k = 2; k <= 10; k++) {
        BacGrid g;
        assert(!bac_grid_new(&g, mesh, 3, k));
        double *c = calloc(g.n, sizeof(double));
        assert(!bac_grid_interpolate(&g, 1, poly, NULL, c));
        for (size_t i = 0; i < 101; i++) {
            double x = -.2 + 1.2 * (double)i / 100, u, ux, uxx;
            assert(!bac_grid_evaluate(&g, c, 1, x, 0, &u));
            assert(!bac_grid_evaluate(&g, c, 1, x, 1, &ux));
            assert(!bac_grid_evaluate(&g, c, 1, x, 2, &uxx));
            assert(fabs(u - (1 + x + x * x)) < 1e-11);
            assert(fabs(ux - (1 + 2 * x)) < 1e-9);
            assert(fabs(uxx - 2) < 1e-7);
        }
        if (k >= 3) {
            double a[] = {1e-6}, r[] = {1e-6};
            for (int est = 0; est < 2; est++) {
                BacSpatialError e;
                assert(!bac_grid_estimate(&g, c, 1, a, r, (BacEstimator)est, &e));
                assert(e.ratio < 1e-5);
                bac_spatial_error_free(&e);
            }
        }
        free(c);
        bac_grid_free(&g);
    }
    puts("spatial kcol2-10 interpolation, derivatives, LOI/SCI passed");
}
int main(void) {
    spatial_test();
    BacProblem p = {0};
    p.components = 1;
    p.initial = initial;
    p.rhs = rhs;
    p.boundary = bc;
    p.rhs_jacobian = jac;
    p.boundary_jacobian = bcjac;
    for (int mode = 0; mode < 2; mode++)
        for (int adapt = 0; adapt < 2; adapt++) {
            BacOptions o = bac_options();
            o.intervals = 5;
            o.absolute = o.relative = 1e-6;
            o.spatial_adaptivity = adapt;
            if (mode)
                o.time = bi_daskr_options(1e-6, 1e-6);
            BacSolver s;
            assert(!bac_solver_new(&s, &p, &o, 0));
            BacAdvanceResult result;
            int rc = bac_solver_advance(&s, .1, &result);
            if (rc)
                fprintf(stderr, "mode%d adapt%d %s\n", mode, adapt, s.error);
            assert(!rc);
            BacView v = bac_solver_view(&s);
            double error = 0;
            for (size_t i = 0; i <= 100; i++) {
                double x = (double)i / 100, u;
                assert(!bac_view_evaluate(&v, x, &u, NULL, NULL, NULL));
                error = fmax(error, fabs(u - exp(-PI * PI * .1) * sin(PI * x)));
            }
            printf("heat mode%d adapt%d error %.9g intervals %zu steps %zu\n", mode, adapt, error,
                   s.grid.nint, bac_solver_statistics(&s).time.accepted_steps);
            assert(error < 5e-5);
            bac_advance_result_free(&result);
            bac_solver_free(&s);
        }
    BacOptions o = bac_options();
    o.absolute = o.relative = 1e-7;
    BacSolver s;
    assert(!bac_solver_new(&s, &p, &o, 0));
    assert(bac_solver_add_event(&s, BI_FALLING, root, NULL) == 0);
    BacAdvanceResult result;
    assert(!bac_solver_advance(&s, .2, &result));
    assert(result.event);
    assert(fabs(result.time - log(2) / (PI * PI)) < 2e-5);
    bac_advance_result_free(&result);
    assert(!bac_solver_advance(&s, .2, &result));
    assert(!result.event);
    bac_advance_result_free(&result);
    bac_solver_free(&s);
    puts("PDE heat, adaptive mesh, source mode, events passed");
    return 0;
}
