#include "bacolikr.h"
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#define PI 3.14159265358979323846
static void initial(void *ctx, double x, double *u) {
    (void)ctx;
    u[0] = sin(PI * x);
}
static void rhs(void *ctx, double t, double x, const double *u, const double *ux, const double *uxx,
                double *out) {
    (void)ctx;
    (void)t;
    (void)x;
    (void)u;
    (void)ux;
    out[0] = uxx[0];
}
static void boundary(void *ctx, int side, double t, const double *u, const double *ux,
                     double *out) {
    (void)ctx;
    (void)side;
    (void)t;
    (void)ux;
    out[0] = u[0];
}
static int jacobian(void *ctx, double t, double x, const double *u, const double *ux,
                    const double *uxx, double *a, double *b, double *c) {
    (void)ctx;
    (void)t;
    (void)x;
    (void)u;
    (void)ux;
    (void)uxx;
    a[0] = b[0] = 0.;
    c[0] = 1.;
    return 1;
}
static int boundary_jacobian(void *ctx, int side, double t, const double *u, const double *ux,
                             double *a, double *b, double *c) {
    (void)ctx;
    (void)side;
    (void)t;
    (void)u;
    (void)ux;
    a[0] = 1.;
    b[0] = c[0] = 0.;
    return 1;
}
static double mass(void *ctx, double x, const double *u, const double *ux, const double *uxx,
                   const double *ut) {
    (void)ctx;
    (void)x;
    (void)ux;
    (void)uxx;
    (void)ut;
    return u[0];
}
static void *read_view(void *ctx) {
    const BacView *v = ctx;
    for (size_t pass = 0; pass < 3; pass++) {
        for (size_t i = 0; i <= 1100; i++) {
            double x = (double)i / 1100., u;
            assert(!bac_view_evaluate(v, x, &u, NULL, NULL, NULL));
            assert(fabs(u - exp(-PI * PI * v->time) * sin(PI * x)) < 4e-6);
        }
        assert(fabs(bac_view_integrate(v, mass, NULL) - 2. / PI * exp(-PI * PI * v->time)) < 3e-6);
    }
    return NULL;
}
int main(void) {
    BacProblem problem = {.components = 1,
                          .initial = initial,
                          .rhs = rhs,
                          .boundary = boundary,
                          .rhs_jacobian = jacobian,
                          .boundary_jacobian = boundary_jacobian};
    for (int source = 0; source < 2; source++) {
        BacOptions options = bac_options();
        options.spatial_adaptivity = 0;
        if (source)
            options.time = bi_daskr_options(1e-6, 1e-6);
        BacSolver solver;
        assert(!bac_solver_new(&solver, &problem, &options, 0.));
        assert(!solver.grid.quadrature);
        double *mesh = solver.grid.mesh;
        BacBasis *basis = solver.grid.collocation_basis;
        for (size_t output = 1; output <= 3; output++) {
            BacAdvanceResult result;
            assert(!bac_solver_advance(&solver, .05 * (double)output, &result));
            bac_advance_result_free(&result);
            assert(solver.grid.mesh == mesh && solver.grid.collocation_basis == basis);
            assert(!solver.grid.quadrature);
        }
        BacView view = bac_solver_view(&solver);
        pthread_t threads[4];
        for (size_t i = 0; i < 4; i++)
            assert(!pthread_create(threads + i, NULL, read_view, &view));
        for (size_t i = 0; i < 4; i++)
            assert(!pthread_join(threads[i], NULL));
        const BacQuadrature *quadrature = solver.grid.quadrature;
        assert(quadrature && bac_grid_quadrature(&solver.grid) == quadrature);
        bac_solver_free(&solver);
    }
    BacOptions options = bac_options();
    options.time = bi_daskr_options(1e-10, 1e-10);
    options.absolute = options.relative = 1e-10;
    options.intervals = 2;
    options.max_remeshes = 1;
    BacSolver solver;
    assert(!bac_solver_new(&solver, &problem, &options, 0.));
    BacAdvanceResult result;
    assert(bac_solver_advance(&solver, .1, &result));
    BacView view = bac_solver_view(&solver);
    double u;
    assert(!bac_view_evaluate(&view, .5, &u, NULL, NULL, NULL));
    assert(isfinite(u) && isfinite(bac_view_integrate(&view, mass, NULL)));
    bac_advance_result_free(&result);
    bac_solver_free(&solver);
    puts("Solver snapshots: grid reuse, lazy/concurrent view caches and failed-remesh ownership "
         "passed");
    return 0;
}
