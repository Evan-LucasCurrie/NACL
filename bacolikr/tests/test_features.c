#include "bacolikr.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#define PI 3.14159265358979323846
typedef struct {
    int robin, conservative, numeric;
} Context;
static void initial(void *ctx, double x, double *out) {
    (void)ctx;
    for (size_t c = 0; c < 2; c++)
        out[c] = (double)c + 2 + ((double)c + 1) * cos(PI * x);
}
static void rhs(void *ctx, double t, double x, const double *u, const double *ux, const double *uxx,
                double *out) {
    (void)ctx;
    (void)t;
    (void)x;
    (void)u;
    (void)ux;
    for (size_t c = 0; c < 2; c++)
        out[c] = uxx[c];
}
static void boundary(void *ctx, int side, double t, const double *u, const double *ux,
                     double *out) {
    Context *c = ctx;
    for (size_t j = 0; j < 2; j++) {
        double val = (double)j + 2 + ((double)j + 1) * exp(-PI * PI * t) * (side ? -1 : 1);
        out[j] = c->robin ? u[j] + ux[j] - val : ux[j];
    }
}
static int rhs_jac(void *ctx, double t, double x, const double *u, const double *ux,
                   const double *uxx, double *a, double *b, double *c) {
    (void)ctx;
    (void)t;
    (void)x;
    (void)u;
    (void)ux;
    (void)uxx;
    memset(a, 0, 4 * sizeof(double));
    memset(b, 0, 4 * sizeof(double));
    memset(c, 0, 4 * sizeof(double));
    c[0] = c[3] = 1;
    return 1;
}
static int boundary_jac(void *ctx, int side, double t, const double *u, const double *ux, double *a,
                        double *b, double *dt) {
    (void)u;
    (void)ux;
    Context *c = ctx;
    memset(a, 0, 4 * sizeof(double));
    memset(b, 0, 4 * sizeof(double));
    a[0] = a[3] = c->robin ? 1 : 0;
    b[0] = b[3] = 1;
    for (size_t j = 0; j < 2; j++)
        dt[j] = c->robin ? ((double)j + 1) * PI * PI * exp(-PI * PI * t) * (side ? -1 : 1) : 0;
    return 1;
}
static int flux(void *ctx, double t, double x, const double *u, const double *ux, double *q,
                double *source) {
    (void)t;
    (void)x;
    (void)u;
    Context *c = ctx;
    if (!c->conservative)
        return 0;
    for (size_t j = 0; j < 2; j++) {
        q[j] = ux[j];
        source[j] = 0;
    }
    return 1;
}
static int fluxjac(void *ctx, double t, double x, const double *u, const double *ux, double *a,
                   double *b, double *c, double *d) {
    (void)t;
    (void)x;
    (void)u;
    (void)ux;
    Context *k = ctx;
    if (k->numeric)
        return 0;
    memset(a, 0, 4 * sizeof(double));
    memset(b, 0, 4 * sizeof(double));
    memset(c, 0, 4 * sizeof(double));
    memset(d, 0, 4 * sizeof(double));
    b[0] = b[3] = 1;
    return 1;
}
static double mass(void *ctx, double x, const double *u, const double *ux, const double *uxx,
                   const double *ut) {
    (void)x;
    (void)ux;
    (void)uxx;
    (void)ut;
    return u[(size_t)(uintptr_t)ctx];
}
static BacProblem problem(Context *ctx) {
    BacProblem p = {0};
    p.components = 2;
    p.ctx = ctx;
    p.initial = initial;
    p.rhs = rhs;
    p.boundary = boundary;
    p.rhs_jacobian = rhs_jac;
    p.boundary_jacobian = boundary_jac;
    p.conservation_form = flux;
    p.conservation_jacobian = fluxjac;
    return p;
}
int main(void) {
    for (int mode = 0; mode < 3; mode++) {
        Context ctx = {mode == 2, mode != 2, mode == 1};
        BacProblem p = problem(&ctx);
        BacOptions o = bac_options();
        o.absolute = o.relative = 1e-7;
        BacSolver s;
        assert(!bac_solver_new(&s, &p, &o, 0));
        BacAdvanceResult r;
        assert(!bac_solver_advance(&s, .1, &r));
        BacView v = bac_solver_view(&s);
        double error = 0;
        for (size_t i = 0; i <= 100; i++) {
            double x = (double)i / 100, u[2];
            assert(!bac_view_evaluate(&v, x, u, NULL, NULL, NULL));
            for (size_t c = 0; c < 2; c++) {
                double exact = (double)c + 2 + ((double)c + 1) * exp(-PI * PI * .1) * cos(PI * x);
                error = fmax(error, fabs(u[c] - exact));
            }
        }
        assert(error < 3e-6);
        if (ctx.conservative)
            for (size_t c = 0; c < 2; c++) {
                double value = bac_view_integrate(&v, mass, (void *)(uintptr_t)c);
                assert(fabs(value - ((double)c + 2)) < 1e-10);
            }
        for (int side = 0; side < 2; side++) {
            double u[2], ux[2], value[2];
            assert(!bac_view_evaluate(&v, (double)side, u, ux, NULL, NULL));
            boundary(&ctx, side, .1, u, ux, value);
            assert(fabs(value[0]) < 1e-10 && fabs(value[1]) < 1e-10);
        }
        bac_advance_result_free(&r);
        bac_solver_free(&s);
    }
    Context ctx = {0, 0, 0};
    BacProblem p = problem(&ctx);
    BacOptions o = bac_options();
    o.time.max_steps = 2;
    BacSolver s;
    assert(!bac_solver_new(&s, &p, &o, 0));
    BacAdvanceResult r;
    assert(bac_solver_advance(&s, .1, &r));
    double t = bi_time(&s.integrator);
    assert(t > 0 && t < .1);
    double state = s.integrator.y[4];
    assert(bac_solver_advance(&s, .1, &r));
    assert(bi_time(&s.integrator) == t && s.integrator.y[4] == state);
    assert(!bac_solver_reinitialize(&s));
    bac_solver_free(&s);
    o.kcol = 2;
    assert(bac_solver_new(&s, &p, &o, 0));
    bac_solver_free(&s);
    puts("PDE systems, Neumann/Robin, analytic/numerical conservative Jacobians, mass, rollback "
         "and invalid options passed");
    return 0;
}
