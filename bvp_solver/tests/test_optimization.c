#include "internal.h"
#include <assert.h>
#include <stdio.h>

static void allocation_tests(void) {
    context *c = calloc(1, sizeof(*c));
    assert(c);
    jmp_buf trap;
    c->trap = &trap;
    if (!setjmp(trap)) {
        double *anchor = bvp_alloc(c, 1, sizeof(double));
        *anchor = 3.;
        allocation *mark = c->memory;
        double *first = bvp_alloc(c, 64, sizeof(double));
        for (size_t i = 0; i < 64; i++)
            first[i] = NAN;
        bvp_release(c, mark);
        assert(c->scratch_cache && c->cached_bytes == 64 * sizeof(double));
        double *second = bvp_alloc(c, 64, sizeof(double));
        assert(first == second && !c->cached_bytes);
        for (size_t i = 0; i < 64; i++)
            assert(second[i] == 0.);
        allocation *inner = c->memory;
        double *third = bvp_alloc(c, 32, sizeof(double));
        third[0] = 17.;
        bvp_release(c, inner);
        assert(*anchor == 3. && second[0] == 0.);
        double *detached = bvp_alloc(c, 32, sizeof(double));
        detached[0] = 23.;
        bvp_detach(c, detached);
        bvp_release(c, mark);
        /* Retaining scratch is bounded even across unusually large trials. */
        for (size_t i = 0; i < 20; i++)
            (void)bvp_alloc(c, 65536, 1);
        bvp_release(c, mark);
        assert(c->cached_bytes <= (1u << 20));
        (void)bvp_alloc(c, (1u << 20) + 1, 1);
        bvp_release(c, mark);
        assert(c->cached_bytes <= (1u << 20));
        bvp_release(c, NULL);
        assert(!c->memory && !c->scratch_cache && !c->cached_bytes);
        assert(detached[0] == 23.);
        free(detached);
        (void)bvp_alloc(c, 1, sizeof(double));
        mark = c->memory;
        (void)bvp_alloc(c, 128, sizeof(double));
        bvp_release(c, mark);
        (void)bvp_alloc(c, 64, sizeof(double));
        bvp_fail(c, BVP_NONFINITE);
    }
    assert(c->error == BVP_NONFINITE);
    bvp_release(c, NULL);
    assert(!c->memory && !c->scratch_cache && !c->cached_bytes);
    free(c);
}

enum { DIM = 17 };
static void rhs(void *ctx, double x, const double *y, const double *p, double *out) {
    (void)p;
    for (size_t i = 0; i < DIM; i++)
        out[i] = ctx && x > .35 ? NAN : y[i] / (double)(i + 2);
}
static int jacobian(void *ctx, double x, const double *y, const double *p, double *out) {
    (void)ctx;
    (void)x;
    (void)y;
    (void)p;
    memset(out, 0, DIM * DIM * sizeof(double));
    for (size_t i = 0; i < DIM; i++)
        out[i * DIM + i] = 1. / (double)(i + 2);
    return 1;
}
static void boundary(void *ctx, const double *a, const double *b, const double *p, double *l,
                     double *r) {
    (void)ctx;
    (void)b;
    (void)p;
    (void)r;
    for (size_t i = 0; i < DIM; i++)
        l[i] = a[i] - 1.;
}
static int boundary_jacobian(void *ctx, const double *a, const double *b, const double *p,
                             double *l, double *r) {
    (void)ctx;
    (void)a;
    (void)b;
    (void)p;
    (void)r;
    memset(l, 0, DIM * DIM * sizeof(double));
    for (size_t i = 0; i < DIM; i++)
        l[i * DIM + i] = 1.;
    return 1;
}
static void guess(void *ctx, double x, double *out) {
    (void)ctx;
    (void)x;
    for (size_t i = 0; i < DIM; i++)
        out[i] = 1.;
}
static void heap_fallback_test(void) {
    bvp_problem p = {.dimension = DIM,
                     .left_conditions = DIM,
                     .is_affine = 1,
                     .rhs = rhs,
                     .boundary = boundary,
                     .jacobian = jacobian,
                     .boundary_jacobian = boundary_jacobian};
    double mesh[5];
    assert(bvp_linspace(0., 1., 5, mesh) == BVP_OK);
    bvp_initial_guess g;
    assert(bvp_guess_from_function(5, DIM, 0, mesh, NULL, guess, NULL, &g) == BVP_OK);
    bvp_options o;
    bvp_default_options(&o);
    o.tolerance = 1e-7;
    o.estimate_condition = 1;
    for (size_t repeat = 0; repeat < 3; repeat++) {
        bvp_solution s;
        assert(bvp_solve(&p, &g, &o, &s) == BVP_OK);
        for (size_t point = 0; point <= 100; point++) {
            double x = (double)point / 100., values[DIM];
            assert(bvp_evaluate(&s, x, values) == BVP_OK);
            for (size_t i = 0; i < DIM; i++)
                assert(fabs(values[i] - exp(x / (double)(i + 2))) < 2e-7);
        }
        assert(isfinite(s.statistics.condition_estimate));
        bvp_solution_free(&s);
    }
    p.user = &p;
    o.max_refinements = 2;
    bvp_solution failed;
    assert(bvp_solve(&p, &g, &o, &failed) == BVP_NONFINITE);
    assert(!failed.mesh && !failed.values);
    bvp_solution_free(&failed);
    bvp_guess_free(&g);
}
int main(void) {
    allocation_tests();
    heap_fallback_test();
    puts("BVP optimization: zeroed scratch reuse, detached lifetime, bounded cache, longjmp "
         "cleanup and 17-state heap fallback passed");
    return 0;
}
