#include "reference/reference_internal.h"
#include <assert.h>
#include <stdio.h>

static void original_sum(cnref_context *c, int n, double h, R base, int stages, R discrete,
                         int total, R extra, R weights, R out) {
    double zero = 0;
    cnref_dcopy(c, n, &zero, 0, out, 1);
    for (int j = 0; j < stages; j++)
        cnref_daxpy(c, n, weights[j], discrete + j * n, 1, out, 1);
    for (int j = 0; j < total - stages; j++)
        cnref_daxpy(c, n, weights[stages + j], extra + j * n, 1, out, 1);
    cnref_dscal(c, n, h, out, 1);
    cnref_daxpy(c, n, 1., base, 1, out, 1);
}
static void stage_tests(void) {
    cnref_context context = {0};
    for (int n = 1; n <= 7; n++) {
        double base[7], stages[63], extra[28], weights[13], actual[7], expected[7];
        for (int i = 0; i < n; i++)
            base[i] = -.2 + .25 * i;
        for (int i = 0; i < 9 * n; i++)
            stages[i] = .125 * (i % 7 - 3);
        for (int i = 0; i < 4 * n; i++)
            extra[i] = .375 * (i % 5 - 2);
        for (int i = 0; i < 13; i++)
            weights[i] = i % 3 == 0 ? -0. : .0625 * (i % 5 - 2);
        for (int s = 4; s <= 9; s++)
            for (int more = 0; more <= 4; more++) {
                original_sum(&context, n, -.3, base, s, stages, s + more, extra, weights, expected);
                cnref_sum_stages(&context, n, -.3, base, s, stages, s + more, extra, weights,
                                 actual);
                assert(!memcmp(actual, expected, (size_t)n * sizeof(double)));
            }
        /* An aliased source must retain the BLAS routine's streaming writes. */
        double aliased[63], old[63];
        memcpy(aliased, stages, sizeof stages);
        memcpy(old, stages, sizeof stages);
        original_sum(&context, n, .75, base, 4, old, 6, extra, weights, old);
        cnref_sum_stages(&context, n, .75, base, 4, aliased, 6, extra, weights, aliased);
        assert(!memcmp(aliased, old, (size_t)n * sizeof(double)));
        memcpy(actual, base, (size_t)n * sizeof(double));
        memcpy(expected, base, (size_t)n * sizeof(double));
        original_sum(&context, n, .75, expected, 4, stages, 6, extra, weights, expected);
        cnref_sum_stages(&context, n, .75, actual, 4, stages, 6, extra, weights, actual);
        assert(!memcmp(actual, expected, (size_t)n * sizeof(double)));
        double aliased_weights[13], original_weights[13];
        memcpy(aliased_weights, weights, sizeof weights);
        memcpy(original_weights, weights, sizeof weights);
        original_sum(&context, n, .75, base, 4, stages, 6, extra, original_weights,
                     original_weights);
        cnref_sum_stages(&context, n, .75, base, 4, stages, 6, extra, aliased_weights,
                         aliased_weights);
        assert(!memcmp(aliased_weights, original_weights, (size_t)n * sizeof(double)));
    }
}
typedef struct {
    double scale;
    size_t calls;
} Callbacks;
static void rhs(void *ctx, double x, const double *y, double *f) {
    (void)x;
    Callbacks *c = ctx;
    c->calls++;
    f[0] = c->scale * y[0];
}
static double boundary(void *ctx, size_t condition, const double *y) {
    (void)ctx;
    (void)condition;
    return y[0] - 1.;
}
static void guess(void *ctx, double x, double *y) {
    (void)ctx;
    y[0] = 1. + x;
}
static void callback_tests(void) {
    for (size_t k = 2; k <= 4; k++) {
        Callbacks callbacks = {1., 0};
        colnewsc_reference_problem problem = {0};
        problem.base = (colnewsc_problem){.dimension = 1,
                                          .left_boundary_conditions = 1,
                                          .a = 0.,
                                          .b = 1.,
                                          .context = &callbacks,
                                          .rhs = rhs,
                                          .boundary = boundary,
                                          .initial_guess = guess};
        colnewsc_options options = colnewsc_default_options();
        options.tolerance = 1e-8;
        options.collocation_points = k;
        colnewsc_reference_solution *solution = NULL;
        assert(!colnewsc_reference_solve(&problem, &options, &solution));
        size_t intervals = colnewsc_reference_mesh_count(solution) - 1;
        size_t extras = k == 3   ? 1
                        : k == 4 ? 3
                                 : 0,
               expected = (intervals + 1) + intervals * (k + extras);
        double a, b, changed;
        callbacks.calls = 0;
        assert(!colnewsc_reference_evaluate(solution, .371, &a, 1));
        assert(callbacks.calls == expected);
        callbacks.calls = 0;
        assert(!colnewsc_reference_evaluate(solution, .371, &b, 1));
        assert(callbacks.calls == expected && !memcmp(&a, &b, sizeof a));
        callbacks.scale = 2.;
        callbacks.calls = 0;
        assert(!colnewsc_reference_evaluate(solution, .371, &changed, 1));
        assert(callbacks.calls == expected && fabs(changed - a) > 1e-13);
        callbacks.scale = 1.;
        assert(!colnewsc_reference_evaluate(solution, .371, &b, 1));
        assert(!memcmp(&a, &b, sizeof a));
        colnewsc_reference_free(solution);
    }
}
int main(void) {
    stage_tests();
    callback_tests();
    puts("Reference optimization: exact stage sums, zero weights, aliases and callback "
         "re-evaluation passed");
    return 0;
}
