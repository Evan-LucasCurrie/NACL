#include "bvp.h"
#include "internal.h"
#include "test_problems.h"
#include <assert.h>
#include <stdio.h>
#define CHECK(e)                                                                                   \
    do {                                                                                           \
        bvp_status status_ = (e);                                                                  \
        if (status_) {                                                                             \
            fprintf(stderr, "%s:%d %s: %s\n", __FILE__, __LINE__, #e, bvp_status_string(status_)); \
            assert(status_ == BVP_OK);                                                             \
        }                                                                                          \
    } while (0)
static bvp_solution run(test_problem *d, bvp_method method, bvp_error_control mode, int estimates) {
    bvp_problem p = tp_make(d);
    double mesh[5], parameter = 2;
    CHECK(bvp_linspace(0, 1, 5, mesh));
    bvp_initial_guess g;
    CHECK(bvp_guess_from_function(5, p.dimension, p.parameter_count, mesh, &parameter, tp_guess, d,
                                  &g));
    bvp_options o;
    bvp_default_options(&o);
    o.method = method;
    o.error_control = mode;
    o.max_intervals = 12000;
    o.estimate_error = estimates;
    o.estimate_condition = estimates;
    bvp_solution s;
    CHECK(bvp_solve(&p, &g, &o, &s));
    bvp_guess_free(&g);
    return s;
}
static void numerical_tests(void) {
    bvp_method methods[] = {BVP_MIRK2, BVP_MIRK4, BVP_MIRK6};
    size_t solves = 0;
    for (size_t m = 0; m < 3; m++)
        for (int analytic = 0; analytic <= 1; analytic++)
            for (int mode = 0; mode < 5; mode++) {
                test_problem p = {0, analytic};
                bvp_solution s = run(&p, methods[m], (bvp_error_control)mode, 1);
                assert(s.statistics.boundary_residual < 1e-10);
                assert(isfinite(s.statistics.richardson_error));
                assert(isfinite(s.statistics.higher_order_error));
                assert(isfinite(s.statistics.deferred_correction_error));
                assert(isfinite(s.statistics.condition_estimate));
                if (mode == BVP_DEFECT || mode == BVP_DEFECT_AND_GLOBAL ||
                    mode == BVP_DEFECT_AND_DERIVATIVE)
                    assert(s.statistics.defect_norm <= 1e-6);
                if (mode == BVP_GLOBAL || mode == BVP_DEFECT_AND_GLOBAL)
                    assert(s.statistics.global_error <= 1e-6);
                if (mode == BVP_HYBRID)
                    assert(s.statistics.global_error + s.statistics.defect_norm <= 1e-6);
                if (mode == BVP_DEFECT_AND_DERIVATIVE)
                    assert(s.statistics.derivative_error <= 1e-6);
                double points[] = {0, 0.0317, 0.1931, 0.4173, 0.6839, 0.9271, 1};
                for (size_t k = 0; k < 7; k++) {
                    double y[2], dy[2];
                    CHECK(bvp_evaluate(&s, points[k], y));
                    CHECK(bvp_derivative(&s, points[k], dy));
                    for (size_t j = 0; j < 2; j++) {
                        assert(fabs(y[j] - exp(points[k])) < 2e-5);
                        assert(fabs(dy[j] - exp(points[k])) < 4e-5);
                    }
                }
                bvp_solution_free(&s);
                solves++;
            }
    for (int kind = 1; kind <= 5; kind++)
        for (size_t m = 0; m < 3; m++) {
            test_problem p = {kind, 1};
            bvp_solution s = run(&p, methods[m], BVP_DEFECT, 0);
            for (size_t k = 0; k <= 100; k++) {
                double y[2];
                CHECK(bvp_evaluate(&s, k / 100.0, y));
                assert(fabs(y[0] - tp_exact(&p, k / 100.0)) < 3e-5);
            }
            if (kind == 3)
                assert(fabs(s.parameters[0] - 3) < 1e-10);
            bvp_solution_free(&s);
            solves++;
        }
    printf("BVP: %zu solver cases passed (all MIRK methods, all error controls, finite "
           "differences, singular endpoints, parameters, nonlinear Bratu and band fallback)\n",
           solves);
}
static void api_tests(void) {
    test_problem p = {0, 1};
    bvp_solution s = run(&p, BVP_MIRK4, BVP_DEFECT, 0);
    double out[2];
    assert(bvp_evaluate(&s, -0.01, out) == BVP_OUT_OF_DOMAIN);
    assert(bvp_evaluate(&s, NAN, out) == BVP_OUT_OF_DOMAIN);
    bvp_initial_guess g;
    CHECK(bvp_extend(&s, -0.2, 1.2, 1, &g));
    assert(g.mesh[0] == -0.2 && g.mesh[g.points - 1] == 1.2);
    assert(fabs(g.values[0] - 0.8) < 1e-6);
    bvp_guess_free(&g);
    CHECK(bvp_solution_initial_guess(&s, &g));
    bvp_problem problem = tp_make(&p);
    bvp_solution continued;
    CHECK(bvp_solve(&problem, &g, NULL, &continued));
    bvp_guess_free(&g);
    bvp_solution_free(&continued);
    CHECK(bvp_save(&s, "bvp_test_solution.tmp"));
    FILE *saved = fopen("bvp_test_solution.tmp", "r");
    assert(saved);
    char signature[32];
    unsigned version;
    assert(fscanf(saved, "%31s %u", signature, &version) == 2);
    assert(strcmp(signature, "BVP_SOLVER_C") == 0 && version == 1);
    assert(fclose(saved) == 0);
    bvp_solution loaded;
    CHECK(bvp_load("bvp_test_solution.tmp", &loaded));
    remove("bvp_test_solution.tmp");
    for (size_t k = 0; k <= 20; k++) {
        double a[2], b[2];
        CHECK(bvp_evaluate(&s, k / 20.0, a));
        CHECK(bvp_evaluate(&loaded, k / 20.0, b));
        assert(memcmp(a, b, sizeof(a)) == 0);
    }
    bvp_solution_free(&loaded);
    double mesh[] = {0, 1}, values[] = {0, 0, 0, 0};
    CHECK(bvp_guess_from_values(2, 2, 0, mesh, values, NULL, &g));
    bvp_options o;
    bvp_default_options(&o);
    o.tolerance = NAN;
    assert(bvp_solve(&problem, &g, &o, &loaded) == BVP_INVALID_INPUT);
    o.tolerance = 1e-16;
    assert(bvp_solve(&problem, &g, &o, &loaded) == BVP_INVALID_INPUT);
    bvp_guess_free(&g);
    double badmesh[] = {0, 0};
    assert(bvp_guess_from_values(2, 2, 0, badmesh, values, NULL, &g) == BVP_INVALID_INPUT);
    double line[3];
    assert(bvp_linspace(1, 1 + DBL_EPSILON, 3, line) == BVP_INVALID_INPUT);
    assert(bvp_linspace(-DBL_MAX, DBL_MAX, 3, line) == BVP_INVALID_INPUT);
    bvp_solution_free(&s);
    puts(
        "BVP API: continuation, extension, persistence, input validation and domain checks passed");
}
static void linear_tests(void) {
    context *c = calloc(1, sizeof(*c));
    jmp_buf top;
    c->trap = &top;
    if (setjmp(top)) {
        fprintf(stderr, "linear test failed: %s\n", bvp_status_string(c->error));
        abort();
    }
    /* Dense known systems with row/column pivoting, every boundary split. */
    for (size_t q = 1; q <= 6; q++)
        for (size_t l = 0; l <= q; l++) {
            matrix m = bvp_matrix(c, q, l, 4);
            size_t n = m.n;
            double *dense = bvp_alloc(c, n * n, sizeof(double));
            for (size_t row = 0; row < n; row++) {
                size_t first = row < l           ? 0
                               : row < l + 3 * q ? ((row - l) / q) * q
                                                 : 3 * q,
                       width = row < l || row >= l + 3 * q ? q : 2 * q;
                double *r = bvp_matrix_row(&m, row, first);
                for (size_t j = 0; j < width; j++) {
                    double v = sin((row + 1) * (first + j + 2));
                    if (first + j == row)
                        v += 3 * q;
                    r[j] = v;
                    dense[row * n + first + j] = v;
                }
            }
            bvp_factor(c, &m);
            double *rhs = bvp_alloc(c, n, sizeof(double));
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    rhs[i] += dense[i * n + j] * (j + 1);
            bvp_linear(c, &m, rhs, 0);
            for (size_t i = 0; i < n; i++)
                assert(fabs(rhs[i] - (i + 1)) < 1e-8);
            memset(rhs, 0, n * sizeof(double));
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    rhs[i] += dense[j * n + i] * (j + 1);
            bvp_linear(c, &m, rhs, 1);
            for (size_t i = 0; i < n; i++)
                assert(fabs(rhs[i] - (i + 1)) < 1e-8);
        }
    /* Force top, interval and bottom COLROW pivots in a well-conditioned system. */
    matrix permutation = bvp_matrix(c, 4, 2, 3);
    for (size_t row = 0; row < 12; row++)
        *bvp_matrix_row(&permutation, row, row) = 2;
    for (size_t first = 0; first < 12; first += 2) {
        size_t start = first < 2 ? 0 : first < 10 ? ((first - 2) / 4) * 4 : 8;
        size_t width = first < 2 || first >= 10 ? 4 : 8;
        double *a = bvp_matrix_row(&permutation, first, start);
        double *b = bvp_matrix_row(&permutation, first + 1, start);
        for (size_t j = 0; j < width; j++)
            SWAP(double, a[j], b[j]);
    }
    double known[12], rhs_perm[12];
    for (size_t i = 0; i < 12; i++)
        known[i] = i * 0.125 - 0.5;
    for (size_t i = 0; i < 12; i++)
        rhs_perm[i] = 2 * known[i ^ 1];
    bvp_factor(c, &permutation);
    assert(permutation.pivot[0] == 1 && permutation.pivot[2] == 3 && permutation.pivot[4] == 5 &&
           permutation.pivot[10] == 11);
    bvp_linear(c, &permutation, rhs_perm, 0);
    assert(memcmp(known, rhs_perm, sizeof(known)) == 0);
    for (size_t i = 0; i < 12; i++)
        rhs_perm[i] = 2 * known[i ^ 1];
    bvp_linear(c, &permutation, rhs_perm, 1);
    assert(memcmp(known, rhs_perm, sizeof(known)) == 0);
    /* Generic band pivots with 240 orders of row scaling. */
    matrix band = bvp_matrix(c, 3, 3, 8);
    double scales[] = {1e-120, 1e120, 1e-60, 1e60, 1};
    double known_band[24], rhs_band[24];
    for (size_t i = 0; i < 24; i++)
        known_band[i] = sin(i + 0.3);
    for (size_t i = 0; i < 24; i++) {
        size_t col = i / 3 * 3 + (i + 1) % 3;
        *bvp_matrix_row(&band, i, col) = scales[i % 5] * (i % 3 + 2);
        rhs_band[i] = scales[i % 5] * (i % 3 + 2) * known_band[col];
    }
    bvp_factor(c, &band);
    bvp_linear(c, &band, rhs_band, 0);
    for (size_t i = 0; i < 24; i++)
        assert(fabs(rhs_band[i] - known_band[i]) < 3e-15);
    memset(rhs_band, 0, sizeof(rhs_band));
    for (size_t i = 0; i < 24; i++) {
        size_t col = i / 3 * 3 + (i + 1) % 3;
        rhs_band[col] += (i % 3 + 2) * known_band[i];
    }
    bvp_linear(c, &band, rhs_band, 1);
    for (size_t i = 0; i < 24; i++)
        assert(fabs(rhs_band[i] * scales[i % 5] - known_band[i]) < 3e-15);
    const double a[] = {1, 2, 2, 4};
    double *inv = bvp_pinv(c, a, 2);
    assert(fabs(inv[0] - 0.04) < 1e-12);
    assert(fabs(inv[1] - 0.08) < 1e-12);
    assert(fabs(inv[2] - 0.08) < 1e-12);
    assert(fabs(inv[3] - 0.16) < 1e-12);
    bvp_release(c, NULL);
    free(c);
    puts("BVP linear algebra: COLROW and band normal/transpose solves and rank-deficient "
         "pseudoinverse passed");
}
int main(void) {
    linear_tests();
    numerical_tests();
    api_tests();
    return 0;
}
