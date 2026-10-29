#include "internal.h"
#include "problems.h"
#include "workloads.h"
#include <assert.h>
static colnewsc_solution *checked_solve(const colnewsc_problem *p, const colnewsc_options *o) {
    colnewsc_solution *s = NULL;
    colnewsc_status e = colnewsc_solve(p, o, &s);
    if (e) {
        fprintf(stderr, "solve failed: %s\n", colnewsc_status_string(e));
        abort();
    }
    assert(s);
    return s;
}
static void all_orders(void) {
    for (size_t k = 1; k <= 7; k++)
        for (int right = 0; right < 2; right++) {
            colnewsc_problem p = workload_problem(0);
            if (right) {
                p.left_boundary_conditions =
                    0; /* Constant boundary shifted by changing interval to [-1,0]. */
                p.a = -1;
                p.b = 0;
            }
            colnewsc_options o = colnewsc_default_options();
            o.collocation_points = k;
            o.tolerance = k == 1 ? 1e-6 : 1e-9;
            colnewsc_solution *s = checked_solve(&p, &o);
            for (int i = 0; i <= 100; i++) {
                double x = p.a + (p.b - p.a) * i / 100., y;
                assert(colnewsc_evaluate(s, x, &y, 1) == 0);
                assert(fabs(y - exp(x)) < o.tolerance * (1 + exp(x)));
            }
            colnewsc_free(s);
        }
}
static void analytic_numerical(void) {
    colnewsc_options o = colnewsc_default_options();
    o.tolerance = 1e-9;
    colnewsc_problem p = workload_problem(0);
    colnewsc_solution *a = checked_solve(&p, &o);
    p.jacobian = NULL;
    p.boundary_gradient = NULL;
    colnewsc_solution *b = checked_solve(&p, &o);
    for (int i = 0; i <= 100; i++) {
        double x = i / 100., u, v;
        assert(!colnewsc_evaluate(a, x, &u, 1));
        assert(!colnewsc_evaluate(b, x, &v, 1));
        assert(fabs(u - v) < 1e-10);
    }
    colnewsc_free(a);
    colnewsc_free(b);
}
static void nonlinear_and_sine(void) {
    for (int id = 1; id <= 2; id++) {
        colnewsc_problem p = workload_problem(id);
        colnewsc_options o = colnewsc_default_options();
        o.tolerance = 1e-9;
        colnewsc_solution *s = checked_solve(&p, &o);
        for (int i = 0; i <= 100; i++) {
            double x = i / 100., y[2];
            assert(!colnewsc_evaluate(s, x, y, 2));
            assert(fabs(y[0] - workload_exact(id, x)) < 1e-8);
        }
        assert(colnewsc_stats(s)->boundary_residual < 1e-10);
        colnewsc_free(s);
    }
}
static void mesh_and_continuation(void) {
    colnewsc_problem p = workload_problem(0);
    colnewsc_options o = colnewsc_default_options();
    o.tolerance = 1e-9;
    double mesh[] = {0, 0.4, 1}, fixed[] = {0.17, 0.67};
    o.initial_mesh = mesh;
    o.initial_mesh_count = 3;
    o.fixed_points = fixed;
    o.fixed_point_count = 2;
    colnewsc_solution *old = checked_solve(&p, &o);
    for (size_t j = 0; j < 2; j++) {
        int found = 0;
        for (size_t i = 0; i < colnewsc_mesh_count(old); i++)
            found |= colnewsc_mesh(old)[i] == fixed[j];
        assert(found);
    }
    colnewsc_solution *s = NULL;
    assert(!colnewsc_solve_with_initial_solution(&p, &o, old, &s));
    for (size_t i = 0; i < colnewsc_mesh_count(s); i++) {
        double y;
        assert(!colnewsc_evaluate(s, colnewsc_mesh(s)[i], &y, 1));
        assert(y == colnewsc_mesh_values(s)[i]);
    }
    colnewsc_free(old);
    colnewsc_free(s);
}
static void bad_rhs(void *c, double x, const double *y, double *f) {
    (void)c;
    (void)x;
    (void)y;
    f[0] = NAN;
}
static void errors(void) {
    colnewsc_problem p = workload_problem(0);
    colnewsc_options o = colnewsc_default_options();
    colnewsc_solution *s = NULL;
    o.tolerance = NAN;
    assert(colnewsc_solve(&p, &o, &s) == COLNEWSC_INVALID_INPUT);
    o = colnewsc_default_options();
    o.collocation_points = 8;
    assert(colnewsc_solve(&p, &o, &s) == COLNEWSC_INVALID_INPUT);
    o = colnewsc_default_options();
    o.max_subintervals = 5;
    assert(colnewsc_solve(&p, &o, &s) == COLNEWSC_MESH_LIMIT);
    o = colnewsc_default_options();
    double bad[] = {0, 0.4, 0.3, 1};
    o.initial_mesh = bad;
    o.initial_mesh_count = 4;
    assert(colnewsc_solve(&p, &o, &s) == COLNEWSC_INVALID_INPUT);
    o = colnewsc_default_options();
    p.rhs = bad_rhs;
    assert(colnewsc_solve(&p, &o, &s) == COLNEWSC_NON_FINITE);
    p = workload_problem(0);
    s = checked_solve(&p, &o);
    double y;
    assert(colnewsc_evaluate(s, NAN, &y, 1) == COLNEWSC_OUTSIDE_INTERVAL);
    assert(colnewsc_evaluate(s, -0.1, &y, 1) == COLNEWSC_OUTSIDE_INTERVAL);
    assert(colnewsc_evaluate(s, 0.5, &y, 0) == COLNEWSC_INVALID_INPUT);
    colnewsc_free(s);
}
static void band_and_dense(void) {
    for (size_t n = 1; n <= 31; n += 3)
        for (size_t l = 0; l <= 4; l++)
            for (size_t u = 0; u <= 4; u++) {
                cn_lu a = {0};
                assert(!cn_band_new(&a, n, l, u));
                double *dense = cn_alloc(n * n, sizeof(double)),
                       *b = cn_alloc(n * 3, sizeof(double)), *x = cn_alloc(n * 3, sizeof(double));
                assert(dense && b && x);
                for (size_t i = 0; i < n; i++)
                    for (size_t j = i > l ? i - l : 0; j < CN_MIN(i + u + 1, n); j++) {
                        double v = i == j ? 8 : sin(i * 1.31 + j * 0.77);
                        dense[i * n + j] = v;
                        cn_band_set(&a, i, j, v);
                    }
                for (size_t i = 0; i < n * 3; i++)
                    x[i] = cos(i * 0.61);
                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j < n; j++)
                        for (size_t q = 0; q < 3; q++)
                            b[i * 3 + q] += dense[i * n + j] * x[j * 3 + q];
                assert(!cn_lu_factor(&a));
                assert(!cn_lu_solve(&a, b, 3));
                for (size_t i = 0; i < n * 3; i++)
                    assert(fabs(b[i] - x[i]) < 2e-12);
                cn_lu_free(&a);
                free(dense);
                free(b);
                free(x);
            }
}
static size_t degree;
static void polynomial_rhs(void *c, double x, const double *y, double *f) {
    (void)c;
    (void)y;
    f[0] = pow(x, (double)degree);
}
static double zero_bc(void *c, size_t i, const double *y) {
    (void)c;
    (void)i;
    return y[0];
}
static void polynomial(void) {
    for (size_t k = 3; k <= 4; k++)
        for (degree = 0; degree < 2 * k - 1; degree++) {
            colnewsc_problem p = workload_problem(0);
            p.rhs = polynomial_rhs;
            p.boundary = zero_bc;
            p.jacobian = NULL;
            p.initial_guess = NULL;
            colnewsc_options o = colnewsc_default_options();
            o.collocation_points = k;
            o.initial_subintervals = 1;
            o.tolerance = 1e-9;
            colnewsc_solution *s = checked_solve(&p, &o);
            for (int i = 0; i <= 100; i++) {
                double x = i / 100., y;
                assert(!colnewsc_evaluate(s, x, &y, 1));
                assert(fabs(y - pow(x, (double)(degree + 1)) / (double)(degree + 1)) < 2e-12);
            }
            colnewsc_free(s);
        }
}
static void upstream(void) {
    for (int id = 0; id < CN_TEST_CASE_COUNT; id++) {
        colnewsc_problem p = colnewsc_test_problem((colnewsc_test_case)id);
        colnewsc_options o = colnewsc_test_options((colnewsc_test_case)id);
        colnewsc_solution *s = checked_solve(&p, &o);
        double error = 0;
        for (int i = 0; i <= 100; i++) {
            double x = p.a + (p.b - p.a) * i / 100., v[20],
                   exact = colnewsc_test_exact((colnewsc_test_case)id, x);
            assert(!colnewsc_evaluate(s, x, v, p.dimension));
            if (isfinite(exact))
                error = fmax(error, fabs(v[0] - exact));
        }
        printf("  %-12s intervals=%zu first-component-error=%.3g\n",
               colnewsc_test_name((colnewsc_test_case)id), colnewsc_mesh_count(s) - 1, error);
        assert(error < 2e-5);
        colnewsc_free(s);
    }
}
int main(void) {
    all_orders();
    analytic_numerical();
    nonlinear_and_sine();
    mesh_and_continuation();
    errors();
    band_and_dense();
    polynomial();
    upstream();
    puts("COLNEWSC C: all tests passed");
    return 0;
}
