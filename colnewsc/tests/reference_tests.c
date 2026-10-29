#include "colnewsc_reference.h"
#include "problems.h"
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "test_data/reference_fixtures.h"
static uint64_t bits(double value) {
    uint64_t result;
    memcpy(&result, &value, sizeof(result));
    return result;
}
static double from_bits(uint64_t value) {
    double result;
    memcpy(&result, &value, sizeof(result));
    return result;
}
static FILE *fixture_open(const char *name) {
    for (size_t i = 0; i < sizeof cn_fixtures / sizeof cn_fixtures[0]; i++) {
        const cn_fixture *fixture = cn_fixtures + i;
        if (strcmp(name, fixture->name))
            continue;
        FILE *stream = tmpfile();
        assert(stream);
        for (size_t line = 0; line < fixture->count; line++)
            assert(fputs(fixture->lines[line], stream) >= 0);
        assert(!ferror(stream));
        rewind(stream);
        return stream;
    }
    fprintf(stderr, "missing embedded fixture: %s\n", name);
    abort();
}
static void fixture_tests(void) {
    size_t samples = 0, values = 0, coefficients = 0, mismatches = 0;
    for (int which = 0; which < CN_TEST_CASE_COUNT; which++) {
        colnewsc_test_case id = (colnewsc_test_case)which;
        colnewsc_reference_problem p = colnewsc_reference_test_problem(id);
        colnewsc_options o = colnewsc_test_options(id);
        colnewsc_reference_solution *s = NULL;
        colnewsc_status status = colnewsc_reference_solve(&p, &o, &s);
        if (status) {
            fprintf(stderr, "reference %s solve: %s\n", colnewsc_test_name(id),
                    colnewsc_status_string(status));
            assert(status == COLNEWSC_OK);
        }
        FILE *f = fixture_open(colnewsc_test_name(id));
        char line[4096];
        size_t work = 0, case_samples = 0, case_mismatches = 0, history_count = 0,
               representation_count = 0;
        const int *history = colnewsc_reference_mesh_history(s, &history_count);
        const double *representation = colnewsc_reference_representation(s, &representation_count);
        while (fgets(line, sizeof(line), f)) {
            if (!strncmp(line, "MESH ", 5)) {
                char *next = line + 5;
                size_t i = 0;
                for (;;) {
                    char *end;
                    long expected = strtol(next, &end, 10);
                    assert(end != next && i < history_count);
                    if (history[i] != expected) {
                        fprintf(stderr, "%s history[%zu]: %d != %ld\n", colnewsc_test_name(id),
                                i, history[i], expected);
                        abort();
                    }
                    i++;
                    if (*end != ',')
                        break;
                    next = end + 1;
                }
                assert(i == history_count);
            } else if (!strncmp(line, "SAMPLE ", 7)) {
                char *end;
                uint64_t raw = strtoull(line + 7, &end, 16);
                double x = from_bits(raw), actual[20];
                status = colnewsc_reference_evaluate(s, x, actual, p.base.dimension);
                if (status) {
                    fprintf(stderr, "%s evaluate %.17g: %s\n", colnewsc_test_name(id), x,
                            colnewsc_status_string(status));
                    abort();
                }
                for (size_t j = 0; j < p.base.dimension; j++) {
                    assert(*end == ',');
                    char *next = end + 1;
                    uint64_t expected = strtoull(next, &end, 16);
                    assert(end != next);
                    if (bits(actual[j]) != expected) {
                        if (case_mismatches < 3)
                            fprintf(stderr,
                                    "%s sample %.17g component %zu: %016" PRIx64 " != %016" PRIx64
                                    " (%.17g vs %.17g)\n",
                                    colnewsc_test_name(id), x, j, bits(actual[j]), expected,
                                    actual[j], from_bits(expected));
                        case_mismatches++;
                    }
                    values++;
                }
                case_samples++;
                samples++;
            } else if (!strncmp(line, "WORK ", 5)) {
                uint64_t expected = strtoull(line + 5, NULL, 16);
                assert(work < representation_count);
                if (bits(representation[work]) != expected) {
                    if (case_mismatches < 3)
                        fprintf(stderr, "%s coefficient %zu: %016" PRIx64 " != %016" PRIx64 "\n",
                                colnewsc_test_name(id), work, bits(representation[work]),
                                expected);
                    case_mismatches++;
                }
                work++;
                coefficients++;
            } else {
                fprintf(stderr, "unexpected fixture line: %s", line);
                abort();
            }
        }
        fclose(f);
        assert(case_samples == 11 && work == representation_count);
        mismatches += case_mismatches;
        printf("Reference %-12s: %zu samples, %zu stored coefficients, %zu bit mismatches\n",
               colnewsc_test_name(id), case_samples, work, case_mismatches);
        colnewsc_reference_free(s);
    }
    assert(samples == 143 && values == 704 && coefficients == 24275);
    printf("COLNEWSC reference: %zu values and %zu coefficients checked; %zu bit mismatches\n",
           values, coefficients, mismatches);
    assert(mismatches == 0);
}
typedef struct {
    int fault, fail_interpolation;
} quality;
static void quality_rhs(void *v, double x, const double *z, double *f) {
    (void)x;
    (void)z;
    quality *q = v;
    if (q->fault == 3) {
        f[0] = 0;
        return;
    }
    for (size_t i = 0; i < 2; i++)
        f[i] = q->fail_interpolation ? NAN
               : q->fault == 1       ? NAN
               : q->fault == 2       ? INFINITY
               : q->fault == 9       ? 1e200
                                     : 0;
}
static double quality_bc(void *v, size_t i, const double *z) {
    quality *q = v;
    return q->fault == 4 ? NAN : q->fault == 5 ? 1 : z[i];
}
static void quality_guess(void *v, double x, double *z, double *dmz) {
    (void)x;
    quality *q = v;
    dmz[0] = dmz[1] = 0;
    if (q->fault == 6) {
        z[0] = 0;
        return;
    }
    double value = q->fault == 7    ? NAN
                   : q->fault == 8  ? from_bits(1)
                   : q->fault == 10 ? from_bits(UINT64_C(0x7fefffffffffffff))
                                    : 0;
    z[0] = z[1] = value;
}
static colnewsc_reference_problem quality_problem(quality *q) {
    colnewsc_reference_problem p = {0};
    p.base.dimension = 2;
    p.base.left_boundary_conditions = 2;
    p.base.a = 0;
    p.base.b = 1;
    p.base.context = q;
    p.base.rhs = quality_rhs;
    p.base.boundary = quality_bc;
    p.reference_guess = quality_guess;
    return p;
}
static void quality_tests(void) {
    colnewsc_options o = colnewsc_default_options();
    o.max_subintervals = 100;
    for (int i = 1; i <= 10; i++) {
        quality q = {i, 0};
        colnewsc_reference_problem p = quality_problem(&q);
        colnewsc_reference_solution *s = NULL;
        colnewsc_status e = colnewsc_reference_solve(&p, &o, &s);
        if (e != (i == 5 ? COLNEWSC_SINGULAR_JACOBIAN : COLNEWSC_NON_FINITE)) {
            fprintf(stderr, "quality fault %d: %s\n", i, colnewsc_status_string(e));
            abort();
        }
        assert(s == NULL);
    }
    quality q = {0, 0};
    colnewsc_reference_problem p = quality_problem(&q);
    colnewsc_reference_solution *s = NULL;
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_OK);
    double result[] = {123, 456};
    assert(colnewsc_reference_evaluate(s, -1, result, 2) == COLNEWSC_OUTSIDE_INTERVAL);
    assert(result[0] == 123 && result[1] == 456);
    q.fail_interpolation = 1;
    assert(colnewsc_reference_evaluate(s, 0.35, result, 2) == COLNEWSC_NON_FINITE);
    assert(result[0] == 123 && result[1] == 456);
    q.fail_interpolation = 0;
    assert(colnewsc_reference_evaluate(s, 0.35, result, 2) == COLNEWSC_OK);
    assert(result[0] == 0 && result[1] == 0);
    colnewsc_reference_free(s);
    o.initial_subintervals = SIZE_MAX;
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_INVALID_INPUT);
    o.initial_subintervals = 5;
    o.max_subintervals = INT32_MAX;
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_INVALID_INPUT);
    puts("COLNEWSC reference: invalid/incomplete callbacks, singular systems, input bounds and "
         "recoverable interpolation faults passed");
}

typedef struct {
    size_t dimension, left;
} linear_problem;
static void linear_rhs(void *v, double x, const double *z, double *f) {
    (void)x;
    (void)z;
    linear_problem *p = v;
    for (size_t i = 0; i < p->dimension; i++)
        f[i] = (i + 1) * 0.25;
}
static double linear_bc(void *v, size_t i, const double *z) {
    linear_problem *p = v;
    return z[i] - 1 - (i < p->left ? 0 : (i + 1) * 0.25);
}
static void linear_guess(void *v, double x, double *z) {
    linear_problem *p = v;
    for (size_t i = 0; i < p->dimension; i++)
        z[i] = 1 + (i + 1) * 0.25 * x;
}
static colnewsc_reference_problem make_linear(linear_problem *d) {
    colnewsc_reference_problem p = {0};
    p.base.dimension = d->dimension;
    p.base.left_boundary_conditions = d->left;
    p.base.a = 0;
    p.base.b = 1;
    p.base.context = d;
    p.base.rhs = linear_rhs;
    p.base.boundary = linear_bc;
    p.base.initial_guess = linear_guess;
    return p;
}
static void partitions_and_limits(void) {
    size_t dimensions[] = {1, 2, 7, 20}, solves = 0;
    for (size_t a = 0; a < 4; a++)
        for (size_t k = 1; k <= 7; k++)
            for (size_t side = 0; side < 3; side++) {
                linear_problem data = {dimensions[a], side == 0   ? 0
                                                      : side == 1 ? dimensions[a] / 2
                                                                  : dimensions[a]};
                colnewsc_reference_problem p = make_linear(&data);
                colnewsc_options o = colnewsc_default_options();
                o.collocation_points = k;
                o.tolerance = 1e-10;
                o.max_subintervals = 100;
                colnewsc_reference_solution *s = NULL;
                assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_OK);
                for (size_t i = 0; i <= 31; i++) {
                    double x = i / 31.0, z[20];
                    assert(colnewsc_reference_evaluate(s, x, z, data.dimension) == COLNEWSC_OK);
                    for (size_t j = 0; j < data.dimension; j++)
                        if (!(fabs(z[j] - (1 + (j + 1) * 0.25 * x)) <=
                              o.tolerance * (1 + fabs(1 + (j + 1) * 0.25 * x)))) {
                            fprintf(stderr,
                                    "linear d=%zu k=%zu left=%zu x=%.17g j=%zu got %.17g expected "
                                    "%.17g\n",
                                    data.dimension, k, data.left, x, j, z[j],
                                    1 + (j + 1) * 0.25 * x);
                            abort();
                        }
                }
                colnewsc_reference_free(s);
                solves++;
            }
    linear_problem data = {1, 1};
    colnewsc_reference_problem p = make_linear(&data);
    double fixed[] = {0.13, 0.37};
    colnewsc_options o = colnewsc_default_options();
    o.max_subintervals = 100;
    o.fixed_points = fixed;
    o.fixed_point_count = 2;
    colnewsc_reference_solution *s = NULL;
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_OK);
    const double *mesh = colnewsc_reference_mesh(s);
    size_t points = colnewsc_reference_mesh_count(s);
    for (size_t i = 0; i < 2; i++) {
        int found = 0;
        for (size_t j = 0; j < points; j++)
            if (mesh[j] == fixed[i])
                found = 1;
        assert(found);
    }
    colnewsc_reference_free(s);
    double supplied[] = {0, 0.5, 1};
    o.initial_mesh = supplied;
    o.initial_mesh_count = 3;
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_INVALID_INPUT);
    quality q = {0, 0};
    p = quality_problem(&q);
    o = colnewsc_default_options();
    o.max_subintervals = 100;
    o.max_refinements = 0;
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_REFINEMENT_LIMIT);
    o.max_refinements = 20;
    o.max_newton_iterations = 0;
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_INVALID_INPUT);
    p = colnewsc_reference_test_problem(CN_TP21);
    o = colnewsc_test_options(CN_TP21);
    o.max_newton_iterations = 1;
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_NEWTON_NOT_CONVERGED);
    p = colnewsc_reference_test_problem(CN_SIWA22_EX1);
    o = colnewsc_test_options(CN_SIWA22_EX1);
    assert(colnewsc_reference_solve(&p, &o, &s) == COLNEWSC_OK);
    colnewsc_reference_diagnostics diagnostics = colnewsc_reference_get_diagnostics(s);
    assert(diagnostics.non_finite_interpolation_evaluations > 0 &&
           diagnostics.has_first_non_finite_interpolation_point &&
           diagnostics.first_non_finite_interpolation_point == 0);
    double z[] = {123, 123, 123};
    mesh = colnewsc_reference_mesh(s);
    assert(colnewsc_reference_evaluate(s, mesh[1] * 0.5, z, 3) == COLNEWSC_NON_FINITE);
    assert(z[0] == 123 && z[1] == 123 && z[2] == 123);
    assert(colnewsc_reference_evaluate(s, mesh[1], z, 3) == COLNEWSC_OK);
    colnewsc_reference_free(s);
    printf("COLNEWSC reference: %zu order/dimension/boundary combinations, fixed points, iteration "
           "bounds and compatibility diagnostics passed\n",
           solves);
}
int main(void) {
    fixture_tests();
    quality_tests();
    partitions_and_limits();
    return 0;
}
