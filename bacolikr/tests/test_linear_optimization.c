#include "linear.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t bits(double v) {
    uint64_t b;
    memcpy(&b, &v, sizeof b);
    return b;
}

static void residual(const double *a, const double *rhs, const double *x, size_t n) {
    for (size_t row = 0; row < n; row++) {
        double value = 0., scale = fabs(rhs[row]);
        for (size_t col = 0; col < n; col++) {
            value += a[row * n + col] * x[col];
            scale += fabs(a[row * n + col] * x[col]);
        }
        assert(fabs(value - rhs[row]) <= 1e-12 * fmax(scale, 1.));
    }
}

static void band_case(size_t n, size_t lower, size_t upper, int pivot) {
    BacBand matrix;
    assert(!bac_band_new(&matrix, n, lower, upper));
    double *a = calloc(n * n, sizeof(double)), *rhs = calloc(n, sizeof(double));
    double *result = calloc(n, sizeof(double)), *weights = calloc(n, sizeof(double));
    double *norms = calloc(n, sizeof(double));
    assert(a && rhs && result && weights && norms);
    for (size_t row = 0; row < n; row++) {
        for (size_t col = 0; col < n; col++) {
            if (row > col + lower || col > row + upper)
                continue;
            double v = row == col ? 4. : .01 * (double)((row + 3 * col) % 11 + 1);
            if (pivot && lower && row == col + 1)
                v = 12.;
            a[row * n + col] = v;
            bac_band_set(&matrix, row, col, v);
        }
        for (size_t col = 0; col < n; col++)
            rhs[row] += a[row * n + col] * (1. + .2 * sin((double)col));
        weights[row] = .3 + (double)row;
    }
    bac_band_row_norms(&matrix, weights, norms);
    for (size_t row = 0; row < n; row++) {
        double expected = 0.;
        for (size_t col = 0; col < n; col++)
            expected += fabs(a[row * n + col]) * weights[col];
        assert(bits(norms[row]) == bits(expected));
    }
    BacLu lu, clone;
    assert(!bac_band_factor(&matrix, &lu));
    assert(!bac_lu_clone(&clone, &lu));
    bac_lu_free(&lu);
    memcpy(result, rhs, n * sizeof(double));
    assert(!bac_lu_solve(&clone, result));
    residual(a, rhs, result, n);
    bac_lu_free(&clone);
    assert(!bac_dense_factor(a, n, &lu));
    assert(!bac_lu_clone(&clone, &lu));
    bac_lu_free(&lu);
    memcpy(result, rhs, n * sizeof(double));
    assert(!bac_lu_solve(&clone, result));
    residual(a, rhs, result, n);
    bac_lu_free(&clone);
    bac_band_free(&matrix);
    free(a);
    free(rhs);
    free(result);
    free(weights);
    free(norms);
}

int main(void) {
    const size_t sizes[] = {1, 2, 3, 17, 63};
    size_t cases = 0;
    for (size_t i = 0; i < sizeof sizes / sizeof sizes[0]; i++)
        for (size_t lower = 0; lower <= 3; lower++)
            for (size_t upper = 0; upper <= 3; upper++)
                for (int pivot = 0; pivot < 2; pivot++) {
                    band_case(sizes[i], lower, upper, pivot);
                    cases++;
                }
    BacBand matrix;
    BacLu lu;
    assert(!bac_band_new(&matrix, 20, 19, 19));
    for (size_t i = 0; i < 20; i++)
        bac_band_set(&matrix, i, i, 1.);
    assert(!bac_band_factor(&matrix, &lu));
    assert(lu.lower_count == 0 && lu.upper_count == 20);
    bac_lu_free(&lu);
    bac_band_set(&matrix, 0, 0, NAN);
    assert(bac_band_factor(&matrix, &lu) == -2);
    bac_band_set(&matrix, 0, 0, 0.);
    assert(bac_band_factor(&matrix, &lu) == -3);
    bac_band_free(&matrix);
    assert(!bac_band_new(&matrix, 2, 1, 1));
    bac_band_set(&matrix, 0, 0, 1.);
    bac_band_set(&matrix, 1, 1, 1.);
    bac_band_set(&matrix, 0, 1, 1e-300);
    assert(!bac_band_factor(&matrix, &lu));
    double rhs[] = {2., 1e300};
    assert(!bac_lu_solve(&lu, rhs) && rhs[0] == 1. && rhs[1] == 1e300);
    rhs[0] = NAN;
    assert(bac_lu_solve(&lu, rhs) == -2);
    bac_lu_free(&lu);
    bac_band_free(&matrix);
    puts("Active-extent linear factors: 160 band/dense/clone cases, exact row norms, zero extents "
         "and tiny coefficients passed");
    assert(cases == 160);
    return 0;
}
