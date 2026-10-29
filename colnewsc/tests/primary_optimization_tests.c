#include "internal.h"
#include <assert.h>

static void check_system(cn_lu *factor, const double *matrix, size_t n, size_t columns) {
    double *rhs = calloc(n * columns, sizeof(double));
    assert(rhs);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            for (size_t q = 0; q < columns; q++)
                rhs[i * columns + q] += matrix[i * n + j] * (1. + (double)(q + j) / 8.);
    assert(!cn_lu_solve(factor, rhs, columns));
    for (size_t i = 0; i < n; i++)
        for (size_t q = 0; q < columns; q++)
            assert(fabs(rhs[i * columns + q] - (1. + (double)(q + i) / 8.)) < 2e-12);
    free(rhs);
}
int main(void) {
    cn_lu factor = {0};
    const size_t sizes[] = {3, 17, 40, 4, 31, 2};
    for (size_t iteration = 0; iteration < 18; iteration++) {
        size_t n = sizes[iteration % 6];
        double *matrix = calloc(n * n, sizeof(double));
        assert(matrix);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                matrix[i * n + j] = i == j ? 4.
                                    : iteration % 2 || i + 1 == j || j + 1 == i
                                        ? .02 * cos((double)(i + j))
                                        : 0.;
        if (iteration % 3 == 1)
            for (size_t j = 0; j < n; j++) {
                double value = matrix[j];
                matrix[j] = matrix[(n - 1) * n + j];
                matrix[(n - 1) * n + j] = value;
            }
        assert(!cn_stage_factor_reuse(matrix, n, &factor));
        for (size_t columns = 1; columns <= 5; columns++)
            check_system(&factor, matrix, n, columns);
        double *data = factor.data;
        assert(!cn_stage_factor_reuse(matrix, n, &factor));
        assert(factor.data == data);
        matrix[0] = NAN;
        assert(cn_stage_factor_reuse(matrix, n, &factor) == COLNEWSC_NON_FINITE);
        free(matrix);
    }
    assert(!cn_band_reset(&factor, 12, 2, 3));
    double *data = factor.data;
    for (size_t i = 0; i < 12; i++)
        cn_band_set(&factor, i, i, 2.);
    assert(!cn_lu_factor(&factor));
    assert(!cn_band_reset(&factor, 12, 2, 3));
    assert(factor.data == data);
    for (size_t i = 0; i < 12 * factor.stride; i++)
        assert(factor.data[i] == 0.);
    cn_lu_free(&factor);

    double slab[16] = {0};
    size_t pivots[4] = {0};
    cn_lu borrowed = {
        .data = slab, .pivots = pivots, .data_capacity = 16, .pivot_capacity = 4, .borrowed = 1};
    double matrix[] = {2., 1., 0., 3.};
    assert(!cn_stage_factor_reuse(matrix, 2, &borrowed));
    assert(borrowed.borrowed && borrowed.data == slab && borrowed.pivots == pivots);
    check_system(&borrowed, matrix, 2, 3);
    cn_lu_free(&borrowed);
    slab[0] = 17.;
    assert(slab[0] == 17.);

    puts("Primary optimization: reusable dense/band factors, all RHS widths, pivoting and borrowed "
         "storage passed");
    return 0;
}
