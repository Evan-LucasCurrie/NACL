#include "abd.h"
#include "test_data/abd_fixture.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static uint64_t bits(double v) {
    uint64_t r;
    memcpy(&r, &v, 8);
    return r;
}
static double *read_values(FILE *f, const char *tag, size_t n) {
    char label[64];
    assert(fscanf(f, "%63s", label) == 1 && !strcmp(label, tag));
    double *a = calloc(n ? n : 1, 8);
    assert(a);
    for (size_t i = 0; i < n; i++) {
        uint64_t v;
        assert(fscanf(f, "%" SCNx64, &v) == 1);
        memcpy(a + i, &v, 8);
    }
    return a;
}
static void same(const double *a, const double *b, size_t n, size_t id) {
    for (size_t i = 0; i < n; i++)
        if (bits(a[i]) != bits(b[i])) {
            fprintf(stderr, "ABD case %zu index %zu: %016" PRIX64 " != %016" PRIX64 "\n", id, i,
                    bits(a[i]), bits(b[i]));
            abort();
        }
}
int main(void) {
    FILE *f = bac_abd_fixture();
    assert(f);
    char label[64];
    size_t successful = 0, singular = 0, matched = 0;
    while (fscanf(f, "%63s", label) == 1) {
        assert(!strcmp(label, "case"));
        size_t id, p, kcol, nint;
        int status;
        assert(fscanf(f, "%zu %zu %zu %zu %d", &id, &p, &kcol, &nint, &status) == 5);
        size_t m = p * kcol, n = m * nint + 2 * p, edge = 2 * p * p, mid = m * (m + 2 * p) * nint;
        double *top = read_values(f, "input_top", edge), *blocks = read_values(f, "input_mid", mid),
               *bottom = read_values(f, "input_bot", edge), *rhs = read_values(f, "rhs", n),
               *ft = read_values(f, "factor_top", edge), *fm = read_values(f, "factor_mid", mid),
               *fb = read_values(f, "factor_bot", edge),
               *sol = read_values(f, "solution", status ? 0 : n);
        assert(fscanf(f, "%63s", label) == 1 && !strcmp(label, "pivots"));
        size_t *pivots = calloc(n, sizeof(size_t));
        assert(pivots);
        for (size_t i = 0; i < n; i++)
            assert(fscanf(f, "%zu", pivots + i) == 1);
        BacBand matrix;
        assert(!bac_band_new(&matrix, n, n - 1, n - 1));
        matrix.npde = p;
        matrix.kcol = kcol;
        matrix.nint = nint;
        for (size_t j = 0; j < 2 * p; j++)
            for (size_t i = 0; i < p; i++) {
                bac_band_set(&matrix, i, j, top[j * p + i]);
                bac_band_set(&matrix, p + m * nint + i, m * nint + j, bottom[j * p + i]);
            }
        for (size_t k = 0; k < nint; k++)
            for (size_t j = 0; j < m + 2 * p; j++)
                for (size_t i = 0; i < m; i++)
                    bac_band_set(&matrix, p + k * m + i, k * m + j,
                                 blocks[(k * (m + 2 * p) + j) * m + i]);
        BacLu lu;
        int rc = bac_band_factor(&matrix, &lu);
        if (status) {
            assert(rc);
            singular++;
        } else {
            assert(!rc && lu.abd);
            same(lu.abd->top, ft, edge, id);
            same(lu.abd->blocks, fm, mid, id);
            same(lu.abd->bottom, fb, edge, id);
            for (size_t i = 0; i < n; i++)
                assert(lu.abd->pivots[i + 1] == pivots[i]);
            BacLu clone;
            assert(!bac_lu_clone(&clone, &lu));
            assert(!bac_lu_solve(&clone, rhs));
            same(rhs, sol, n, id);
            bac_lu_free(&clone);
            matched += edge * 2 + mid + n;
            successful++;
        }
        bac_lu_free(&lu);
        bac_band_free(&matrix);
        free(top);
        free(blocks);
        free(bottom);
        free(rhs);
        free(ft);
        free(fm);
        free(fb);
        free(sol);
        free(pivots);
    }
    fclose(f);
    assert(successful == 48 && singular == 12);
    printf("ABD: %zu successful, %zu singular cases, %zu exact binary64 fields plus every pivot\n",
           successful, singular, matched);
    return 0;
}
