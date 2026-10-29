/* Compact band and dense partial-pivot elimination. See LICENSE. */
#include "linear.h"
#include "abd.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
int bac_band_new(BacBand *a, size_t n, size_t lower, size_t upper) {
    memset(a, 0, sizeof *a);
    if (!n)
        return -1;
    a->n = n;
    a->lower = MIN(lower, n - 1);
    a->upper = MIN(upper, n - 1);
    a->stride = a->lower + a->upper + 1;
    a->values = calloc(n * a->stride, sizeof(double));
    return a->values ? 0 : -1;
}
void bac_band_free(BacBand *a) {
    free(a->values);
    memset(a, 0, sizeof *a);
}
double bac_band_get(const BacBand *a, size_t r, size_t c) {
    if (r >= a->n || c >= a->n || r > c + a->lower || c > r + a->upper)
        return 0;
    return a->values[c * a->stride + a->upper + r - c];
}
void bac_band_set(BacBand *a, size_t r, size_t c, double v) {
    if (r < a->n && c < a->n && r <= c + a->lower && c <= r + a->upper)
        a->values[c * a->stride + a->upper + r - c] = v;
}
void bac_band_add(BacBand *a, size_t r, size_t c, double v) {
    if (r < a->n && c < a->n && r <= c + a->lower && c <= r + a->upper)
        a->values[c * a->stride + a->upper + r - c] += v;
}
void bac_band_row_norms(const BacBand *a, const double *w, double *o) {
    memset(o, 0, a->n * sizeof(double));
    for (size_t c = 0; c < a->n; c++) {
        size_t s = c > a->upper ? c - a->upper : 0, e = MIN(c + a->lower + 1, a->n);
        const double *entries = a->values + c * a->stride + a->upper + s - c;
        double weight = w[c];
        for (size_t r = s; r < e; r++)
            o[r] += fabs(entries[r - s]) * weight;
    }
    for (size_t r = 0; r < a->n; r++)
        o[r] = fmax(o[r], DBL_MIN);
}
typedef struct {
    size_t upper, stride;
    double *values;
} FactorWorkspace;
static size_t ix(const FactorWorkspace *a, size_t r, size_t c) {
    return c * a->stride + a->upper + r - c;
}
static double get(const FactorWorkspace *a, size_t r, size_t c) {
    return a->values[ix(a, r, c)];
}
int bac_band_factor(const BacBand *m, BacLu *a) {
    memset(a, 0, sizeof *a);
    size_t n = m->n;
    if (!n)
        return -1;
    if (m->npde) {
        a->n = n;
        return bac_abd_factor(m, &a->abd);
    }
    a->n = n;
    a->lower = m->lower;
    a->upper = MIN(m->upper + m->lower, n - 1);
    a->stride = a->lower + a->upper + 1;
    if (n > SIZE_MAX / (3 * sizeof(size_t)) || n > SIZE_MAX / a->stride / sizeof(double))
        return -1;
    FactorWorkspace work = {a->upper, a->stride, calloc(n * a->stride, sizeof(double))};
    a->pivots = malloc(3 * n * sizeof(size_t));
    if (!work.values || !a->pivots) {
        free(work.values);
        bac_lu_free(a);
        return -1;
    }
    for (size_t c = 0; c < n; c++) {
        size_t s = c > m->upper ? c - m->upper : 0, e = MIN(c + m->lower + 1, n);
        for (size_t r = s; r < e; r++) {
            double v = m->values[c * m->stride + m->upper + r - c];
            if (!isfinite(v)) {
                free(work.values);
                bac_lu_free(a);
                return -2;
            }
            work.values[ix(&work, r, c)] = v;
        }
    }
    for (size_t c = 0; c < n; c++) {
        size_t re = MIN(c + a->lower + 1, n), ce = MIN(c + a->upper + 1, n), p = c;
        for (size_t r = c + 1; r < re; r++)
            if (fabs(get(&work, r, c)) > fabs(get(&work, p, c)))
                p = r;
        double d = get(&work, p, c);
        if (d == 0 || !isfinite(d)) {
            free(work.values);
            bac_lu_free(a);
            return -3;
        }
        a->pivots[c] = p;
        if (p != c)
            for (size_t j = c; j < ce; j++) {
                size_t ai = ix(&work, c, j), bi = ix(&work, p, j);
                double t = work.values[ai];
                work.values[ai] = work.values[bi];
                work.values[bi] = t;
            }
        d = get(&work, c, c);
        for (size_t r = c + 1; r < re; r++) {
            size_t ri = ix(&work, r, c);
            double factor = work.values[ri] / d;
            work.values[ri] = factor;
            if (factor == 0.)
                continue;
            for (size_t j = c + 1; j < ce; j++)
                work.values[ix(&work, r, j)] -= factor * get(&work, c, j);
        }
    }
    a->lower_lengths = a->pivots + n;
    a->upper_lengths = a->pivots + 2 * n;
    for (size_t c = 0; c < n; c++) {
        size_t count = MIN(a->lower, n - c - 1);
        while (count && get(&work, c + count, c) == 0.)
            count--;
        a->lower_lengths[c] = count;
        a->lower_count += count;
    }
    for (size_t r = 0; r < n; r++) {
        size_t count = MIN(a->upper + 1, n - r);
        while (count > 1 && get(&work, r, r + count - 1) == 0.)
            count--;
        a->upper_lengths[r] = count;
        a->upper_count += count;
    }
    a->values = work.values;
    return 0;
}
int bac_dense_factor(const double *v, size_t n, BacLu *a) {
    memset(a, 0, sizeof *a);
    if (!n)
        return -1;
    if (n > SIZE_MAX / n / sizeof(double))
        return -1;
    a->n = n;
    a->dense = 1;
    a->values = malloc(n * n * sizeof(double));
    a->pivots = malloc(n * sizeof(size_t));
    if (!a->values || !a->pivots) {
        bac_lu_free(a);
        return -1;
    }
    memcpy(a->values, v, n * n * sizeof(double));
    for (size_t i = 0; i < n * n; i++)
        if (!isfinite(a->values[i])) {
            bac_lu_free(a);
            return -2;
        }
    for (size_t c = 0; c < n; c++) {
        size_t pivot = c;
        for (size_t r = c + 1; r < n; r++)
            if (fabs(a->values[r * n + c]) > fabs(a->values[pivot * n + c]))
                pivot = r;
        double diagonal = a->values[pivot * n + c];
        if (diagonal == 0. || !isfinite(diagonal)) {
            bac_lu_free(a);
            return -3;
        }
        a->pivots[c] = pivot;
        if (pivot != c)
            for (size_t j = c; j < n; j++) {
                double swap = a->values[c * n + j];
                a->values[c * n + j] = a->values[pivot * n + j];
                a->values[pivot * n + j] = swap;
            }
        diagonal = a->values[c * n + c];
        for (size_t r = c + 1; r < n; r++) {
            double factor = a->values[r * n + c] / diagonal;
            a->values[r * n + c] = factor;
            if (factor == 0.)
                continue;
            for (size_t j = c + 1; j < n; j++)
                a->values[r * n + j] -= factor * a->values[c * n + j];
        }
    }
    return 0;
}
int bac_lu_solve(const BacLu *a, double *r) {
    if (a->abd)
        return bac_abd_solve(a->abd, r);
    size_t n = a->n;
    for (size_t i = 0; i < n; i++)
        if (!isfinite(r[i]))
            return -2;
    if (a->dense) {
        for (size_t c = 0; c < n; c++) {
            size_t p = a->pivots[c];
            double swap = r[c];
            r[c] = r[p];
            r[p] = swap;
            for (size_t row = c + 1; row < n; row++)
                r[row] -= a->values[row * n + c] * r[c];
        }
        for (size_t row = n; row-- > 0;) {
            for (size_t c = row + 1; c < n; c++)
                r[row] -= a->values[row * n + c] * r[c];
            r[row] /= a->values[row * n + row];
            if (!isfinite(r[row]))
                return -2;
        }
        return 0;
    }
    for (size_t c = 0; c < n; c++) {
        size_t p = a->pivots[c];
        double t = r[c];
        r[c] = r[p];
        r[p] = t;
        double v = r[c];
        size_t count = a->lower_lengths[c];
        const double *factors = a->values + c * a->stride + a->upper + 1;
        for (size_t j = 0; j < count; j++)
            r[c + 1 + j] -= factors[j] * v;
    }
    for (size_t row = n; row-- > 0;) {
        double v = r[row];
        size_t count = a->upper_lengths[row];
        const double *factors = a->values + row * a->stride + a->upper;
        for (size_t j = 1; j < count; j++)
            v -= factors[j * (a->stride - 1)] * r[row + j];
        r[row] = v / factors[0];
        if (!isfinite(r[row]))
            return -2;
    }
    return 0;
}
void bac_lu_free(BacLu *a) {
    bac_abd_free(a->abd);
    free(a->values);
    free(a->pivots);
    memset(a, 0, sizeof *a);
}
int bac_lu_clone(BacLu *a, const BacLu *b) {
    *a = *b;
    if (b->abd) {
        a->abd = NULL;
        return bac_abd_clone(&a->abd, b->abd);
    }
    size_t count = b->dense ? b->n * b->n : b->n * b->stride;
    a->values = malloc(count * sizeof(double));
    a->pivots = malloc((b->dense ? 1 : 3) * b->n * sizeof(size_t));
    if (!a->values || !a->pivots) {
        bac_lu_free(a);
        return -1;
    }
    memcpy(a->values, b->values, count * sizeof(double));
    memcpy(a->pivots, b->pivots, (b->dense ? 1 : 3) * b->n * sizeof(size_t));
    if (!b->dense) {
        a->lower_lengths = a->pivots + b->n;
        a->upper_lengths = a->pivots + 2 * b->n;
    }
    return 0;
}
