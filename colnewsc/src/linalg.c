/* Partial-pivoting dense and band LU factorizations. */
#include "internal.h"
void cn_lu_free(cn_lu *m) {
    if (!m->borrowed) {
        if (m->storage)
            free(m->storage);
        else {
            free(m->data);
            free(m->pivots);
        }
    }
    memset(m, 0, sizeof *m);
}
static colnewsc_status reserve(cn_lu *m, size_t values, size_t pivots) {
    if (m->data && m->pivots && values <= m->data_capacity && pivots <= m->pivot_capacity)
        return COLNEWSC_OK;
    if (values > (SIZE_MAX - _Alignof(size_t)) / sizeof(double))
        return COLNEWSC_ALLOCATION_FAILED;
    size_t bytes = values * sizeof(double);
    size_t offset = (bytes + _Alignof(size_t) - 1) & ~((size_t)_Alignof(size_t) - 1);
    if (pivots > (SIZE_MAX - offset) / sizeof(size_t))
        return COLNEWSC_ALLOCATION_FAILED;
    void *storage = malloc(offset + pivots * sizeof(size_t));
    if (!storage)
        return COLNEWSC_ALLOCATION_FAILED;
    if (!m->borrowed) {
        if (m->storage)
            free(m->storage);
        else {
            free(m->data);
            free(m->pivots);
        }
    }
    m->storage = storage;
    m->data = storage;
    m->pivots = (size_t *)((unsigned char *)storage + offset);
    m->data_capacity = values;
    m->pivot_capacity = pivots;
    m->borrowed = 0;
    return COLNEWSC_OK;
}
colnewsc_status cn_band_reset(cn_lu *m, size_t n, size_t lower, size_t upper) {
    if (!n)
        return COLNEWSC_INVALID_INPUT;
    m->n = n;
    m->band = 1;
    m->lower = CN_MIN(lower, n - 1);
    m->upper = CN_MIN(upper, n - 1);
    m->stride = 2 * m->lower + m->upper + 1;
    if (n > SIZE_MAX / m->stride)
        return COLNEWSC_ALLOCATION_FAILED;
    CN_TRY(reserve(m, n * m->stride, n));
    memset(m->data, 0, n * m->stride * sizeof(double));
    return COLNEWSC_OK;
}
colnewsc_status cn_band_new(cn_lu *m, size_t n, size_t lower, size_t upper) {
    memset(m, 0, sizeof *m);
    return cn_band_reset(m, n, lower, upper);
}
void cn_band_set(cn_lu *m, size_t i, size_t j, double v) {
    m->data[i * m->stride + j + m->lower - i] = v;
}
static double get(const cn_lu *m, size_t i, size_t j) {
    return m->data[i * m->stride + (m->band ? j + m->lower - i : j)];
}
static void set(cn_lu *m, size_t i, size_t j, double v) {
    m->data[i * m->stride + (m->band ? j + m->lower - i : j)] = v;
}
static colnewsc_status factor_finite(cn_lu *m) {
    size_t n = m->n;
    if (!m->band) {
        for (size_t k = 0; k < n; k++) {
            size_t pivot = k;
            for (size_t i = k + 1; i < n; i++)
                if (fabs(m->data[i * n + k]) >= fabs(m->data[pivot * n + k]))
                    pivot = i;
            double *row = m->data + k * n, diag = m->data[pivot * n + k];
            if (!isfinite(diag))
                return COLNEWSC_NON_FINITE;
            if (fabs(diag) <= DBL_MIN)
                return COLNEWSC_SINGULAR_JACOBIAN;
            m->pivots[k] = pivot;
            if (pivot != k) {
                double *other = m->data + pivot * n;
                for (size_t j = k; j < n; j++) {
                    double t = row[j];
                    row[j] = other[j];
                    other[j] = t;
                }
            }
            for (size_t i = k + 1; i < n; i++) {
                double *target = m->data + i * n;
                double multiplier = target[k] / diag;
                target[k] = multiplier;
                if (multiplier != 0)
                    for (size_t j = k + 1; j < n; j++)
                        target[j] -= multiplier * row[j];
            }
        }
        return COLNEWSC_OK;
    }
    for (size_t k = 0; k < n; k++) {
        size_t end = m->band ? CN_MIN(k + m->lower + 1, n) : n, pivot = k;
        for (size_t i = k + 1; i < end; i++)
            if (fabs(get(m, i, k)) >= fabs(get(m, pivot, k)))
                pivot = i;
        double diag = get(m, pivot, k);
        if (!isfinite(diag))
            return COLNEWSC_NON_FINITE;
        if (fabs(diag) <= DBL_MIN)
            return COLNEWSC_SINGULAR_JACOBIAN;
        m->pivots[k] = pivot;
        size_t last = m->band ? CN_MIN(k + m->upper + m->lower + 1, n) : n;
        if (pivot != k)
            for (size_t j = k; j < last; j++) {
                double t = get(m, k, j);
                set(m, k, j, get(m, pivot, j));
                set(m, pivot, j, t);
            }
        for (size_t i = k + 1; i < end; i++) {
            double multiplier = get(m, i, k) / diag;
            set(m, i, k, multiplier);
            if (multiplier != 0)
                for (size_t j = k + 1; j < last; j++)
                    set(m, i, j, get(m, i, j) - multiplier * get(m, k, j));
        }
    }
    return COLNEWSC_OK;
}
colnewsc_status cn_lu_factor(cn_lu *m) {
    if (!cn_finite(m->data, m->n * m->stride))
        return COLNEWSC_NON_FINITE;
    return factor_finite(m);
}
colnewsc_status cn_stage_factor_reuse(const double *data, size_t n, cn_lu *m) {
    size_t lower = 0, upper = 0;
    if (!cn_finite(data, n * n)) {
        return COLNEWSC_NON_FINITE;
    }
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (data[i * n + j] != 0) {
                if (i > j)
                    lower = CN_MAX(lower, i - j);
                else
                    upper = CN_MAX(upper, j - i);
                if (2 * lower + upper + 1 >= n)
                    goto dense;
            }
    if (2 * lower + upper + 1 < n) {
        CN_TRY(cn_band_reset(m, n, lower, upper));
        for (size_t i = 0; i < n; i++)
            for (size_t j = i > lower ? i - lower : 0; j < CN_MIN(i + upper + 1, n); j++)
                cn_band_set(m, i, j, data[i * n + j]);
    } else {
    dense:
        m->n = n;
        m->stride = n;
        m->band = 0;
        m->lower = m->upper = 0;
        CN_TRY(reserve(m, n * n, n));
        memcpy(m->data, data, n * n * sizeof(double));
    }
    /* Input was checked before packing; the copied entries and zero fill
     * cannot introduce non-finite values before elimination starts. */
    return factor_finite(m);
}
colnewsc_status cn_stage_factor(double *data, size_t n, cn_lu *m) {
    memset(m, 0, sizeof *m);
    colnewsc_status e = cn_stage_factor_reuse(data, n, m);
    free(data);
    if (e)
        cn_lu_free(m);
    return e;
}
static colnewsc_status solve_many(const cn_lu *m, double *rhs, size_t columns) {
    size_t n = m->n;
    for (size_t k = 0; k < n; k++) {
        size_t p = m->pivots[k];
        if (p != k)
            for (size_t q = 0; q < columns; q++) {
                double t = rhs[k * columns + q];
                rhs[k * columns + q] = rhs[p * columns + q];
                rhs[p * columns + q] = t;
            }
        size_t end = m->band ? CN_MIN(k + m->lower + 1, n) : n;
        for (size_t i = k + 1; i < end; i++) {
            double a = get(m, i, k);
            if (a != 0)
                for (size_t q = 0; q < columns; q++)
                    rhs[i * columns + q] -= a * rhs[k * columns + q];
        }
    }
    for (size_t i = n; i-- > 0;) {
        size_t end = m->band ? CN_MIN(i + m->upper + m->lower + 1, n) : n;
        for (size_t j = i + 1; j < end; j++) {
            double a = get(m, i, j);
            if (a != 0)
                for (size_t q = 0; q < columns; q++)
                    rhs[i * columns + q] -= a * rhs[j * columns + q];
        }
        double diag = get(m, i, i);
        for (size_t q = 0; q < columns; q++)
            rhs[i * columns + q] /= diag;
    }
    return cn_finite(rhs, n * columns) ? COLNEWSC_OK : COLNEWSC_NON_FINITE;
}
/* The most common collocation systems have one or two state components.
 * Constant RHS counts let the compiler remove short inner loops without
 * changing the pivot or accumulation order used by the general kernel. */
colnewsc_status cn_lu_solve(const cn_lu *m, double *rhs, size_t columns) {
    if (columns == 1) {
        size_t n = m->n;
        for (size_t k = 0; k < n; k++) {
            size_t p = m->pivots[k];
            if (p != k) {
                double value = rhs[k];
                rhs[k] = rhs[p];
                rhs[p] = value;
            }
            size_t end = m->band ? CN_MIN(k + m->lower + 1, n) : n;
            for (size_t i = k + 1; i < end; i++) {
                double a = get(m, i, k);
                if (a != 0)
                    rhs[i] -= a * rhs[k];
            }
        }
        for (size_t i = n; i-- > 0;) {
            size_t end = m->band ? CN_MIN(i + m->upper + m->lower + 1, n) : n;
            for (size_t j = i + 1; j < end; j++) {
                double a = get(m, i, j);
                if (a != 0)
                    rhs[i] -= a * rhs[j];
            }
            rhs[i] /= get(m, i, i);
        }
        return cn_finite(rhs, n) ? COLNEWSC_OK : COLNEWSC_NON_FINITE;
    }
    if (columns == 2)
        return solve_many(m, rhs, 2);
    if (columns == 3)
        return solve_many(m, rhs, 3);
    return solve_many(m, rhs, columns);
}
