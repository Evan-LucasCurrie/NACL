/* Pivoted band LU, COLROW/CRSLVE/CRSLVE2 and one-sided Jacobi SVD. */
#include "internal.h"
matrix bvp_matrix(context *ctx, size_t q, size_t left, size_t points) {
    matrix m = {0};
    m.q = q;
    m.left = left;
    m.n = q * points;
    m.intervals = points - 1;
    m.colrow = left < q;
    if (m.colrow)
        m.count = (2 * m.intervals + 1) * q * q;
    else {
        m.lower = MIN(left + q - 1, m.n - 1);
        m.upper = MIN(2 * q - left - 1, m.n - 1);
        m.width = m.lower + MIN(m.upper + m.lower, m.n - 1) + 1;
        m.count = m.n * m.width;
    }
    m.a = bvp_alloc(ctx, m.count, sizeof(double));
    m.pivot = bvp_alloc(ctx, m.n, sizeof(size_t));
    m.ends = bvp_alloc(ctx, m.n, sizeof(size_t));
    return m;
}
double *bvp_matrix_row(matrix *m, size_t row, size_t first) {
    if (!m->colrow)
        return m->a + row * m->width + first + m->lower - row;
    size_t q = m->q, l = m->left, k = m->intervals;
    if (row < l)
        return m->a + row * q + first;
    if (row < l + k * q) {
        size_t local = row - l;
        return m->a + l * q + local * 2 * q + first - (local / q) * q;
    }
    return m->a + l * q + k * 2 * q * q + (row - l - k * q) * q + first - k * q;
}
static void accept(context *c, double *largest, double magnitude) {
    if (!isfinite(magnitude) || *largest + magnitude == *largest)
        bvp_fail(c, BVP_SINGULAR);
    *largest = MAX(*largest, magnitude);
}
static void factor_band(context *c, matrix *m) {
    size_t n = m->n, lo = m->lower, w = m->width, last = 0, cs = w - 1;
    double *a = m->a;
    for (size_t k = 0; k < n; k++) {
        size_t lr = MIN(k + lo, n - 1), di = k * w + lo, p = k;
        double magnitude = fabs(a[di]);
        for (size_t r = k + 1; r <= lr; r++) {
            double v = fabs(a[di + (r - k) * cs]);
            if (v > magnitude) {
                magnitude = v;
                p = r;
            }
        }
        if (magnitude == 0 || !isfinite(magnitude))
            bvp_fail(c, BVP_SINGULAR);
        m->pivot[k] = p;
        last = MAX(last, MIN(p + m->upper, n - 1));
        m->ends[k] = last;
        size_t len = last - k + 1;
        if (p != k)
            for (size_t j = 0; j < len; j++)
                SWAP(double, a[di + j], a[p * w + lo - (p - k) + j]);
        double d = a[di];
        for (size_t r = k + 1; r <= lr; r++) {
            size_t rs = di + (r - k) * cs;
            double v = a[rs] / d;
            a[rs] = v;
            if (v != 0)
                for (size_t j = 1; j < len; j++)
                    a[rs + j] -= v * a[di + j];
        }
    }
    bvp_finite(c, a, m->count);
}
void bvp_factor(context *c, matrix *m) {
    bvp_finite(c, m->a, m->count);
    if (!m->colrow) {
        factor_band(c, m);
        return;
    }
    size_t q = m->q, l = m->left, r = q - l, k = m->intervals, w = 2 * q, bs = q * w, start = l * q,
           bottom = start + k * bs;
    double *a = m->a, largest = 0;
    size_t *p = m->pivot;
    for (size_t row = 0; row < l; row++) {
        size_t pivot = row;
        double mag = fabs(a[row * q + row]);
        for (size_t col = row + 1; col < q; col++)
            if (fabs(a[row * q + col]) > mag) {
                mag = fabs(a[row * q + col]);
                pivot = col;
            }
        accept(c, &largest, mag);
        p[row] = pivot;
        if (pivot != row) {
            for (size_t j = row; j < l; j++)
                SWAP(double, a[j * q + row], a[j * q + pivot]);
            for (size_t j = 0; j < q; j++)
                SWAP(double, a[start + j * w + row], a[start + j * w + pivot]);
        }
        double diagonal = a[row * q + row];
        for (size_t col = row + 1; col < q; col++) {
            double mult = a[row * q + col] / diagonal;
            a[row * q + col] = mult;
            for (size_t j = row + 1; j < l; j++)
                a[j * q + col] -= mult * a[j * q + row];
            for (size_t j = 0; j < q; j++)
                a[start + j * w + col] -= mult * a[start + j * w + row];
        }
    }
    for (size_t block = 0; block < k; block++) {
        size_t base = start + block * bs, global = block * q;
        for (size_t col = l; col < q; col++) {
            size_t row = col - l, pivot = row;
            double mag = fabs(a[base + row * w + col]);
            for (size_t j = row + 1; j < q; j++)
                if (fabs(a[base + j * w + col]) > mag) {
                    mag = fabs(a[base + j * w + col]);
                    pivot = j;
                }
            accept(c, &largest, mag);
            p[global + col] = global + l + pivot;
            if (pivot != row)
                for (size_t j = col; j < w; j++)
                    SWAP(double, a[base + row * w + j], a[base + pivot * w + j]);
            double diagonal = a[base + row * w + col];
            for (size_t j = row + 1; j < q; j++)
                a[base + j * w + col] /= diagonal;
            for (size_t j = col + 1; j < w; j++) {
                double v = a[base + row * w + j];
                for (size_t i = row + 1; i < q; i++)
                    a[base + i * w + j] -= v * a[base + i * w + col];
            }
        }
        int isbottom = block + 1 == k;
        size_t next = isbottom ? bottom : base + bs, nw = isbottom ? q : w, nr = isbottom ? r : q;
        for (size_t row = r; row < q; row++) {
            size_t col = row + l, pivot = col;
            double mag = fabs(a[base + row * w + col]);
            for (size_t j = col + 1; j < w; j++)
                if (fabs(a[base + row * w + j]) > mag) {
                    mag = fabs(a[base + row * w + j]);
                    pivot = j;
                }
            accept(c, &largest, mag);
            p[global + col] = global + pivot;
            if (pivot != col) {
                for (size_t i = row; i < q; i++)
                    SWAP(double, a[base + i * w + col], a[base + i * w + pivot]);
                for (size_t i = 0; i < nr; i++)
                    SWAP(double, a[next + i * nw + col - q], a[next + i * nw + pivot - q]);
            }
            double diagonal = a[base + row * w + col];
            for (size_t j = col + 1; j < w; j++) {
                double mult = a[base + row * w + j] / diagonal;
                a[base + row * w + j] = mult;
                for (size_t i = row + 1; i < q; i++)
                    a[base + i * w + j] -= mult * a[base + i * w + col];
                for (size_t i = 0; i < nr; i++)
                    a[next + i * nw + j - q] -= mult * a[next + i * nw + col - q];
            }
        }
    }
    for (size_t col = l; col + 1 < q; col++) {
        size_t row = col - l, pivot = row;
        double mag = fabs(a[bottom + row * q + col]);
        for (size_t j = row + 1; j < r; j++)
            if (fabs(a[bottom + j * q + col]) > mag) {
                mag = fabs(a[bottom + j * q + col]);
                pivot = j;
            }
        accept(c, &largest, mag);
        p[k * q + col] = k * q + l + pivot;
        if (pivot != row)
            for (size_t j = col; j < q; j++)
                SWAP(double, a[bottom + row * q + j], a[bottom + pivot * q + j]);
        double diagonal = a[bottom + row * q + col];
        for (size_t j = row + 1; j < r; j++)
            a[bottom + j * q + col] /= diagonal;
        for (size_t j = col + 1; j < q; j++) {
            double v = a[bottom + row * q + j];
            for (size_t i = row + 1; i < r; i++)
                a[bottom + i * q + j] -= v * a[bottom + i * q + col];
        }
    }
    accept(c, &largest, fabs(a[bottom + (r - 1) * q + q - 1]));
    bvp_finite(c, a, m->count);
}
static void solve_band(context *c, const matrix *m, double *x, int transpose) {
    size_t n = m->n, lo = m->lower, w = m->width, cs = w - 1;
    const double *a = m->a;
    if (!transpose) {
        for (size_t k = 0; k < n; k++) {
            SWAP(double, x[k], x[m->pivot[k]]);
            double v = x[k];
            for (size_t j = k + 1; j < MIN(k + 1 + lo, n); j++)
                x[j] -= a[k * w + lo + (j - k) * cs] * v;
        }
        for (size_t row = n; row-- > 0;) {
            double v = x[row];
            for (size_t j = row + 1; j <= m->ends[row]; j++)
                v -= a[row * w + lo + j - row] * x[j];
            x[row] = v / a[row * w + lo];
        }
    } else {
        for (size_t row = 0; row < n; row++) {
            x[row] /= a[row * w + lo];
            double v = x[row];
            for (size_t j = row + 1; j <= m->ends[row]; j++)
                x[j] -= a[row * w + lo + j - row] * v;
        }
        for (size_t k = n; k-- > 0;) {
            double v = x[k];
            for (size_t j = k + 1; j < MIN(k + 1 + lo, n); j++)
                v -= a[k * w + lo + (j - k) * cs] * x[j];
            x[k] = v;
            SWAP(double, x[k], x[m->pivot[k]]);
        }
    }
    bvp_finite(c, x, n);
}
void bvp_linear(context *c, const matrix *m, double *x, int tr) {
    bvp_finite(c, x, m->n);
    if (!m->colrow) {
        solve_band(c, m, x, tr);
        return;
    }
    size_t q = m->q, l = m->left, r = q - l, k = m->intervals, w = 2 * q, bs = q * w, start = l * q,
           bottom = start + k * bs;
    const double *a = m->a;
    const size_t *p = m->pivot;
    if (!tr) {
        for (size_t col = 0; col < l; col++) {
            x[col] /= a[col * q + col];
            double v = -x[col];
            for (size_t row = col + 1; row < l; row++)
                x[row] += a[row * q + col] * v;
        }
        for (size_t b = 0; b < k; b++) {
            size_t base = start + b * bs, g = b * q;
            for (size_t col = 0; col < l; col++) {
                double v = -x[g + col];
                for (size_t row = 0; row < q; row++)
                    x[g + l + row] += a[base + row * w + col] * v;
            }
            for (size_t col = l; col < q; col++) {
                size_t index = g + col;
                SWAP(double, x[index], x[p[index]]);
                double v = -x[index];
                for (size_t row = col - l + 1; row < q; row++)
                    x[g + l + row] += a[base + row * w + col] * v;
            }
            for (size_t col = q; col < q + l; col++) {
                size_t index = g + col;
                x[index] /= a[base + (col - l) * w + col];
                double v = -x[index];
                for (size_t row = col - l + 1; row < q; row++)
                    x[g + l + row] += a[base + row * w + col] * v;
            }
        }
        size_t g = k * q;
        for (size_t col = 0; col < l; col++) {
            double v = -x[g + col];
            for (size_t row = 0; row < r; row++)
                x[g + l + row] += a[bottom + row * q + col] * v;
        }
        for (size_t col = l; col + 1 < q; col++) {
            size_t index = g + col;
            SWAP(double, x[index], x[p[index]]);
            double v = -x[index];
            for (size_t row = col - l + 1; row < r; row++)
                x[g + l + row] += a[bottom + row * q + col] * v;
        }
        for (size_t row = r; row-- > 0;) {
            size_t col = l + row, index = g + col;
            x[index] /= a[bottom + row * q + col];
            double v = -x[index];
            for (size_t j = 0; j < row; j++)
                x[g + l + j] += a[bottom + j * q + col] * v;
        }
        for (size_t b = k; b-- > 0;) {
            size_t base = start + b * bs, g = b * q;
            for (size_t row = q; row-- > r;) {
                size_t col = row + l, index = g + col;
                double v = x[index];
                for (size_t j = col + 1; j < w; j++)
                    v -= a[base + row * w + j] * x[g + j];
                x[index] = v;
                SWAP(double, x[index], x[p[index]]);
            }
            for (size_t col = q; col < w; col++) {
                double v = -x[g + col];
                for (size_t row = 0; row < r; row++)
                    x[g + l + row] += a[base + row * w + col] * v;
            }
            for (size_t row = r; row-- > 0;) {
                size_t col = row + l, index = g + col;
                x[index] /= a[base + row * w + col];
                double v = -x[index];
                for (size_t j = 0; j < row; j++)
                    x[g + l + j] += a[base + j * w + col] * v;
            }
        }
        for (size_t row = l; row-- > 0;) {
            double v = x[row];
            for (size_t col = row + 1; col < q; col++)
                v -= a[row * q + col] * x[col];
            x[row] = v;
            SWAP(double, x[row], x[p[row]]);
        }
    } else {
        for (size_t row = 0; row < l; row++) {
            SWAP(double, x[row], x[p[row]]);
            double v = -x[row];
            for (size_t col = row + 1; col < q; col++)
                x[col] += v * a[row * q + col];
        }
        for (size_t b = 0; b < k; b++) {
            size_t base = start + b * bs, e = b * q + l, next = (b + 1) * q;
            for (size_t row = 0; row < r; row++) {
                size_t index = e + row;
                for (size_t j = 0; j < row; j++)
                    x[index] -= a[base + j * w + l + row] * x[e + j];
                x[index] /= a[base + row * w + l + row];
            }
            for (size_t col = 0; col < q; col++)
                for (size_t row = 0; row < r; row++)
                    x[next + col] -= a[base + row * w + q + col] * x[e + row];
            for (size_t row = 0; row < l; row++) {
                size_t index = next + row;
                SWAP(double, x[index], x[p[index]]);
                double v = -x[index];
                for (size_t col = row + 1; col < q; col++)
                    x[next + col] += v * a[base + (r + row) * w + q + col];
            }
        }
        size_t node = k * q, e = node + l;
        for (size_t row = 0; row < r; row++) {
            size_t index = e + row;
            for (size_t j = 0; j < row; j++)
                x[index] -= a[bottom + j * q + l + row] * x[e + j];
            x[index] /= a[bottom + row * q + l + row];
        }
        for (size_t row = r - 1; row-- > 0;) {
            size_t index = e + row;
            for (size_t j = row + 1; j < r; j++)
                x[index] -= a[bottom + j * q + l + row] * x[e + j];
            SWAP(double, x[index], x[p[index]]);
        }
        for (size_t col = l; col-- > 0;)
            for (size_t row = 0; row < r; row++)
                x[node + col] -= a[bottom + row * q + col] * x[e + row];
        for (size_t b = k; b-- > 0;) {
            size_t base = start + b * bs, node = b * q, next = node + q, e = node + l;
            for (size_t col = l; col-- > 0;) {
                size_t index = next + col;
                for (size_t j = l; j-- > col + 1;)
                    x[index] -= a[base + (r + j) * w + q + col] * x[next + j];
                x[index] /= a[base + (r + col) * w + q + col];
            }
            for (size_t row = r; row-- > 0;) {
                size_t index = e + row;
                for (size_t j = row + 1; j < q; j++)
                    x[index] -= a[base + j * w + l + row] * x[e + j];
                SWAP(double, x[index], x[p[index]]);
            }
            for (size_t col = 0; col < l; col++)
                for (size_t row = 0; row < q; row++)
                    x[node + col] -= a[base + row * w + col] * x[e + row];
        }
        for (size_t col = l; col-- > 0;) {
            for (size_t row = col + 1; row < l; row++)
                x[col] -= a[row * q + col] * x[row];
            x[col] /= a[col * q + col];
        }
    }
    bvp_finite(c, x, m->n);
}
static double column_norm(const double *a, size_t n, size_t col) {
    double v = 0;
    for (size_t row = 0; row < n; row++)
        v = hypot(v, a[row * n + col]);
    return v;
}
double *bvp_pinv(context *c, const double *a, size_t n) {
    double *inverse = bvp_alloc(c, n * n, sizeof(double));
    allocation *mark = c->memory;
    double scale = bvp_maximum(a, n * n);
    if (scale == 0)
        return inverse;
    double *b = bvp_alloc(c, n * n, sizeof(double)), *v = bvp_alloc(c, n * n, sizeof(double));
    for (size_t i = 0; i < n * n; i++)
        b[i] = a[i] / scale;
    for (size_t i = 0; i < n; i++)
        v[i * n + i] = 1;
    int converged = 0;
    for (size_t iter = 0; iter < 30 + 6 * n; iter++) {
        double largest = 0;
        for (size_t col = 0; col < n; col++)
            largest = MAX(largest, column_norm(b, n, col));
        double negligible = n * DBL_EPSILON * largest;
        int rotated = 0;
        for (size_t p = 0; p < n; p++)
            for (size_t q = p + 1; q < n; q++) {
                double np = column_norm(b, n, p), nq = column_norm(b, n, q);
                if (np <= negligible || nq <= negligible)
                    continue;
                double corr = 0;
                for (size_t row = 0; row < n; row++)
                    corr += (b[row * n + p] / np) * (b[row * n + q] / nq);
                if (fabs(corr) <= 4 * DBL_EPSILON)
                    continue;
                double pair = MAX(np, nq), ps = np / pair, qs = nq / pair, cross = corr * ps * qs,
                       diff = qs * qs - ps * ps, tau = diff / (2 * cross);
                double tangent = isfinite(tau) ? copysign(1.0, tau) / (fabs(tau) + hypot(tau, 1.0))
                                               : cross / diff,
                       co = 1 / hypot(tangent, 1.0), si = co * tangent;
                for (size_t row = 0; row < n; row++) {
                    double bp = b[row * n + p], bq = b[row * n + q];
                    b[row * n + p] = co * bp - si * bq;
                    b[row * n + q] = si * bp + co * bq;
                    double vp = v[row * n + p], vq = v[row * n + q];
                    v[row * n + p] = co * vp - si * vq;
                    v[row * n + q] = si * vp + co * vq;
                }
                rotated = 1;
            }
        if (!rotated) {
            converged = 1;
            break;
        }
    }
    if (!converged)
        bvp_fail(c, BVP_INVALID_INPUT);
    double largest = 0;
    for (size_t k = 0; k < n; k++)
        largest = MAX(largest, column_norm(b, n, k));
    double cutoff = n * DBL_EPSILON * largest;
    for (size_t k = 0; k < n; k++) {
        double sigma = column_norm(b, n, k);
        if (sigma <= cutoff)
            continue;
        for (size_t row = 0; row < n; row++) {
            double weight = (v[row * n + k] / sigma) / scale;
            for (size_t col = 0; col < n; col++)
                inverse[row * n + col] += weight * (b[col * n + k] / sigma);
        }
    }
    bvp_finite(c, inverse, n * n);
    bvp_release(c, mark);
    return inverse;
}
