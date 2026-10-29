/* Source-order LAMDEC/LAMSOL block elimination and substitution. */
#include "abd.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define T(i, j) (a->top[((j) - 1) * p + (i) - 1])
#define B(i, j) (a->bottom[((j) - 1) * p + (i) - 1])
#define A(i, j, k) (a->blocks[(((k) - 1) * (m + 2 * p) + (j) - 1) * m + (i) - 1])
#define SWAP(x, y)                                                                                 \
    do {                                                                                           \
        double swap_value = (x);                                                                   \
        (x) = (y);                                                                                 \
        (y) = swap_value;                                                                          \
    } while (0)
void bac_abd_free(BacAbd *a) {
    if (a) {
        free(a->top);
        free(a->blocks);
        free(a->bottom);
        free(a->pivots);
        free(a);
    }
}
int bac_abd_factor(const BacBand *matrix, BacAbd **out) {
    *out = NULL;
    size_t p = matrix->npde, m = p * matrix->kcol, nblocks = matrix->nint, overlap = 2 * p,
           columns = m + overlap, n = m * nblocks + overlap, nrowel = m - p;
    if (!p || matrix->kcol < 2 || !nblocks || matrix->n != n)
        return -1;
    for (size_t i = 0; i < matrix->n * matrix->stride; i++)
        if (!isfinite(matrix->values[i]))
            return -2;
    BacAbd *a = calloc(1, sizeof *a);
    if (!a)
        return -1;
    a->p = p;
    a->m = m;
    a->nblocks = nblocks;
    a->top = calloc(p * overlap, sizeof(double));
    a->blocks = calloc(m * columns * nblocks, sizeof(double));
    a->bottom = calloc(p * overlap, sizeof(double));
    a->pivots = calloc(n + 1, sizeof(size_t));
    if (!a->top || !a->blocks || !a->bottom || !a->pivots)
        goto error;
    for (size_t j = 1; j <= overlap; j++)
        for (size_t i = 1; i <= p; i++) {
            T(i, j) = bac_band_get(matrix, i - 1, j - 1);
            B(i, j) = bac_band_get(matrix, p + m * nblocks + i - 1, m * nblocks + j - 1);
        }
    for (size_t k = 1; k <= nblocks; k++)
        for (size_t j = 1; j <= columns; j++)
            for (size_t i = 1; i <= m; i++)
                A(i, j, k) = bac_band_get(matrix, p + (k - 1) * m + i - 1, (k - 1) * m + j - 1);
    double pivmax = 0.;
    for (size_t i = 1; i <= p; i++) {
        size_t ipvt = i;
        double colmax = fabs(T(i, i));
        for (size_t j = i + 1; j <= overlap; j++) {
            double temp = fabs(T(i, j));
            if (temp > colmax) {
                ipvt = j;
                colmax = temp;
            }
        }
        if (pivmax + colmax == pivmax)
            goto singular;
        pivmax = fmax(colmax, pivmax);
        a->pivots[i] = ipvt;
        if (ipvt != i) {
            for (size_t l = i; l <= p; l++)
                SWAP(T(l, ipvt), T(l, i));
            for (size_t l = 1; l <= m; l++)
                SWAP(A(l, ipvt, 1), A(l, i, 1));
        }
        double colpiv = T(i, i);
        for (size_t l = i + 1; l <= p; l++)
            T(l, i) /= colpiv;
        for (size_t l = 1; l <= m; l++)
            A(l, i, 1) /= colpiv;
        for (size_t j = i + 1; j <= overlap; j++) {
            double colmlt = T(i, j);
            for (size_t l = i + 1; l <= p; l++)
                T(l, j) -= colmlt * T(l, i);
            for (size_t l = 1; l <= m; l++)
                A(l, j, 1) -= colmlt * A(l, i, 1);
        }
    }
    size_t incr = 0;
    for (size_t k = 1; k <= nblocks; k++) {
        for (size_t j = p + 1; j <= m; j++) {
            size_t jminn = j - p, ipvt = jminn;
            double rowmax = fabs(A(jminn, j, k));
            for (size_t i = jminn + 1; i <= m; i++) {
                double temp = fabs(A(i, j, k));
                if (temp > rowmax) {
                    ipvt = i;
                    rowmax = temp;
                }
            }
            if (pivmax + rowmax == pivmax)
                goto singular;
            pivmax = fmax(rowmax, pivmax);
            a->pivots[incr + j] = incr + ipvt + p;
            if (ipvt != jminn)
                for (size_t l = j; l <= columns; l++)
                    SWAP(A(ipvt, l, k), A(jminn, l, k));
            double rowpiv = A(jminn, j, k);
            for (size_t i = jminn + 1; i <= m; i++)
                A(i, j, k) /= rowpiv;
            for (size_t l = j + 1; l <= columns; l++) {
                double rowmlt = A(jminn, l, k);
                for (size_t i = jminn + 1; i <= m; i++)
                    A(i, l, k) -= rowmlt * A(i, j, k);
            }
        }
        for (size_t i = nrowel + 1; i <= m; i++) {
            size_t iplusn = i + p, ipvt = iplusn;
            double colmax = fabs(A(i, ipvt, k));
            for (size_t j = iplusn + 1; j <= columns; j++) {
                double temp = fabs(A(i, j, k));
                if (temp > colmax) {
                    ipvt = j;
                    colmax = temp;
                }
            }
            if (pivmax + colmax == pivmax)
                goto singular;
            pivmax = fmax(colmax, pivmax);
            a->pivots[incr + iplusn] = incr + ipvt;
            size_t irwblk = iplusn - m;
            if (ipvt != iplusn) {
                for (size_t l = i; l <= m; l++)
                    SWAP(A(l, ipvt, k), A(l, iplusn, k));
                size_t ipvblk = ipvt - m;
                if (k != nblocks) {
                    for (size_t l = 1; l <= m; l++)
                        SWAP(A(l, ipvblk, k + 1), A(l, irwblk, k + 1));
                } else
                    for (size_t l = 1; l <= p; l++)
                        SWAP(B(l, ipvblk), B(l, irwblk));
            }
            double colpiv = A(i, iplusn, k);
            for (size_t l = i + 1; l <= m; l++)
                A(l, iplusn, k) /= colpiv;
            if (k != nblocks) {
                for (size_t l = 1; l <= m; l++)
                    A(l, irwblk, k + 1) /= colpiv;
            } else
                for (size_t l = 1; l <= p; l++)
                    B(l, irwblk) /= colpiv;
            for (size_t j = iplusn + 1; j <= columns; j++) {
                double colmlt = A(i, j, k);
                for (size_t l = i + 1; l <= m; l++)
                    A(l, j, k) -= colmlt * A(l, iplusn, k);
                size_t jrwblk = j - m;
                if (k != nblocks) {
                    for (size_t l = 1; l <= m; l++)
                        A(l, jrwblk, k + 1) -= colmlt * A(l, irwblk, k + 1);
                } else
                    for (size_t l = 1; l <= p; l++)
                        B(l, jrwblk) -= colmlt * B(l, irwblk);
            }
        }
        incr += m;
    }
    for (size_t j = p + 1; j < overlap; j++) {
        size_t jminn = j - p, ipvt = jminn;
        double rowmax = fabs(B(jminn, j));
        for (size_t i = jminn + 1; i <= p; i++) {
            double temp = fabs(B(i, j));
            if (temp > rowmax) {
                ipvt = i;
                rowmax = temp;
            }
        }
        if (pivmax + rowmax == pivmax)
            goto singular;
        pivmax = fmax(rowmax, pivmax);
        a->pivots[incr + j] = incr + ipvt + p;
        if (ipvt != jminn)
            for (size_t l = j; l <= overlap; l++)
                SWAP(B(ipvt, l), B(jminn, l));
        double rowpiv = B(jminn, j);
        for (size_t i = jminn + 1; i <= p; i++)
            B(i, j) /= rowpiv;
        for (size_t l = j + 1; l <= overlap; l++) {
            double rowmlt = B(jminn, l);
            for (size_t i = jminn + 1; i <= p; i++)
                B(i, l) -= rowmlt * B(i, j);
        }
    }
    if (pivmax + fabs(B(p, overlap)) == pivmax)
        goto singular;
    *out = a;
    return 0;
singular:
    bac_abd_free(a);
    return -3;
error:
    bac_abd_free(a);
    return -1;
}
int bac_abd_solve(const BacAbd *a, double *rhs) {
    size_t p = a->p, m = a->m, overlap = 2 * p, columns = m + overlap, nrowel = m - p,
           n = m * a->nblocks + overlap;
    for (size_t i = 0; i < n; i++)
        if (!isfinite(rhs[i]))
            return -2;
    /* Both temporary vectors have the same lifetime. Keep the source's zero
     * initialization, but use a single allocation for each substitution. */
    if (n >= SIZE_MAX / (2 * sizeof(double)))
        return -1;
    double *b = calloc(2 * (n + 1), sizeof(double));
    if (!b) {
        return -1;
    }
    double *x = b + n + 1;
    memcpy(b + 1, rhs, n * sizeof(double));
    for (size_t j = 1; j < p; j++) {
        double bj = -b[j];
        for (size_t i = j + 1; i <= p; i++)
            b[i] += T(i, j) * bj;
    }
    size_t incr = 0;
    for (size_t k = 1; k <= a->nblocks; k++) {
        size_t incrtp = incr + p;
        for (size_t j = 1; j <= p; j++) {
            double bj = -b[incr + j];
            for (size_t i = 1; i <= m; i++)
                b[incrtp + i] += A(i, j, k) * bj;
        }
        for (size_t j = p + 1; j <= m; j++) {
            size_t incrj = incr + j, jpivot = a->pivots[incrj];
            if (jpivot != incrj)
                SWAP(b[incrj], b[jpivot]);
            double bincrj = -b[incrj];
            for (size_t i = j - p + 1; i <= m; i++)
                b[incrtp + i] += A(i, j, k) * bincrj;
        }
        for (size_t j = m + 1; j < m + p; j++) {
            double bj = -b[incr + j];
            for (size_t i = j - p + 1; i <= m; i++)
                b[incrtp + i] += A(i, j, k) * bj;
        }
        incr += m;
    }
    size_t incrtp = incr + p;
    for (size_t j = 1; j <= p; j++) {
        double bj = -b[incr + j];
        for (size_t i = 1; i <= p; i++)
            b[incrtp + i] += B(i, j) * bj;
    }
    for (size_t j = p + 1; j < overlap; j++) {
        size_t incrj = incr + j, jpivot = a->pivots[incrj];
        if (jpivot != incrj)
            SWAP(b[incrj], b[jpivot]);
        double bincrj = -b[incrj];
        for (size_t i = j - p + 1; i <= p; i++)
            b[incrtp + i] += B(i, j) * bincrj;
    }
    for (size_t ll = 1; ll <= p; ll++) {
        size_t j = overlap + 1 - ll, incrj = incr + j;
        x[incrj] = b[incrj] / B(p + 1 - ll, j);
        if (ll != p) {
            double xincrj = -x[incrj];
            for (size_t i = 1; i <= p - ll; i++)
                b[incrtp + i] += B(i, j) * xincrj;
        }
    }
    for (size_t l = 1; l <= a->nblocks; l++) {
        size_t k = a->nblocks + 1 - l;
        incr -= m;
        for (size_t l1 = nrowel + 1; l1 <= m; l1++) {
            size_t i = m + nrowel + 1 - l1, iplusn = i + p, incrn = incr + iplusn;
            double dotprd = b[incrn];
            for (size_t j = iplusn + 1; j <= columns; j++)
                dotprd -= A(i, j, k) * x[incr + j];
            x[incrn] = dotprd / A(i, iplusn, k);
            size_t ipvtn = a->pivots[incrn];
            if (ipvtn != incrn)
                SWAP(x[incrn], x[ipvtn]);
        }
        size_t ip = incr + p;
        for (size_t j = m + 1; j <= columns; j++) {
            double xincrj = -x[incr + j];
            for (size_t i = 1; i <= nrowel; i++)
                b[ip + i] += A(i, j, k) * xincrj;
        }
        for (size_t ll = 1; ll <= nrowel; ll++) {
            size_t j = m + 1 - ll, incrj = incr + j;
            x[incrj] = b[incrj] / A(nrowel + 1 - ll, j, k);
            if (ll != nrowel) {
                double xincrj = -x[incrj];
                for (size_t i = 1; i <= nrowel - ll; i++)
                    b[ip + i] += A(i, j, k) * xincrj;
            }
        }
    }
    for (size_t l = 1; l <= p; l++) {
        size_t i = p + 1 - l;
        double dotprd = b[i];
        for (size_t j = i + 1; j <= overlap; j++)
            dotprd -= T(i, j) * x[j];
        x[i] = dotprd / T(i, i);
        size_t ipvti = a->pivots[i];
        if (i != ipvti)
            SWAP(x[i], x[ipvti]);
    }
    int result = 0;
    for (size_t i = 1; i <= n; i++)
        if (!isfinite(x[i])) {
            result = -2;
            break;
        }
    if (!result)
        memcpy(rhs, x + 1, n * sizeof(double));
    free(b);
    return result;
}
int bac_abd_clone(BacAbd **out, const BacAbd *a) {
    *out = NULL;
    BacAbd *d = calloc(1, sizeof *d);
    if (!d)
        return -1;
    *d = *a;
    size_t edge = 2 * a->p * a->p, mid = a->m * (a->m + 2 * a->p) * a->nblocks,
           n = a->m * a->nblocks + 2 * a->p;
    d->top = malloc(edge * sizeof(double));
    d->bottom = malloc(edge * sizeof(double));
    d->blocks = malloc(mid * sizeof(double));
    d->pivots = malloc((n + 1) * sizeof(size_t));
    if (!d->top || !d->bottom || !d->blocks || !d->pivots) {
        bac_abd_free(d);
        return -1;
    }
    memcpy(d->top, a->top, edge * sizeof(double));
    memcpy(d->bottom, a->bottom, edge * sizeof(double));
    memcpy(d->blocks, a->blocks, mid * sizeof(double));
    memcpy(d->pivots, a->pivots, (n + 1) * sizeof(size_t));
    *out = d;
    return 0;
}
size_t bac_abd_workspace(const BacLu *lu, double *out) {
    const BacAbd *a = lu->abd;
    if (!a)
        return 0;
    size_t edge = 2 * a->p * a->p, mid = a->m * (a->m + 2 * a->p) * a->nblocks;
    if (out) {
        memcpy(out, a->top, edge * sizeof(double));
        memcpy(out + edge, a->blocks, mid * sizeof(double));
        memcpy(out + edge + mid, a->bottom, edge * sizeof(double));
    }
    return 2 * edge + mid;
}
