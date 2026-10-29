/* C1 spline collocation, interpolation and spatial error estimation. See LICENSE. */
#include "spatial.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdatomic.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#include "coefficients.inc"
static double pythag(double a, double b) {
    double p = fmax(fabs(a), fabs(b)), q = fmin(fabs(a), fabs(b));
    while (q != 0) {
        double z = q / p, r = z * z, t = 4 + r;
        if (t == 4)
            break;
        double s = r / t;
        p = p + (2 * s) * p;
        q *= s;
    }
    return p;
}
int bac_gauss(size_t n, int weighted, double *d, double *weights) {
    if (!n || n > 128)
        return -1;
    double *e = calloc(n, sizeof(double)), *z = weighted ? calloc(n * n, sizeof(double)) : NULL;
    if (!e || (weighted && !z)) {
        free(e);
        free(z);
        return -1;
    }
    memset(d, 0, n * sizeof(double));
    for (size_t j = 1; j < n; j++)
        e[j - 1] = (double)j / sqrt((double)(4 * j * j) - 1);
    if (weighted)
        for (size_t j = 0; j < n; j++)
            z[j * n + j] = 1;
    for (size_t l = 0; l < n; l++) {
        size_t iterations = 0;
        for (;;) {
            size_t m = l;
            while (m < n - 1) {
                double tst = fabs(d[m]) + fabs(d[m + 1]);
                if (tst + fabs(e[m]) == tst)
                    break;
                m++;
            }
            double p = d[l];
            if (m == l) {
                if (!weighted) {
                    size_t i = l;
                    while (i > 0 && p < d[i - 1]) {
                        d[i] = d[i - 1];
                        i--;
                    }
                    d[i] = p;
                }
                break;
            }
            if (iterations++ == 30) {
                free(e);
                free(z);
                return -2;
            }
            double g = (d[l + 1] - p) / (2 * e[l]), r = pythag(g, 1);
            g = d[m] - p + e[l] / (g + copysign(r, g));
            double s = 1, c = 1;
            p = 0;
            for (size_t ii = m; ii-- > l;) {
                size_t i = ii;
                double f = s * e[i], b = c * e[i];
                r = pythag(f, g);
                e[i + 1] = r;
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                if (weighted)
                    for (size_t k = 0; k < n; k++) {
                        double v = z[(i + 1) * n + k];
                        z[(i + 1) * n + k] = s * z[i * n + k] + c * v;
                        z[i * n + k] = c * z[i * n + k] - s * v;
                    }
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0;
        }
    }
    if (weighted)
        for (size_t i = 0; i + 1 < n; i++) {
            size_t k = i;
            double p = d[i];
            for (size_t j = i + 1; j < n; j++)
                if (d[j] < p) {
                    k = j;
                    p = d[j];
                }
            if (k != i) {
                d[k] = d[i];
                d[i] = p;
                for (size_t j = 0; j < n; j++) {
                    double t = z[i * n + j];
                    z[i * n + j] = z[k * n + j];
                    z[k * n + j] = t;
                }
            }
        }
    for (size_t j = 0; j < n; j++) {
        d[j] = (d[j] + 1) / 2;
        if (weights)
            weights[j] = weighted ? (2 * z[j * n] * z[j * n]) / 2 : 0;
    }
    free(e);
    free(z);
    return 0;
}
static void extend_basis(const BacGrid *g, double x, size_t il, size_t hi, size_t *order, double *v,
                         double *dp, double *dm) {
    while (*order < hi) {
        size_t j = *order;
        dp[j - 1] = g->knots[il + j - 1] - x;
        dm[j - 1] = x - g->knots[il - j];
        double prev = 0;
        for (size_t l = 1; l <= j; l++) {
            size_t op = j + 1 - l;
            double value = v[l - 1] / (dp[l - 1] + dm[op - 1]);
            v[l - 1] = value * dp[l - 1] + prev;
            prev = value * dm[op - 1];
        }
        v[j] = prev;
        (*order)++;
    }
}
static void source_basis(const BacGrid *g, double x, size_t start, size_t nd, double v[13][13]) {
    size_t k = g->kcol + 2, il = start + k, order = 1;
    double dp[13] = {0}, dm[13] = {0};
    memset(v, 0, 169 * sizeof(double));
    v[0][0] = 1;
    extend_basis(g, x, il, k + 1 - nd, &order, v[0], dp, dm);
    if (nd == 1)
        return;
    size_t derivative = nd;
    for (size_t a = 2; a <= nd; a++) {
        for (size_t j = derivative; j <= k; j++)
            v[derivative - 1][j - 1] = v[0][j - derivative];
        derivative--;
        extend_basis(g, x, il, k + 1 - derivative, &order, v[0], dp, dm);
    }
    double work[92] = {0};
    for (size_t i = 1; i <= k; i++)
        work[i * (i + 3) / 2] = 1;
    for (size_t m = 2; m <= nd; m++) {
        size_t kmd = k + 1 - m, i = il, j = k, jj = j * (j + 1) / 2, jm = jj - j;
        for (size_t a = 0; a < kmd; a++) {
            double factor = (double)kmd / (g->knots[i + kmd - 1] - g->knots[i - 1]);
            for (size_t l = 1; l <= j; l++)
                work[l + jj] = (work[l + jj] - work[l + jm]) * factor;
            i--;
            j--;
            jj = jm;
            jm -= j;
        }
        for (size_t q = 1; q <= k; q++) {
            double value = 0;
            size_t lo = MAX(q, m), idx = lo * (lo + 1) / 2;
            for (size_t a = lo; a <= k; a++) {
                value = work[q + idx] * v[m - 1][a - 1] + value;
                idx += a + 1;
            }
            v[m - 1][q - 1] = value;
        }
    }
}
int bac_grid_basis(const BacGrid *g, double x, BacBasis *b) {
    if (!isfinite(x) || x < g->mesh[0] || x > g->mesh[g->nint])
        return -1;
    size_t i;
    if (x == g->mesh[g->nint])
        i = g->nint - 1;
    else {
        size_t lo = 0, hi = g->nint + 1;
        while (lo < hi) {
            size_t m = (lo + hi) / 2;
            if (g->mesh[m] <= x)
                lo = m + 1;
            else
                hi = m;
        }
        i = lo - 1;
    }
    b->start = i * g->kcol;
    b->len = g->kcol + 2;
    double v[13][13];
    source_basis(g, x, b->start, 3, v);
    memcpy(b->b, v[0], 13 * sizeof(double));
    memcpy(b->bx, v[1], 13 * sizeof(double));
    memcpy(b->bxx, v[2], 13 * sizeof(double));
    return 0;
}
void bac_basis_sum(const BacBasis *b, const double *coeff, size_t npde, int derivative,
                   double *out) {
    const double *v = derivative == 0 ? b->b : derivative == 1 ? b->bx : b->bxx;
    if (npde == 1) {
        double value = 0;
        for (size_t j = 0; j < b->len; j++)
            value += v[j] * coeff[b->start + j];
        out[0] = value;
        return;
    }
    memset(out, 0, npde * sizeof(double));
    for (size_t j = 0; j < b->len; j++)
        for (size_t m = 0; m < npde; m++)
            out[m] += v[j] * coeff[(b->start + j) * npde + m];
}
int bac_grid_evaluate(const BacGrid *g, const double *coeff, size_t npde, double x, int derivative,
                      double *out) {
    if (!npde || derivative < 0 || derivative > (int)g->kcol + 1)
        return -1;
    BacBasis b;
    if (bac_grid_basis(g, x, &b)) {
        memset(out, 0, npde * sizeof(double));
        return -1;
    }
    if (derivative <= 2)
        bac_basis_sum(&b, coeff, npde, derivative, out);
    else {
        double v[13][13];
        source_basis(g, x, b.start, (size_t)derivative + 1, v);
        memset(out, 0, npde * sizeof(double));
        for (size_t j = 0; j < b.len; j++)
            for (size_t m = 0; m < npde; m++)
                out[m] += v[derivative][j] * coeff[(b.start + j) * npde + m];
    }
    return 0;
}
struct BacViewCache {
    atomic_flag lock;
    size_t count, capacity;
    double *x;
    BacBasis *basis;
};
int bac_grid_cached_basis(const BacGrid *g, double x, BacBasis *out) {
    BacViewCache *c = g->view_cache;
    if (!c)
        return bac_grid_basis(g, x, out);
    while (atomic_flag_test_and_set_explicit(&c->lock, memory_order_acquire)) {
    };
    for (size_t i = 0; i < c->count; i++)
        if (memcmp(c->x + i, &x, sizeof(double)) == 0) {
            *out = c->basis[i];
            atomic_flag_clear_explicit(&c->lock, memory_order_release);
            return 0;
        }
    atomic_flag_clear_explicit(&c->lock, memory_order_release);
    BacBasis b;
    if (bac_grid_basis(g, x, &b))
        return -1;
    while (atomic_flag_test_and_set_explicit(&c->lock, memory_order_acquire)) {
    };
    if (c->count == c->capacity && c->capacity < 1024) {
        size_t capacity = c->capacity ? 2 * c->capacity : 8;
        double *xs = malloc(capacity * sizeof(double));
        BacBasis *basis = malloc(capacity * sizeof(BacBasis));
        if (xs && basis) {
            if (c->count) {
                memcpy(xs, c->x, c->count * sizeof(double));
                memcpy(basis, c->basis, c->count * sizeof(BacBasis));
            }
            free(c->x);
            free(c->basis);
            c->x = xs;
            c->basis = basis;
            c->capacity = capacity;
        } else {
            free(xs);
            free(basis);
        }
    }
    if (c->count < c->capacity) {
        c->x[c->count] = x;
        c->basis[c->count++] = b;
    }
    atomic_flag_clear_explicit(&c->lock, memory_order_release);
    *out = b;
    return 0;
}
int bac_grid_flux_geometry(BacGrid *g) {
    if (g->flux_basis)
        return 0;
    size_t k = g->kcol;
    g->flux_basis = malloc(g->nint * (k + 2) * sizeof(BacBasis));
    g->flux_derivatives = calloc(k * (k + 2), sizeof(double));
    if (!g->flux_basis || !g->flux_derivatives)
        return -1;
    double nodes[13], bary[13];
    nodes[0] = 0;
    nodes[k + 1] = 1;
    bac_gauss(k, 1, nodes + 1, NULL);
    for (size_t i = 0; i < k + 2; i++) {
        double product = 1;
        for (size_t j = 0; j < k + 2; j++)
            if (i != j)
                product *= nodes[i] - nodes[j];
        bary[i] = 1 / product;
    }
    for (size_t j = 1; j <= k; j++) {
        double sum = 0;
        for (size_t i = 0; i < k + 2; i++) {
            double v = i == j ? 0 : bary[i] / (bary[j] * (nodes[j] - nodes[i]));
            g->flux_derivatives[(j - 1) * (k + 2) + i] = v;
            sum += v;
        }
        g->flux_derivatives[(j - 1) * (k + 2) + j] = -sum;
    }
    for (size_t i = 0; i < g->nint; i++)
        for (size_t j = 0; j < k + 2; j++) {
            double x = j == 0       ? g->mesh[i]
                       : j == k + 1 ? g->mesh[i + 1]
                                    : g->collocation[i * k + j];
            bac_grid_basis(g, x, g->flux_basis + i * (k + 2) + j);
        }
    return 0;
}
static int cache_geometry(BacGrid *g) {
    g->view_cache = calloc(1, sizeof(BacViewCache));
    if (!g->view_cache)
        return -1;
    atomic_flag_clear(&g->view_cache->lock);
    return 0;
}
const BacQuadrature *bac_grid_quadrature(const BacGrid *grid) {
    BacGrid *g = (BacGrid *)grid;
    BacViewCache *cache = g->view_cache;
    if (!cache)
        return NULL;
    while (atomic_flag_test_and_set_explicit(&cache->lock, memory_order_acquire)) {
    };
    if (g->quadrature) {
        const BacQuadrature *result = g->quadrature;
        atomic_flag_clear_explicit(&cache->lock, memory_order_release);
        return result;
    }
    size_t k = g->kcol;
    BacQuadrature *quadrature = malloc(g->nint * (k + 3) * sizeof(BacQuadrature));
    if (!quadrature) {
        atomic_flag_clear_explicit(&cache->lock, memory_order_release);
        return NULL;
    }
    double nodes[13], weights[13];
    bac_gauss(k + 3, 1, nodes, weights);
    for (size_t i = 0; i < g->nint; i++)
        for (size_t j = 0; j < k + 3; j++) {
            BacQuadrature *p = quadrature + i * (k + 3) + j;
            double h = g->mesh[i + 1] - g->mesh[i];
            p->x = g->mesh[i] + h * nodes[j];
            p->weight = h * weights[j];
            bac_grid_basis(g, p->x, &p->basis);
        }
    g->quadrature = quadrature;
    atomic_flag_clear_explicit(&cache->lock, memory_order_release);
    return quadrature;
}
int bac_grid_new(BacGrid *g, const double *mesh, size_t ni, size_t k) {
    memset(g, 0, sizeof *g);
    if (!ni || k < 2 || k > 10)
        return -1;
    for (size_t i = 0; i <= ni; i++)
        if (!isfinite(mesh[i]) ||
            (i && (!(mesh[i] > mesh[i - 1]) || !isfinite(mesh[i] - mesh[i - 1]))))
            return -1;
    g->nint = ni;
    g->kcol = k;
    g->n = k * ni + 2;
    g->nknots = k * (ni + 1) + 4;
    g->mesh = malloc((ni + 1) * sizeof(double));
    g->knots = malloc(g->nknots * sizeof(double));
    g->collocation = malloc(g->n * sizeof(double));
    g->collocation_basis = malloc(g->n * sizeof(BacBasis));
    if (!g->mesh || !g->knots || !g->collocation || !g->collocation_basis) {
        bac_grid_free(g);
        return -1;
    }
    memcpy(g->mesh, mesh, (ni + 1) * sizeof(double));
    size_t j = 0;
    for (size_t q = 0; q < k + 2; q++)
        g->knots[j++] = mesh[0];
    for (size_t i = 1; i < ni; i++)
        for (size_t q = 0; q < k; q++)
            g->knots[j++] = mesh[i];
    for (size_t q = 0; q < k + 2; q++)
        g->knots[j++] = mesh[ni];
    double gp[13];
    if (bac_gauss(k, 0, gp, NULL)) {
        bac_grid_free(g);
        return -1;
    }
    g->collocation[0] = mesh[0];
    for (size_t i = 0; i < ni; i++)
        for (size_t q = 0; q < k; q++)
            g->collocation[1 + i * k + q] = mesh[i] + (mesh[i + 1] - mesh[i]) * gp[q];
    g->collocation[g->n - 1] = mesh[ni];
    BacBand a;
    if (bac_band_new(&a, g->n, k + 1, k + 1)) {
        bac_grid_free(g);
        return -1;
    }
    for (size_t i = 0; i < g->n; i++) {
        BacBasis *b = g->collocation_basis + i;
        bac_grid_basis(g, g->collocation[i], b);
        for (size_t l = 0; l < b->len; l++)
            bac_band_set(&a, i, b->start + l, b->b[l]);
    }
    int rc = bac_band_factor(&a, &g->interpolation);
    bac_band_free(&a);
    if (rc || cache_geometry(g)) {
        bac_grid_free(g);
        return rc ? rc : -1;
    }
    return 0;
}
static int solve_interpolation(const BacGrid *g, size_t npde, double *out) {
    double *r = malloc(g->n * sizeof(double));
    if (!r)
        return -1;
    for (size_t m = 0; m < npde; m++) {
        for (size_t i = 0; i < g->n; i++)
            r[i] = out[i * npde + m];
        if (bac_lu_solve(&g->interpolation, r)) {
            free(r);
            return -2;
        }
        for (size_t i = 0; i < g->n; i++)
            out[i * npde + m] = r[i];
    }
    free(r);
    return 0;
}
int bac_grid_interpolate(const BacGrid *g, size_t npde, void (*f)(void *, double, double *),
                         void *ctx, double *out) {
    if (!npde || !f)
        return -1;
    for (size_t i = 0; i < g->n; i++) {
        f(ctx, g->collocation[i], out + i * npde);
        for (size_t m = 0; m < npde; m++)
            if (!isfinite(out[i * npde + m]))
                return -2;
    }
    return solve_interpolation(g, npde, out);
}
int bac_grid_project(const BacGrid *target, const BacGrid *source, const double *coeff, size_t npde,
                     double *out) {
    if (!npde || target->mesh[0] != source->mesh[0] ||
        target->mesh[target->nint] != source->mesh[source->nint])
        return -1;
    for (size_t i = 0; i < target->n; i++) {
        if (bac_grid_evaluate(source, coeff, npde, target->collocation[i], 0, out + i * npde))
            return -1;
        for (size_t m = 0; m < npde; m++)
            if (!isfinite(out[i * npde + m]))
                return -2;
    }
    return solve_interpolation(target, npde, out);
}
static const double scpoints[11][8] = {
    {0},
    {0},
    {0},
    {.5},
    {.3110177634953864, .6889822365046136},
    {.2113248654051871, .5, .7886751345948129},
    {.1526267046965671, .3747185964571342, .6252814035428658, .8473732953034329},
    {.1152723378341063, .2895425974880943, .5, .7104574025119057, .8847276621658937},
    {.09007700226825652, .2296976813063206, .405661288754607, .594338711245393, .7703023186936794,
     .9099229977317435},
    {.07229898685756272, .1863109301186906, .3341852231986051, .5, .6658147768013949,
     .8136890698813094, .9277010131424373},
    {.05929571219129399, .1539696908715823, .2792835119457421, .4241841678533667, .5758158321466333,
     .7207164880542579, .8460303091284177, .940704287808706}};
typedef struct {
    BacBasis basis;
    double weight, h[2], hd[2], g[10];
} EPoint;
typedef struct {
    double width, correction[2];
    size_t nn;
    BacBasis nodes[10];
    EPoint quadrature[13];
} EInterval;
struct BacEstimatorGeometry {
    size_t quad;
    double power;
    BacBasis *endpoints;
    EInterval *intervals;
};
static BacEstimatorGeometry *geometry_new(const BacGrid *g, BacEstimator est) {
    int sci = est == BAC_SCI;
    size_t q = g->kcol + (sci ? 3 : 2), ni = g->nint;
    if (g->kcol < 3 || (sci && ni < 2))
        return NULL;
    BacEstimatorGeometry *e = calloc(1, sizeof *e);
    if (!e)
        return NULL;
    e->quad = q;
    e->power = 1.0 / (double)(g->kcol + (sci ? 2 : 1));
    e->endpoints = malloc((ni + 1) * sizeof(BacBasis));
    e->intervals = calloc(ni, sizeof(EInterval));
    if (!e->endpoints || !e->intervals) {
        free(e->endpoints);
        free(e->intervals);
        free(e);
        return NULL;
    }
    double gp[13], ip[13], w[13];
    bac_gauss(q, 1, gp, w);
    bac_gauss(q, 0, ip, NULL);
    double (*h)[2] = NULL, (*hd)[2] = NULL, (*gg)[10] = NULL;
    if (sci) {
        h = calloc(ni * q, sizeof *h);
        hd = calloc(ni * q, sizeof *hd);
        gg = calloc(ni * q, sizeof *gg);
        if (!h || !hd || !gg) {
            free(h);
            free(hd);
            free(gg);
            free(e->endpoints);
            free(e->intervals);
            free(e);
            return NULL;
        }
        sci_coeff(g->kcol, g->mesh, ni, ip, h, hd, gg);
    }
    for (size_t i = 0; i <= ni; i++)
        bac_grid_basis(g, g->mesh[i], e->endpoints + i);
    for (size_t i = 0; i < ni; i++) {
        EInterval *v = e->intervals + i;
        double a = g->mesh[i], width = g->mesh[i + 1] - a;
        v->width = width;
        size_t k = sci ? g->kcol : g->kcol - 1, nn = k - 2;
        double nodes[10];
        for (size_t j = 0; j < nn; j++)
            nodes[j] = a + width * scpoints[k][j];
        if (sci) {
            if (i == 0) {
                double z = g->mesh[1], hh = g->mesh[2] - z;
                nodes[nn++] = z + hh * scpoints[k][0];
                nodes[nn++] = k == 3 ? g->mesh[2] : z + hh * scpoints[k][1];
            } else if (i == ni - 1) {
                double z = g->mesh[i - 1], hh = a - z;
                memmove(nodes + 2, nodes, nn * sizeof(double));
                nodes[1] = z + hh * scpoints[k][nn - 1];
                nodes[0] = k == 3 ? z : z + hh * scpoints[k][nn - 2];
                nn += 2;
            } else {
                double z = g->mesh[i - 1], hh = a - z;
                memmove(nodes + 1, nodes, nn * sizeof(double));
                nodes[0] = z + hh * scpoints[k][nn - 1];
                nn++;
                z = g->mesh[i + 1];
                hh = g->mesh[i + 2] - z;
                nodes[nn++] = z + hh * scpoints[k][0];
            }
        }
        v->nn = nn;
        for (size_t j = 0; j < nn; j++)
            bac_grid_basis(g, nodes[j], v->nodes + j);
        for (size_t j = 0; j < q; j++) {
            EPoint *p = v->quadrature + j;
            bac_grid_basis(g, a + width * gp[j], &p->basis);
            p->weight = width * w[j];
            if (sci) {
                memcpy(p->h, h[i * q + j], sizeof p->h);
                memcpy(p->hd, hd[i * q + j], sizeof p->hd);
                memcpy(p->g, gg[i * q + j], sizeof p->g);
            } else
                low_coeff(ip[j], g->kcol - 1, p->h, p->hd, p->g);
        }
        v->correction[0] = v->correction[1] = 1;
        if (sci) {
            if (i == 0) {
                double r = (g->mesh[2] - g->mesh[1]) / width;
                v->correction[0] = v->correction[1] = r;
            } else if (i == ni - 1) {
                double r = (g->mesh[i] - g->mesh[i - 1]) / width;
                v->correction[0] = v->correction[1] = r;
            } else {
                v->correction[0] = (g->mesh[i] - g->mesh[i - 1]) / width;
                v->correction[1] = (g->mesh[i + 2] - g->mesh[i + 1]) / width;
            }
        }
    }
    free(h);
    free(hd);
    free(gg);
    return e;
}
void bac_grid_free(BacGrid *g) {
    free(g->flux_basis);
    free(g->flux_derivatives);
    free(g->quadrature);
    if (g->view_cache) {
        free(g->view_cache->x);
        free(g->view_cache->basis);
        free(g->view_cache);
    }
    free(g->boundary_lift);
    free(g->mesh);
    free(g->knots);
    free(g->collocation);
    free(g->collocation_basis);
    bac_lu_free(&g->interpolation);
    for (size_t i = 0; i < 2; i++)
        if (g->geometry[i]) {
            free(g->geometry[i]->endpoints);
            free(g->geometry[i]->intervals);
            free(g->geometry[i]);
        }
    memset(g, 0, sizeof *g);
}
void bac_spatial_error_free(BacSpatialError *e) {
    free(e->indicators);
    free(e->component_ratios);
    memset(e, 0, sizeof *e);
}
int bac_grid_estimate(BacGrid *g, const double *c, size_t npde, const double *atol,
                      const double *rtol, BacEstimator est, BacSpatialError *out) {
    memset(out, 0, sizeof *out);
    if (!npde || g->kcol < 3 || (est == BAC_SCI && g->nint < 2))
        return -1;
    for (size_t m = 0; m < npde; m++)
        if (!isfinite(atol[m]) || !isfinite(rtol[m]) || atol[m] < 0 || rtol[m] < 0 ||
            atol[m] + rtol[m] == 0)
            return -1;
    for (size_t i = 0; i < g->n * npde; i++)
        if (!isfinite(c[i]))
            return -2;
    if (!g->geometry[est])
        g->geometry[est] = geometry_new(g, est);
    BacEstimatorGeometry *e = g->geometry[est];
    if (!e)
        return -1;
    size_t ni = g->nint;
    double *sq = calloc(ni * npde, sizeof(double)), *work = malloc(15 * npde * sizeof(double));
    out->indicators = calloc(ni, sizeof(double));
    out->component_ratios = calloc(npde, sizeof(double));
    if (!sq || !work || !out->indicators || !out->component_ratios) {
        free(sq);
        free(work);
        bac_spatial_error_free(out);
        return -1;
    }
    double *nv = work, *ul = nv + 10 * npde, *ur = ul + npde, *dl = ur + npde, *dr = dl + npde,
           *exact = dr + npde;
    for (size_t i = 0; i < ni; i++) {
        EInterval *v = e->intervals + i;
        for (size_t j = 0; j < v->nn; j++)
            bac_basis_sum(v->nodes + j, c, npde, 0, nv + j * npde);
        bac_basis_sum(e->endpoints + i, c, npde, 0, ul);
        bac_basis_sum(e->endpoints + i + 1, c, npde, 0, ur);
        bac_basis_sum(e->endpoints + i, c, npde, 1, dl);
        bac_basis_sum(e->endpoints + i + 1, c, npde, 1, dr);
        for (size_t j = 0; j < e->quad; j++) {
            EPoint *p = v->quadrature + j;
            bac_basis_sum(&p->basis, c, npde, 0, exact);
            for (size_t m = 0; m < npde; m++) {
                double interp;
                if (est == BAC_LOI)
                    interp = p->h[0] * ul[m] + p->h[1] * ur[m] + p->hd[0] * dl[m] * v->width +
                             p->hd[1] * dr[m] * v->width;
                else if (g->kcol == 3)
                    interp = p->g[0] * nv[m] + p->g[1] * nv[npde + m] + p->g[2] * nv[2 * npde + m] +
                             p->h[0] * ul[m] + p->h[1] * ur[m] + p->hd[0] * dl[m] +
                             p->hd[1] * dr[m];
                else
                    interp =
                        p->h[0] * ul[m] + p->h[1] * ur[m] + p->hd[0] * dl[m] + p->hd[1] * dr[m];
                if (est == BAC_LOI || g->kcol != 3)
                    for (size_t l = 0; l < v->nn; l++)
                        interp += p->g[l] * nv[l * npde + m];
                double tol = atol[m] + rtol[m] * fabs(exact[m]);
                if (tol == 0) {
                    free(sq);
                    free(work);
                    bac_spatial_error_free(out);
                    return -2;
                }
                double scaled = (exact[m] - interp) / tol;
                sq[i * npde + m] += scaled * scaled * p->weight;
            }
        }
        if (est == BAC_SCI)
            for (size_t m = 0; m < npde; m++)
                sq[i * npde + m] = sq[i * npde + m] * v->correction[0] * v->correction[1];
    }
    for (size_t m = 0; m < npde; m++) {
        double sum = 0;
        for (size_t i = 0; i < ni; i++)
            sum += sq[i * npde + m];
        out->component_ratios[m] = sqrt(sum);
        out->ratio = fmax(out->ratio, out->component_ratios[m]);
    }
    double sum = 0, mx = 0;
    for (size_t i = 0; i < ni; i++) {
        double v = 0;
        for (size_t m = 0; m < npde; m++)
            v += sq[i * npde + m];
        out->indicators[i] = pow(sqrt(v), e->power);
        sum += out->indicators[i];
        mx = fmax(mx, out->indicators[i]);
    }
    out->distribution_ratio = sum == 0 ? 1 : mx / (sum / (double)ni);
    out->intervals = ni;
    out->npde = npde;
    out->estimator = est;
    free(sq);
    free(work);
    if (!isfinite(out->ratio) || !isfinite(sum)) {
        bac_spatial_error_free(out);
        return -2;
    }
    return 0;
}
int bac_grid_refined_mesh(const BacGrid *g, const BacSpatialError *e, size_t max, size_t prescribed,
                          int capacity, double **out, size_t *count) {
    size_t old = g->nint, min = e->estimator == BAC_SCI ? 2 : 1;
    if (e->intervals != old || !isfinite(e->ratio) || e->ratio < 0 || max < min)
        return -1;
    double factor = fmax(
        .5,
        fmin(2, pow(e->ratio / .2, 1.0 / (double)(g->kcol + (e->estimator == BAC_LOI ? 1 : 2)))));
    size_t n = prescribed ? prescribed : MAX((size_t)((double)old * factor), min);
    if (n < min)
        return -1;
    if (n == old && !prescribed)
        n++;
    if (n > max) {
        if (capacity && !prescribed)
            n = max;
        else
            return -3;
    }
    double *mesh = malloc((n + 1) * sizeof(double));
    if (!mesh)
        return -1;
    double total = 0;
    for (size_t i = 0; i < old; i++) {
        if (!isfinite(e->indicators[i]) || e->indicators[i] < 0) {
            free(mesh);
            return -1;
        }
        total += e->indicators[i];
    }
    mesh[0] = g->mesh[0];
    mesh[n] = g->mesh[old];
    if (total == 0) {
        for (size_t j = 1; j < n; j++)
            mesh[j] = mesh[0] + (mesh[n] - mesh[0]) * ((double)j / (double)n);
    } else {
        double *cum = malloc(old * sizeof(double));
        if (!cum) {
            free(mesh);
            return -1;
        }
        double sum = e->indicators[0];
        cum[0] = sum;
        for (size_t i = 1; i < old; i++) {
            sum = e->indicators[i] + sum;
            cum[i] = sum;
        }
        double mean = cum[old - 1] / (double)n, target = mean;
        size_t i = 0;
        for (size_t j = 1; j < n; j++) {
            while (i < old - 1 && target > cum[i])
                i++;
            double density = e->indicators[i];
            if (density == 0) {
                free(cum);
                free(mesh);
                return -2;
            }
            double accumulated = i == 0 ? 0 : cum[i - 1];
            mesh[j] = g->mesh[i] + (g->mesh[i + 1] - g->mesh[i]) * (target - accumulated) / density;
            target += mean;
        }
        free(cum);
    }
    for (size_t i = 1; i <= n; i++)
        if (!(mesh[i] > mesh[i - 1])) {
            free(mesh);
            return -2;
        }
    *out = mesh;
    *count = n;
    return 0;
}
int bac_grid_projection_error(const BacGrid *target, const BacGrid *source, const double *original,
                              const double *projected, size_t npde, const double *atol,
                              const double *rtol, double *ratio, double *maximum) {
    double nodes[13], weights[13];
    size_t q = MAX(source->kcol, target->kcol) + 2;
    bac_gauss(q, 1, nodes, weights);
    double *v = calloc(3 * npde, sizeof(double));
    if (!v)
        return -1;
    double *a = v, *b = v + npde, *sq = v + 2 * npde;
    *maximum = 0;
    double x = source->mesh[0];
    size_t si = 0, ti = 0;
    while (si < source->nint && ti < target->nint) {
        double end = fmin(source->mesh[si + 1], target->mesh[ti + 1]);
        for (size_t j = 0; j < q; j++) {
            double z = x + (end - x) * nodes[j];
            bac_grid_evaluate(source, original, npde, z, 0, a);
            bac_grid_evaluate(target, projected, npde, z, 0, b);
            for (size_t m = 0; m < npde; m++) {
                double tol = atol[m] + rtol[m] * fmax(fabs(a[m]), fabs(b[m]));
                double e = fabs(a[m] - b[m]) / tol;
                sq[m] += e * e * (end - x) * weights[j];
                *maximum = fmax(*maximum, e);
            }
        }
        x = end;
        if (source->mesh[si + 1] == end)
            si++;
        if (target->mesh[ti + 1] == end)
            ti++;
    }
    *ratio = 0;
    for (size_t m = 0; m < npde; m++)
        *ratio = fmax(*ratio, sqrt(sq[m]));
    free(v);
    return 0;
}
void bac_projection_error_free(BacProjectionError *e) {
    free(e->interval_maximum);
    free(e->interval_squared);
    memset(e, 0, sizeof *e);
}
int bac_grid_projection_measure(const BacGrid *target, const BacGrid *source,
                                const double *original, const double *projected, size_t npde,
                                const double *atol, const double *rtol, BacProjectionError *error) {
    memset(error, 0, sizeof *error);
    error->intervals = target->nint;
    error->interval_maximum = calloc(target->nint, sizeof(double));
    error->interval_squared = calloc(target->nint, sizeof(double));
    double nodes[13], weights[13];
    size_t q = MAX(source->kcol, target->kcol) + 2;
    bac_gauss(q, 1, nodes, weights);
    double *v = calloc(3 * npde, sizeof(double));
    if (!v || !error->interval_maximum || !error->interval_squared) {
        free(v);
        bac_projection_error_free(error);
        return -1;
    }
    double *a = v, *b = v + npde, *sq = v + 2 * npde;
    double x = source->mesh[0];
    size_t si = 0, ti = 0;
    while (si < source->nint && ti < target->nint) {
        double end = fmin(source->mesh[si + 1], target->mesh[ti + 1]), width = end - x;
        for (size_t j = 0; j < q; j++) {
            double z = x + width * nodes[j];
            bac_grid_evaluate(source, original, npde, z, 0, a);
            bac_grid_evaluate(target, projected, npde, z, 0, b);
            for (size_t m = 0; m < npde; m++) {
                double tol = atol[m] + rtol[m] * fmax(fabs(a[m]), fabs(b[m]));
                if (tol == 0) {
                    free(v);
                    bac_projection_error_free(error);
                    return -2;
                }
                double e = fabs((a[m] - b[m]) / tol), contribution = e * e * width * weights[j];
                sq[m] += contribution;
                error->interval_squared[ti] += contribution;
                error->interval_maximum[ti] = fmax(error->interval_maximum[ti], e);
            }
        }
        x = end;
        if (source->mesh[si + 1] == end)
            si++;
        if (target->mesh[ti + 1] == end)
            ti++;
    }
    double max = 0;
    for (size_t m = 0; m < npde; m++)
        max = fmax(max, sq[m]);
    error->ratio = sqrt(max);
    for (size_t i = 0; i < target->nint; i++)
        error->maximum = fmax(error->maximum, error->interval_maximum[i]);
    free(v);
    return isfinite(error->ratio) && isfinite(error->maximum) ? 0 : -2;
}
int bac_grid_refine_projection(const BacGrid *g, const BacGrid *source,
                               const BacProjectionError *error, size_t max, double **out,
                               size_t *intervals) {
    size_t worst = 0;
    for (size_t i = 1; i < g->nint; i++)
        if (error->interval_squared[i] >= error->interval_squared[worst])
            worst = i;
    double *mesh = malloc((2 * g->nint + source->nint + 2) * sizeof(double));
    if (!mesh)
        return -1;
    size_t n = 0;
    for (size_t i = 0; i < g->nint; i++) {
        mesh[n++] = g->mesh[i];
        if (error->interval_maximum[i] > 1 || (error->ratio > .2 && i == worst)) {
            size_t prev = n;
            for (size_t j = 0; j <= source->nint; j++)
                if (source->mesh[j] > g->mesh[i] && source->mesh[j] < g->mesh[i + 1])
                    mesh[n++] = source->mesh[j];
            if (n == prev) {
                double midpoint = g->mesh[i] + .5 * (g->mesh[i + 1] - g->mesh[i]);
                if (!(midpoint > g->mesh[i] && midpoint < g->mesh[i + 1])) {
                    free(mesh);
                    return -2;
                }
                mesh[n++] = midpoint;
            }
        }
    }
    mesh[n++] = g->mesh[g->nint];
    if (n - 1 > max) {
        free(mesh);
        return -3;
    }
    *out = mesh;
    *intervals = n - 1;
    return 0;
}
int bac_source_projection_new(BacSourceProjection *p, const BacGrid *target, const BacGrid *source,
                              size_t npde) {
    memset(p, 0, sizeof *p);
    if (!npde || target->mesh[0] != source->mesh[0] ||
        target->mesh[target->nint] != source->mesh[source->nint])
        return -1;
    size_t n = target->n * npde, w = (target->kcol + 2) * npde - 1;
    BacBand a;
    if (bac_band_new(&a, n, w, w))
        return -1;
    a.npde = npde;
    a.kcol = target->kcol;
    a.nint = target->nint;
    p->basis = malloc(target->n * sizeof(BacBasis));
    if (!p->basis) {
        bac_band_free(&a);
        return -1;
    }
    for (size_t row = 0; row < target->n; row++) {
        const BacBasis *b = target->collocation_basis + row;
        for (size_t j = 0; j < b->len; j++)
            for (size_t m = 0; m < npde; m++)
                bac_band_set(&a, row * npde + m, (b->start + j) * npde + m, b->b[j]);
        bac_grid_basis(source, target->collocation[row], p->basis + row);
    }
    int rc = bac_band_factor(&a, &p->factor);
    bac_band_free(&a);
    if (rc) {
        bac_source_projection_free(p);
        return rc;
    }
    p->source = source;
    p->target = target;
    p->npde = npde;
    return 0;
}
int bac_source_projection_apply(const BacSourceProjection *p, const double *in, double *out) {
    for (size_t i = 0; i < p->source->n * p->npde; i++)
        if (!isfinite(in[i]))
            return -1;
    for (size_t row = 0; row < p->target->n; row++)
        bac_basis_sum(p->basis + row, in, p->npde, 0, out + row * p->npde);
    return bac_lu_solve(&p->factor, out);
}
void bac_source_projection_free(BacSourceProjection *p) {
    free(p->basis);
    bac_lu_free(&p->factor);
    memset(p, 0, sizeof *p);
}
