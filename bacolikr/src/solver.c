/* Strong and conservative PDE collocation with adaptive spatial control. */
#include "bacolikr.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
static int fail(BacSolver *s, const char *msg) {
    snprintf(s->error, sizeof s->error, "%s", msg);
    return -1;
}
static int compat(const BacSolver *s) {
    return s->options.time.method == BI_DASKR;
}
BacOptions bac_options(void) {
    BacOptions o = {0};
    o.intervals = 10;
    o.kcol = 4;
    o.absolute = o.relative = 1e-6;
    o.estimator = BAC_LOI;
    o.max_intervals = 500;
    o.max_remeshes = 25;
    o.spatial_adaptivity = 1;
    o.time = bi_options(1e-8, 1e-6);
    return o;
}
static void state(const BacBasis *b, const double *y, size_t p, int side, int original, double *u,
                  double *ux, double *uxx) {
    if (original && side >= 0) {
        if (side == 0)
            for (size_t c = 0; c < p; c++) {
                u[c] = b->b[0] * y[c];
                ux[c] = b->bx[0] * y[c] + b->bx[1] * y[p + c];
                uxx[c] = b->bxx[0] * y[c] + b->bxx[1] * y[p + c] + b->bxx[2] * y[2 * p + c];
            }
        else {
            size_t last = b->len - 1, end = (b->start + last) * p;
            for (size_t c = 0; c < p; c++) {
                u[c] = b->b[last] * y[end + c];
                ux[c] = b->bx[last] * y[end + c] + b->bx[last - 1] * y[end - p + c];
                uxx[c] = b->bxx[last] * y[end + c] + b->bxx[last - 1] * y[end - p + c] +
                         b->bxx[last - 2] * y[end - 2 * p + c];
            }
        }
        return;
    }
    if (p == 1) {
        double value = 0, first = 0, second = 0;
        for (size_t j = 0; j < b->len; j++) {
            double v = y[b->start + j];
            value += b->b[j] * v;
            first += b->bx[j] * v;
            second += b->bxx[j] * v;
        }
        u[0] = value;
        ux[0] = first;
        uxx[0] = second;
        return;
    }
    memset(u, 0, p * sizeof(double));
    memset(ux, 0, p * sizeof(double));
    memset(uxx, 0, p * sizeof(double));
    for (size_t j = 0; j < b->len; j++)
        for (size_t m = 0; m < p; m++) {
            double v = y[(b->start + j) * p + m];
            u[m] += b->b[j] * v;
            ux[m] += b->bx[j] * v;
            uxx[m] += b->bxx[j] * v;
        }
}
/* Flux callbacks consume values and first derivatives only. */
static void flux_state(const BacBasis *b, const double *y, size_t p, double *u, double *ux) {
    if (p == 1) {
        double value = 0, first = 0;
        for (size_t j = 0; j < b->len; j++) {
            double v = y[b->start + j];
            value += b->b[j] * v;
            first += b->bx[j] * v;
        }
        u[0] = value;
        ux[0] = first;
        return;
    }
    memset(u, 0, p * sizeof(double));
    memset(ux, 0, p * sizeof(double));
    for (size_t j = 0; j < b->len; j++)
        for (size_t m = 0; m < p; m++) {
            double v = y[(b->start + j) * p + m];
            u[m] += b->b[j] * v;
            ux[m] += b->bx[j] * v;
        }
}
static int boundary_derivatives(const BacSolver *s, int side, double t, const double *u,
                                const double *ux, double *du, double *dx, double *dt) {
    size_t p = s->problem.components;
    int original = compat(s);
    int (*jac)(void *, int, double, const double *, const double *, double *, double *, double *) =
        original && s->problem.source_boundary_jacobian ? s->problem.source_boundary_jacobian
                                                        : s->problem.boundary_jacobian;
    if (jac && jac(s->problem.ctx, side, t, u, ux, du, dx, dt))
        return 0;
    double stack[48], *work = p <= 8 ? stack : malloc(6 * p * sizeof(double));
    if (!work)
        return -1;
    double *uu = work, *xx = uu + p, *base = xx + p, *value = base + p, *minus = value + p;
    memcpy(uu, u, p * sizeof(double));
    memcpy(xx, ux, p * sizeof(double));
    s->problem.boundary(s->problem.ctx, side, t, u, ux, base);
    double uround = 2 * DBL_EPSILON;
    for (size_t group = 0; group < 2; group++)
        for (size_t d = 0; d < p; d++) {
            double *v = group ? xx : uu, old = v[d],
                   delta = sqrt(original ? uround : DBL_EPSILON) *
                           fmax(fabs(old), original ? uround : 1),
                   actual = (old + delta) - old;
            v[d] = old + actual;
            s->problem.boundary(s->problem.ctx, side, t, uu, xx, value);
            double *mat = group ? dx : du;
            if (original) {
                double d2 = (old - delta) - old;
                v[d] = old + d2;
                s->problem.boundary(s->problem.ctx, side, t, uu, xx, minus);
                double inv = 1 / (actual - d2);
                for (size_t c = 0; c < p; c++)
                    mat[c * p + d] = (value[c] - minus[c]) * inv;
            } else
                for (size_t c = 0; c < p; c++)
                    mat[c * p + d] = (value[c] - base[c]) / actual;
            v[d] = old;
        }
    double delta = sqrt(original ? uround : DBL_EPSILON) * fmax(fabs(t), original ? uround : 1);
    if (original) {
        double d1 = (t + delta) - t, d2 = (t - delta) - t;
        s->problem.boundary(s->problem.ctx, side, t + d1, uu, xx, value);
        s->problem.boundary(s->problem.ctx, side, t + d2, uu, xx, minus);
        for (size_t c = 0; c < p; c++)
            dt[c] = (value[c] - minus[c]) * (1 / (d1 - d2));
    } else {
        s->problem.boundary(s->problem.ctx, side, t + delta, u, ux, value);
        for (size_t c = 0; c < p; c++)
            dt[c] = (value[c] - base[c]) / delta;
    }
    if (work != stack)
        free(work);
    return 0;
}
static int rhs_derivatives(const BacSolver *s, double t, double x, const double *u,
                           const double *ux, const double *uxx, double *du, double *dx,
                           double *dxx) {
    size_t p = s->problem.components;
    int (*jac)(void *, double, double, const double *, const double *, const double *, double *,
               double *, double *) = compat(s) && s->problem.source_rhs_jacobian
                                         ? s->problem.source_rhs_jacobian
                                         : s->problem.rhs_jacobian;
    if (jac && jac(s->problem.ctx, t, x, u, ux, uxx, du, dx, dxx))
        return 1;
    if (compat(s))
        return 0;
    double stack[40], *w = p <= 8 ? stack : malloc(5 * p * sizeof(double));
    if (!w)
        return -1;
    double *uu = w, *xx = w + p, *xxx = w + 2 * p, *base = w + 3 * p, *value = w + 4 * p;
    memcpy(uu, u, p * sizeof(double));
    memcpy(xx, ux, p * sizeof(double));
    memcpy(xxx, uxx, p * sizeof(double));
    s->problem.rhs(s->problem.ctx, t, x, u, ux, uxx, base);
    for (size_t group = 0; group < 3; group++)
        for (size_t d = 0; d < p; d++) {
            double *v = group == 0   ? uu
                        : group == 1 ? xx
                                     : xxx,
                   *mat = group == 0   ? du
                          : group == 1 ? dx
                                       : dxx,
                   old = v[d], delta = sqrt(DBL_EPSILON) * fmax(fabs(old), 1);
            v[d] += delta;
            double actual = v[d] - old;
            s->problem.rhs(s->problem.ctx, t, x, uu, xx, xxx, value);
            for (size_t c = 0; c < p; c++)
                mat[c * p + d] = (value[c] - base[c]) / actual;
            v[d] = old;
        }
    if (w != stack)
        free(w);
    return 1;
}
static void source_boundary_matrix(BacBand *a, size_t row, const BacBasis *b, size_t p,
                                   const double *du, const double *dx, double scale) {
    for (size_t c = 0; c < p; c++)
        for (size_t d = 0; d < p; d++) {
            double first, second;
            if (row == 0) {
                second = b->bx[1] * dx[c * p + d];
                first = du[c * p + d] - second;
            } else {
                first = b->bx[b->len - 2] * dx[c * p + d];
                second = du[c * p + d] - first;
            }
            size_t start = row == 0 ? 0 : row - 1;
            bac_band_set(a, row * p + c, start * p + d, scale * first);
            bac_band_set(a, row * p + c, (start + 1) * p + d, scale * second);
        }
}
/* Degree kcol+1 flux interpolation through endpoints and Gaussian nodes. */
static int uses_conservation_form(const BacSolver *s, double t, const double *y) {
    if (compat(s) || !s->problem.conservation_form)
        return 0;
    size_t p = s->problem.components;
    double stack[32], *w = p <= 8 ? stack : malloc(4 * p * sizeof(double));
    if (!w)
        return -1;
    flux_state(s->grid.collocation_basis, y, p, w, w + p);
    int use = s->problem.conservation_form(s->problem.ctx, t, s->grid.mesh[0], w, w + p, w + 2 * p,
                                           w + 3 * p);
    if (w != stack)
        free(w);
    return use != 0;
}
static int conservative_rhs(const BacSolver *s, double t, const double *y, double *out) {
    if (compat(s) || !s->problem.conservation_form)
        return 0;
    const BacGrid *g = &s->grid;
    size_t p = s->problem.components, k = g->kcol;
    double stack[136], *w = p <= 8 ? stack : malloc((5 * p + (k + 2) * p) * sizeof(double));
    if (!w)
        return -1;
    double *u = w, *ux = u + p, *unused = ux + p, *src = unused + p, *flux = src + p,
           *fluxes = flux + p;
    flux_state(g->collocation_basis, y, p, u, ux);
    if (!s->problem.conservation_form(s->problem.ctx, t, g->mesh[0], u, ux, flux, src)) {
        if (w != stack)
            free(w);
        return 0;
    }
    if (bac_grid_flux_geometry((BacGrid *)g)) {
        if (w != stack)
            free(w);
        return -1;
    }
    memset(out, 0, g->n * p * sizeof(double));
    for (size_t i = 0; i < g->nint; i++) {
        double h = g->mesh[i + 1] - g->mesh[i];
        for (size_t j = 0; j < k + 2; j++) {
            double x = j == 0       ? g->mesh[i]
                       : j == k + 1 ? g->mesh[i + 1]
                                    : g->collocation[i * k + j];
            const BacBasis *basis = g->flux_basis + i * (k + 2) + j;
            flux_state(basis, y, p, u, ux);
            if (!s->problem.conservation_form(s->problem.ctx, t, x, u, ux, fluxes + j * p, src)) {
                if (w != stack)
                    free(w);
                return -1;
            }
            for (size_t m = 0; m < p; m++)
                if (!isfinite(fluxes[j * p + m]) || !isfinite(src[m])) {
                    if (w != stack)
                        free(w);
                    return -1;
                }
            if (j > 0 && j <= k)
                memcpy(out + (i * k + j) * p, src, p * sizeof(double));
        }
        for (size_t j = 0; j < k; j++)
            for (size_t m = 0; m < p; m++) {
                double value = 0, center = fluxes[(j + 1) * p + m];
                for (size_t l = 0; l < k + 2; l++)
                    if (l != j + 1)
                        value +=
                            g->flux_derivatives[j * (k + 2) + l] * (fluxes[l * p + m] - center);
                out[(i * k + j + 1) * p + m] += value / h;
            }
    }
    if (w != stack)
        free(w);
    return 1;
}
static void residual(void *ctx, double t, const double *y, const double *yp, double cj,
                     double *out) {
    BacSolver *s = ctx;
    const BacGrid *g = &s->grid;
    size_t p = s->problem.components, n = g->n * p;
    double stack[32], *w = p <= 8 ? stack : malloc(4 * p * sizeof(double));
    if (!w) {
        for (size_t i = 0; i < n; i++)
            out[i] = NAN;
        return;
    }
    double *u = w, *ux = w + p, *uxx = w + 2 * p, *f = w + 3 * p;
    /* Flux values can occupy the residual buffer until each row is consumed. */
    int use = conservative_rhs(s, t, y, out);
    for (size_t row = 0; row < g->n; row++) {
        const BacBasis *b = g->collocation_basis + row;
        int side = row == 0 ? 0 : row == g->n - 1 ? 1 : -1;
        if (side >= 0 || use == 0)
            state(b, y, p, side, compat(s), u, ux, uxx);
        if (side >= 0) {
            s->problem.boundary(s->problem.ctx, side, t, u, ux, f);
            for (size_t c = 0; c < p; c++)
                out[row * p + c] = f[c] * (compat(s) && isfinite(cj) ? cj : 1);
        } else {
            if (use > 0)
                memcpy(f, out + row * p, p * sizeof(double));
            else if (use < 0)
                for (size_t c = 0; c < p; c++)
                    f[c] = NAN;
            else
                s->problem.rhs(s->problem.ctx, t, g->collocation[row], u, ux, uxx, f);
            for (size_t c = 0; c < p; c++)
                out[row * p + c] = -f[c];
            for (size_t j = 0; j < b->len; j++)
                if (b->b[j] != 0)
                    for (size_t c = 0; c < p; c++)
                        out[row * p + c] += b->b[j] * yp[(b->start + j) * p + c];
        }
    }
    if (w != stack)
        free(w);
}
static int conservative_jacobian(BacSolver *s, double t, const double *y, double cj,
                                 BacBand *matrix) {
    if (!s->problem.conservation_jacobian)
        return 0;
    BacGrid *g = &s->grid;
    size_t p = s->problem.components, k = g->kcol, pp = p * p;
    if (bac_grid_flux_geometry(g))
        return -1;
    size_t count = 4 * p + 4 * pp * (k + 2);
    double stack[256], *w = p <= 2 ? stack : malloc(count * sizeof(double));
    if (!w)
        return -1;
    memset(w, 0, count * sizeof(double));
    double *u = w, *ux = u + p, *uxx = ux + p, *dt = uxx + p, *jac = dt + p;
    for (size_t i = 0; i < g->nint; i++) {
        for (size_t j = 0; j < k + 2; j++) {
            const BacBasis *b = g->flux_basis + i * (k + 2) + j;
            flux_state(b, y, p, u, ux);
            double x = j == 0       ? g->mesh[i]
                       : j == k + 1 ? g->mesh[i + 1]
                                    : g->collocation[i * k + j],
                   *a = jac + j * 4 * pp;
            if (!s->problem.conservation_jacobian(s->problem.ctx, t, x, u, ux, a, a + pp,
                                                  a + 2 * pp, a + 3 * pp)) {
                if (w != stack)
                    free(w);
                return 0;
            }
            for (size_t l = 0; l < 4 * pp; l++)
                if (!isfinite(a[l])) {
                    bac_band_set(matrix, 0, 0, NAN);
                    if (w != stack)
                        free(w);
                    return 1;
                }
        }
        double h = g->mesh[i + 1] - g->mesh[i];
        for (size_t j = 0; j < k; j++) {
            size_t row = i * k + j + 1;
            for (size_t q = 0; q < k + 2; q++) {
                const BacBasis *b = g->flux_basis + i * (k + 2) + q;
                double weight = g->flux_derivatives[j * (k + 2) + q], *du = jac + q * 4 * pp,
                       *dx = du + pp;
                for (size_t l = 0; l < b->len; l++) {
                    double bb = b->b[l], bx = b->bx[l];
                    if (bb == 0 && bx == 0)
                        continue;
                    for (size_t c = 0; c < p; c++)
                        for (size_t d = 0; d < p; d++)
                            bac_band_add(matrix, row * p + c, (b->start + l) * p + d,
                                         -weight * (du[c * p + d] * bb + dx[c * p + d] * bx) / h);
                }
            }
            const BacBasis *b = g->flux_basis + i * (k + 2) + j + 1;
            double *du = jac + (j + 1) * 4 * pp + 2 * pp, *dx = du + pp;
            for (size_t l = 0; l < b->len; l++)
                for (size_t c = 0; c < p; c++)
                    for (size_t d = 0; d < p; d++)
                        bac_band_add(matrix, row * p + c, (b->start + l) * p + d,
                                     -du[c * p + d] * b->b[l] - dx[c * p + d] * b->bx[l] +
                                         (c == d ? cj * b->b[l] : 0));
        }
    }
    for (size_t side = 0; side < 2; side++) {
        size_t row = side ? g->n - 1 : 0;
        const BacBasis *b = g->collocation_basis + row;
        state(b, y, p, -1, 0, u, ux, uxx);
        double *du = jac, *dx = du + pp;
        boundary_derivatives(s, (int)side, t, u, ux, du, dx, dt);
        for (size_t l = 0; l < b->len; l++)
            for (size_t c = 0; c < p; c++)
                for (size_t d = 0; d < p; d++)
                    bac_band_set(matrix, row * p + c, (b->start + l) * p + d,
                                 du[c * p + d] * b->b[l] + dx[c * p + d] * b->bx[l]);
    }
    if (w != stack)
        free(w);
    return 1;
}
static int jacobian(void *ctx, double t, const double *y, const double *yp, double cj, BacBand *a,
                    double *dense) {
    (void)yp;
    (void)dense;
    BacSolver *s = ctx;
    const BacGrid *g = &s->grid;
    size_t p = s->problem.components, pp = p * p;
    int conservative = uses_conservation_form(s, t, y);
    if (conservative)
        return conservative < 0 ? -1 : conservative_jacobian(s, t, y, cj, a);
    size_t count = 4 * p + 3 * pp;
    double stack[224], *w = p <= 8 ? stack : malloc(count * sizeof(double));
    if (!w)
        return -1;
    memset(w, 0, count * sizeof(double));
    double *u = w, *ux = u + p, *uxx = ux + p, *dt = uxx + p, *du = dt + p, *dx = du + pp,
           *dxx = dx + pp;
    for (size_t row = 0; row < g->n; row++) {
        const BacBasis *b = g->collocation_basis + row;
        int side = row == 0 ? 0 : row == g->n - 1 ? 1 : -1;
        state(b, y, p, side, compat(s), u, ux, uxx);
        if (side >= 0) {
            if (boundary_derivatives(s, side, t, u, ux, du, dx, dt)) {
                if (w != stack)
                    free(w);
                return -1;
            }
            memset(dxx, 0, pp * sizeof(double));
            if (compat(s)) {
                source_boundary_matrix(a, row, b, p, du, dx, cj);
                continue;
            }
        } else {
            int rc = rhs_derivatives(s, t, g->collocation[row], u, ux, uxx, du, dx, dxx);
            if (rc <= 0) {
                if (w != stack)
                    free(w);
                return rc;
            }
        }
        for (size_t j = 0; j < b->len; j++) {
            size_t col = b->start + j;
            double bb = b->b[j], bx = b->bx[j], bxx = b->bxx[j];
            if (bb == 0 && bx == 0 && bxx == 0)
                continue;
            for (size_t c = 0; c < p; c++)
                for (size_t d = 0; d < p; d++) {
                    double v;
                    if (compat(s)) {
                        v = -du[c * p + d] * bb - dx[c * p + d] * bx - dxx[c * p + d] * bxx;
                        v += cj * (c == d ? bb : 0);
                    } else {
                        v = du[c * p + d] * bb + dx[c * p + d] * bx + dxx[c * p + d] * bxx;
                        if (side < 0)
                            v = -v + (c == d ? cj * bb : 0);
                    }
                    bac_band_set(a, row * p + c, col * p + d, v);
                }
        }
    }
    a->npde = compat(s) ? p : 0;
    a->kcol = g->kcol;
    a->nint = g->nint;
    if (w != stack)
        free(w);
    return 1;
}
static DaeSystem system_for(BacSolver *s) {
    DaeSystem d = {0};
    d.n = s->grid.n * s->problem.components;
    d.ctx = s;
    d.residual = residual;
    d.jacobian = jacobian;
    d.banded = 1;
    d.lower = d.upper = (s->grid.kcol + 1) * s->problem.components;
    d.npde = compat(s) ? s->problem.components : 0;
    d.kcol = s->grid.kcol;
    d.nint = s->grid.nint;
    return d;
}
static double roundoff(size_t n, double scale) {
    double e = (double)n * DBL_EPSILON, b = e / (1 - e) * scale;
    return isfinite(scale) && scale >= 0 && e < 1 && isfinite(b) ? b : 0;
}
static int enforce_boundary(BacSolver *s, double t, double *y, const double *initial) {
    BacGrid *g = &s->grid;
    size_t p = s->problem.components, n = g->n * p, pp = p * p, wid = (g->kcol + 1) * p;
    double *work = calloc(6 * p + 2 * pp + 3 * n, sizeof(double));
    if (!work)
        return -1;
    double *u = work, *ux = u + p, *uxx = ux + p, *dt = uxx + p, *av = dt + p, *ad = av + p,
           *du = ad + p, *dx = du + pp, *target = dx + pp, *r = target + n, *ro = r + n;
    if (initial)
        memcpy(target, initial, n * sizeof(double));
    else
        for (size_t row = 0; row < g->n; row++)
            bac_basis_sum(g->collocation_basis + row, y, p, 0, target + row * p);
    int rc = -1;
    for (size_t iter = 0; iter < (compat(s) ? 1 : 15); iter++) {
        BacBand a;
        BacLu lu = {0};
        if (bac_band_new(&a, n, wid, wid))
            break;
        a.npde = compat(s) ? p : 0;
        a.kcol = g->kcol;
        a.nint = g->nint;
        double maxres = 0;
        int resolved = 1;
        for (size_t row = 0; row < g->n; row++) {
            const BacBasis *b = g->collocation_basis + row;
            int side = row == 0 ? 0 : row == g->n - 1 ? 1 : -1;
            state(b, y, p, side, compat(s), u, ux, uxx);
            if (side >= 0) {
                s->problem.boundary(s->problem.ctx, side, t, u, ux, r + row * p);
                boundary_derivatives(s, side, t, u, ux, du, dx, dt);
                memset(av, 0, p * sizeof(double));
                memset(ad, 0, p * sizeof(double));
                for (size_t j = 0; j < b->len; j++) {
                    for (size_t d = 0; d < p; d++) {
                        av[d] += fabs(b->b[j] * y[(b->start + j) * p + d]);
                        ad[d] += fabs(b->bx[j] * y[(b->start + j) * p + d]);
                    }
                    if (!compat(s))
                        for (size_t c = 0; c < p; c++)
                            for (size_t d = 0; d < p; d++)
                                bac_band_set(&a, row * p + c, (b->start + j) * p + d,
                                             du[c * p + d] * b->b[j] + dx[c * p + d] * b->bx[j]);
                }
                if (compat(s))
                    source_boundary_matrix(&a, row, b, p, du, dx, 1);
                for (size_t c = 0; c < p; c++) {
                    double scale = 0;
                    for (size_t d = 0; d < p; d++)
                        scale += fabs(du[c * p + d]) * av[d] + fabs(dx[c * p + d]) * ad[d];
                    ro[row * p + c] = roundoff(2 * b->len + 2 * p + 4, scale);
                }
            } else {
                for (size_t c = 0; c < p; c++) {
                    r[row * p + c] = compat(s) ? -target[row * p + c] : u[c] - target[row * p + c];
                    if (compat(s))
                        for (size_t j = 0; j < b->len; j++)
                            r[row * p + c] += b->b[j] * y[(b->start + j) * p + c];
                    double scale = fabs(target[row * p + c]);
                    for (size_t j = 0; j < b->len; j++) {
                        bac_band_set(&a, row * p + c, (b->start + j) * p + c, b->b[j]);
                        scale += fabs(b->b[j] * y[(b->start + j) * p + c]);
                    }
                    ro[row * p + c] =
                        fabs(r[row * p + c]) >= 1e-11 ? roundoff(2 * b->len + 1, scale) : 0;
                }
            }
        }
        for (size_t i = 0; i < n; i++) {
            if (!isfinite(r[i])) {
                bac_band_free(&a);
                goto done;
            }
            maxres = fmax(maxres, fabs(r[i]));
            if (fabs(r[i]) >= fmax(1e-11, ro[i]))
                resolved = 0;
        }
        if (!compat(s) && maxres < 1e-11) {
            bac_band_free(&a);
            rc = 0;
            break;
        }
        if (bac_band_factor(&a, &lu) || bac_lu_solve(&lu, r)) {
            bac_lu_free(&lu);
            bac_band_free(&a);
            break;
        }
        bac_lu_free(&lu);
        bac_band_free(&a);
        double mc = 0;
        for (size_t i = 0; i < n; i++)
            mc = fmax(mc, fabs(r[i]) / (1 + fabs(y[i])));
        if (!compat(s) && resolved && mc <= roundoff(2 * (g->kcol + 2) + 2 * p + 4, 1)) {
            rc = 0;
            break;
        }
        for (size_t i = 0; i < n; i++)
            y[i] -= r[i];
        if (compat(s))
            rc = 0;
    }
done:
    free(work);
    return rc;
}
static int consistent_derivative(BacSolver *s, double t, const double *y, double *yp) {
    BacGrid *g = &s->grid;
    size_t p = s->problem.components, n = g->n * p, pp = p * p, wid = (g->kcol + 1) * p;
    BacBand a;
    BacLu lu = {0};
    if (bac_band_new(&a, n, wid, wid))
        return -1;
    a.npde = compat(s) ? p : 0;
    a.kcol = g->kcol;
    a.nint = g->nint;
    double *work = calloc(4 * p + 2 * pp + n, sizeof(double));
    if (!work) {
        bac_band_free(&a);
        return -1;
    }
    double *u = work, *ux = u + p, *uxx = ux + p, *dt = uxx + p, *du = dt + p, *dx = du + pp,
           *cons = dx + pp;
    int use = conservative_rhs(s, t, y, cons);
    for (size_t row = 0; row < g->n; row++) {
        const BacBasis *b = g->collocation_basis + row;
        int side = row == 0 ? 0 : row == g->n - 1 ? 1 : -1;
        state(b, y, p, side, compat(s), u, ux, uxx);
        if (side >= 0) {
            boundary_derivatives(s, side, t, u, ux, du, dx, dt);
            for (size_t c = 0; c < p; c++)
                yp[row * p + c] = -dt[c];
            if (compat(s))
                source_boundary_matrix(&a, row, b, p, du, dx, 1);
            else
                for (size_t j = 0; j < b->len; j++)
                    for (size_t c = 0; c < p; c++)
                        for (size_t d = 0; d < p; d++)
                            bac_band_set(&a, row * p + c, (b->start + j) * p + d,
                                         du[c * p + d] * b->b[j] + dx[c * p + d] * b->bx[j]);
        } else {
            if (use > 0)
                memcpy(yp + row * p, cons + row * p, p * sizeof(double));
            else
                s->problem.rhs(s->problem.ctx, t, g->collocation[row], u, ux, uxx, yp + row * p);
            for (size_t j = 0; j < b->len; j++)
                for (size_t c = 0; c < p; c++)
                    bac_band_set(&a, row * p + c, (b->start + j) * p + c, b->b[j]);
        }
    }
    int rc = use < 0 ? -1 : bac_band_factor(&a, &lu);
    if (!rc)
        rc = bac_lu_solve(&lu, yp);
    bac_lu_free(&lu);
    bac_band_free(&a);
    free(work);
    return rc;
}
static void initial_cb(void *ctx, double x, double *out) {
    BacSolver *s = ctx;
    s->problem.initial(s->problem.ctx, x, out);
}
static int initial_coefficients(BacSolver *s, double t, double *y) {
    size_t p = s->problem.components, n = s->grid.n * p;
    double *initial = malloc(n * sizeof(double));
    if (!initial)
        return -1;
    for (size_t i = 0; i < s->grid.n; i++)
        initial_cb(s, s->grid.collocation[i], initial + i * p);
    int rc;
    if (compat(s)) {
        BacBand mat;
        BacLu lu = {0};
        size_t wid = (s->grid.kcol + 1) * p;
        rc = bac_band_new(&mat, n, wid, wid);
        if (!rc) {
            mat.npde = p;
            mat.kcol = s->grid.kcol;
            mat.nint = s->grid.nint;
            for (size_t row = 0; row < s->grid.n; row++) {
                BacBasis *b = s->grid.collocation_basis + row;
                for (size_t j = 0; j < b->len; j++)
                    for (size_t c = 0; c < p; c++)
                        bac_band_set(&mat, row * p + c, (b->start + j) * p + c, b->b[j]);
            }
            memcpy(y, initial, n * sizeof(double));
            rc = bac_band_factor(&mat, &lu);
            if (!rc)
                rc = bac_lu_solve(&lu, y);
            bac_lu_free(&lu);
            bac_band_free(&mat);
        }
    } else
        rc = bac_grid_interpolate(&s->grid, p, initial_cb, s, y);
    if (!rc)
        rc = enforce_boundary(s, t, y, initial);
    free(initial);
    return rc;
}
static TimeOptions coefficient_options(const BacSolver *s, const BacGrid *g, double **atol,
                                       double **rtol) {
    TimeOptions o = s->options.time;
    size_t p = s->problem.components, n = g->n * p;
    *atol = malloc(n * sizeof(double));
    *rtol = malloc(n * sizeof(double));
    if (*atol && *rtol)
        for (size_t i = 0; i < n; i++) {
            (*atol)[i] = s->atol[i % p];
            (*rtol)[i] = s->rtol[i % p];
        }
    o.atol = *atol;
    o.rtol = *rtol;
    o.atol_len = o.rtol_len = n;
    return o;
}
static int estimate(BacSolver *s, const double *y, BacSpatialError *e) {
    int rc = bac_grid_estimate(&s->grid, y, s->problem.components, s->atol, s->rtol,
                               s->options.estimator, e);
    if (rc || compat(s) || s->options.estimator != BAC_LOI || s->grid.nint < 2)
        return rc;
    BacSpatialError sci;
    if ((rc = bac_grid_estimate(&s->grid, y, s->problem.components, s->atol, s->rtol, BAC_SCI,
                                &sci))) {
        bac_spatial_error_free(e);
        return rc;
    }
    e->ratio = 0;
    for (size_t m = 0; m < s->problem.components; m++) {
        e->component_ratios[m] = hypot(e->component_ratios[m], sci.component_ratios[m]);
        e->ratio = fmax(e->ratio, e->component_ratios[m]);
    }
    double sum = 0, mx = 0;
    for (size_t i = 0; i < s->grid.nint; i++) {
        e->indicators[i] = pow(hypot(pow(e->indicators[i], (double)s->grid.kcol + 1),
                                     pow(sci.indicators[i], (double)s->grid.kcol + 2)),
                               1.0 / ((double)s->grid.kcol + 1));
        sum += e->indicators[i];
        mx = fmax(mx, e->indicators[i]);
    }
    e->distribution_ratio = sum > 0 ? mx * (double)s->grid.nint / sum : 1;
    bac_spatial_error_free(&sci);
    return 0;
}
static int spatial_passes(const BacSolver *s, const BacSpatialError *e, int correction) {
    if (e->distribution_ratio > 2)
        return 0;
    if (!compat(s) && s->options.estimator == BAC_LOI && s->grid.nint >= 2)
        return e->ratio < (correction ? .4 : 1);
    return correction ? (e->ratio < .4 && (s->grid.nint <= 15 || e->ratio > .1)) : e->ratio < 1;
}
int bac_solver_new(BacSolver *s, const BacProblem *p, const BacOptions *options, double t) {
    memset(s, 0, sizeof *s);
    if (!p || !p->components || !p->initial || !p->rhs || !p->boundary || !isfinite(t))
        return fail(s, "invalid PDE callbacks, dimension, or initial time");
    s->problem = *p;
    s->options = options ? *options : bac_options();
    BacOptions *o = &s->options;
    if (o->kcol < 3 || o->kcol > 10 || !o->max_remeshes || !o->intervals ||
        o->intervals > o->max_intervals)
        return fail(s, "invalid spatial order or mesh limits");
    if (o->estimator == BAC_SCI && o->intervals < 2)
        return fail(s, "SCI requires at least two intervals");
    size_t np = p->components;
    s->atol = malloc(np * sizeof(double));
    s->rtol = malloc(np * sizeof(double));
    if (!s->atol || !s->rtol)
        return fail(s, "out of memory");
    if ((o->atol && o->atol_len != 1 && o->atol_len != np) ||
        (o->rtol && o->rtol_len != 1 && o->rtol_len != np))
        return fail(s, "tolerance length must be one or PDE components");
    for (size_t m = 0; m < np; m++) {
        s->atol[m] = o->atol ? o->atol[o->atol_len == 1 ? 0 : m] : o->absolute;
        s->rtol[m] = o->rtol ? o->rtol[o->rtol_len == 1 ? 0 : m] : o->relative;
        if (!isfinite(s->atol[m]) || !isfinite(s->rtol[m]) || s->atol[m] < 0 || s->rtol[m] < 0 ||
            s->atol[m] + s->rtol[m] == 0)
            return fail(s, "invalid PDE tolerances");
    }
    double *mesh = NULL;
    if (!o->mesh) {
        mesh = malloc((o->intervals + 1) * sizeof(double));
        if (!mesh)
            return fail(s, "out of memory");
        for (size_t i = 0; i <= o->intervals; i++)
            mesh[i] = (double)i / (double)o->intervals;
    }
    int rc = bac_grid_new(&s->grid, o->mesh ? o->mesh : mesh, o->intervals, o->kcol);
    free(mesh);
    if (rc)
        return fail(s, "invalid spatial grid");
    double *y = NULL, *yp = NULL;
    size_t remeshes = 0;
    for (;;) {
        free(y);
        y = calloc(s->grid.n * np, sizeof(double));
        if (!y || initial_coefficients(s, t, y)) {
            free(y);
            return fail(s, "initial or boundary coefficients failed");
        }
        if (!o->spatial_adaptivity || compat(s))
            break;
        BacSpatialError e;
        if (estimate(s, y, &e)) {
            free(y);
            return fail(s, "initial spatial estimate failed");
        }
        s->statistics.spatial_error_ratio = e.ratio;
        if (spatial_passes(s, &e, 1)) {
            bac_spatial_error_free(&e);
            break;
        }
        if (remeshes++ >= o->max_remeshes) {
            bac_spatial_error_free(&e);
            free(y);
            return fail(s, "initial spatial mesh failed within max_remeshes");
        }
        double *newmesh = NULL;
        size_t ni = 0;
        rc = bac_grid_refined_mesh(&s->grid, &e, o->max_intervals, 0, !compat(s), &newmesh, &ni);
        bac_spatial_error_free(&e);
        if (rc) {
            free(y);
            return fail(s, "initial spatial mesh limit");
        }
        BacGrid newgrid;
        rc = bac_grid_new(&newgrid, newmesh, ni, o->kcol);
        free(newmesh);
        if (rc) {
            free(y);
            return fail(s, "initial refinement grid failed");
        }
        bac_grid_free(&s->grid);
        s->grid = newgrid;
    }
    yp = calloc(s->grid.n * np, sizeof(double));
    if (!yp || consistent_derivative(s, t, y, yp)) {
        free(y);
        free(yp);
        return fail(s, "consistent initial derivative failed");
    }
    double *a = NULL, *r = NULL;
    TimeOptions to = coefficient_options(s, &s->grid, &a, &r);
    rc = (!a || !r) ? -1 : bi_new(&s->integrator, t, y, yp, s->grid.n * np, &to);
    free(a);
    free(r);
    free(y);
    free(yp);
    if (rc)
        return fail(s, "time integrator initialization failed");
    s->statistics.remeshes = remeshes;
    s->statistics.intervals = s->grid.nint;
    return 0;
}
void bac_solver_clear_events(BacSolver *s) {
    free(s->events);
    free(s->event_contexts);
    free(s->event_directions);
    s->events = NULL;
    s->event_contexts = NULL;
    s->event_directions = NULL;
    s->event_count = 0;
    bi_clear_event_history(&s->integrator);
}
void bac_solver_free(BacSolver *s) {
    bac_solver_clear_events(s);
    bac_grid_free(&s->grid);
    bi_free(&s->integrator);
    free(s->atol);
    free(s->rtol);
    memset(s, 0, sizeof *s);
}
BacView bac_solver_view(const BacSolver *s) {
    BacView v = {bi_time(&s->integrator), &s->grid, s->problem.components, bi_state(&s->integrator),
                 bi_derivative(&s->integrator)};
    return v;
}
int bac_view_evaluate(const BacView *v, double x, double *u, double *ux, double *uxx, double *ut) {
    if (!isfinite(x) || x < v->grid->mesh[0] || x > v->grid->mesh[v->grid->nint])
        return -1;
    BacBasis b;
    if (bac_grid_cached_basis(v->grid, x, &b))
        return -1;
    if (u)
        bac_basis_sum(&b, v->y, v->npde, 0, u);
    if (ux)
        bac_basis_sum(&b, v->y, v->npde, 1, ux);
    if (uxx)
        bac_basis_sum(&b, v->y, v->npde, 2, uxx);
    if (ut)
        bac_basis_sum(&b, v->yp, v->npde, 0, ut);
    return 0;
}
double bac_view_integrate(const BacView *v,
                          double (*fn)(void *, double, const double *, const double *,
                                       const double *, const double *),
                          void *ctx) {
    size_t p = v->npde, q = v->grid->kcol + 3;
    if (!fn)
        return NAN;
    const BacQuadrature *quadrature = bac_grid_quadrature(v->grid);
    if (!quadrature)
        return NAN;
    double stack[32], *w = p <= 8 ? stack : malloc(4 * p * sizeof(double));
    if (!w)
        return NAN;
    double sum = 0;
    for (size_t i = 0; i < v->grid->nint * q; i++) {
        const BacQuadrature *point = quadrature + i;
        bac_basis_sum(&point->basis, v->y, p, 0, w);
        bac_basis_sum(&point->basis, v->y, p, 1, w + p);
        bac_basis_sum(&point->basis, v->y, p, 2, w + 2 * p);
        bac_basis_sum(&point->basis, v->yp, p, 0, w + 3 * p);
        sum += point->weight * fn(ctx, point->x, w, w + p, w + 2 * p, w + 3 * p);
    }
    if (w != stack)
        free(w);
    return sum;
}
static IntegrationStatistics add_stats(IntegrationStatistics a, IntegrationStatistics b) {
    a.accepted_steps += b.accepted_steps;
    a.rejected_steps += b.rejected_steps;
    a.residual_evaluations += b.residual_evaluations;
    a.jacobian_evaluations += b.jacobian_evaluations;
    a.newton_iterations += b.newton_iterations;
    a.linear_solves += b.linear_solves;
    a.error_test_failures += b.error_test_failures;
    a.convergence_failures += b.convergence_failures;
    a.maximum_order = MAX(a.maximum_order, b.maximum_order);
    return a;
}
BacStatistics bac_solver_statistics(const BacSolver *s) {
    BacStatistics st = s->statistics;
    st.time = add_stats(st.time, s->integrator.statistics);
    return st;
}
int bac_solver_add_event(BacSolver *s, int direction, BacEventFunction fn, void *ctx) {
    size_t n = s->event_count;
    if (!fn || direction < BI_FALLING || direction > BI_RISING)
        return -1;
    BacEventFunction *f = malloc((n + 1) * sizeof *f);
    void **contexts = malloc((n + 1) * sizeof *contexts);
    int *dirs = malloc((n + 1) * sizeof *dirs);
    if (!f || !contexts || !dirs) {
        free(f);
        free(contexts);
        free(dirs);
        return -1;
    }
    if (n) {
        memcpy(f, s->events, n * sizeof *f);
        memcpy(contexts, s->event_contexts, n * sizeof *contexts);
        memcpy(dirs, s->event_directions, n * sizeof *dirs);
    }
    f[n] = fn;
    contexts[n] = ctx;
    dirs[n] = direction;
    free(s->events);
    free(s->event_contexts);
    free(s->event_directions);
    s->events = f;
    s->event_contexts = contexts;
    s->event_directions = dirs;
    s->event_count++;
    return (int)n;
}
void bac_advance_result_free(BacAdvanceResult *r) {
    free(r->indices);
    free(r->directions);
    memset(r, 0, sizeof *r);
}
/* Two endpoint cardinal functions preserve every interior collocation value. */
static int project_dense_boundary(BacSolver *s, double t, double *y, double *yp) {
    BacGrid *g = &s->grid;
    size_t p = s->problem.components, size = 2 * p, pp = p * p;
    double *w = calloc(5 * p + 2 * pp + size * size + size, sizeof(double));
    if (!w)
        return -1;
    double *r = w, *u = r + size, *ux = u + p, *dtx = ux + p, *du = dtx + p, *dx = du + pp,
           *mat = dx + pp;
    if (!g->boundary_lift) {
        g->boundary_lift = calloc(2 * g->n, sizeof(double));
        double *lr = calloc(g->n, sizeof(double));
        if (!g->boundary_lift || !lr) {
            free(lr);
            free(w);
            free(g->boundary_lift);
            g->boundary_lift = NULL;
            return -1;
        }
        for (size_t col = 0; col < 2; col++) {
            memset(lr, 0, g->n * sizeof(double));
            lr[col == 0 ? 0 : g->n - 1] = 1;
            if (bac_lu_solve(&g->interpolation, lr)) {
                free(lr);
                free(w);
                free(g->boundary_lift);
                g->boundary_lift = NULL;
                return -1;
            }
            memcpy(g->boundary_lift + col * g->n, lr, g->n * sizeof(double));
        }
        free(lr);
        for (size_t side = 0; side < 2; side++) {
            const BacBasis *b = g->collocation_basis + (side ? g->n - 1 : 0);
            for (size_t col = 0; col < 2; col++) {
                g->lift_values[side][col] = g->lift_derivatives[side][col] = 0;
                for (size_t j = 0; j < b->len; j++) {
                    g->lift_values[side][col] +=
                        b->b[j] * g->boundary_lift[col * g->n + b->start + j];
                    g->lift_derivatives[side][col] +=
                        b->bx[j] * g->boundary_lift[col * g->n + b->start + j];
                }
            }
        }
    }
    double *lift = g->boundary_lift;
    double (*values)[2] = g->lift_values, (*derivatives)[2] = g->lift_derivatives;
    int rc = -1;
    for (size_t iter = 0; iter < 16; iter++) {
        double scale = 1;
        for (size_t side = 0; side < 2; side++) {
            const BacBasis *b = g->collocation_basis + (side ? g->n - 1 : 0);
            bac_basis_sum(b, y, p, 0, u);
            bac_basis_sum(b, y, p, 1, ux);
            s->problem.boundary(s->problem.ctx, (int)side, t, u, ux, r + side * p);
            boundary_derivatives(s, (int)side, t, u, ux, du, dx, dtx);
            for (size_t c = 0; c < p; c++) {
                scale = fmax(scale, fabs(u[c]));
                for (size_t col = 0; col < 2; col++)
                    for (size_t d = 0; d < p; d++)
                        mat[(side * p + c) * size + col * p + d] =
                            du[c * p + d] * values[side][col] +
                            dx[c * p + d] * derivatives[side][col];
            }
        }
        BacLu lu = {0};
        if (bac_dense_factor(mat, size, &lu) || bac_lu_solve(&lu, r)) {
            bac_lu_free(&lu);
            break;
        }
        bac_lu_free(&lu);
        double max = 0;
        for (size_t j = 0; j < size; j++)
            max = fmax(max, fabs(r[j]));
        for (size_t i = 0; i < g->n; i++)
            for (size_t c = 0; c < p; c++)
                for (size_t col = 0; col < 2; col++)
                    y[i * p + c] -= lift[col * g->n + i] * r[col * p + c];
        if (max <= 8 * DBL_EPSILON * scale) {
            double *ut = calloc(2 * p, sizeof(double));
            if (!ut)
                break;
            for (size_t side = 0; side < 2; side++) {
                const BacBasis *b = g->collocation_basis + (side ? g->n - 1 : 0);
                bac_basis_sum(b, y, p, 0, u);
                bac_basis_sum(b, y, p, 1, ux);
                boundary_derivatives(s, (int)side, t, u, ux, du, dx, dtx);
                bac_basis_sum(b, yp, p, 0, ut);
                bac_basis_sum(b, yp, p, 1, ut + p);
                for (size_t c = 0; c < p; c++) {
                    double val = dtx[c];
                    for (size_t d = 0; d < p; d++)
                        val += du[c * p + d] * ut[d] + dx[c * p + d] * ut[p + d];
                    r[side * p + c] = val;
                    for (size_t col = 0; col < 2; col++)
                        for (size_t d = 0; d < p; d++)
                            mat[(side * p + c) * size + col * p + d] =
                                du[c * p + d] * values[side][col] +
                                dx[c * p + d] * derivatives[side][col];
                }
            }
            free(ut);
            if (bac_dense_factor(mat, size, &lu) || bac_lu_solve(&lu, r)) {
                bac_lu_free(&lu);
                break;
            }
            bac_lu_free(&lu);
            for (size_t i = 0; i < g->n; i++)
                for (size_t c = 0; c < p; c++)
                    for (size_t col = 0; col < 2; col++)
                        yp[i * p + c] -= lift[col * g->n + i] * r[col * p + c];
            rc = 0;
            break;
        }
    }
    free(w);
    return rc;
}
static int root_function(void *ctx, double t, const double *y, const double *yp, double *out,
                         size_t count) {
    BacSolver *s = ctx;
    size_t n = s->grid.n * s->problem.components;
    double *yy = NULL, *yd = NULL;
    if (!compat(s)) {
        yy = malloc(n * sizeof(double));
        yd = malloc(n * sizeof(double));
        if (!yy || !yd) {
            free(yy);
            free(yd);
            return -1;
        }
        memcpy(yy, y, n * sizeof(double));
        memcpy(yd, yp, n * sizeof(double));
        if (project_dense_boundary(s, t, yy, yd)) {
            free(yy);
            free(yd);
            return -1;
        }
        y = yy;
        yp = yd;
    }
    BacView v = {t, &s->grid, s->problem.components, y, yp};
    for (size_t i = 0; i < count; i++)
        out[i] = s->events[i](s->event_contexts[i], &v);
    free(yy);
    free(yd);
    return 0;
}
static int roots_to(BacSolver *s, double output, BacAdvanceResult *out) {
    if (!s->event_count)
        return 0;
    RootEvent root = {0};
    int rc = bi_locate_roots(&s->integrator, root_function, s, s->event_count, s->event_directions,
                             output, &root);
    if (rc <= 0)
        return rc;
    if (compat(s))
        rc = bi_publish_dense_output(&s->integrator, root.time);
    else {
        rc = project_dense_boundary(s, root.time, root.state, root.derivative);
        if (!rc)
            rc = bi_reset(&s->integrator, root.time, root.state, root.derivative, s->integrator.n);
    }
    if (rc) {
        bi_root_free(&root);
        return -1;
    }
    out->time = root.time;
    out->event = 1;
    out->count = root.count;
    out->indices = root.indices;
    out->directions = root.directions;
    root.indices = NULL;
    root.directions = NULL;
    bi_root_free(&root);
    return 1;
}
int bac_solver_reinitialize(BacSolver *s) {
    BdfIntegrator *b = &s->integrator;
    size_t n = b->n;
    double *y = malloc(n * sizeof(double)), *yp = malloc(n * sizeof(double));
    if (!y || !yp) {
        free(y);
        free(yp);
        return fail(s, "out of memory");
    }
    memcpy(y, bi_state(b), n * sizeof(double));
    memcpy(yp, bi_derivative(b), n * sizeof(double));
    double t = compat(s) ? b->t : bi_time(b);
    int rc = 0;
    if (!compat(s)) {
        rc = enforce_boundary(s, t, y, NULL);
        if (!rc)
            rc = consistent_derivative(s, t, y, yp);
    }
    if (!rc)
        rc = bi_reset(b, t, y, yp, n);
    free(y);
    free(yp);
    return rc ? fail(s, "reinitialization failed") : 0;
}
typedef struct {
    const BacGrid *source, *target;
    size_t p;
    const BacSourceProjection *original;
} ProjectionContext;
static int project_history(void *ctx, const double *in, size_t oldn, double *out, size_t newn) {
    ProjectionContext *p = ctx;
    if (oldn != p->source->n * p->p || newn != p->target->n * p->p)
        return -1;
    return p->original ? bac_source_projection_apply(p->original, in, out)
                       : bac_grid_project(p->target, p->source, in, p->p, out);
}
static int replace_snapshot(BdfIntegrator *accepted, const BdfIntegrator *b) {
    return bi_clone_from(accepted, b);
}
/* The accepted grid normally is the active grid. Retain its ownership only
 * when a spatial trial actually replaces it; a state snapshot needs no new
 * knots, interpolation factorization or quadrature geometry. */
static void retire_trial_grid(BacGrid *old, BacGrid *accepted, int *owns_accepted) {
    if (!*owns_accepted && old->mesh == accepted->mesh) {
        *accepted = *old;
        *owns_accepted = 1;
        memset(old, 0, sizeof *old);
    } else {
        bac_grid_free(old);
    }
}
int bac_solver_advance(BacSolver *s, double tout, BacAdvanceResult *out) {
    memset(out, 0, sizeof *out);
    if (!isfinite(tout) || tout < bi_time(&s->integrator))
        return fail(s, "output time must be finite and no earlier than current time");
    BdfIntegrator accepted = {0};
    BacGrid accepted_grid = s->grid;
    int owns_accepted_grid = 0;
    BacStatistics astats = s->statistics;
    if (bi_clone(&accepted, &s->integrator)) {
        bi_free(&accepted);
        return fail(s, "snapshot allocation failed");
    }
    int result = -1;
    size_t consecutive = 0;
    int projection_refined = 0, initial_step = compat(s) && !accepted.source_count &&
                                               bac_solver_statistics(s).time.accepted_steps == 0;
    double *checkpoint = NULL;
    size_t p = s->problem.components;
    if (compat(s) && s->integrator.t > bi_time(&s->integrator)) {
        double output = fmin(tout, s->integrator.t);
        int rc = roots_to(s, output, out);
        if (rc < 0)
            goto rollback;
        if (rc > 0) {
            result = 0;
            goto done;
        }
        if (tout <= s->integrator.t) {
            if (bi_publish_dense_output(&s->integrator, tout))
                goto rollback;
            out->time = tout;
            result = 0;
            goto done;
        }
        bi_clear_published_output(&s->integrator);
        if (replace_snapshot(&accepted, &s->integrator))
            goto rollback;
    }
    while (bi_time(&s->integrator) < tout) {
        if (bac_solver_statistics(s).time.accepted_steps >= s->options.time.max_steps) {
            fail(s, "maximum cumulative time steps reached");
            goto rollback;
        }
        double old_t = bi_time(&s->integrator);
        DaeSystem sys = system_for(s);
        StepOutcome step;
        if (bi_step(&s->integrator, &sys, tout, !compat(s), &step)) {
            fail(s, s->integrator.error);
            goto rollback;
        }
        if (s->options.spatial_adaptivity) {
            const double *spatial_state =
                compat(s) ? s->integrator.source_states : bi_state(&s->integrator);
            BacSpatialError error;
            if (estimate(s, spatial_state, &error)) {
                fail(s, "spatial error estimation failed");
                goto rollback;
            }
            s->statistics.spatial_error_ratio = error.ratio;
            int preserves = projection_refined && consecutive > 0 && error.ratio < .4 &&
                            error.distribution_ratio <= 2;
            if (!spatial_passes(s, &error, initial_step || consecutive > 0) && !preserves) {
                consecutive++;
                s->statistics.spatial_rejections++;
                if (consecutive > s->options.max_remeshes) {
                    bac_spatial_error_free(&error);
                    fail(s, "spatial tolerance failed after max_remeshes");
                    goto rollback;
                }
                size_t prescribed = !initial_step && (consecutive == 1 || consecutive == 5)
                                        ? accepted_grid.nint
                                        : 0;
                double *mesh = NULL;
                size_t ni = 0;
                int rc = bac_grid_refined_mesh(&s->grid, &error, s->options.max_intervals,
                                               prescribed, !compat(s), &mesh, &ni);
                bac_spatial_error_free(&error);
                if (rc) {
                    fail(s, "spatial adaptation exceeds mesh capacity");
                    goto rollback;
                }
                BacGrid ng = {0};
                rc = bac_grid_new(&ng, mesh, ni, s->options.kcol);
                free(mesh);
                if (rc) {
                    fail(s, "replacement grid failed");
                    goto rollback;
                }
                if (initial_step) {
                    BacGrid old = s->grid;
                    s->grid = ng;
                    double *y = calloc(s->grid.n * p, sizeof(double)),
                           *yp = calloc(s->grid.n * p, sizeof(double));
                    rc = (!y || !yp) ? -1 : initial_coefficients(s, bi_time(&accepted), y);
                    if (!rc)
                        rc = consistent_derivative(s, bi_time(&accepted), y, yp);
                    double *a = NULL, *r = NULL;
                    TimeOptions to = coefficient_options(s, &s->grid, &a, &r);
                    BdfIntegrator next = {0};
                    if (!rc)
                        rc = (!a || !r)
                                 ? -1
                                 : bi_new(&next, bi_time(&accepted), y, yp, s->grid.n * p, &to);
                    free(a);
                    free(r);
                    free(y);
                    free(yp);
                    if (rc) {
                        bac_grid_free(&s->grid);
                        s->grid = old;
                        fail(s, "initial spatial restart failed");
                        goto rollback;
                    }
                    s->statistics.time = add_stats(s->statistics.time, s->integrator.statistics);
                    bi_free(&s->integrator);
                    s->integrator = next;
                    retire_trial_grid(&old, &accepted_grid, &owns_accepted_grid);
                    s->statistics.remeshes++;
                    s->statistics.intervals = s->grid.nint;
                    continue;
                }
                if (!checkpoint) {
                    checkpoint = malloc(accepted.n * sizeof(double));
                    if (!checkpoint) {
                        bac_grid_free(&ng);
                        goto rollback;
                    }
                    memcpy(checkpoint, compat(s) ? accepted.source_states : bi_state(&accepted),
                           accepted.n * sizeof(double));
                }
                projection_refined = 0;
                double *projected = NULL;
                for (;;) {
                    free(projected);
                    projected = calloc(ng.n * p, sizeof(double));
                    if (!projected ||
                        bac_grid_project(&ng, &accepted_grid, checkpoint, p, projected)) {
                        free(projected);
                        bac_grid_free(&ng);
                        fail(s, "mesh projection failed");
                        goto rollback;
                    }
                    if (compat(s))
                        break;
                    BacProjectionError transfer;
                    if (bac_grid_projection_measure(&ng, &accepted_grid, checkpoint, projected, p,
                                                    s->atol, s->rtol, &transfer)) {
                        free(projected);
                        bac_grid_free(&ng);
                        goto rollback;
                    }
                    if (transfer.ratio <= .2 && transfer.maximum <= 1) {
                        bac_projection_error_free(&transfer);
                        break;
                    }
                    projection_refined = 1;
                    double *refined = NULL;
                    size_t nn = 0;
                    rc = bac_grid_refine_projection(&ng, &accepted_grid, &transfer,
                                                    s->options.max_intervals, &refined, &nn);
                    bac_projection_error_free(&transfer);
                    if (rc) {
                        free(projected);
                        bac_grid_free(&ng);
                        fail(s, "mesh transfer refinement exceeds capacity");
                        goto rollback;
                    }
                    size_t nm = nn + 1;
                    BacGrid replacement = {0};
                    rc = bac_grid_new(&replacement, refined, nm - 1, s->options.kcol);
                    free(refined);
                    if (rc) {
                        free(projected);
                        bac_grid_free(&ng);
                        goto rollback;
                    }
                    bac_grid_free(&ng);
                    ng = replacement;
                }
                BacSourceProjection source_map = {0};
                if (compat(s)) {
                    if (bac_source_projection_new(&source_map, &ng, &accepted_grid, p) ||
                        bac_source_projection_apply(&source_map, checkpoint, projected)) {
                        bac_source_projection_free(&source_map);
                        free(projected);
                        bac_grid_free(&ng);
                        fail(s, "source mesh projection failed");
                        goto rollback;
                    }
                }
                ProjectionContext projection = {&accepted_grid, &ng, p,
                                                compat(s) ? &source_map : NULL};
                double *a = NULL, *r = NULL;
                TimeOptions to = coefficient_options(s, &ng, &a, &r);
                int warm = consecutive < 5 && a && r
                               ? bi_project_history(&s->integrator, &accepted, ng.n * p, &to,
                                                    project_history, &projection)
                               : 0;
                if (warm < 0) {
                    bac_source_projection_free(&source_map);
                    free(a);
                    free(r);
                    free(projected);
                    bac_grid_free(&ng);
                    fail(s, s->integrator.error);
                    goto rollback;
                }
                if (warm) {
                    s->statistics.warm_remeshes++;
                    if (compat(s) && accepted.last_order >= 2) {
                        size_t size = bac_abd_workspace(&source_map.factor, NULL), n = ng.n * p,
                               offset = (accepted.last_order - 2) * n;
                        if (offset + n <= size) {
                            double *workspace = malloc(size * sizeof(double));
                            if (workspace) {
                                bac_abd_workspace(&source_map.factor, workspace);
                                bi_restore_source_spare(&s->integrator, workspace + offset, n);
                                free(workspace);
                            }
                        }
                    }
                } else {
                    BacGrid old = s->grid;
                    s->grid = ng;
                    double *yp = calloc(ng.n * p, sizeof(double));
                    rc = yp ? 0 : -1;
                    if (!rc && !compat(s))
                        rc = enforce_boundary(s, old_t, projected, NULL);
                    if (!rc)
                        rc = consistent_derivative(s, old_t, projected, yp);
                    BdfIntegrator next = {0};
                    if (!rc)
                        rc = (!a || !r) ? -1 : bi_new(&next, old_t, projected, yp, ng.n * p, &to);
                    s->grid = old;
                    free(yp);
                    if (rc) {
                        bac_source_projection_free(&source_map);
                        free(a);
                        free(r);
                        free(projected);
                        bac_grid_free(&ng);
                        fail(s, "cold spatial restart failed");
                        goto rollback;
                    }
                    s->statistics.time = add_stats(s->statistics.time, s->integrator.statistics);
                    bi_free(&s->integrator);
                    s->integrator = next;
                    s->statistics.cold_remeshes++;
                }
                bac_source_projection_free(&source_map);
                free(a);
                free(r);
                free(projected);
                retire_trial_grid(&s->grid, &accepted_grid, &owns_accepted_grid);
                s->grid = ng;
                s->statistics.remeshes++;
                s->statistics.intervals = s->grid.nint;
                continue;
            }
            bac_spatial_error_free(&error);
        }
        consecutive = 0;
        projection_refined = 0;
        initial_step = 0;
        free(checkpoint);
        checkpoint = NULL;
        if (replace_snapshot(&accepted, &s->integrator))
            goto rollback;
        if (owns_accepted_grid)
            bac_grid_free(&accepted_grid);
        accepted_grid = s->grid;
        owns_accepted_grid = 0;
        astats = s->statistics;
        double output = compat(s) ? fmin(tout, s->integrator.t) : bi_time(&s->integrator);
        int rc = roots_to(s, output, out);
        if (rc < 0) {
            fail(s, s->integrator.error);
            goto rollback;
        }
        if (rc > 0) {
            result = 0;
            goto done;
        }
        if (compat(s) && tout <= s->integrator.t) {
            if (bi_publish_dense_output(&s->integrator, tout))
                goto rollback;
            out->time = tout;
            result = 0;
            goto done;
        }
    }
    out->time = bi_time(&s->integrator);
    result = 0;
    goto done;
rollback:
    if (owns_accepted_grid) {
        bac_grid_free(&s->grid);
        s->grid = accepted_grid;
        memset(&accepted_grid, 0, sizeof accepted_grid);
        owns_accepted_grid = 0;
    }
    bi_free(&s->integrator);
    s->integrator = accepted;
    memset(&accepted, 0, sizeof accepted);
    s->statistics = astats;
    if (!s->error[0])
        fail(s, "PDE advancement failed");
done:
    free(checkpoint);
    if (owns_accepted_grid)
        bac_grid_free(&accepted_grid);
    bi_free(&accepted);
    return result;
}
