/* Adaptive Gaussian collocation in C11. No global mutable state. */
#include "internal.h"
colnewsc_options colnewsc_default_options(void) {
    colnewsc_options o = {0};
    o.tolerance = 1e-5;
    o.collocation_points = 3;
    o.initial_subintervals = 5;
    o.max_subintervals = 10000;
    o.max_newton_iterations = 100;
    o.max_refinements = 20;
    return o;
}
const char *colnewsc_status_string(colnewsc_status e) {
    static const char *const names[] = {"success",
                                        "invalid input",
                                        "non-finite evaluation",
                                        "singular collocation Jacobian",
                                        "Newton iteration did not converge",
                                        "mesh limit",
                                        "refinement limit",
                                        "evaluation outside interval",
                                        "allocation failed"};
    return (unsigned)e < sizeof names / sizeof *names ? names[e] : "unknown error";
}
void colnewsc_free(colnewsc_solution *s) {
    if (s) {
        free(s->mesh);
        free(s->values);
        free(s->stages);
        free(s->extra);
        free(s->stats.estimated_errors);
        free(s);
    }
}
size_t colnewsc_dimension(const colnewsc_solution *s) {
    return s ? s->dimension : 0;
}
size_t colnewsc_mesh_count(const colnewsc_solution *s) {
    return s ? s->count : 0;
}
const double *colnewsc_mesh(const colnewsc_solution *s) {
    return s ? s->mesh : NULL;
}
const double *colnewsc_mesh_values(const colnewsc_solution *s) {
    return s ? s->values : NULL;
}
const colnewsc_statistics *colnewsc_stats(const colnewsc_solution *s) {
    return s ? &s->stats : NULL;
}
colnewsc_status colnewsc_evaluate(const colnewsc_solution *s, double x, double *out, size_t count) {
    if (!s || !out || count != s->dimension)
        return COLNEWSC_INVALID_INPUT;
    if (!isfinite(x) || x < s->mesh[0] || x > s->mesh[s->count - 1])
        return COLNEWSC_OUTSIDE_INTERVAL;
    size_t lo = 0, hi = s->count;
    while (lo < hi) {
        size_t m = lo + (hi - lo) / 2;
        if (s->mesh[m] <= x)
            lo = m + 1;
        else
            hi = m;
    }
    size_t i = CN_MIN(lo ? lo - 1 : 0, s->count - 2);
    if (x == s->mesh[i])
        memcpy(out, s->values + i * count, count * sizeof(double));
    else if (x == s->mesh[s->count - 1])
        memcpy(out, s->values + (s->count - 1) * count, count * sizeof(double));
    else
        cn_evaluate_interval(s, i, x, out);
    return cn_finite(out, count) ? COLNEWSC_OK : COLNEWSC_NON_FINITE;
}
static int compare_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static colnewsc_status validate(const colnewsc_problem *p, const colnewsc_options *o, double **mesh,
                                size_t *count) {
    if (!p || !o || !p->rhs || !p->boundary || !p->dimension || p->dimension > 160 ||
        p->left_boundary_conditions > p->dimension)
        return COLNEWSC_INVALID_INPUT;
    if (!isfinite(p->a) || !isfinite(p->b) || p->a >= p->b || !isfinite(p->b - p->a))
        return COLNEWSC_INVALID_INPUT;
    if (o->collocation_points < 1 || o->collocation_points > 7 || !isfinite(o->tolerance) ||
        o->tolerance <= 0)
        return COLNEWSC_INVALID_INPUT;
    if (o->component_tolerances) {
        if (o->component_tolerance_count != p->dimension)
            return COLNEWSC_INVALID_INPUT;
        for (size_t q = 0; q < p->dimension; q++)
            if (!isfinite(o->component_tolerances[q]) || o->component_tolerances[q] <= 0)
                return COLNEWSC_INVALID_INPUT;
    } else if (o->component_tolerance_count)
        return COLNEWSC_INVALID_INPUT;
    if (!o->initial_subintervals || !o->max_subintervals || !o->max_newton_iterations ||
        !o->max_refinements)
        return COLNEWSC_INVALID_INPUT;
    if (o->max_subintervals >
        SIZE_MAX / (p->dimension * o->collocation_points * (p->dimension + 1) * sizeof(double)) - 1)
        return COLNEWSC_INVALID_INPUT;
    if ((o->fixed_point_count && !o->fixed_points) || (o->initial_mesh_count && !o->initial_mesh))
        return COLNEWSC_INVALID_INPUT;
    size_t n;
    if (o->initial_mesh) {
        n = o->initial_mesh_count;
        if (n < 2 || o->initial_mesh[0] != p->a || o->initial_mesh[n - 1] != p->b ||
            !cn_finite(o->initial_mesh, n))
            return COLNEWSC_INVALID_INPUT;
        for (size_t i = 1; i < n; i++)
            if (o->initial_mesh[i - 1] >= o->initial_mesh[i])
                return COLNEWSC_INVALID_INPUT;
    } else {
        if (o->initial_subintervals > o->max_subintervals)
            return COLNEWSC_MESH_LIMIT;
        n = o->initial_subintervals + 1;
    }
    if (o->fixed_point_count > SIZE_MAX - n)
        return COLNEWSC_ALLOCATION_FAILED;
    double *x = cn_alloc(n + o->fixed_point_count, sizeof(double));
    if (!x)
        return COLNEWSC_ALLOCATION_FAILED;
    for (size_t i = 0; i < n; i++)
        x[i] = o->initial_mesh ? o->initial_mesh[i]
               : i == 0        ? p->a
               : i + 1 == n    ? p->b
                               : p->a + (p->b - p->a) * ((double)i / (double)(n - 1));
    for (size_t i = 1; i < n; i++)
        if (x[i - 1] >= x[i]) {
            free(x);
            return COLNEWSC_INVALID_INPUT;
        }
    for (size_t i = 0; i < o->fixed_point_count; i++) {
        double v = o->fixed_points[i];
        if (!isfinite(v) || v <= p->a || v >= p->b) {
            free(x);
            return COLNEWSC_INVALID_INPUT;
        }
        x[n + i] = v;
    }
    n += o->fixed_point_count;
    qsort(x, n, sizeof(double), compare_double);
    size_t unique = 1;
    for (size_t i = 1; i < n; i++)
        if (x[i] != x[unique - 1])
            x[unique++] = x[i];
    if (unique - 1 > o->max_subintervals) {
        free(x);
        return COLNEWSC_MESH_LIMIT;
    }
    *mesh = x;
    *count = unique;
    return COLNEWSC_OK;
}
colnewsc_status cn_checked_rhs(const colnewsc_problem *p, double x, const double *y, double *f) {
    size_t d = p->dimension;
    if (!isfinite(x) || !cn_finite(y, d))
        return COLNEWSC_NON_FINITE;
    for (size_t q = 0; q < d; q++)
        f[q] = NAN;
    p->rhs(p->context, x, y, f);
    return cn_finite(f, d) ? COLNEWSC_OK : COLNEWSC_NON_FINITE;
}
static colnewsc_status initialize(const colnewsc_problem *p, const double *mesh, size_t count,
                                  const cn_basis *basis, const colnewsc_solution *initial,
                                  colnewsc_solution **out) {
    size_t d = p->dimension, k = basis->k;
    colnewsc_solution *s = cn_alloc(1, sizeof *s);
    if (!s)
        return COLNEWSC_ALLOCATION_FAILED;
    s->count = count;
    s->dimension = d;
    s->basis = *basis;
    s->mesh = cn_alloc(count, sizeof(double));
    s->values = cn_alloc(count * d, sizeof(double));
    s->stages = cn_alloc((count - 1) * k * d, sizeof(double));
    double *y = cn_alloc(d, sizeof(double));
    colnewsc_status e = COLNEWSC_OK;
    if (!s->mesh || !s->values || !s->stages || !y) {
        e = COLNEWSC_ALLOCATION_FAILED;
        goto done;
    }
    memcpy(s->mesh, mesh, count * sizeof(double));
    for (size_t i = 0; i < count; i++) {
        double *v = s->values + i * d;
        if (initial)
            e = colnewsc_evaluate(initial, mesh[i], v, d);
        else if (p->initial_guess)
            p->initial_guess(p->context, mesh[i], v);
        if (e)
            goto done;
        if (!cn_finite(v, d)) {
            e = COLNEWSC_NON_FINITE;
            goto done;
        }
    }
    for (size_t i = 0; i < count - 1; i++)
        for (size_t j = 0; j < k; j++) {
            double h = mesh[i + 1] - mesh[i], *v = s->stages + (i * k + j) * d;
            if (initial) {
                double x = mesh[i] + h * basis->nodes[j];
                e = colnewsc_evaluate(initial, x, y, d);
                if (e)
                    goto done;
                e = cn_checked_rhs(p, x, y, v);
                if (e)
                    goto done;
            } else
                for (size_t q = 0; q < d; q++)
                    v[q] = (s->values[(i + 1) * d + q] - s->values[i * d + q]) / h;
        }
    if (!cn_finite(s->stages, (count - 1) * k * d))
        e = COLNEWSC_NON_FINITE;
done:
    free(y);
    if (e)
        colnewsc_free(s);
    else
        *out = s;
    return e;
}
static colnewsc_status gradient(const colnewsc_problem *p, size_t q, const double *y, double *g,
                                double *z) {
    size_t d = p->dimension;
    memset(g, 0, d * sizeof(double));
    if (!p->boundary_gradient || !p->boundary_gradient(p->context, q, y, g)) {
        memcpy(z, y, d * sizeof(double));
        for (size_t j = 0; j < d; j++) {
            double h = cbrt(DBL_EPSILON) * (1 + fabs(y[j]));
            z[j] = y[j] + h;
            double plus = p->boundary(p->context, q, z);
            z[j] = y[j] - h;
            double minus = p->boundary(p->context, q, z);
            z[j] = y[j];
            g[j] = (plus - minus) / (2 * h);
        }
    }
    return cn_finite(g, d) ? COLNEWSC_OK : COLNEWSC_NON_FINITE;
}
static colnewsc_status jacobian(const colnewsc_problem *p, double x, const double *y, double *jac,
                                double *z, double *plus, double *minus) {
    size_t d = p->dimension;
    memset(jac, 0, d * d * sizeof(double));
    if (!p->jacobian || !p->jacobian(p->context, x, y, jac)) {
        memcpy(z, y, d * sizeof(double));
        for (size_t j = 0; j < d; j++) {
            double h = cbrt(DBL_EPSILON) * (1 + fabs(y[j]));
            z[j] = y[j] + h;
            CN_TRY(cn_checked_rhs(p, x, z, plus));
            z[j] = y[j] - h;
            CN_TRY(cn_checked_rhs(p, x, z, minus));
            z[j] = y[j];
            for (size_t i = 0; i < d; i++)
                jac[i * d + j] = (plus[i] - minus[i]) / (2 * h);
        }
    }
    return cn_finite(jac, d * d) ? COLNEWSC_OK : COLNEWSC_NON_FINITE;
}
static void stage_state(const colnewsc_solution *s, size_t i, size_t j, double *y) {
    size_t d = s->dimension, k = s->basis.k;
    double h = s->mesh[i + 1] - s->mesh[i];
    if (d == 1) {
        double value = s->values[i];
        for (size_t l = 0; l < k; l++) {
            double w = h * s->basis.a[j * k + l];
            value += w * s->stages[i * k + l];
        }
        y[0] = value;
        return;
    }
    if (d == 2) {
        double first = s->values[2 * i], second = s->values[2 * i + 1];
        for (size_t l = 0; l < k; l++) {
            double w = h * s->basis.a[j * k + l];
            first += w * s->stages[(i * k + l) * 2];
            second += w * s->stages[(i * k + l) * 2 + 1];
        }
        y[0] = first;
        y[1] = second;
        return;
    }
    memcpy(y, s->values + i * d, d * sizeof(double));
    for (size_t l = 0; l < k; l++) {
        double w = h * s->basis.a[j * k + l];
        for (size_t q = 0; q < d; q++)
            y[q] += w * s->stages[(i * k + l) * d + q];
    }
}
typedef struct {
    double *state, *rhs, *stages, *boundary, *allocation;
} residual_work;
static int residual_new(residual_work *w, size_t d, size_t stages) {
    memset(w, 0, sizeof *w);
    w->allocation = cn_alloc(3 * d + stages, sizeof(double));
    if (!w->allocation)
        return 0;
    w->state = w->allocation;
    w->rhs = w->state + d;
    w->boundary = w->rhs + d;
    w->stages = w->boundary + d;
    return 1;
}
static colnewsc_status residual(const colnewsc_problem *p, const colnewsc_solution *s,
                                residual_work *w, double *norm_out) {
    size_t d = s->dimension, k = s->basis.k, n = s->count - 1;
    double norm = 0;
    for (size_t i = 0; i < n; i++) {
        double h = s->mesh[i + 1] - s->mesh[i];
        for (size_t j = 0; j < k; j++) {
            stage_state(s, i, j, w->state);
            CN_TRY(cn_checked_rhs(p, s->mesh[i] + h * s->basis.nodes[j], w->state, w->rhs));
            for (size_t q = 0; q < d; q++) {
                size_t off = (i * k + j) * d + q;
                double v = s->stages[off], fq = w->rhs[q], diff = v - fq;
                if (!isfinite(diff))
                    return COLNEWSC_NON_FINITE;
                w->stages[off] = fq - v;
                norm = fmax(norm, fabs(diff) / (1 + fmax(fabs(v), fabs(fq))));
            }
        }
        for (size_t q = 0; q < d; q++) {
            double r = s->values[(i + 1) * d + q] - s->values[i * d + q];
            for (size_t j = 0; j < k; j++)
                r -= h * s->basis.b[j] * s->stages[(i * k + j) * d + q];
            if (!isfinite(r))
                return COLNEWSC_NON_FINITE;
            norm = fmax(norm, fabs(r) / (1 + fmax(fabs(s->values[(i + 1) * d + q]),
                                                  fabs(s->values[i * d + q]))));
        }
    }
    for (size_t q = 0; q < d; q++) {
        size_t i = q < p->left_boundary_conditions ? 0 : n;
        const double *y = s->values + i * d;
        double v = p->boundary(p->context, q, y);
        if (!isfinite(v))
            return COLNEWSC_NON_FINITE;
        w->boundary[q] = -v;
        double ymax = 0;
        for (size_t j = 0; j < d; j++)
            ymax = fmax(ymax, fabs(y[j]));
        norm = fmax(norm, fabs(v) / (1 + ymax));
    }
    if (!isfinite(norm))
        return COLNEWSC_NON_FINITE;
    *norm_out = norm;
    return COLNEWSC_OK;
}
static double scaled_correction(const colnewsc_solution *s, const double *dy, const double *dk) {
    size_t d = s->dimension, k = s->basis.k;
    double norm = 0;
    for (size_t i = 0; i < s->count * d; i++) {
        double v = dy[i] / (1 + fabs(s->values[i]));
        norm += v * v;
    }
    for (size_t i = 0; i < s->count - 1; i++) {
        double h = s->mesh[i + 1] - s->mesh[i];
        for (size_t j = 0; j < k; j++)
            for (size_t q = 0; q < d; q++) {
                double v = h * dk[(i * k + j) * d + q] /
                           (1 + fabs(s->values[i * d + q]) + fabs(s->values[(i + 1) * d + q]));
                norm += v * v;
            }
    }
    return sqrt(norm);
}
static colnewsc_status simplified_correction(const colnewsc_solution *s, size_t left,
                                             const double *maps, const cn_lu *factors,
                                             const cn_lu *global, residual_work *res, double *dy,
                                             double *a) {
    size_t d = s->dimension, k = s->basis.k, n = s->count - 1, cols = d + 1, mi = k * d * cols;
    for (size_t i = 0; i < n; i++) {
        double h = s->mesh[i + 1] - s->mesh[i];
        for (size_t j = 0; j < k; j++)
            for (size_t q = 0; q < d; q++)
                a[q * k + j] = res->stages[(i * k + j) * d + q];
        CN_TRY(cn_lu_solve(factors + i, a, 1));
        for (size_t q = 0; q < d; q++) {
            size_t row = left + i * d + q;
            dy[row] = s->values[i * d + q] - s->values[(i + 1) * d + q];
            for (size_t j = 0; j < k; j++) {
                dy[row] += h * s->basis.b[j] * (s->stages[(i * k + j) * d + q] + a[q * k + j]);
                res->stages[(i * k + j) * d + q] = a[q * k + j];
            }
        }
    }
    for (size_t q = 0; q < d; q++)
        dy[q < left ? q : n * d + q] = res->boundary[q];
    CN_TRY(cn_lu_solve(global, dy, 1));
    for (size_t i = 0; i < n; i++) {
        const double *map = maps + i * mi;
        for (size_t j = 0; j < k; j++)
            for (size_t q = 0; q < d; q++) {
                size_t r = q * k + j, stage = (i * k + j) * d + q;
                for (size_t c = 0; c < d; c++)
                    res->stages[stage] += map[r * cols + c] * dy[i * d + c];
            }
    }
    return COLNEWSC_OK;
}
static colnewsc_status newton(const colnewsc_problem *p, colnewsc_solution *s,
                              const colnewsc_options *o, size_t *iterations) {
    size_t d = s->dimension, k = s->basis.k, n = s->count - 1, kd = k * d, cols = d + 1,
           mi = kd * cols, left = p->left_boundary_conditions;
    size_t nv = s->count * d, ns = n * kd;
    double tol = o->tolerance;
    if (o->component_tolerances) {
        tol = INFINITY;
        for (size_t q = 0; q < d; q++)
            tol = fmin(tol, o->component_tolerances[q]);
    }
    double target = fmin(fmax(tol * 0.002, 2e-13), 1e-10);
    residual_work work = {0}, tw = {0};
    colnewsc_solution trial = *s;
    trial.values = NULL;
    trial.stages = NULL;
    trial.extra = NULL;
    double *jac = NULL, *z = NULL, *plus = NULL, *minus = NULL, *bg = NULL, *rhs = NULL, *dk = NULL,
           *maps = NULL, *tn = NULL, *a = NULL, *matrix_work = NULL, *factor_values = NULL;
    size_t *factor_pivots = NULL;
    cn_lu *factors = NULL, global = {0};
    colnewsc_status e = COLNEWSC_ALLOCATION_FAILED;
    if (!residual_new(&work, d, ns) || !residual_new(&tw, d, ns))
        goto done;
    jac = cn_alloc(d * d, sizeof(double));
    z = cn_alloc(d, sizeof(double));
    plus = cn_alloc(d, sizeof(double));
    minus = cn_alloc(d, sizeof(double));
    bg = cn_alloc(d, sizeof(double));
    rhs = cn_alloc(nv, sizeof(double));
    dk = cn_alloc(ns, sizeof(double));
    maps = cn_alloc(n * mi, sizeof(double));
    tn = cn_alloc(nv, sizeof(double));
    a = cn_alloc(kd, sizeof(double));
    factors = cn_alloc(n, sizeof(cn_lu));
    trial.values = cn_alloc(nv, sizeof(double));
    trial.stages = cn_alloc(ns, sizeof(double));
    matrix_work = cn_alloc(kd * kd, sizeof(double));
    if (!jac || !z || !plus || !minus || !bg || !rhs || !dk || !maps || !tn || !a || !factors ||
        !trial.values || !trial.stages || !matrix_work)
        goto done;
    /* Small stage systems share one slab. Larger systems retain compact
     * per-interval factors so uncoupled large systems keep band-sized memory. */
    if (d <= 2 && kd <= 16) {
        if (n > SIZE_MAX / (kd * kd))
            goto done;
        factor_values = cn_alloc(n * kd * kd, sizeof(double));
        factor_pivots = cn_alloc(n * kd, sizeof(size_t));
        if (!factor_values || !factor_pivots)
            goto done;
        for (size_t i = 0; i < n; i++) {
            factors[i].data = factor_values + i * kd * kd;
            factors[i].pivots = factor_pivots + i * kd;
            factors[i].data_capacity = kd * kd;
            factors[i].pivot_capacity = kd;
            factors[i].borrowed = 1;
        }
    }
#define NT(x)                                                                                      \
    do {                                                                                           \
        e = (x);                                                                                   \
        if (e)                                                                                     \
            goto done;                                                                             \
    } while (0)
    for (size_t it = 0; it < o->max_newton_iterations; it++) {
        double norm;
        NT(residual(p, s, &work, &norm));
        if (o->trace)
            fprintf(stderr, "Newton: %zu intervals, residual %.4e\n", n, norm);
        if (norm <= target) {
            e = COLNEWSC_OK;
            goto done;
        }
        (*iterations)++;
        NT(cn_band_reset(&global, nv, d + left - 1, CN_MAX(d - 1, d - left)));
        for (size_t i = 0; i < n; i++) {
            double h = s->mesh[i + 1] - s->mesh[i], *m = matrix_work, *map = maps + i * mi;
            for (size_t j = 0; j < k; j++) {
                stage_state(s, i, j, work.state);
                e = jacobian(p, s->mesh[i] + h * s->basis.nodes[j], work.state, jac, z, plus,
                             minus);
                if (e) {
                    goto done;
                }
                for (size_t q = 0; q < d; q++) {
                    size_t r = q * k + j;
                    map[r * cols + d] = work.stages[(i * k + j) * d + q];
                    for (size_t c = 0; c < d; c++)
                        map[r * cols + c] = jac[q * d + c];
                    for (size_t l = 0; l < k; l++)
                        for (size_t c = 0; c < d; c++) {
                            size_t col = c * k + l;
                            m[r * kd + col] =
                                (col == r ? 1.0 : 0.0) - h * s->basis.a[j * k + l] * jac[q * d + c];
                        }
                }
            }
            NT(cn_stage_factor_reuse(m, kd, factors + i));
            NT(cn_lu_solve(factors + i, map, cols));
            for (size_t q = 0; q < d; q++) {
                size_t row = left + i * d + q;
                double r = s->values[i * d + q] - s->values[(i + 1) * d + q];
                for (size_t j = 0; j < k; j++)
                    r += h * s->basis.b[j] *
                         (s->stages[(i * k + j) * d + q] + map[(q * k + j) * cols + d]);
                rhs[row] = r;
                for (size_t c = 0; c < d; c++) {
                    double v = c == q ? -1.0 : 0.0;
                    for (size_t j = 0; j < k; j++)
                        v -= h * s->basis.b[j] * map[(q * k + j) * cols + c];
                    cn_band_set(&global, row, i * d + c, v);
                }
                cn_band_set(&global, row, (i + 1) * d + q, 1);
            }
        }
        for (size_t q = 0; q < d; q++) {
            size_t i = q < left ? 0 : n, row = q < left ? q : n * d + q;
            rhs[row] = work.boundary[q];
            NT(gradient(p, q, s->values + i * d, bg, z));
            for (size_t c = 0; c < d; c++)
                cn_band_set(&global, row, i * d + c, bg[c]);
        }
        NT(cn_lu_factor(&global));
        NT(cn_lu_solve(&global, rhs, 1));
        for (size_t i = 0; i < n; i++) {
            const double *map = maps + i * mi;
            for (size_t j = 0; j < k; j++)
                for (size_t q = 0; q < d; q++) {
                    size_t r = q * k + j, stage = (i * k + j) * d + q;
                    dk[stage] = map[r * cols + d];
                    for (size_t c = 0; c < d; c++)
                        dk[stage] += map[r * cols + c] * rhs[i * d + c];
                }
        }
        double alpha = 1, correction_norm = scaled_correction(s, rhs, dk);
        int accepted = 0;
        for (size_t damping = 0; damping < 16; damping++) {
            for (size_t i = 0; i < nv; i++)
                trial.values[i] = s->values[i] + alpha * rhs[i];
            for (size_t i = 0; i < ns; i++)
                trial.stages[i] = s->stages[i] + alpha * dk[i];
            double next;
            e = residual(p, &trial, &tw, &next);
            if (!e &&
                (next <= target ||
                 (!simplified_correction(&trial, left, maps, factors, &global, &tw, tn, a) &&
                  scaled_correction(s, tn, tw.stages) < correction_norm * (1 - 0.05 * alpha)))) {
                accepted = 1;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            e = COLNEWSC_NEWTON_NOT_CONVERGED;
            goto done;
        }
        double *swap = s->values;
        s->values = trial.values;
        trial.values = swap;
        swap = s->stages;
        s->stages = trial.stages;
        trial.stages = swap;
        free(s->extra);
        s->extra = NULL;
    }
    e = COLNEWSC_NEWTON_NOT_CONVERGED;
done:
    if (factors)
        for (size_t i = 0; i < n; i++)
            cn_lu_free(factors + i);
    cn_lu_free(&global);
    free(factors);
    free(work.allocation);
    free(tw.allocation);
    free(jac);
    free(z);
    free(plus);
    free(minus);
    free(bg);
    free(rhs);
    free(dk);
    free(maps);
    free(tn);
    free(a);
    free(matrix_work);
    free(factor_values);
    free(factor_pivots);
    free(trial.values);
    free(trial.stages);
    return e;
#undef NT
}
static colnewsc_status double_mesh(const colnewsc_solution *s, size_t maximum, double **out,
                                   size_t *count) {
    size_t required = 2 * (s->count - 1);
    if (required > maximum)
        return COLNEWSC_MESH_LIMIT;
    double *x = cn_alloc(required + 1, sizeof(double));
    if (!x)
        return COLNEWSC_ALLOCATION_FAILED;
    for (size_t i = 0; i < s->count - 1; i++) {
        double mid = s->mesh[i] + 0.5 * (s->mesh[i + 1] - s->mesh[i]);
        if (mid <= s->mesh[i] || mid >= s->mesh[i + 1]) {
            free(x);
            return COLNEWSC_INVALID_INPUT;
        }
        x[2 * i] = s->mesh[i];
        x[2 * i + 1] = mid;
    }
    x[required] = s->mesh[s->count - 1];
    *out = x;
    *count = required + 1;
    return COLNEWSC_OK;
}
static colnewsc_status estimates(const colnewsc_solution *coarse, const colnewsc_solution *fine,
                                 const colnewsc_options *o, double *intervals, double *errors) {
    size_t d = fine->dimension, k = fine->basis.k, order = (k >= 2 && k <= 4) ? 2 * k : k + 1;
    double divisor = pow(2, (double)order) - 1;
    static const double samples[] = {0, 1.0 / 6, 1.0 / 3, 0.5, 2.0 / 3, 5.0 / 6, 1};
    double a[160], b[160];
    memset(errors, 0, d * sizeof(double));
    memset(intervals, 0, (coarse->count - 1) * sizeof(double));
    for (size_t i = 0; i < coarse->count - 1; i++)
        for (size_t sample = 0; sample < 7; sample++) {
            double x = coarse->mesh[i] + (coarse->mesh[i + 1] - coarse->mesh[i]) * samples[sample];
            if (sample == 0 || sample == 6) {
                size_t node = i + (sample == 6);
                memcpy(a, coarse->values + node * d, d * sizeof(double));
            } else
                cn_evaluate_interval(coarse, i, x, a);
            if (sample == 0 || sample == 3 || sample == 6) {
                size_t node = 2 * i + sample / 3;
                memcpy(b, fine->values + node * d, d * sizeof(double));
            } else
                cn_evaluate_interval(fine, 2 * i + (sample > 3), x, b);
            for (size_t q = 0; q < d; q++) {
                double e = fabs(a[q] - b[q]) / divisor;
                if (!isfinite(e))
                    return COLNEWSC_NON_FINITE;
                double tol = o->component_tolerances ? o->component_tolerances[q] : o->tolerance;
                errors[q] = fmax(errors[q], e);
                intervals[i] = fmax(intervals[i], e / (tol * (1 + fabs(b[q]))));
            }
        }
    return COLNEWSC_OK;
}
static colnewsc_status solve_impl(const colnewsc_problem *p, const colnewsc_options *o,
                                  const colnewsc_solution *initial, colnewsc_solution **out) {
    if (!out)
        return COLNEWSC_INVALID_INPUT;
    *out = NULL;
    double *mesh = NULL, *ratios = NULL, *errors = NULL;
    size_t count = 0, iterations = 0;
    cn_basis basis;
    colnewsc_solution *current = NULL, *fine = NULL, *next = NULL;
    colnewsc_status e = validate(p, o, &mesh, &count);
    if (e)
        return e;
    cn_basis_new(&basis, o->collocation_points);
#define ST(x)                                                                                      \
    do {                                                                                           \
        e = (x);                                                                                   \
        if (e)                                                                                     \
            goto done;                                                                             \
    } while (0)
    ST(initialize(p, mesh, count, &basis, initial, &current));
    free(mesh);
    mesh = NULL;
    for (;;) {
        e = newton(p, current, o, &iterations);
        if (!e)
            break;
        if (e != COLNEWSC_NEWTON_NOT_CONVERGED && e != COLNEWSC_SINGULAR_JACOBIAN)
            goto done;
        ST(double_mesh(current, o->max_subintervals, &mesh, &count));
        ST(initialize(p, mesh, count, &basis, initial, &next));
        free(mesh);
        mesh = NULL;
        colnewsc_free(current);
        current = next;
        next = NULL;
    }
    ST(cn_prepare(p, current));
    double previous_ratio = INFINITY;
    for (size_t refinement = 1; refinement <= o->max_refinements; refinement++) {
        ST(double_mesh(current, o->max_subintervals, &mesh, &count));
        ST(initialize(p, mesh, count, &basis, current, &fine));
        free(mesh);
        mesh = NULL;
        ST(newton(p, fine, o, &iterations));
        ST(cn_prepare(p, fine));
        ratios = cn_alloc(current->count - 1, sizeof(double));
        errors = cn_alloc(p->dimension, sizeof(double));
        if (!ratios || !errors) {
            e = COLNEWSC_ALLOCATION_FAILED;
            goto done;
        }
        ST(estimates(current, fine, o, ratios, errors));
        double ratio = 0;
        for (size_t i = 0; i < current->count - 1; i++)
            ratio = fmax(ratio, ratios[i]);
        if (o->trace)
            fprintf(stderr, "Refinement %zu: %zu intervals, error ratio %.4e\n", refinement,
                    fine->count - 1, ratio);
        if (ratio <= 0.2) {
            size_t d = p->dimension, n = fine->count - 1;
            double br = 0;
            for (size_t q = 0; q < d; q++) {
                size_t i = q < p->left_boundary_conditions ? 0 : n;
                double v = p->boundary(p->context, q, fine->values + i * d);
                if (!isfinite(v)) {
                    e = COLNEWSC_NON_FINITE;
                    goto done;
                }
                br = fmax(br, fabs(v));
            }
            fine->stats.newton_iterations = iterations;
            fine->stats.refinements = refinement;
            fine->stats.subintervals = n;
            fine->stats.estimated_errors = errors;
            errors = NULL;
            fine->stats.estimated_error_ratio = ratio;
            fine->stats.boundary_residual = br;
            *out = fine;
            fine = NULL;
            e = COLNEWSC_OK;
            goto done;
        }
        size_t order = o->collocation_points >= 2 && o->collocation_points <= 4
                           ? 2 * o->collocation_points
                           : o->collocation_points + 1;
        double minimum = ratio >= 0.5 * previous_ratio ? 2 : 1;
        previous_ratio = ratio;
        /* At most eight pieces per coarse interval. */
        mesh = cn_alloc(8 * (current->count - 1) + 1, sizeof(double));
        if (!mesh) {
            e = COLNEWSC_ALLOCATION_FAILED;
            goto done;
        }
        count = 1;
        mesh[0] = current->mesh[0];
        for (size_t i = 0; i < current->count - 1; i++) {
            size_t parts =
                (size_t)fmin(fmax(ceil(pow(ratios[i] / 0.1, 1.0 / (double)order)), minimum), 8);
            for (size_t j = 1; j <= parts; j++)
                mesh[count++] = j == parts
                                    ? current->mesh[i + 1]
                                    : current->mesh[i] + (current->mesh[i + 1] - current->mesh[i]) *
                                                             (double)j / (double)parts;
        }
        free(ratios);
        ratios = NULL;
        free(errors);
        errors = NULL;
        if (count - 1 > o->max_subintervals / 2) {
            e = COLNEWSC_MESH_LIMIT;
            goto done;
        }
        for (size_t i = 1; i < count; i++)
            if (mesh[i - 1] >= mesh[i]) {
                e = COLNEWSC_INVALID_INPUT;
                goto done;
            }
        int same = count == fine->count && !memcmp(mesh, fine->mesh, count * sizeof(double));
        if (same) {
            colnewsc_free(current);
            current = fine;
            fine = NULL;
        } else {
            ST(initialize(p, mesh, count, &basis, fine, &next));
            ST(newton(p, next, o, &iterations));
            ST(cn_prepare(p, next));
            colnewsc_free(current);
            current = next;
            next = NULL;
        }
        free(mesh);
        mesh = NULL;
        colnewsc_free(fine);
        fine = NULL;
    }
    e = COLNEWSC_REFINEMENT_LIMIT;
done:
    free(mesh);
    free(ratios);
    free(errors);
    colnewsc_free(current);
    colnewsc_free(fine);
    colnewsc_free(next);
    return e;
#undef ST
}
colnewsc_status colnewsc_solve(const colnewsc_problem *p, const colnewsc_options *o,
                               colnewsc_solution **out) {
    return solve_impl(p, o, NULL, out);
}
colnewsc_status colnewsc_solve_with_initial_solution(const colnewsc_problem *p,
                                                     const colnewsc_options *o,
                                                     const colnewsc_solution *initial,
                                                     colnewsc_solution **out) {
    if (!p || !initial || p->dimension != initial->dimension || p->a != initial->mesh[0] ||
        p->b != initial->mesh[initial->count - 1]) {
        if (out)
            *out = NULL;
        return COLNEWSC_INVALID_INPUT;
    }
    return solve_impl(p, o, initial, out);
}
