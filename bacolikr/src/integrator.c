/* Native C implicit DAE integration (BSD-2-Clause).
 * Predictor, correctors and DDSTP order controller retain source arithmetic. */
#include "integrator.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdatomic.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define PH(b, j) ((b)->phi + (j) * (b)->n)
#define UROUND (2.0 * DBL_EPSILON)
/* A factorization and its matrix are immutable after construction. Snapshots
 * keep the cache alive without duplicating the O(n * bandwidth) storage. */
struct BiNewtonCache {
    atomic_size_t references;
    BacBand matrix;
    BacLu factor;
};
static void cache_retain(struct BiNewtonCache *cache) {
    if (cache)
        atomic_fetch_add_explicit(&cache->references, 1, memory_order_relaxed);
}
static void cache_release(struct BiNewtonCache *cache) {
    if (cache && atomic_fetch_sub_explicit(&cache->references, 1, memory_order_acq_rel) == 1) {
        bac_band_free(&cache->matrix);
        bac_lu_free(&cache->factor);
        free(cache);
    }
}
/* Workspaces are not numerical state. A snapshot owns independent work
 * arrays, allocated only when that snapshot is subsequently advanced. */
static double *workspace(double **data, size_t *capacity, size_t length) {
    if (length > *capacity) {
        if (length > SIZE_MAX / sizeof(double))
            return NULL;
        double *next = realloc(*data, length * sizeof(double));
        if (!next)
            return NULL;
        *data = next;
        *capacity = length;
    }
    return *data;
}
static void move_workspaces(BdfIntegrator *destination, BdfIntegrator *source) {
    if (destination->n != source->n)
        return;
    destination->step_work = source->step_work;
    destination->step_work_len = source->step_work_len;
    destination->newton_work = source->newton_work;
    destination->newton_work_len = source->newton_work_len;
    source->step_work = source->newton_work = NULL;
    source->step_work_len = source->newton_work_len = 0;
}
static int fail(BdfIntegrator *b, const char *s) {
    snprintf(b->error, sizeof b->error, "%s", s);
    return -1;
}
static double signum(double a) {
    return copysign(1., a);
}
static double clamp(double x, double a, double z) {
    return MAX(a, MIN(z, x));
}
static double *zeros(size_t n) {
    if (n > SIZE_MAX / sizeof(double))
        return NULL;
    return calloc(n ? n : 1, sizeof(double));
}
static double *copy(const double *p, size_t n) {
    double *a = zeros(n);
    if (a && n)
        memcpy(a, p, n * sizeof(double));
    return a;
}
static int finite_vector(const double *a, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (!isfinite(a[i]))
            return 0;
    return 1;
}
static double root_tol(double a, double z) {
    return 100. * DBL_EPSILON * MAX(fabs(a) + fabs(z) + fabs(z - a), DBL_MIN);
}
static double maxnorm(const double *a, const double *w, size_t n) {
    double m = 0.;
    for (size_t i = 0; i < n; i++) {
        double v = fabs(a[i] / w[i]);
        if (isnan(v))
            return NAN;
        m = MAX(m, v);
    }
    return m;
}
static double rmsnorm(const double *a, const double *w, size_t n) {
    double m = maxnorm(a, w, n), s = 0.;
    if (m == 0. || !isfinite(m))
        return m;
    for (size_t i = 0; i < n; i++) {
        double v = (a[i] / w[i]) / m;
        s += v * v;
    }
    return m * sqrt(s / (double)n);
}
static double daskrnorm(const double *a, const double *w, size_t n) {
    double m = 0., s = 0.;
    for (size_t i = 0; i < n; i++) {
        double v = fabs(a[i] * (1. / w[i]));
        if (v > m)
            m = v;
    }
    if (m <= 0.)
        return 0.;
    for (size_t i = 0; i < n; i++) {
        double v = (a[i] * (1. / w[i])) / m;
        s += v * v;
    }
    return m * sqrt(s / (double)n);
}
static double norm(const BdfIntegrator *b, const double *a, const double *w) {
    return b->options.method == BI_DASKR     ? daskrnorm(a, w, b->n)
           : b->options.error_norm == BI_RMS ? rmsnorm(a, w, b->n)
                                             : maxnorm(a, w, b->n);
}
TimeOptions bi_options(double a, double r) {
    TimeOptions o = {0};
    o.absolute = a;
    o.relative = r;
    o.max_step = INFINITY;
    o.max_steps = 100000;
    o.max_order = 5;
    return o;
}
TimeOptions bi_daskr_options(double a, double r) {
    TimeOptions o = bi_options(a, r);
    o.method = BI_DASKR;
    o.error_norm = BI_RMS;
    return o;
}
static double absolute(const TimeOptions *o, size_t i) {
    return o->atol ? o->atol[o->atol_len == 1 ? 0 : i] : o->absolute;
}
static double relative(const TimeOptions *o, size_t i) {
    return o->rtol ? o->rtol[o->rtol_len == 1 ? 0 : i] : o->relative;
}
int bi_validate_options(const TimeOptions *o, size_t n, char *e, size_t cap) {
    const char *s = NULL;
    if (!o || !n)
        s = "a DAE must have at least one unknown";
    else if ((o->atol && o->atol_len != 1 && o->atol_len != n) ||
             (o->rtol && o->rtol_len != 1 && o->rtol_len != n))
        s = "tolerance must have length one or dimension";
    else if (o->method != BI_VERIFIED && o->method != BI_DASKR)
        s = "invalid time method";
    else if (o->error_norm != BI_MAXIMUM && o->error_norm != BI_RMS)
        s = "invalid time error norm";
    else if (o->initial_step < 0. || !isfinite(o->initial_step))
        s = "initial_step must be finite and positive, or zero for automatic";
    else if (isnan(o->max_step) || o->max_step <= 0.)
        s = "max_step must be positive";
    else if (!isfinite(o->min_step) || o->min_step < 0. || o->min_step > o->max_step)
        s = "invalid min_step";
    else if (!o->max_steps || o->max_order < 1 || o->max_order > 5)
        s = "max_steps must be positive and max_order between one and five";
    if (!s)
        for (size_t i = 0; i < n; i++) {
            double a = absolute(o, i), r = relative(o, i);
            if (!isfinite(a) || !isfinite(r) || a < 0. || r < 0.) {
                s = "tolerances must be finite and nonnegative";
                break;
            }
            if (a == 0. && r == 0.) {
                s = "both error tolerances are zero";
                break;
            }
        }
    if (s) {
        if (e && cap)
            snprintf(e, cap, "%s", s);
        return -1;
    }
    return 0;
}
static int weights(BdfIntegrator *b, const double *y, const double *z, double *w) {
    for (size_t i = 0; i < b->n; i++) {
        double m = b->options.method == BI_DASKR ? fabs(y[i]) : MAX(fabs(y[i]), fabs(z[i]));
        w[i] = relative(&b->options, i) * m + absolute(&b->options, i);
        if (!isfinite(w[i]) || w[i] <= 0.)
            return fail(b, "error weights are zero or nonfinite");
    }
    return 0;
}
void bi_invalidate_jacobian(BdfIntegrator *b) {
    cache_release(b->newton_cache);
    b->newton_cache = NULL;
    b->matrix = (BacBand){0};
    b->factor = (BacLu){0};
    b->has_matrix = 0;
}
void bi_free(BdfIntegrator *b) {
    if (!b)
        return;
    bi_invalidate_jacobian(b);
    free(b->y);
    free(b->yp);
    free(b->phi);
    free(b->source_states);
    free(b->atol_owned);
    free(b->rtol_owned);
    free(b->published_y);
    free(b->published_yp);
    free(b->initial_root_y);
    free(b->initial_root_yp);
    free(b->last_event_indices);
    free(b->step_work);
    free(b->newton_work);
    memset(b, 0, sizeof *b);
}
int bi_new(BdfIntegrator *b, double t, const double *y, const double *yp, size_t n,
           const TimeOptions *o) {
    memset(b, 0, sizeof *b);
    if (bi_validate_options(o, n, b->error, sizeof b->error))
        return -1;
    if (!isfinite(t) || !y || !yp || !finite_vector(y, n) || !finite_vector(yp, n))
        return fail(b, "initial time, state, and derivative must be finite");
    if (n > SIZE_MAX / 7)
        return fail(b, "dimension overflow");
    b->n = n;
    b->t = t;
    b->options = *o;
    b->y = copy(y, n);
    b->yp = copy(yp, n);
    b->phi = zeros(7 * n);
    b->source_states = zeros(6 * n);
    b->published_y = zeros(n);
    b->published_yp = zeros(n);
    b->initial_root_y = zeros(n);
    b->initial_root_yp = zeros(n);
    if (o->atol) {
        b->atol_owned = copy(o->atol, o->atol_len);
        b->options.atol = b->atol_owned;
    }
    if (o->rtol) {
        b->rtol_owned = copy(o->rtol, o->rtol_len);
        b->options.rtol = b->rtol_owned;
    }
    if (!b->y || !b->yp || !b->phi || !b->source_states || !b->published_y || !b->published_yp ||
        !b->initial_root_y || !b->initial_root_yp || (o->atol && !b->atol_owned) ||
        (o->rtol && !b->rtol_owned)) {
        bi_free(b);
        return fail(b, "out of memory");
    }
    b->order = 1;
    b->startup = 1;
    b->newton_convergence_factor = 100.;
    return 0;
}
int bi_clone(BdfIntegrator *d, const BdfIntegrator *s) {
    BdfIntegrator q;
    if (bi_new(&q, s->t, s->y, s->yp, s->n, &s->options))
        return -1;
    double *y = q.y, *yp = q.yp, *phi = q.phi, *states = q.source_states, *py = q.published_y,
           *pp = q.published_yp, *iy = q.initial_root_y, *ip = q.initial_root_yp, *a = q.atol_owned,
           *r = q.rtol_owned;
    TimeOptions o = q.options;
    q = *s;
    q.y = y;
    q.yp = yp;
    q.phi = phi;
    q.source_states = states;
    q.published_y = py;
    q.published_yp = pp;
    q.initial_root_y = iy;
    q.initial_root_yp = ip;
    q.atol_owned = a;
    q.rtol_owned = r;
    q.options = o;
    cache_retain(q.newton_cache);
    q.step_work = q.newton_work = NULL;
    q.step_work_len = q.newton_work_len = 0;
    q.last_event_indices = NULL;
    memcpy(q.phi, s->phi, 7 * s->n * sizeof(double));
    memcpy(q.source_states, s->source_states, 6 * s->n * sizeof(double));
    memcpy(q.published_y, s->published_y, s->n * sizeof(double));
    memcpy(q.published_yp, s->published_yp, s->n * sizeof(double));
    memcpy(q.initial_root_y, s->initial_root_y, s->n * sizeof(double));
    memcpy(q.initial_root_yp, s->initial_root_yp, s->n * sizeof(double));
    if (s->last_event_count) {
        q.last_event_indices = malloc(s->last_event_count * sizeof(size_t));
        if (!q.last_event_indices)
            goto oom;
        memcpy(q.last_event_indices, s->last_event_indices, s->last_event_count * sizeof(size_t));
    }
    *d = q;
    return 0;
oom:
    bi_free(&q);
    return -1;
}
int bi_clone_from(BdfIntegrator *d, const BdfIntegrator *s) {
    if (d == s)
        return 0;
    /* Tolerances and event storage are normally fixed between accepted steps.
     * The uncommon dimension/registry changes use the transactional slow path. */
    if (d->n != s->n || !d->y || !d->yp || !d->phi || !!d->options.atol != !!s->options.atol ||
        !!d->options.rtol != !!s->options.rtol ||
        (s->options.atol && d->options.atol_len != s->options.atol_len) ||
        (s->options.rtol && d->options.rtol_len != s->options.rtol_len) ||
        d->last_event_count < s->last_event_count) {
        BdfIntegrator q = {0};
        if (bi_clone(&q, s))
            return -1;
        move_workspaces(&q, d);
        bi_free(d);
        *d = q;
        return 0;
    }
    BdfIntegrator q = *s;
    q.y = d->y;
    q.yp = d->yp;
    q.phi = d->phi;
    q.source_states = d->source_states;
    q.published_y = d->published_y;
    q.published_yp = d->published_yp;
    q.initial_root_y = d->initial_root_y;
    q.initial_root_yp = d->initial_root_yp;
    q.atol_owned = d->atol_owned;
    q.rtol_owned = d->rtol_owned;
    q.last_event_indices = d->last_event_indices;
    q.step_work = d->step_work;
    q.newton_work = d->newton_work;
    q.step_work_len = d->step_work_len;
    q.newton_work_len = d->newton_work_len;
    if (s->options.atol) {
        memcpy(q.atol_owned, s->options.atol, s->options.atol_len * sizeof(double));
        q.options.atol = q.atol_owned;
    }
    if (s->options.rtol) {
        memcpy(q.rtol_owned, s->options.rtol, s->options.rtol_len * sizeof(double));
        q.options.rtol = q.rtol_owned;
    }
    memcpy(q.y, s->y, s->n * sizeof(double));
    memcpy(q.yp, s->yp, s->n * sizeof(double));
    memcpy(q.phi, s->phi, 7 * s->n * sizeof(double));
    if (s->source_count)
        memcpy(q.source_states, s->source_states, s->source_count * s->n * sizeof(double));
    if (s->published) {
        memcpy(q.published_y, s->published_y, s->n * sizeof(double));
        memcpy(q.published_yp, s->published_yp, s->n * sizeof(double));
    }
    if (s->has_initial_root) {
        memcpy(q.initial_root_y, s->initial_root_y, s->n * sizeof(double));
        memcpy(q.initial_root_yp, s->initial_root_yp, s->n * sizeof(double));
    }
    if (s->last_event_count)
        memcpy(q.last_event_indices, s->last_event_indices, s->last_event_count * sizeof(size_t));
    cache_retain(q.newton_cache);
    bi_invalidate_jacobian(d);
    *d = q;
    return 0;
}
int bi_reset(BdfIntegrator *b, double t, const double *y, const double *yp, size_t n) {
    BdfIntegrator q;
    if (bi_new(&q, t, y, yp, n, &b->options))
        return fail(b, q.error);
    q.statistics = b->statistics;
    if (b->has_last_event && fabs(t - b->last_event_time) <= root_tol(t, b->last_event_time)) {
        q.has_last_event = 1;
        q.last_event_time = b->last_event_time;
        q.last_event_count = b->last_event_count;
        q.last_event_indices = malloc(q.last_event_count * sizeof(size_t));
        if (q.last_event_count && !q.last_event_indices) {
            bi_free(&q);
            return fail(b, "out of memory");
        }
        if (q.last_event_count)
            memcpy(q.last_event_indices, b->last_event_indices,
                   q.last_event_count * sizeof(size_t));
    }
    move_workspaces(&q, b);
    bi_free(b);
    *b = q;
    return 0;
}
void bi_clear_published_output(BdfIntegrator *b) {
    b->published = 0;
}
double bi_time(const BdfIntegrator *b) {
    return b->published ? b->published_time : b->t;
}
const double *bi_state(const BdfIntegrator *b) {
    return b->published ? b->published_y : b->y;
}
const double *bi_derivative(const BdfIntegrator *b) {
    return b->published ? b->published_yp : b->yp;
}
void bi_clear_event_history(BdfIntegrator *b) {
    b->has_last_event = 0;
    b->last_event_count = 0;
    b->has_root_cursor = b->published;
    b->root_cursor = b->published_time;
}
static int evaluate(BdfIntegrator *b, const DaeSystem *s, double t, const double *y,
                    const double *yp, double cj, double *out) {
    for (size_t i = 0; i < b->n; i++)
        out[i] = NAN;
    s->residual(s->ctx, t, y, yp, cj, out);
    b->statistics.residual_evaluations++;
    return finite_vector(out, b->n) ? 0 : fail(b, "DAE residual returned a nonfinite value");
}
static int make_matrix(BdfIntegrator *b, const DaeSystem *s, double t, double h, double cj,
                       const double *y, const double *yp, const double *res, const double *w,
                       BacBand *m) {
    size_t n = b->n, l = s->banded ? MIN(s->lower, n - 1) : n - 1,
           u = s->banded ? MIN(s->upper, n - 1) : n - 1;
    int supplied = 0;
    if (bac_band_new(m, n, l, u))
        return fail(b, "out of memory");
    b->statistics.jacobian_evaluations++;
    if (s->jacobian) {
        double *dense = s->banded ? NULL : zeros(n * n);
        if (!s->banded && !dense) {
            bac_band_free(m);
            return fail(b, "out of memory");
        }
        supplied = s->jacobian(s->ctx, t, y, yp, cj, s->banded ? m : NULL, dense);
        if (supplied > 0 && !s->banded)
            for (size_t r = 0; r < n; r++)
                for (size_t c = 0; c < n; c++)
                    bac_band_set(m, r, c, dense[r * n + c]);
        free(dense);
        if (supplied < 0) {
            bac_band_free(m);
            return fail(b, "Jacobian callback failed");
        }
        if (supplied > 0)
            return 0;
        memset(m->values, 0, n * m->stride * sizeof(double));
    }
    int source = b->options.method == BI_DASKR, abd = source && s->npde > 0;
    size_t groups = MIN(l + u + 1, n);
    double *work = zeros(4 * n);
    if (!work) {
        bac_band_free(m);
        return fail(b, "out of memory");
    }
    double *wy = work, *wp = work + n, *rr = work + 2 * n, *inc = work + 3 * n;
    memcpy(wy, y, n * sizeof(double));
    memcpy(wp, yp, n * sizeof(double));
    if (abd) {
        m->npde = s->npde;
        m->kcol = s->kcol;
        m->nint = s->nint;
        groups = (s->kcol + 2) * s->npde;
    }
    for (size_t g = 0; g < groups; g++) {
        size_t first = g, step = groups;
        if (abd) {
            size_t top = 2 * s->npde, height = s->kcol * s->npde, width = (s->kcol + 2) * s->npde;
            if (g < top) {
                first = g;
                step = 2 * height;
            } else if (g < width - top) {
                first = g;
                step = height;
            } else {
                first = g;
                step = 2 * height;
            }
        }
        for (size_t c = first; c < n; c += step) {
            double state = MAX(fabs(y[c]), fabs(h * yp[c])),
                   weight = source ? 1. / (1. / w[c]) : w[c],
                   sq = sqrt(source ? UROUND : DBL_EPSILON),
                   scale = (abd || !source) ? sq * MAX(state, weight) : MAX(sq * state, weight),
                   sgn = source           ? copysign(1., h * yp[c])
                         : h * yp[c] < 0. ? -1.
                                          : 1.;
            inc[c] = (y[c] + sgn * scale) - y[c];
            if (inc[c] == 0.)
                inc[c] = sgn * sqrt(DBL_EPSILON);
            wy[c] = y[c] + inc[c];
            wp[c] = yp[c] + cj * inc[c];
        }
        if (evaluate(b, s, t, wy, wp, source ? cj : NAN, rr)) {
            free(work);
            bac_band_free(m);
            return -1;
        }
        for (size_t c = first; c < n; c += step) {
            double recip = 1. / inc[c];
            if (abd) {
                size_t height = s->kcol * s->npde;
                if (c < 2 * s->npde)
                    for (size_t r = 0; r < s->npde; r++)
                        bac_band_set(m, r, c, (rr[r] - res[r]) * recip);
                for (size_t k = 0; k < s->nint; k++)
                    if (c >= k * height && c < k * height + (s->kcol + 2) * s->npde)
                        for (size_t r = s->npde + k * height; r < s->npde + (k + 1) * height; r++)
                            bac_band_set(m, r, c, (rr[r] - res[r]) * recip);
                if (c >= s->nint * height)
                    for (size_t r = s->npde + s->nint * height; r < n; r++)
                        bac_band_set(m, r, c, (rr[r] - res[r]) * recip);
            } else
                for (size_t r = c > u ? c - u : 0; r < MIN(n, c + l + 1); r++)
                    bac_band_set(m, r, c, (rr[r] - res[r]) * recip);
            wy[c] = y[c];
            wp[c] = yp[c];
        }
    }
    free(work);
    return 0;
}
static int refresh(BdfIntegrator *b, const DaeSystem *s, double t, double h, double cj,
                   const double *y, const double *yp, const double *r, const double *w) {
    BacBand m = {0};
    BacLu f = {0};
    if (make_matrix(b, s, t, h, cj, y, yp, r, w, &m))
        return -1;
    if (bac_band_factor(&m, &f)) {
        bac_band_free(&m);
        return fail(b, "singular Newton Jacobian");
    }
    struct BiNewtonCache *cache = malloc(sizeof *cache);
    if (!cache) {
        bac_band_free(&m);
        bac_lu_free(&f);
        return fail(b, "out of memory");
    }
    atomic_init(&cache->references, 1);
    cache->matrix = m;
    cache->factor = f;
    bi_invalidate_jacobian(b);
    b->newton_cache = cache;
    b->matrix = m;
    b->factor = f;
    b->has_matrix = 1;
    b->matrix_banded = s->banded;
    b->matrix_lower = s->lower;
    b->matrix_upper = s->upper;
    b->matrix_cj = cj;
    b->matrix_step = b->statistics.accepted_steps;
    return 0;
}
static int cache_valid(const BdfIntegrator *b, const DaeSystem *s, double cj, int age) {
    double ratio = cj / b->matrix_cj;
    return b->has_matrix && ratio >= .6 && ratio <= 5. / 3. && b->matrix_banded == s->banded &&
           (!s->banded || (b->matrix_lower == s->lower && b->matrix_upper == s->upper)) &&
           (!age || b->statistics.accepted_steps - b->matrix_step < 20);
}
static int newton(BdfIntegrator *b, const DaeSystem *s, double t, double h, double cj,
                  const double *prediction, const double *pd, const double *w, double *y,
                  double *yp, double *err) {
    size_t n = b->n;
    double *work = workspace(&b->newton_work, &b->newton_work_len, 6 * n);
    if (!work)
        return fail(b, "out of memory");
    double *r = work, *sc = work + n, *ny = work + 2 * n, *np = work + 3 * n, *nr = work + 4 * n,
           *corr = work + 5 * n;
    int result = -1, source = b->options.method == BI_DASKR;
    if (source) {
        if (cj != b->previous_cj)
            b->newton_convergence_factor = 100.;
        b->previous_cj = cj;
        int cached = cache_valid(b, s, cj, 0);
        double tolnew = 100. * UROUND * daskrnorm(prediction, w, n);
        for (size_t attempt = 0; attempt < 2; attempt++) {
            memcpy(y, prediction, n * sizeof(double));
            memcpy(yp, pd, n * sizeof(double));
            memset(err, 0, n * sizeof(double));
            if (evaluate(b, s, t, y, yp, cj, r))
                goto done;
            int old = cached;
            if (!cached) {
                b->newton_convergence_factor = 100.;
                if (refresh(b, s, t, h, cj, y, yp, r, w))
                    goto done;
            }
            cached = 0;
            double scale = 2. / (1. + cj / b->matrix_cj), first = 0.;
            for (size_t it = 0; it < 4; it++) {
                for (size_t i = 0; i < n; i++)
                    r[i] *= scale;
                if (bac_lu_solve(&b->factor, r)) {
                    fail(b, "linear solve failed");
                    goto done;
                }
                b->statistics.linear_solves++;
                b->statistics.newton_iterations++;
                for (size_t i = 0; i < n; i++) {
                    y[i] -= r[i];
                    err[i] -= r[i];
                    yp[i] -= cj * r[i];
                }
                double v = daskrnorm(r, w, n);
                if (!isfinite(v)) {
                    fail(b, "DASKR Newton correction is nonfinite");
                    goto done;
                }
                if (!it) {
                    first = v;
                    if (v <= tolnew) {
                        result = 0;
                        goto done;
                    }
                } else {
                    double rate = pow(v / first, 1. / (double)it);
                    if (rate > .9)
                        break;
                    b->newton_convergence_factor = rate / (1. - rate);
                }
                if (b->newton_convergence_factor * v <= .33) {
                    result = 0;
                    goto done;
                }
                if (it == 3)
                    break;
                if (evaluate(b, s, t, y, yp, cj, r))
                    goto done;
            }
            if (!old) {
                fail(b, "DASKR Newton iteration did not converge with a fresh matrix");
                goto done;
            }
        }
        fail(b, "DASKR Newton iteration did not converge");
        goto done;
    }
    memcpy(y, prediction, n * sizeof(double));
    memcpy(yp, pd, n * sizeof(double));
    if (evaluate(b, s, t, y, yp, NAN, r))
        goto done;
    if (!cache_valid(b, s, cj, 1) && refresh(b, s, t, h, cj, y, yp, r, w))
        goto done;
    bac_band_row_norms(&b->matrix, w, sc);
    double previous = INFINITY;
    size_t mit = 0, refreshes = 0;
    double roundoff = MIN(100. * DBL_EPSILON * maxnorm(prediction, w, n), .0001);
    for (size_t iteration = 0; iteration < 20; iteration++) {
        memcpy(corr, r, n * sizeof(double));
        if (bac_lu_solve(&b->factor, corr)) {
            fail(b, "linear solve failed");
            goto done;
        }
        double scale = 2. / (1. + cj / b->matrix_cj);
        for (size_t i = 0; i < n; i++)
            corr[i] *= scale;
        b->statistics.linear_solves++;
        b->statistics.newton_iterations++;
        double v = maxnorm(corr, w, n);
        if (!isfinite(v)) {
            fail(b, "Newton correction is nonfinite");
            goto done;
        }
        double old = rmsnorm(r, sc, n), lambda = 1.;
        int failed = 0;
        for (;;) {
            for (size_t i = 0; i < n; i++) {
                ny[i] = y[i] - lambda * corr[i];
                np[i] = yp[i] - lambda * cj * corr[i];
            }
            int valid = !evaluate(b, s, t, ny, np, NAN, nr);
            double nv = valid ? rmsnorm(nr, sc, n) : INFINITY;
            if (valid && (nv <= old * 1.1 || v * lambda <= .0001))
                break;
            lambda *= .5;
            if (lambda < 1. / 64.) {
                failed = 1;
                break;
            }
        }
        double scaled = lambda * v, rate = mit ? scaled / MAX(previous, DBL_MIN) : 0.;
        int at_roundoff = scaled <= roundoff;
        if (failed || (rate >= .95 && !at_roundoff) || mit >= 9) {
            if (refreshes >= 2) {
                fail(b, "Newton convergence failed after refreshing the Jacobian");
                goto done;
            }
            if (refresh(b, s, t, h, cj, y, yp, r, w))
                goto done;
            bac_band_row_norms(&b->matrix, w, sc);
            previous = INFINITY;
            mit = 0;
            refreshes++;
            continue;
        }
        memcpy(y, ny, n * sizeof(double));
        memcpy(yp, np, n * sizeof(double));
        memcpy(r, nr, n * sizeof(double));
        double factor = rate / MAX(1. - rate, .05);
        if (at_roundoff || (mit > 0 && scaled * factor <= .0001)) {
            if (refresh(b, s, t, h, cj, y, yp, r, w))
                goto done;
            memcpy(corr, r, n * sizeof(double));
            if (bac_lu_solve(&b->factor, corr)) {
                fail(b, "linear solve failed");
                goto done;
            }
            b->statistics.linear_solves++;
            if (maxnorm(corr, w, n) <= .0001) {
                for (size_t i = 0; i < n; i++) {
                    ny[i] = y[i] - corr[i];
                    np[i] = yp[i] - cj * corr[i];
                }
                if (evaluate(b, s, t, ny, np, NAN, nr))
                    goto done;
                memcpy(y, ny, n * sizeof(double));
                memcpy(yp, np, n * sizeof(double));
                for (size_t i = 0; i < n; i++)
                    err[i] = ny[i] - prediction[i];
                result = 0;
                goto done;
            }
            bac_band_row_norms(&b->matrix, w, sc);
            previous = INFINITY;
            mit = 0;
            continue;
        }
        previous = scaled;
        mit++;
    }
    fail(b, "Newton iteration limit reached");
done:
    if (result && !source)
        bi_invalidate_jacobian(b);
    return result;
}
static void restore_attempt(BdfIntegrator *b, double *phi, double psi[7], const double beta[7],
                            size_t ns, size_t order, double h, double reduction) {
    for (size_t j = ns; j <= order; j++)
        for (size_t i = 0; i < b->n; i++)
            phi[j * b->n + i] *= 1. / beta[j];
    for (size_t j = 1; j <= order; j++)
        psi[j - 1] = psi[j] - h;
    if (!b->last_order) {
        psi[0] = h * reduction;
        for (size_t i = 0; i < b->n; i++)
            phi[b->n + i] *= reduction;
    }
    memcpy(b->phi, phi, 7 * b->n * sizeof(double));
    memcpy(b->psi, psi, 7 * sizeof(double));
    b->stable_steps = ns;
}
int bi_step(BdfIntegrator *b, const DaeSystem *s, double bound, int clip, StepOutcome *out) {
    size_t n = b->n;
    if (!s || !s->residual || s->n != n)
        return fail(b, "DAE system dimension does not match integrator state");
    if (!isfinite(bound) || bound == b->t)
        return fail(b, "step bound must be finite and different from the current time");
    if (b->statistics.accepted_steps >= b->options.max_steps)
        return fail(b, "maximum number of accepted steps exceeded");
    int source = b->options.method == BI_DASKR;
    double direction = signum(bound - b->t);
    if (b->initialized && signum(b->h) != direction) {
        double *y = copy(b->y, n), *yp = copy(b->yp, n);
        if (!y || !yp) {
            free(y);
            free(yp);
            return fail(b, "out of memory");
        }
        int status = bi_reset(b, b->t, y, yp, n);
        free(y);
        free(yp);
        if (status)
            return -1;
    }
    double *work = workspace(&b->step_work, &b->step_work_len, 14 * n);
    if (!work)
        return fail(b, "out of memory");
    double *phi = work, *prediction = work + 7 * n, *pd = work + 8 * n, *w = work + 9 * n,
           *y = work + 10 * n, *yp = work + 11 * n, *error = work + 12 * n, *delta = work + 13 * n;
    int result = -1;
    if (!b->initialized) {
        if (weights(b, b->y, b->y, w))
            goto done;
        double dn = norm(b, b->yp, w), startup = source ? .001 * fabs(bound - b->t) : .001,
               def = dn > .5 / startup ? .5 / dn : startup;
        double mag = b->options.initial_step > 0. ? b->options.initial_step : def;
        b->h = direction * MIN(MIN(mag, b->options.max_step), clip ? fabs(bound - b->t) : INFINITY);
        if (source) {
            b->has_initial_root = 1;
            b->initial_root_time = b->t;
            b->initial_root_h = b->h;
            memcpy(b->initial_root_y, b->y, n * sizeof(double));
            memcpy(b->initial_root_yp, b->yp, n * sizeof(double));
            memcpy(b->source_states, b->y, n * sizeof(double));
            b->source_count = 1;
        }
        memcpy(PH(b, 0), b->y, n * sizeof(double));
        for (size_t i = 0; i < n; i++)
            PH(b, 1)[i] = b->h * b->yp[i];
        b->psi[0] = b->h;
        b->initialized = 1;
    }
    double h = direction *
               MIN(MIN(fabs(b->h), b->options.max_step), clip ? fabs(bound - b->t) : INFINITY),
           roundmin =
               source ? 4. * UROUND * MAX(fabs(b->t), fabs(bound)) : 16. * DBL_EPSILON * fabs(b->t),
           hmin = MAX(b->options.min_step, roundmin);
    size_t order = b->order, ef = 0, nf = 0;
    if (clip && !source && h == bound - b->t && (fabs(h) < roundmin || b->t + h == b->t)) {
        *out = (StepOutcome){b->t, bound, 0., b->last_order, 0.};
        b->t = bound;
        result = 0;
        goto done;
    }
    for (size_t attempt = 0; attempt < 25; attempt++) {
        if (fabs(h) < hmin || b->t + h == b->t) {
            fail(b, "step size underflow");
            goto done;
        }
        memcpy(phi, b->phi, 7 * n * sizeof(double));
        double psi[7], alpha[7], beta[7], gamma[7], sigma[7];
        memcpy(psi, b->psi, sizeof psi);
        memcpy(alpha, b->alpha, sizeof alpha);
        memcpy(beta, b->beta, sizeof beta);
        memcpy(gamma, b->gamma, sizeof gamma);
        memcpy(sigma, b->sigma, sizeof sigma);
        if (!b->last_order) {
            psi[0] = h;
            if (!source)
                for (size_t i = 0; i < n; i++)
                    phi[n + i] = h * b->yp[i];
        }
        size_t ns = (h != b->hold || order != b->last_order) ? 0 : b->stable_steps;
        ns = MIN(ns + 1, b->last_order + 2);
        if (!source || order + 1 >= ns) {
            alpha[0] = 1.;
            beta[0] = 1.;
            gamma[0] = 0.;
            sigma[0] = 1.;
            double temp = h;
            for (size_t j = 1; j <= order; j++) {
                double old = psi[j - 1];
                psi[j - 1] = temp;
                beta[j] = beta[j - 1] * psi[j - 1] / old;
                temp = old + h;
                alpha[j] = h / temp;
                sigma[j] = (double)j * sigma[j - 1] * alpha[j];
                gamma[j] = gamma[j - 1] + alpha[j - 1] / h;
            }
            psi[order] = temp;
        }
        double asum = 0., azero = 0.;
        for (size_t j = 1; j <= order; j++) {
            asum -= 1. / (double)j;
            azero -= alpha[j - 1];
        }
        double cj = -asum / h, ck = MAX(fabs(alpha[order] + asum - azero), alpha[order]);
        for (size_t j = ns; j <= order; j++)
            for (size_t i = 0; i < n; i++)
                phi[j * n + i] *= beta[j];
        memset(prediction, 0, n * sizeof(double));
        memset(pd, 0, n * sizeof(double));
        if (source)
            memcpy(prediction, phi, n * sizeof(double));
        for (size_t j = source ? 1 : 0; j <= order; j++)
            for (size_t i = 0; i < n; i++) {
                prediction[i] += phi[j * n + i];
                pd[i] += gamma[j] * phi[j * n + i];
            }
        double nt = (!source && h == bound - b->t) ? bound : b->t + h;
        if (weights(b, source ? b->phi : prediction, source ? b->phi : b->y, w))
            goto done;
        if (newton(b, s, nt, h, cj, prediction, pd, w, y, yp, error)) {
            nf++;
            b->statistics.rejected_steps++;
            b->statistics.convergence_failures++;
            b->startup = 0;
            if (nf >= 10) {
                fail(b, "Newton convergence failed repeatedly");
                goto done;
            }
            if (source) {
                restore_attempt(b, phi, psi, beta, ns, order, h, .25);
                memcpy(b->alpha, alpha, sizeof alpha);
                memcpy(b->beta, beta, sizeof beta);
                memcpy(b->gamma, gamma, sizeof gamma);
                memcpy(b->sigma, sigma, sizeof sigma);
            }
            h *= .25;
            if (nf >= 3 && !source)
                order = 1;
            continue;
        }
        if (!source && weights(b, y, b->y, w))
            goto done;
        double enorm = norm(b, error, w), local = ck * enorm, erk = sigma[order] * enorm,
               terk = (order + 1) * erk, estimated = erk, erkm1 = 0., terkm1 = INFINITY;
        size_t next = order;
        if (order > 1) {
            for (size_t i = 0; i < n; i++)
                delta[i] = phi[order * n + i] + error[i];
            erkm1 = sigma[order - 1] * norm(b, delta, w);
            terkm1 = order * erkm1;
            int lower;
            if (order == 2)
                lower = terkm1 <= .5 * terk;
            else {
                for (size_t i = 0; i < n; i++)
                    delta[i] += phi[(order - 1) * n + i];
                double erkm2 = sigma[order - 2] * norm(b, delta, w);
                lower = MAX(terkm1, (order - 1) * erkm2) <= terk;
            }
            if (lower) {
                next = order - 1;
                estimated = erkm1;
            }
        }
        if (!isfinite(local) || local > 1.) {
            ef++;
            b->statistics.rejected_steps++;
            b->statistics.error_test_failures++;
            b->startup = 0;
            size_t used = order;
            double reduction;
            if (ef == 1) {
                order = next;
                reduction =
                    clamp(.9 * pow(2. * estimated + .0001, -1. / (double)(order + 1)), .25, .9);
            } else {
                order = ef >= 3 ? 1 : next;
                reduction = .25;
            }
            if (source) {
                restore_attempt(b, phi, psi, beta, ns, used, h, reduction);
                memcpy(b->alpha, alpha, sizeof alpha);
                memcpy(b->beta, beta, sizeof beta);
                memcpy(b->gamma, gamma, sizeof gamma);
                memcpy(b->sigma, sigma, sizeof sigma);
            }
            h *= reduction;
            continue;
        }
        if (next == order - 1 || order == b->options.max_order)
            b->startup = 0;
        double next_h;
        if (b->startup) {
            next = MIN(order + 1, b->options.max_order);
            next_h = 2. * h;
        } else {
            if (next != order - 1 && order < b->options.max_order && order + 1 < ns &&
                (long)order - (long)b->last_order != 1) {
                for (size_t i = 0; i < n; i++)
                    delta[i] = error[i] - phi[(order + 1) * n + i];
                double upper = norm(b, delta, w),
                       erkp1 = source ? (1. / (double)(order + 2)) * upper
                                      : upper / (double)(order + 2),
                       terkp1 = (order + 2) * erkp1;
                if (order == 1) {
                    if (terkp1 < .5 * terk) {
                        next = order + 1;
                        estimated = erkp1;
                    }
                } else if (terkm1 <= MIN(terk, terkp1)) {
                    next = order - 1;
                    estimated = erkm1;
                } else if (terkp1 < terk) {
                    next = order + 1;
                    estimated = erkp1;
                }
            }
            double factor = pow(2. * estimated + .0001, -1. / (double)(next + 1));
            next_h = factor >= 2. ? 2. * h : factor <= 1. ? h * clamp(factor, .5, .9) : h;
        }
        if (order < b->options.max_order)
            memcpy(phi + (order + 1) * n, error, n * sizeof(double));
        for (size_t i = 0; i < n; i++)
            phi[order * n + i] += error[i];
        for (size_t j = order; j-- > 0;)
            for (size_t i = 0; i < n; i++)
                phi[j * n + i] += phi[(j + 1) * n + i];
        double old = b->t;
        b->last_start = old;
        b->has_last_start = 1;
        b->published = 0;
        b->t = nt;
        memcpy(b->y, y, n * sizeof(double));
        memcpy(b->yp, yp, n * sizeof(double));
        memcpy(b->phi, phi, 7 * n * sizeof(double));
        memcpy(b->psi, psi, sizeof psi);
        memcpy(b->alpha, alpha, sizeof alpha);
        memcpy(b->beta, beta, sizeof beta);
        memcpy(b->gamma, gamma, sizeof gamma);
        memcpy(b->sigma, sigma, sizeof sigma);
        if (source) {
            size_t count = MIN(b->source_count, 5);
            memmove(b->source_states + n, b->source_states, count * n * sizeof(double));
            memcpy(b->source_states, b->phi, n * sizeof(double));
            b->source_count = count + 1;
        }
        b->hold = h;
        b->h = direction * MIN(fabs(next_h), b->options.max_step);
        b->stable_steps = ns;
        b->last_order = order;
        b->order = next;
        b->statistics.accepted_steps++;
        b->statistics.maximum_order = MAX(b->statistics.maximum_order, order);
        *out = (StepOutcome){old, nt, h, order, local};
        b->error[0] = '\0';
        result = 0;
        goto done;
    }
    fail(b, "too many rejected steps");
done:
    return result;
}
int bi_interpolate(const BdfIntegrator *b, double t, double *y, double *yp) {
    if (!b->has_last_start)
        return -1;
    double tol = root_tol(b->last_start, b->t), lo = MIN(b->last_start, b->t),
           hi = MAX(b->last_start, b->t);
    if (!isfinite(t) || t < lo - tol || t > hi + tol)
        return -1;
    t = clamp(t, lo, hi);
    double distance = t - b->t, c = 1., d = 0., gamma = distance / b->psi[0];
    memcpy(y, b->phi, b->n * sizeof(double));
    memset(yp, 0, b->n * sizeof(double));
    for (size_t j = 1; j <= b->last_order; j++) {
        d = d * gamma + c / b->psi[j - 1];
        c *= gamma;
        gamma = (distance + b->psi[j - 1]) / b->psi[j];
        for (size_t i = 0; i < b->n; i++) {
            y[i] += c * PH(b, j)[i];
            yp[i] += d * PH(b, j)[i];
        }
    }
    return 0;
}
int bi_publish_dense_output(BdfIntegrator *b, double t) {
    if (bi_interpolate(b, t, b->published_y, b->published_yp))
        return fail(b, "interpolation time is outside the last accepted step");
    b->published_time = t;
    b->published = 1;
    return 0;
}
static int initial_matrix(BdfIntegrator *b, const DaeSystem *s, const double *y, const double *yp,
                          const double *r, const unsigned char *d, BacBand *m) {
    size_t n = b->n, l = s->banded ? MIN(s->lower, n - 1) : n - 1,
           u = s->banded ? MIN(s->upper, n - 1) : n - 1, groups = MIN(l + u + 1, n);
    if (bac_band_new(m, n, l, u))
        return fail(b, "out of memory");
    double *a = zeros(4 * n);
    if (!a) {
        bac_band_free(m);
        return fail(b, "out of memory");
    }
    double *wy = a, *wp = a + n, *rr = a + 2 * n, *inc = a + 3 * n;
    memcpy(wy, y, n * sizeof(double));
    memcpy(wp, yp, n * sizeof(double));
    for (size_t g = 0; g < groups; g++) {
        for (size_t c = g; c < n; c += groups) {
            double v = d[c] ? yp[c] : y[c], delta = sqrt(DBL_EPSILON) * MAX(fabs(v), 1.);
            inc[c] = (v + delta) - v;
            if (d[c])
                wp[c] = v + inc[c];
            else
                wy[c] = v + inc[c];
        }
        if (evaluate(b, s, b->t, wy, wp, NAN, rr)) {
            free(a);
            bac_band_free(m);
            return -1;
        }
        for (size_t c = g; c < n; c += groups) {
            for (size_t row = c > u ? c - u : 0; row < MIN(n, c + l + 1); row++)
                bac_band_set(m, row, c, (rr[row] - r[row]) / inc[c]);
            wy[c] = y[c];
            wp[c] = yp[c];
        }
    }
    b->statistics.jacobian_evaluations++;
    free(a);
    return 0;
}
int bi_initialize_consistent(BdfIntegrator *b, const DaeSystem *s, const unsigned char *d) {
    size_t n = b->n;
    if (s->n != n || !d)
        return fail(b, "consistent initialization dimension mismatch");
    double *a = zeros(6 * n);
    if (!a)
        return fail(b, "out of memory");
    double *y = a, *yp = a + n, *r = a + 2 * n, *cor = a + 3 * n, *ahead = a + 4 * n,
           *tangent = a + 5 * n;
    memcpy(y, b->y, n * sizeof(double));
    memcpy(yp, b->yp, n * sizeof(double));
    int converged = 0, result = -1;
    BacBand m = {0};
    BacLu f = {0};
    for (size_t it = 0; it < 15; it++) {
        if (evaluate(b, s, b->t, y, yp, NAN, r) || initial_matrix(b, s, y, yp, r, d, &m))
            goto done;
        if (bac_band_factor(&m, &f)) {
            fail(b, "singular initial Jacobian");
            goto done;
        }
        memcpy(cor, r, n * sizeof(double));
        if (bac_lu_solve(&f, cor)) {
            fail(b, "linear solve failed");
            goto done;
        }
        bac_band_free(&m);
        bac_lu_free(&f);
        b->statistics.linear_solves++;
        b->statistics.newton_iterations++;
        double v = 0.;
        for (size_t i = 0; i < n; i++) {
            v = MAX(v, fabs(cor[i]) / (1. + fabs(d[i] ? yp[i] : y[i])));
            if (d[i])
                yp[i] -= cor[i];
            else
                y[i] -= cor[i];
        }
        if (v <= 1e-10) {
            converged = 1;
            break;
        }
    }
    if (!converged) {
        fail(b, "consistent initialization Newton iteration did not converge");
        goto done;
    }
    if (evaluate(b, s, b->t, y, yp, NAN, r) || initial_matrix(b, s, y, yp, r, d, &m))
        goto done;
    double ds = 1.;
    for (size_t i = 0; i < n; i++)
        ds = MAX(ds, fabs(yp[i]));
    double eps = sqrt(DBL_EPSILON) * (1. + fabs(b->t)) / ds;
    for (size_t i = 0; i < n; i++)
        ahead[i] = y[i] + eps * yp[i];
    if (evaluate(b, s, b->t + eps, ahead, yp, NAN, tangent))
        goto done;
    for (size_t i = 0; i < n; i++)
        tangent[i] = (tangent[i] - r[i]) / eps;
    if (bac_band_factor(&m, &f) || bac_lu_solve(&f, tangent)) {
        fail(b, "singular initial Jacobian");
        goto done;
    }
    b->statistics.linear_solves++;
    for (size_t i = 0; i < n; i++)
        if (!d[i])
            yp[i] -= tangent[i];
    result = bi_reset(b, b->t, y, yp, n);
done:
    bac_band_free(&m);
    bac_lu_free(&f);
    free(a);
    return result;
}
int bi_restore_source_spare(BdfIntegrator *b, const double *difference, size_t n) {
    if (b->options.method != BI_DASKR || !b->last_order)
        return fail(b, "source spare history requires accepted DASKR history");
    if (n != b->n || !finite_vector(difference, n))
        return fail(b, "invalid source spare history");
    memcpy(PH(b, b->last_order + 1), difference, n * sizeof(double));
    return 0;
}
int bi_project_history(BdfIntegrator *b, const BdfIntegrator *s, size_t n, const TimeOptions *o,
                       BiProjection project, void *ctx) {
    if (!s->initialized || !s->last_order)
        return 0;
    if (!project || bi_validate_options(o, n, b->error, sizeof b->error))
        return -1;
    double *y = zeros(n), *yp = zeros(n);
    if (!y || !yp) {
        free(y);
        free(yp);
        return fail(b, "out of memory");
    }
    BdfIntegrator q = {0};
    int result = -1;
    if (s->options.method == BI_DASKR && o->method == BI_DASKR) {
        size_t count = s->last_order + 1;
        if (s->source_count < count) {
            fail(b, "DASKR warm transfer is missing accepted solution history");
            goto done;
        }
        double *states = zeros(7 * n);
        if (!states) {
            fail(b, "out of memory");
            goto done;
        }
        for (size_t j = 0; j < count; j++)
            if (project(ctx, s->source_states + j * s->n, s->n, states + j * n, n) ||
                !finite_vector(states + j * n, n)) {
                free(states);
                fail(b, "invalid projected DASKR states");
                goto done;
            }
        /* DIVDIF: ordinary divided differences, rescaled into PHI. */
        double times[7] = {0.};
        for (size_t j = 1; j < count; j++)
            times[j] = -s->psi[j - 1];
        double *dd = copy(states, 7 * n);
        if (!dd) {
            free(states);
            fail(b, "out of memory");
            goto done;
        }
        for (size_t k = 1; k < count; k++)
            for (size_t j = count; j-- > k;)
                for (size_t i = 0; i < n; i++)
                    dd[j * n + i] =
                        (dd[(j - 1) * n + i] - dd[j * n + i]) / (times[j - k] - times[j]);
        double prod = 1.;
        memcpy(y, states, n * sizeof(double));
        for (size_t j = 1; j < count; j++) {
            prod *= s->psi[j - 1];
            for (size_t i = 0; i < n; i++)
                dd[j * n + i] *= prod;
        }
        double c = 1., d = 0., g = 0.;
        for (size_t j = 1; j <= s->last_order; j++) {
            d = d * g + c / s->psi[j - 1];
            c *= g;
            g = s->psi[j - 1] / s->psi[j];
            for (size_t i = 0; i < n; i++)
                yp[i] += d * dd[j * n + i];
        }
        if (bi_new(&q, s->t, y, yp, n, o)) {
            free(dd);
            free(states);
            fail(b, q.error);
            goto done;
        }
        memcpy(q.phi, dd, 7 * n * sizeof(double));
        if (s->last_order == 1 && n == b->n)
            memcpy(PH(&q, 2), PH(b, 2), n * sizeof(double));
        memcpy(q.source_states, states, count * n * sizeof(double));
        q.source_count = count;
        free(dd);
        free(states);
        q.stable_steps = s->stable_steps;
        q.startup = s->startup;
    } else {
        if (project(ctx, s->y, s->n, y, n) || project(ctx, s->yp, s->n, yp, n) ||
            bi_new(&q, s->t, y, yp, n, o)) {
            fail(b, "invalid projected BDF state");
            goto done;
        }
        for (size_t j = 0; j < 7; j++)
            if (project(ctx, PH(s, j), s->n, PH(&q, j), n) || !finite_vector(PH(&q, j), n)) {
                fail(b, "invalid projected BDF history");
                goto done;
            }
        q.stable_steps = 0;
        q.startup = 0;
    }
    q.statistics = b->statistics;
    q.order = MIN(MIN(s->order, s->last_order), o->max_order);
    q.h = signum(s->h) * MIN(MIN(fabs(s->h), fabs(s->hold)), o->max_step);
    memcpy(q.psi, s->psi, sizeof q.psi);
    memcpy(q.alpha, s->alpha, sizeof q.alpha);
    memcpy(q.beta, s->beta, sizeof q.beta);
    memcpy(q.gamma, s->gamma, sizeof q.gamma);
    memcpy(q.sigma, s->sigma, sizeof q.sigma);
    q.previous_cj = s->previous_cj;
    q.newton_convergence_factor = s->newton_convergence_factor;
    q.hold = s->hold;
    q.last_order = s->last_order;
    q.initialized = 1;
    q.has_last_event = s->has_last_event;
    q.last_event_time = s->last_event_time;
    q.last_event_count = s->last_event_count;
    if (q.last_event_count) {
        q.last_event_indices = malloc(q.last_event_count * sizeof(size_t));
        if (!q.last_event_indices) {
            fail(b, "out of memory");
            goto done;
        }
        memcpy(q.last_event_indices, s->last_event_indices, q.last_event_count * sizeof(size_t));
    }
    move_workspaces(&q, b);
    bi_free(b);
    *b = q;
    memset(&q, 0, sizeof q);
    result = 1;
done:
    bi_free(&q);
    free(y);
    free(yp);
    return result;
}
void bi_root_free(RootEvent *r) {
    if (!r)
        return;
    free(r->state);
    free(r->derivative);
    free(r->indices);
    free(r->directions);
    memset(r, 0, sizeof *r);
}
static int root_eval(BdfIntegrator *b, BiRootFunction fn, void *ctx, double t, size_t nr,
                     double *values, double *y, double *yp, int source) {
    if (source && (t - b->t) * b->psi[0] > 0.) {
        if (bi_interpolate(b, b->t, y, yp))
            return fail(b, "root interpolation failed");
        double factor = (t - b->t) / b->psi[0];
        for (size_t i = 0; i < b->n; i++)
            y[i] += factor * PH(b, 1)[i];
    } else if (bi_interpolate(b, t, y, yp))
        return fail(b, "root interpolation failed");
    for (size_t i = 0; i < nr; i++)
        values[i] = NAN;
    if (fn(ctx, t, y, yp, values, nr) || !finite_vector(values, nr))
        return fail(b, "root function failed or returned a nonfinite value");
    return 0;
}
static int remember_event(BdfIntegrator *b, RootEvent *event, double t, const size_t *indices,
                          const int *dirs, size_t count) {
    RootEvent r = {0};
    r.time = t;
    r.count = count;
    r.state = zeros(b->n);
    r.derivative = zeros(b->n);
    r.indices = malloc(count * sizeof(size_t));
    r.directions = malloc(count * sizeof(int));
    size_t *saved = malloc(count * sizeof(size_t));
    if (!r.state || !r.derivative || (count && (!r.indices || !r.directions || !saved))) {
        bi_root_free(&r);
        free(saved);
        return fail(b, "out of memory");
    }
    if (bi_interpolate(b, t, r.state, r.derivative)) {
        bi_root_free(&r);
        free(saved);
        return fail(b, "root interpolation failed");
    }
    if (count) {
        memcpy(r.indices, indices, count * sizeof(size_t));
        memcpy(r.directions, dirs, count * sizeof(int));
        memcpy(saved, indices, count * sizeof(size_t));
    }
    free(b->last_event_indices);
    b->last_event_indices = saved;
    b->last_event_count = count;
    b->has_last_event = 1;
    b->last_event_time = t;
    *event = r;
    return 1;
}
static int select_root(const double *a, const double *z, size_t n, size_t *selected, int *zero) {
    int found = 0;
    double maximum = 0.;
    *zero = 0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(z[i]) == 0.) {
            *zero = 1;
            continue;
        }
        if (signum(a[i]) == signum(z[i]))
            continue;
        double fraction = fabs(z[i] / (z[i] - a[i]));
        if (fraction > maximum) {
            maximum = fraction;
            *selected = i;
            found = 1;
        }
    }
    return found;
}
static int daskr_roots(BdfIntegrator *b, BiRootFunction fn, void *ctx, size_t nr,
                       const int *directions, double output, RootEvent *event) {
    double h = b->psi[0], forward = signum(h), resolution = (fabs(b->t) + fabs(h)) * UROUND * 100.,
           left = b->has_root_cursor ? b->root_cursor : b->last_start,
           right = (output - b->t) * h >= 0. ? b->t : output;
    if ((left - b->last_start) * forward < 0.)
        left = b->last_start;
    if ((right - left) * h <= 0.)
        return 0;
    double *a = zeros(5 * nr + 2 * b->n);
    size_t *indices = malloc((nr ? nr : 1) * sizeof(size_t));
    int *dirs = malloc((nr ? nr : 1) * sizeof(int));
    unsigned char *enabled = calloc(nr ? nr : 1, 1);
    if (!a || !indices || !dirs || !enabled) {
        free(a);
        free(indices);
        free(dirs);
        free(enabled);
        return fail(b, "out of memory");
    }
    double *vl = a, *vr = a + nr, *vm = a + 2 * nr, *initial = a + 3 * nr, *probe = a + 4 * nr,
           *y = a + 5 * nr, *yp = y + b->n;
    int result = -1;
    if (root_eval(b, fn, ctx, left, nr, vl, y, yp, 1))
        goto done;
    int continuation = b->has_last_event && b->last_event_time == left;
    if (!b->has_root_cursor && b->has_initial_root && b->initial_root_time == left) {
        if (fn(ctx, left, b->initial_root_y, b->initial_root_yp, initial, nr) ||
            !finite_vector(initial, nr)) {
            fail(b, "initial root function failed");
            goto done;
        }
        int zero = 0;
        for (size_t i = 0; i < nr; i++)
            zero |= initial[i] == 0.;
        if (zero) {
            double res = (fabs(left) + fabs(b->initial_root_h)) * UROUND * 100.,
                   factor = MAX(res / fabs(b->initial_root_h), .1),
                   time = left + factor * b->initial_root_h;
            for (size_t i = 0; i < b->n; i++)
                y[i] = b->initial_root_y[i] + factor * (b->initial_root_h * b->initial_root_yp[i]);
            if (fn(ctx, time, y, b->initial_root_yp, probe, nr) || !finite_vector(probe, nr)) {
                fail(b, "initial root probe failed");
                goto done;
            }
            for (size_t i = 0; i < nr; i++)
                if (probe[i] == 0.) {
                    fail(b, "DASKR root function is zero at initial time and nearby probe");
                    goto done;
                }
            left = time;
            if (root_eval(b, fn, ctx, left, nr, vl, y, yp, 1))
                goto done;
        }
    } else if (continuation) {
        int zero = 0;
        for (size_t i = 0; i < nr; i++) {
            initial[i] = vl[i] == 0.;
            zero |= (int)initial[i];
        }
        if (zero) {
            double previous = left, increment = copysign(resolution, h);
            left += increment;
            if ((left - b->t) * h >= 0.) {
                if (bi_interpolate(b, previous, y, yp)) {
                    fail(b, "root interpolation failed");
                    goto done;
                }
                for (size_t i = 0; i < b->n; i++)
                    y[i] += (increment / h) * PH(b, 1)[i];
                if (fn(ctx, left, y, yp, vl, nr) || !finite_vector(vl, nr)) {
                    fail(b, "root probe failed");
                    goto done;
                }
            } else if (root_eval(b, fn, ctx, left, nr, vl, y, yp, 1))
                goto done;
            size_t count = 0;
            for (size_t i = 0; i < nr; i++)
                if (vl[i] == 0.) {
                    if (initial[i]) {
                        fail(b, "DASKR root function is zero at returned root and nearby probe");
                        goto done;
                    }
                    indices[count] = i;
                    dirs[count++] = BI_ANY;
                }
            if (count) {
                b->root_cursor = left;
                b->has_root_cursor = 1;
                result = remember_event(b, event, left, indices, dirs, count);
                goto done;
            }
        }
    }
    if ((right - left) * h <= 0.) {
        b->has_root_cursor = 1;
        b->root_cursor = left;
        result = 0;
        goto done;
    }
    if (root_eval(b, fn, ctx, right, nr, vr, y, yp, 1))
        goto done;
    for (size_t i = 0; i < nr; i++) {
        int crossing = (vr[i] - vl[i]) * h > 0. ? BI_RISING : BI_FALLING;
        enabled[i] = !directions || directions[i] == BI_ANY || directions[i] == crossing;
        if (!enabled[i])
            vl[i] = vr[i] = 1.;
    }
    size_t selected = 0;
    int zero = 0, found = select_root(vl, vr, nr, &selected, &zero);
    if (!found && !zero) {
        b->has_root_cursor = 1;
        b->root_cursor = right;
        result = 0;
        goto done;
    }
    if (found) {
        int previous_side = 0, last_side = 1, finished = 0;
        double alpha = 1.;
        for (size_t it = 0; it < 2000; it++) {
            if (previous_side != last_side)
                alpha = 1.;
            else if (last_side == 0)
                alpha *= 2.;
            else
                alpha *= .5;
            double middle =
                right - (right - left) * vr[selected] / (vr[selected] - alpha * vl[selected]);
            if (fabs(middle - left) < .5 * resolution) {
                double fraction = fabs(right - left) / resolution,
                       inset = fraction > 5. ? .1 : .5 / fraction;
                middle = left + inset * (right - left);
            }
            if (fabs(right - middle) < .5 * resolution) {
                double fraction = fabs(right - left) / resolution,
                       inset = fraction > 5. ? .1 : .5 / fraction;
                middle = right - inset * (right - left);
            }
            if (root_eval(b, fn, ctx, middle, nr, vm, y, yp, 1))
                goto done;
            for (size_t i = 0; i < nr; i++)
                if (!enabled[i])
                    vm[i] = 1.;
            size_t next = 0;
            int newzero = 0, newfound = select_root(vl, vm, nr, &next, &newzero);
            previous_side = last_side;
            int exact = 0;
            if (newfound) {
                selected = next;
                right = middle;
                memcpy(vr, vm, nr * sizeof(double));
                last_side = 1;
            } else if (newzero) {
                right = middle;
                memcpy(vr, vm, nr * sizeof(double));
                exact = 1;
            } else {
                left = middle;
                memcpy(vl, vm, nr * sizeof(double));
                last_side = 0;
            }
            if (exact || fabs(right - left) <= resolution) {
                finished = 1;
                break;
            }
        }
        if (!finished) {
            fail(b, "DASKR Illinois root search iteration limit exceeded");
            goto done;
        }
    }
    {
        size_t count = 0;
        for (size_t i = 0; i < nr; i++) {
            double sg = 0.;
            if (fabs(vr[i]) == 0.)
                sg = -signum(vl[i]);
            else if (signum(vl[i]) != signum(vr[i]))
                sg = signum(vr[i] - vl[i]);
            if (sg != 0.) {
                indices[count] = i;
                dirs[count++] = sg * forward > 0. ? BI_RISING : BI_FALLING;
            }
        }
        b->root_cursor = right;
        b->has_root_cursor = 1;
        result = remember_event(b, event, right, indices, dirs, count);
    }
done:
    free(a);
    free(indices);
    free(dirs);
    free(enabled);
    return result;
}
int bi_locate_roots(BdfIntegrator *b, BiRootFunction fn, void *ctx, size_t nr,
                    const int *directions, double output, RootEvent *event) {
    memset(event, 0, sizeof *event);
    if (!b->has_last_start)
        return fail(b, "root location requires an accepted step");
    if (!fn)
        return fail(b, "missing root callback");
    if (b->options.method == BI_DASKR)
        return daskr_roots(b, fn, ctx, nr, directions, output, event);
    size_t n = b->n;
    double *work = zeros(4 * nr + 2 * n);
    size_t *indices = malloc((nr ? nr : 1) * sizeof(size_t));
    int *dirs = malloc((nr ? nr : 1) * sizeof(int));
    if (!work || !indices || !dirs) {
        free(work);
        free(indices);
        free(dirs);
        return fail(b, "out of memory");
    }
    double *vs = work, *ve = work + nr, *vals = work + 2 * nr, *times = work + 3 * nr,
           *y = work + 4 * nr, *yp = y + n;
    int result = -1;
    double start = b->last_start, tol = root_tol(start, b->t);
    if (root_eval(b, fn, ctx, start, nr, vs, y, yp, 0) ||
        root_eval(b, fn, ctx, b->t, nr, ve, y, yp, 0))
        goto done;
    size_t count = 0;
    for (size_t k = 0; k < nr; k++) {
        double a = vs[k], z = ve[k], left = start;
        int restarted = 0;
        if (b->has_last_event &&
            fabs(start - b->last_event_time) <= root_tol(start, b->last_event_time))
            for (size_t j = 0; j < b->last_event_count; j++)
                if (b->last_event_indices[j] == k)
                    restarted = 1;
        if (restarted || a == 0.) {
            double increment = MIN(4. * tol, .1 * fabs(b->t - start));
            left = start + signum(b->t - start) * increment;
            if (left == start || left == b->t)
                continue;
            if (root_eval(b, fn, ctx, left, nr, vals, y, yp, 0))
                goto done;
            a = vals[k];
        }
        if (a == 0. || (z != 0. && signum(a) == signum(z)))
            continue;
        int crossing = (z - a) * (b->t - start) > 0. ? BI_RISING : BI_FALLING;
        if (directions && directions[k] != BI_ANY && directions[k] != crossing)
            continue;
        double right = b->t, lv = a, rv = z;
        if (z != 0.)
            for (size_t it = 0; it < 100; it++) {
                if (fabs(right - left) <= tol)
                    break;
                double fraction = clamp(lv / (lv - rv), .1, .9),
                       middle = left + fraction * (right - left);
                if (middle == left || middle == right)
                    break;
                if (root_eval(b, fn, ctx, middle, nr, vals, y, yp, 0))
                    goto done;
                double v = vals[k];
                if (v == 0.) {
                    left = right = middle;
                    break;
                }
                if (signum(lv) != signum(v)) {
                    right = middle;
                    rv = v;
                } else {
                    left = middle;
                    lv = v;
                }
            }
        indices[count] = k;
        dirs[count] = crossing;
        times[count++] = z == 0. ? b->t : .5 * (left + right);
    }
    if (!count) {
        result = 0;
        goto done;
    }
    double direction = signum(b->t - start), first = times[0];
    for (size_t j = 1; j < count; j++)
        if (direction * (times[j] - start) < direction * (first - start))
            first = times[j];
    size_t kept = 0;
    for (size_t j = 0; j < count; j++)
        if (fabs(times[j] - first) <= 4. * tol) {
            indices[kept] = indices[j];
            dirs[kept++] = dirs[j];
        }
    result = remember_event(b, event, first, indices, dirs, kept);
done:
    free(work);
    free(indices);
    free(dirs);
    return result;
}
