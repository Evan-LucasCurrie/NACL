/* Adaptive MIRK boundary-value solver in C11. */
#include "internal.h"
void bvp_fail(context *c, bvp_status status) {
    c->error = status;
    longjmp(*c->trap, 1);
}
void *bvp_alloc(context *c, size_t n, size_t size) {
    if (size && n > SIZE_MAX / size)
        bvp_fail(c, BVP_INVALID_INPUT);
    size_t bytes = (n ? n : 1) * (size ? size : 1);
    allocation **cached = &c->scratch_cache;
    allocation *a = NULL;
    while (*cached) {
        if ((*cached)->bytes >= bytes && (*cached)->bytes - bytes <= bytes) {
            a = *cached;
            *cached = a->next;
            c->cached_bytes -= a->bytes;
            memset(a->ptr, 0, bytes);
            break;
        }
        cached = &(*cached)->next;
    }
    if (!a) {
        a = malloc(sizeof(*a));
        if (!a)
            bvp_fail(c, BVP_NO_MEMORY);
        a->ptr = calloc(n ? n : 1, size ? size : 1);
        if (!a->ptr) {
            free(a);
            bvp_fail(c, BVP_NO_MEMORY);
        }
        a->bytes = bytes;
    }
    a->next = c->memory;
    c->memory = a;
    return a->ptr;
}
void bvp_release(context *c, allocation *mark) {
    while (c->memory != mark) {
        allocation *a = c->memory;
        c->memory = a->next;
        const size_t cache_limit = 1u << 20;
        if (mark && a->bytes && a->bytes <= cache_limit &&
            c->cached_bytes <= cache_limit - a->bytes) {
            a->next = c->scratch_cache;
            c->scratch_cache = a;
            c->cached_bytes += a->bytes;
        } else {
            free(a->ptr);
            free(a);
        }
    }
    if (!mark) {
        while (c->scratch_cache) {
            allocation *a = c->scratch_cache;
            c->scratch_cache = a->next;
            free(a->ptr);
            free(a);
        }
        c->cached_bytes = 0;
    }
}
void bvp_detach(context *c, void *ptr) {
    allocation **p = &c->memory;
    while (*p) {
        if ((*p)->ptr == ptr) {
            allocation *a = *p;
            *p = a->next;
            free(a);
            return;
        }
        p = &(*p)->next;
    }
}
void bvp_finite(context *c, const double *a, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (!isfinite(a[i]))
            bvp_fail(c, BVP_NONFINITE);
}
double bvp_maximum(const double *a, size_t n) {
    double v = 0;
    for (size_t i = 0; i < n; i++)
        v = MAX(v, fabs(a[i]));
    return v;
}
static int valid_mesh(const double *x, size_t n) {
    if (!x || n < 2)
        return 0;
    for (size_t i = 0; i < n; i++)
        if (!isfinite(x[i]) || (i && (!(x[i] > x[i - 1]) || !isfinite(x[i] - x[i - 1]))))
            return 0;
    return 1;
}
static double dot(size_t rows, size_t count, const double *a, size_t astride, const double *b,
                  size_t bstride) {
    double sum = 0;
    if (rows == 1) {
        double even = 0, odd = 0;
        size_t p = count / 2 * 2;
        for (size_t j = 0; j < p; j += 2) {
            even = fma(a[j * astride], b[j * bstride], even);
            odd = fma(a[(j + 1) * astride], b[(j + 1) * bstride], odd);
        }
        sum = even + odd;
        if (p < count)
            sum = fma(a[p * astride], b[p * bstride], sum);
    } else
        for (size_t j = 0; j < count; j++)
            sum = fma(a[j * astride], b[j * bstride], sum);
    return sum;
}
static void rhs(context *c, double x, const double *y, double *f, double *raw) {
    size_t n = c->n, q = c->q;
    c->stats.rhs_evaluations++;
    memset(f, 0, q * sizeof(double));
    c->p->rhs(c->p->user, x, y, y + n, f);
    if (c->p->singular_term) {
        if (x <= c->a) {
            memcpy(raw, f, n * sizeof(double));
            for (size_t i = 0; i < n; i++)
                f[i] = dot(n, n, c->endpoint_inverse + i * n, 1, raw, 1);
        } else
            for (size_t i = 0; i < n; i++)
                f[i] += dot(n, n, c->p->singular_term + i * n, 1, y, 1) / (x - c->a);
    }
    bvp_finite(c, f, q);
}
static void jacobian(context *c, double x, const double *y, double *jac, double *work) {
    size_t n = c->n, q = c->q, qq = q * q;
    c->stats.jacobian_evaluations++;
    memset(jac, 0, qq * sizeof(double));
    double *raw = work, *pert = raw + qq, *plus = pert + q, *minus = plus + q;
    if (c->p->jacobian && c->p->jacobian(c->p->user, x, y, y + n, jac)) {
        if (c->p->singular_term) {
            if (x <= c->a) {
                memcpy(raw, jac, qq * sizeof(double));
                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j < q; j++) {
                        double sum = 0;
                        for (size_t k = 0; k < n; k++)
                            sum += c->endpoint_inverse[i * n + k] * raw[k * q + j];
                        jac[i * q + j] = sum;
                    }
            } else
                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j < n; j++)
                        jac[i * q + j] += c->p->singular_term[i * n + j] / (x - c->a);
        }
    } else {
        memcpy(pert, y, q * sizeof(double));
        rhs(c, x, y, minus, raw);
        for (size_t j = 0; j < q; j++) {
            double step = sqrt(DBL_EPSILON) * (y[j] != 0 ? fabs(y[j]) : 1);
            pert[j] = y[j] + step;
            rhs(c, x, pert, plus, raw);
            pert[j] = y[j];
            for (size_t i = 0; i < n; i++)
                jac[i * q + j] = (plus[i] - minus[i]) / step;
        }
    }
    bvp_finite(c, jac, qq);
}
static void project(context *c, double *y) {
    if (c->projection) {
        allocation *mark = c->memory;
        double stack[16], *raw = c->n <= 16 ? stack : bvp_alloc(c, c->n, sizeof(double));
        memcpy(raw, y, c->n * sizeof(double));
        for (size_t i = 0; i < c->n; i++)
            y[i] = dot(c->n, c->n, c->projection + i * c->n, 1, raw, 1);
        bvp_release(c, mark);
    }
}
static void boundary(context *c, const double *ya, const double *yb, double *a, double *b) {
    c->p->boundary(c->p->user, ya, yb, ya + c->n, a, b);
    bvp_finite(c, a, c->left);
    bvp_finite(c, b, c->q - c->left);
}
static void boundary_jac(context *c, const double *ya, const double *yb, double *a, double *b) {
    size_t q = c->q, l = c->left;
    memset(a, 0, l * q * sizeof(double));
    memset(b, 0, (q - l) * q * sizeof(double));
    if (!(c->p->boundary_jacobian &&
          c->p->boundary_jacobian(c->p->user, ya, yb, ya + c->n, a, b))) {
        allocation *mark = c->memory;
        double stack[64], *work = q <= 16 ? stack : bvp_alloc(c, 4 * q, sizeof(double));
        double *aa = work, *bb = aa + q, *base = bb + q, *plus = base + q;
        if (q <= 16)
            memset(work, 0, 4 * q * sizeof(double));
        memcpy(aa, ya, q * sizeof(double));
        memcpy(bb, yb, q * sizeof(double));
        boundary(c, ya, yb, base, base + l);
        for (size_t j = 0; j < q; j++) {
            double sa = sqrt(DBL_EPSILON) * (ya[j] != 0 ? fabs(ya[j]) : 1),
                   sb = sqrt(DBL_EPSILON) * (yb[j] != 0 ? fabs(yb[j]) : 1);
            aa[j] = ya[j] + sa;
            bb[j] = yb[j] + sb;
            boundary(c, aa, bb, plus, plus + l);
            aa[j] = ya[j];
            bb[j] = yb[j];
            for (size_t i = 0; i < l; i++)
                a[i * q + j] = (plus[i] - base[i]) / sa;
            for (size_t i = 0; i < q - l; i++)
                b[i * q + j] = (plus[l + i] - base[l + i]) / sb;
        }
        bvp_release(c, mark);
    }
    bvp_finite(c, a, l * q);
    bvp_finite(c, b, (q - l) * q);
}
static double boundary_norm(context *c, const double *ya, const double *yb) {
    allocation *mark = c->memory;
    size_t q = c->q, l = c->left;
    double stack[72], *work = q <= 8 ? stack : bvp_alloc(c, q + q * q, sizeof(double));
    double *bc = work, *jac = bc + q;
    if (q <= 8)
        memset(work, 0, (q + q * q) * sizeof(double));
    boundary(c, ya, yb, bc, bc + l);
    boundary_jac(c, ya, yb, jac, jac + l * q);
    double norm = 0;
    for (size_t i = 0; i < q; i++) {
        const double *y = i < l ? ya : yb;
        double scale = 0;
        for (size_t j = 0; j < q; j++)
            scale += fabs(jac[i * q + j]) * (1 + fabs(y[j]));
        norm = MAX(norm, fabs(bc[i]) / (scale > 0 ? scale : 1));
    }
    bvp_release(c, mark);
    return norm;
}
static void residual(context *c, const double *mesh, size_t nodes, const double *y,
                     const tableau *t, double *phi, double *stages) {
    size_t q = c->q, sc = t->s + t->extra;
    allocation *mark = c->memory;
    double stack[32], *work = q <= 16 ? stack : bvp_alloc(c, 2 * q, sizeof(double));
    double *arg = work, *raw = arg + q;
    if (q <= 16)
        memset(work, 0, 2 * q * sizeof(double));
    for (size_t i = 0; i < nodes - 1; i++) {
        double h = mesh[i + 1] - mesh[i], *ks = stages + i * sc * q;
        memset(ks, 0, t->s * q * sizeof(double));
        for (size_t r = 0; r < t->s; r++) {
            for (size_t j = 0; j < q; j++) {
                double sum = 0;
                for (size_t k = 0; k < t->s; k++)
                    sum = fma(ks[k * q + j], t->x[r][k], sum);
                arg[j] = (1 - t->v[r]) * y[i * q + j] + t->v[r] * y[(i + 1) * q + j] + h * sum;
            }
            rhs(c, mesh[i] + t->c[r] * h, arg, ks + r * q, raw);
        }
        for (size_t j = 0; j < q; j++)
            phi[c->left + i * q + j] =
                y[(i + 1) * q + j] - y[i * q + j] - h * dot(q, t->s, ks + j, q, t->b, 1);
    }
    boundary(c, y, y + (nodes - 1) * q, phi, phi + c->left + (nodes - 1) * q);
    bvp_finite(c, phi, nodes * q);
    bvp_release(c, mark);
}
static void assemble(context *c, const double *mesh, size_t nodes, const double *y,
                     const tableau *t, const double *stages, matrix *m) {
    size_t q = c->q, qq = q * q, sc = t->s + t->extra;
    memset(m->a, 0, m->count * sizeof(double));
    allocation *mark = c->memory;
    size_t scratch_count = (2 * t->s + 6) * qq + 4 * q;
    double stack[104], *scratch = q <= 2 ? stack : bvp_alloc(c, scratch_count, sizeof(double));
    if (q <= 2)
        memset(scratch, 0, scratch_count * sizeof(double));
    double *dkl = scratch, *dkr = dkl + t->s * qq;
    double *arg = dkr + t->s * qq, *jac = arg + q, *endpoint = jac + qq;
    double endpoint_x = NAN;
    double *ml = endpoint + qq, *mr = ml + qq, *work = mr + qq;
    double *bc = work + qq + 3 * q;
    for (size_t i = 0; i < nodes - 1; i++) {
        double h = mesh[i + 1] - mesh[i];
        const double *ks = stages + i * sc * q;
        memset(dkl, 0, t->s * qq * sizeof(double));
        memset(dkr, 0, t->s * qq * sizeof(double));
        int lzero[9], rzero[9];
        for (size_t r = 0; r < 9; r++)
            lzero[r] = rzero[r] = 1;
        for (size_t r = 0; r < t->s; r++) {
            for (size_t j = 0; j < q; j++) {
                double sum = 0;
                for (size_t k = 0; k < r; k++)
                    if (t->x[r][k] != 0)
                        sum += t->x[r][k] * ks[k * q + j];
                arg[j] = (1 - t->v[r]) * y[i * q + j] + t->v[r] * y[(i + 1) * q + j] + h * sum;
            }
            double sx = mesh[i] + t->c[r] * h;
            if (r == 0 && c->p->is_affine && memcmp(&sx, &endpoint_x, sizeof(double)) == 0)
                memcpy(jac, endpoint, qq * sizeof(double));
            else
                jacobian(c, sx, arg, jac, work);
            if (r == 1) {
                memcpy(endpoint, jac, qq * sizeof(double));
                endpoint_x = sx;
            }
            for (size_t side = 0; side < 2; side++) {
                double *dk = side ? dkr : dkl, *mid = side ? mr : ml,
                       prefactor = side ? t->v[r] : 1 - t->v[r];
                int *zero = side ? rzero : lzero;
                memset(mid, 0, qq * sizeof(double));
                for (size_t k = 0; k < r; k++)
                    if (t->x[r][k] != 0 && !zero[k]) {
                        zero[r] = 0;
                        for (size_t e = 0; e < qq; e++)
                            mid[e] += t->x[r][k] * dk[k * qq + e];
                    }
                if (!zero[r])
                    for (size_t e = 0; e < qq; e++)
                        mid[e] *= h;
                if (prefactor != 0 && zero[r]) {
                    for (size_t e = 0; e < qq; e++)
                        dk[r * qq + e] = prefactor * jac[e];
                    zero[r] = 0;
                    continue;
                }
                if (prefactor != 0) {
                    for (size_t j = 0; j < q; j++)
                        mid[j * q + j] += prefactor;
                    zero[r] = 0;
                }
                if (!zero[r])
                    for (size_t row = 0; row < q; row++)
                        for (size_t col = 0; col < q; col++) {
                            double sum = 0;
                            for (size_t k = 0; k < q; k++)
                                sum += jac[row * q + k] * mid[k * q + col];
                            dk[r * qq + row * q + col] = sum;
                        }
            }
        }
        for (size_t row = 0; row < q; row++) {
            double *entries = bvp_matrix_row(m, c->left + i * q + row, i * q);
            for (size_t col = 0; col < q; col++) {
                double dl = 0, dr = 0;
                for (size_t r = 0; r < t->s; r++) {
                    dl -= t->b[r] * dkl[r * qq + row * q + col];
                    dr -= t->b[r] * dkr[r * qq + row * q + col];
                }
                dl *= h;
                dr *= h;
                if (row == col) {
                    dl += -1;
                    dr += 1;
                }
                entries[col] = dl;
                entries[q + col] = dr;
            }
        }
    }
    boundary_jac(c, y, y + (nodes - 1) * q, bc, bc + c->left * q);
    for (size_t row = 0; row < c->left; row++)
        memcpy(bvp_matrix_row(m, row, 0), bc + row * q, q * sizeof(double));
    for (size_t row = 0; row < q - c->left; row++)
        memcpy(bvp_matrix_row(m, c->left + (nodes - 1) * q + row, (nodes - 1) * q),
               bc + (c->left + row) * q, q * sizeof(double));
    bvp_release(c, mark);
}
static double criterion(context *c, const double *d, size_t count) {
    double norm = bvp_maximum(d, count);
    if (!isfinite(norm) || norm > sqrt(DBL_MAX))
        bvp_fail(c, BVP_NEWTON_FAILED);
    return 0.5 * (norm * norm);
}
static void offset_add(double *phi, const double *offset, size_t n) {
    if (offset)
        for (size_t i = 0; i < n; i++)
            phi[i] += offset[i];
}
static double *newton_inner(context *c, const double *mesh, size_t nodes, double *y,
                            const tableau *t, const bvp_options *o, const double *offset) {
    size_t q = c->q, count = nodes * q, sc = t->s + t->extra, stage_len = (nodes - 1) * q * sc;
    double *result = bvp_alloc(c, stage_len, sizeof(double));
    allocation *mark = c->memory;
    double *stages = bvp_alloc(c, stage_len, sizeof(double)),
           *trial_stages = bvp_alloc(c, stage_len, sizeof(double)),
           *phi = bvp_alloc(c, count, sizeof(double)),
           *trial_phi = bvp_alloc(c, count, sizeof(double)),
           *delta = bvp_alloc(c, count, sizeof(double)),
           *natural = bvp_alloc(c, count, sizeof(double)),
           *trial = bvp_alloc(c, count, sizeof(double));
    residual(c, mesh, nodes, y, t, phi, stages);
    offset_add(phi, offset, count);
    matrix factors = bvp_matrix(c, q, c->left, nodes);
    int fixed = 0;
    size_t full_steps = 0;
    double lambda = 0, previous_norm = 0, crit = 0, tol = o->tolerance / 100;
    for (size_t iteration = 0; iteration < o->max_newton_iterations * 100; iteration++) {
        double initial = crit;
        if (!fixed) {
            if (full_steps > o->max_newton_iterations)
                break;
            assemble(c, mesh, nodes, y, t, stages, &factors);
            bvp_factor(c, &factors);
            memcpy(delta, phi, count * sizeof(double));
            bvp_linear(c, &factors, delta, 0);
            full_steps++;
            initial = criterion(c, delta, count);
        }
        c->stats.newton_iterations++;
        int small = 1;
        for (size_t j = 0; j < count; j++)
            if (!(fabs(delta[j]) <= tol * (fabs(y[j]) + 1))) {
                small = 0;
                break;
            }
        if (small) {
            if (c->left == q && fixed && !c->p->is_affine) {
                fixed = 0;
                lambda = 0;
                continue;
            }
            for (size_t j = 0; j < count; j++)
                y[j] -= delta[j];
            residual(c, mesh, nodes, y, t, phi, stages);
            if (boundary_norm(c, y, y + (nodes - 1) * q) <= tol) {
                memcpy(result, stages, stage_len * sizeof(double));
                bvp_release(c, mark);
                return result;
            }
            offset_add(phi, offset, count);
            fixed = 0;
            lambda = 0;
            continue;
        }
        if (fixed) {
            for (size_t j = 0; j < count; j++)
                trial[j] = y[j] - delta[j];
            residual(c, mesh, nodes, trial, t, trial_phi, trial_stages);
            offset_add(trial_phi, offset, count);
            memcpy(natural, trial_phi, count * sizeof(double));
            bvp_linear(c, &factors, natural, 0);
            double tc = criterion(c, natural, count);
            int contracting = tc <= 0.5 * crit;
            if (!contracting) {
                fixed = 0;
                lambda = 0;
            }
            if (contracting || tc < crit) {
                memcpy(y, trial, count * sizeof(double));
                SWAP(double *, phi, trial_phi);
                SWAP(double *, stages, trial_stages);
                SWAP(double *, delta, natural);
                if (contracting)
                    crit = tc;
            }
            continue;
        }
        if (lambda == 0)
            lambda = 1;
        else {
            double change = 0;
            for (size_t j = 0; j < count; j++)
                change = MAX(change, fabs(natural[j] - delta[j]));
            lambda = CLAMP(lambda * previous_norm / change, 0.01, 1);
        }
        previous_norm = sqrt(2 * initial);
        int accepted = 0;
        for (size_t damping = 0; damping < 100; damping++) {
            for (size_t j = 0; j < count; j++)
                trial[j] = y[j] - lambda * delta[j];
            project(c, trial);
            residual(c, mesh, nodes, trial, t, trial_phi, trial_stages);
            offset_add(trial_phi, offset, count);
            memcpy(natural, trial_phi, count * sizeof(double));
            bvp_linear(c, &factors, natural, 0);
            double tc = criterion(c, natural, count);
            if (tc <= (1 - 2 * lambda * 0.01) * initial) {
                crit = tc;
                fixed = lambda == 1;
                accepted = 1;
                break;
            }
            lambda =
                MAX(0.1 * lambda, (lambda * lambda * initial) / ((2 * lambda - 1) * initial + tc));
            if (lambda < 0.01)
                break;
        }
        if (!accepted)
            bvp_fail(c, BVP_NEWTON_FAILED);
        memcpy(y, trial, count * sizeof(double));
        SWAP(double *, phi, trial_phi);
        SWAP(double *, stages, trial_stages);
        memcpy(delta, natural, count * sizeof(double));
    }
    bvp_fail(c, BVP_NEWTON_FAILED);
    return NULL;
}
static bvp_status newton(context *c, const double *mesh, size_t nodes, double *y, const tableau *t,
                         const bvp_options *o, const double *offset, double **stages) {
    jmp_buf local;
    jmp_buf *parent = c->trap;
    allocation *mark = c->memory;
    c->trap = &local;
    if (setjmp(local)) {
        bvp_status error = c->error;
        bvp_release(c, mark);
        c->trap = parent;
        return error;
    }
    *stages = newton_inner(c, mesh, nodes, y, t, o, offset);
    c->trap = parent;
    return BVP_OK;
}
static void continuous(context *c, const double *mesh, size_t nodes, const double *y,
                       const tableau *t, double *stages) {
    size_t q = c->q, sc = t->s + t->extra;
    allocation *mark = c->memory;
    double stack[32], *work = q <= 16 ? stack : bvp_alloc(c, 2 * q, sizeof(double));
    double *arg = work, *raw = arg + q;
    if (q <= 16)
        memset(work, 0, 2 * q * sizeof(double));
    for (size_t i = 0; i < nodes - 1; i++) {
        double h = mesh[i + 1] - mesh[i], *ks = stages + i * sc * q;
        for (size_t r = 0; r < t->extra; r++) {
            size_t stage = t->s + r;
            for (size_t j = 0; j < q; j++) {
                double sum = 0;
                for (size_t k = 0; k < stage; k++)
                    sum += t->xe[r][k] * ks[k * q + j];
                arg[j] = t->ve[r] * y[(i + 1) * q + j] + (1 - t->ve[r]) * y[i * q + j] + h * sum;
            }
            rhs(c, mesh[i] + t->ce[r] * h, arg, ks + stage * q, raw);
        }
    }
    bvp_release(c, mark);
}
static void interp_weights(const double *w, size_t sc, size_t q, double h, const double *left,
                           const double *ks, int derivative, double *out) {
    for (size_t j = 0; j < q; j++) {
        double sum = 0;
        for (size_t r = 0; r < sc; r++)
            sum += w[r] * ks[r * q + j];
        out[j] = derivative ? sum : left[j] + h * sum;
    }
}
static void interp(bvp_method method, double tau, double h, size_t q, const double *left,
                   const double *ks, int derivative, double *out) {
    double w[8], wp[8];
    size_t sc = bvp_weights(method, tau, w, wp);
    interp_weights(derivative ? wp : w, sc, q, h, left, ks, derivative, out);
}
static bvp_status find_interval(const double *mesh, size_t nodes, double x, size_t *index,
                                double *tau) {
    double a = mesh[0], b = mesh[nodes - 1], margin = 4 * DBL_EPSILON * MAX(fabs(a), fabs(b));
    if (!isfinite(x) || (x < a && a - x > margin) || (x > b && x - b > margin))
        return BVP_OUT_OF_DOMAIN;
    size_t lo = 0, hi = nodes;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (mesh[mid] <= x)
            lo = mid + 1;
        else
            hi = mid;
    }
    *index = MIN(lo ? lo - 1 : 0, nodes - 2);
    *tau = (x - mesh[*index]) / (mesh[*index + 1] - mesh[*index]);
    return BVP_OK;
}
static double *interpolate_mesh(context *c, bvp_method method, const double *old, size_t oldn,
                                const double *y, const double *ks, const double *mesh,
                                size_t nodes) {
    double *out = bvp_alloc(c, nodes * c->q, sizeof(double));
    size_t sc = method == 2 ? 2 : method == 4 ? 4 : 8;
    size_t interval = 0;
    double a = old[0], b = old[oldn - 1], margin = 4 * DBL_EPSILON * MAX(fabs(a), fabs(b));
    for (size_t j = 0; j < nodes; j++) {
        double x = mesh[j];
        if (!isfinite(x) || (x < a && a - x > margin) || (x > b && x - b > margin))
            bvp_fail(c, BVP_OUT_OF_DOMAIN);
        /* Both meshes are strictly increasing. Retain the same right-sided
         * interval selection without restarting a binary search at each node. */
        while (interval + 1 < oldn && old[interval + 1] <= x)
            interval++;
        size_t i = MIN(interval, oldn - 2);
        double tau = (x - old[i]) / (old[i + 1] - old[i]);
        interp(method, tau, old[i + 1] - old[i], c->q, y + i * c->q, ks + i * sc * c->q, 0,
               out + j * c->q);
    }
    return out;
}
static double *defects(context *c, bvp_method method, const double *mesh, size_t nodes,
                       const double *y, const double *ks, const tableau *t) {
    size_t q = c->q, sc = t->s + t->extra;
    double *out = bvp_alloc(c, nodes - 1, sizeof(double));
    allocation *mark = c->memory;
    double stack[64], *work = q <= 16 ? stack : bvp_alloc(c, 4 * q, sizeof(double));
    double *value = work, *deriv = value + q, *f = deriv + q, *raw = f + q;
    if (q <= 16)
        memset(work, 0, 4 * q * sizeof(double));
    double w[2][8], wp[2][8];
    for (size_t r = 0; r < 2; r++)
        bvp_weights(method, t->defect[r], w[r], wp[r]);
    for (size_t i = 0; i < nodes - 1; i++) {
        double h = mesh[i + 1] - mesh[i];
        for (size_t r = 0; r < 2; r++) {
            double x = mesh[i] + t->defect[r] * h;
            interp_weights(w[r], sc, q, h, y + i * q, ks + i * sc * q, 0, value);
            interp_weights(wp[r], sc, q, h, y + i * q, ks + i * sc * q, 1, deriv);
            bvp_finite(c, value, q);
            bvp_finite(c, deriv, q);
            rhs(c, x, value, f, raw);
            for (size_t j = 0; j < q; j++) {
                double e = fabs(deriv[j] - f[j]) / (1 + fabs(f[j]));
                if (!isfinite(e))
                    bvp_fail(c, BVP_NONFINITE);
                out[i] = MAX(out[i], e);
            }
        }
    }
    bvp_release(c, mark);
    return out;
}
static double *half_mesh(context *c, const double *mesh, size_t nodes, size_t limit) {
    if (nodes - 1 > limit / 2)
        bvp_fail(c, BVP_MESH_LIMIT);
    double *out = bvp_alloc(c, 2 * nodes - 1, sizeof(double));
    for (size_t i = 0; i < nodes - 1; i++) {
        out[2 * i] = mesh[i];
        out[2 * i + 1] = (mesh[i] + mesh[i + 1]) / 2;
    }
    out[2 * nodes - 2] = mesh[nodes - 1];
    if (!valid_mesh(out, 2 * nodes - 1))
        bvp_fail(c, BVP_INVALID_INPUT);
    return out;
}
static double *half_values(context *c, const double *y, size_t nodes) {
    size_t q = c->q;
    double *out = bvp_alloc(c, (2 * nodes - 1) * q, sizeof(double));
    for (size_t i = 0; i < nodes - 1; i++)
        for (size_t j = 0; j < q; j++) {
            out[2 * i * q + j] = y[i * q + j];
            out[(2 * i + 1) * q + j] = (y[i * q + j] + y[(i + 1) * q + j]) / 2;
        }
    memcpy(out + (2 * nodes - 2) * q, y + (nodes - 1) * q, q * sizeof(double));
    return out;
}
static double global_difference(const double *a, const double *b, size_t count, double *error) {
    double norm = 0;
    for (size_t i = 0; i < count; i++) {
        double v = (a[i] - b[i]) / (1 + fabs(a[i]));
        if (error)
            error[i] = v;
        norm = MAX(norm, fabs(v));
    }
    return norm;
}
static void error_intervals(const double *err, size_t q, size_t nodes, double *out) {
    for (size_t i = 0; i < nodes - 1; i++)
        out[i] = MAX(bvp_maximum(err + i * q, q), bvp_maximum(err + (i + 1) * q, q));
}
static double correction(context *c, const double *mesh, size_t nodes, const double *y,
                         const bvp_options *o, const tableau *t, double *errors, int higher) {
    allocation *mark = c->memory;
    size_t count = nodes * c->q;
    tableau high;
    bvp_tableau(o->method == 2 ? 4 : o->method == 4 ? 6 : 8, &high);
    double *copy = bvp_alloc(c, count, sizeof(double)), *offset = NULL, *ks = NULL;
    memcpy(copy, y, count * sizeof(double));
    if (!higher) {
        offset = bvp_alloc(c, count, sizeof(double));
        double *hs = bvp_alloc(c, (nodes - 1) * c->q * (high.s + high.extra), sizeof(double));
        residual(c, mesh, nodes, y, &high, offset, hs);
        memset(offset, 0, c->left * sizeof(double));
        memset(offset + c->left + (nodes - 1) * c->q, 0, (c->q - c->left) * sizeof(double));
    }
    bvp_status e = newton(c, mesh, nodes, copy, higher ? &high : t, o, offset, &ks);
    if (e)
        bvp_fail(c, e);
    double norm = global_difference(y, copy, count, errors);
    bvp_release(c, mark);
    return norm;
}
static double richardson(context *c, const double *mesh, size_t nodes, const double *y,
                         const double *ks, const bvp_options *o, const tableau *t, int derivative,
                         double *errors) {
    allocation *mark = c->memory;
    size_t q = c->q, finen = 2 * nodes - 1, sc = t->s + t->extra;
    double *fm = half_mesh(c, mesh, nodes, o->max_intervals * 2),
           *fy = interpolate_mesh(c, o->method, mesh, nodes, y, ks, fm, finen), *fks = NULL;
    bvp_status e = newton(c, fm, finen, fy, t, o, NULL, &fks);
    if (e)
        bvp_fail(c, e);
    double factor = pow(2.0, o->method), rf = factor / (factor - 1), norm = 0;
    if (!derivative) {
        for (size_t i = 0; i < nodes; i++)
            norm = MAX(norm, global_difference(y + i * q, fy + 2 * i * q, q, NULL));
        norm *= rf;
    } else {
        continuous(c, fm, finen, fy, t, fks);
        double *cd = bvp_alloc(c, q, sizeof(double)), *fd = bvp_alloc(c, q, sizeof(double));
        for (size_t i = 0; i < nodes - 1; i++) {
            double h = mesh[i + 1] - mesh[i];
            errors[i] = 0;
            for (size_t k = 0; k <= 16; k++) {
                double tau = k / 16.0;
                size_t right = tau >= 0.5, j = 2 * i + right;
                interp(o->method, tau, h, q, y + i * q, ks + i * sc * q, 1, cd);
                interp(o->method, right ? 2 * tau - 1 : 2 * tau, h * 0.5, q, fy + j * q,
                       fks + j * sc * q, 1, fd);
                bvp_finite(c, cd, q);
                bvp_finite(c, fd, q);
                for (size_t component = 0; component < c->n; component++) {
                    double er =
                        rf * fabs(cd[component] - fd[component]) / (1 + fabs(fd[component]));
                    if (!isfinite(er))
                        bvp_fail(c, BVP_NONFINITE);
                    errors[i] = MAX(errors[i], er);
                }
            }
            norm = MAX(norm, errors[i]);
        }
    }
    bvp_release(c, mark);
    return norm;
}
static double conditioning(context *c, const double *mesh, size_t nodes, const double *y,
                           const tableau *t) {
    allocation *mark = c->memory;
    size_t q = c->q, count = q * nodes;
    double *phi = bvp_alloc(c, count, sizeof(double)),
           *ks = bvp_alloc(c, (nodes - 1) * q * (t->s + t->extra), sizeof(double));
    residual(c, mesh, nodes, y, t, phi, ks);
    matrix m = bvp_matrix(c, q, c->left, nodes);
    assemble(c, mesh, nodes, y, t, ks, &m);
    bvp_factor(c, &m);
    double *state = bvp_alloc(c, count, sizeof(double)),
           *force = bvp_alloc(c, count, sizeof(double)),
           *vector = bvp_alloc(c, count, sizeof(double)),
           *result = bvp_alloc(c, count, sizeof(double)),
           *dual = bvp_alloc(c, count, sizeof(double)), *middle = bvp_alloc(c, q, sizeof(double)),
           *f = bvp_alloc(c, q, sizeof(double)), *raw = bvp_alloc(c, q, sizeof(double));
    for (size_t j = 0; j < count; j++) {
        state[j] = 1 / (1 + fabs(y[j]));
        force[j] = 1;
        vector[j] = 1.0 / count;
    }
    for (size_t i = 0; i < nodes - 1; i++) {
        double h = mesh[i + 1] - mesh[i];
        for (size_t j = 0; j < q; j++)
            middle[j] = 0.5 * (y[i * q + j] + y[(i + 1) * q + j]);
        rhs(c, mesh[i] + 0.5 * h, middle, f, raw);
        for (size_t j = 0; j < q; j++)
            force[c->left + i * q + j] = h * (1 + fabs(f[j]));
    }
    double estimate = 0;
    size_t previous = count;
    for (size_t it = 0; it < 8; it++) {
        for (size_t j = 0; j < count; j++)
            result[j] = vector[j] * state[j];
        bvp_linear(c, &m, result, 1);
        double next = 0;
        for (size_t j = 0; j < count; j++) {
            result[j] *= force[j];
            next += fabs(result[j]);
        }
        if (next <= estimate)
            break;
        estimate = next;
        for (size_t j = 0; j < count; j++)
            dual[j] = (result[j] < 0 ? -1 : 1) * force[j];
        bvp_linear(c, &m, dual, 0);
        size_t index = 0;
        double product = 0;
        for (size_t j = 0; j < count; j++) {
            dual[j] *= state[j];
            if (fabs(dual[j]) >= fabs(dual[index]))
                index = j;
            product += dual[j] * vector[j];
        }
        if (index == previous || fabs(dual[index]) <= product)
            break;
        previous = index;
        memset(vector, 0, count * sizeof(double));
        vector[index] = 1;
    }
    for (size_t j = 0; j < count; j++)
        result[j] = (j % 2 ? -1 : 1) * (1 + (double)j / (count - 1)) * state[j];
    bvp_linear(c, &m, result, 1);
    double sum = 0;
    for (size_t j = 0; j < count; j++)
        sum += fabs(result[j] * force[j]);
    estimate = MAX(estimate, 2 * sum / (3 * count));
    bvp_release(c, mark);
    return estimate;
}
static double *select_mesh(context *c, const double *mesh, size_t nodes, const double *defect,
                           const double *global, bvp_error_control mode, const bvp_options *o,
                           size_t minimum, size_t *newn) {
    size_t intervals = nodes - 1;
    double p = o->method, *monitor = bvp_alloc(c, intervals, sizeof(double)), total = 0, peak = 0;
    for (size_t i = 0; i < intervals; i++) {
        double de = pow(defect[i] / o->tolerance, 1 / (p + 1)),
               ge = pow(global[i] / o->tolerance, 1 / p),
               num = mode == BVP_DEFECT || mode == BVP_DEFECT_AND_GLOBAL ? de
                     : mode == BVP_GLOBAL                                ? ge
                                                                         : de + ge;
        monitor[i] = num / (mesh[i + 1] - mesh[i]);
        double mass = monitor[i] * (mesh[i + 1] - mesh[i]);
        peak = MAX(peak, mass);
        total += mass;
    }
    double rho = mode == BVP_DEFECT ? 1 : 2;
    if (!isfinite(total) || total == 0 || peak <= rho * (total / intervals)) {
        *newn = 2 * nodes - 1;
        return half_mesh(c, mesh, nodes, o->max_intervals);
    }
    double requested = 1.3 * total + 1;
    if (requested >= (double)SIZE_MAX)
        bvp_fail(c, BVP_MESH_LIMIT);
    size_t count = (size_t)requested;
    if (fabs(((double)count - intervals) / intervals) < 0.1)
        count = (size_t)(1.1 * intervals);
    count = CLAMP(count, MAX(intervals / 2, 1), 4 * intervals);
    count = MAX(count, minimum);
    if (count > o->max_intervals)
        bvp_fail(c, BVP_MESH_LIMIT);
    double *out = bvp_alloc(c, count + 1, sizeof(double));
    out[0] = mesh[0];
    double zeta = total / count, marker = mesh[0], integral = 0;
    size_t k = 0, point = 0;
    while (k < intervals && point < count) {
        double piece = monitor[k] * (mesh[k + 1] - marker);
        if (integral + piece > zeta) {
            out[point + 1] = (zeta - integral) / monitor[k] + marker;
            marker = out[point + 1];
            integral = 0;
            point++;
        } else {
            integral += piece;
            marker = mesh[k + 1];
            k++;
        }
    }
    out[count] = mesh[intervals];
    if (!valid_mesh(out, count + 1)) {
        *newn = 2 * nodes - 1;
        return half_mesh(c, mesh, nodes, o->max_intervals);
    }
    *newn = count + 1;
    return out;
}
static void adopt(context *c, void *ptr) {
    allocation *a = malloc(sizeof(*a));
    if (!a) {
        free(ptr);
        bvp_fail(c, BVP_NO_MEMORY);
    }
    a->ptr = ptr;
    a->bytes = 0;
    a->next = c->memory;
    c->memory = a;
}
static void free_ptr(context *c, void *ptr) {
    bvp_detach(c, ptr);
    free(ptr);
}
static bvp_status try_correction(context *c, const double *mesh, size_t nodes, const double *y,
                                 const bvp_options *o, const tableau *t, double *errors,
                                 double *norm) {
    jmp_buf local;
    jmp_buf *parent = c->trap;
    allocation *mark = c->memory;
    c->trap = &local;
    if (setjmp(local)) {
        bvp_status e = c->error;
        bvp_release(c, mark);
        c->trap = parent;
        return e;
    }
    *norm = correction(c, mesh, nodes, y, o, t, errors, 0);
    c->trap = parent;
    return BVP_OK;
}
static bvp_status try_derivative(context *c, const double *mesh, size_t nodes, const double *y,
                                 const double *ks, const bvp_options *o, const tableau *t,
                                 double *errors, double *norm) {
    jmp_buf local;
    jmp_buf *parent = c->trap;
    allocation *mark = c->memory;
    c->trap = &local;
    if (setjmp(local)) {
        bvp_status e = c->error;
        bvp_release(c, mark);
        c->trap = parent;
        return e;
    }
    *norm = richardson(c, mesh, nodes, y, ks, o, t, 1, errors);
    c->trap = parent;
    return BVP_OK;
}
static void init_statistics(bvp_statistics *s) {
    memset(s, 0, sizeof(*s));
    s->derivative_error = s->global_error = s->richardson_error = s->higher_order_error =
        s->deferred_correction_error = s->condition_estimate = s->conditioning_error = NAN;
}
void bvp_default_options(bvp_options *o) {
    *o = (bvp_options){BVP_MIRK4, 1e-6, BVP_DEFECT, 3000, 40, 1000, 0, 0, 0};
}
const char *bvp_status_string(bvp_status e) {
    static const char *messages[] = {"success",
                                     "invalid input",
                                     "non-finite problem evaluation",
                                     "singular Newton matrix",
                                     "Newton iteration failed",
                                     "mesh limit exceeded",
                                     "refinement limit exceeded",
                                     "evaluation outside domain",
                                     "solution I/O error",
                                     "out of memory"};
    return e >= 0 && e <= BVP_NO_MEMORY ? messages[e] : "unknown status";
}
bvp_status bvp_solve(const bvp_problem *p, const bvp_initial_guess *initial,
                     const bvp_options *options, bvp_solution *out) {
    if (!out)
        return BVP_INVALID_INPUT;
    memset(out, 0, sizeof(*out));
    bvp_options defaults;
    if (!options) {
        bvp_default_options(&defaults);
        options = &defaults;
    }
    if (!p || !initial || !p->rhs || !p->boundary || !valid_mesh(initial->mesh, initial->points) ||
        !initial->values || (!initial->parameters && initial->parameter_count) ||
        p->dimension == 0 || p->dimension > SIZE_MAX - p->parameter_count)
        return BVP_INVALID_INPUT;
    size_t q = p->dimension + p->parameter_count;
    if (p->left_conditions > q || initial->dimension != p->dimension ||
        initial->parameter_count != p->parameter_count || !isfinite(options->tolerance) ||
        options->tolerance < 100 * DBL_EPSILON || !options->max_intervals ||
        !options->max_newton_iterations || !options->max_refinements ||
        (options->method != BVP_MIRK2 && options->method != BVP_MIRK4 &&
         options->method != BVP_MIRK6) ||
        options->error_control < BVP_DEFECT || options->error_control > BVP_DEFECT_AND_DERIVATIVE)
        return BVP_INVALID_INPUT;
    if (initial->points - 1 > options->max_intervals)
        return BVP_MESH_LIMIT;
    if (q > SIZE_MAX / q / sizeof(double) / 20 ||
        options->max_intervals > (SIZE_MAX / sizeof(double) / q / q - 1) / 4 ||
        options->max_newton_iterations > SIZE_MAX / 100)
        return BVP_INVALID_INPUT;
    context *c = calloc(1, sizeof(*c));
    if (!c)
        return BVP_NO_MEMORY;
    jmp_buf top;
    c->trap = &top;
    if (setjmp(top)) {
        bvp_status e = c->error;
        bvp_release(c, NULL);
        free(c);
        return e;
    }
    c->p = p;
    c->n = p->dimension;
    c->q = q;
    c->left = p->left_conditions;
    c->a = initial->mesh[0];
    init_statistics(&c->stats);
    bvp_finite(c, initial->values, initial->points * p->dimension);
    bvp_finite(c, initial->parameters, p->parameter_count);
    if (p->singular_term) {
        size_t n = c->n;
        bvp_finite(c, p->singular_term, n * n);
        double *eye = bvp_alloc(c, n * n, sizeof(double));
        for (size_t i = 0; i < n * n; i++)
            eye[i] = -p->singular_term[i];
        for (size_t i = 0; i < n; i++)
            eye[i * n + i] += 1;
        c->endpoint_inverse = bvp_pinv(c, eye, n);
        double *si = bvp_pinv(c, p->singular_term, n);
        c->projection = bvp_alloc(c, n * n, sizeof(double));
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) {
                double sum = 0;
                for (size_t k = 0; k < n; k++)
                    sum += si[i * n + k] * p->singular_term[k * n + j];
                c->projection[i * n + j] = (i == j ? 1 : 0) - sum;
            }
        free_ptr(c, eye);
        free_ptr(c, si);
    }
    size_t nodes = initial->points;
    double *mesh = bvp_alloc(c, nodes, sizeof(double)),
           *y = bvp_alloc(c, nodes * q, sizeof(double));
    memcpy(mesh, initial->mesh, nodes * sizeof(double));
    for (size_t i = 0; i < nodes; i++) {
        memcpy(y + i * q, initial->values + i * c->n, c->n * sizeof(double));
        if (p->parameter_count)
            memcpy(y + i * q + c->n, initial->parameters, p->parameter_count * sizeof(double));
    }
    project(c, y);
    tableau t;
    bvp_tableau(options->method, &t);
    int global_phase = options->error_control != BVP_DEFECT_AND_GLOBAL, failed_coarsening = 0;
    double last_error = INFINITY;
    size_t stagnant = 0, previous = nodes - 1;
    for (size_t refinement = 0; refinement <= options->max_refinements; refinement++) {
        allocation *iteration_mark = c->memory;
        double *old = bvp_alloc(c, nodes * q, sizeof(double));
        memcpy(old, y, nodes * q * sizeof(double));
        double *ks = NULL, *newmesh = NULL, *newy = NULL;
        size_t newnodes = 0;
        bvp_status e = newton(c, mesh, nodes, y, &t, options, NULL, &ks);
        if (e == BVP_OK) {
            continuous(c, mesh, nodes, y, &t, ks);
            double *defect = defects(c, options->method, mesh, nodes, y, ks, &t);
            double dn = bvp_maximum(defect, nodes - 1);
            if ((options->error_control == BVP_DEFECT ||
                 options->error_control == BVP_DEFECT_AND_DERIVATIVE) &&
                dn > 0.1) {
                if (refinement == options->max_refinements)
                    bvp_fail(c, BVP_REFINEMENT_LIMIT);
                newmesh = half_mesh(c, mesh, nodes, options->max_intervals);
                newy = half_values(c, old, nodes);
                newnodes = 2 * nodes - 1;
            } else {
                double global = NAN, derivative = NAN;
                double *global_sub = bvp_alloc(c, nodes - 1, sizeof(double));
                int retry = 0;
                if (options->error_control == BVP_DEFECT_AND_DERIVATIVE) {
                    e = try_derivative(c, mesh, nodes, y, ks, options, &t, global_sub, &derivative);
                    if (e == BVP_NEWTON_FAILED || e == BVP_SINGULAR || e == BVP_NONFINITE)
                        retry = 1;
                    else if (e)
                        bvp_fail(c, e);
                }
                if (options->error_control == BVP_DEFECT_AND_GLOBAL && dn <= options->tolerance)
                    global_phase = 1;
                if (!retry && global_phase && options->error_control != BVP_DEFECT &&
                    options->error_control != BVP_DEFECT_AND_DERIVATIVE) {
                    double *errors = bvp_alloc(c, nodes * q, sizeof(double));
                    e = try_correction(c, mesh, nodes, y, options, &t, errors, &global);
                    if (e == BVP_NEWTON_FAILED || e == BVP_SINGULAR)
                        retry = 1;
                    else if (e)
                        bvp_fail(c, e);
                    else
                        error_intervals(errors, q, nodes, global_sub);
                }
                if (retry) {
                    if (refinement == options->max_refinements)
                        bvp_fail(c, BVP_REFINEMENT_LIMIT);
                    newmesh = half_mesh(c, mesh, nodes, options->max_intervals);
                    newnodes = 2 * nodes - 1;
                    newy =
                        interpolate_mesh(c, options->method, mesh, nodes, y, ks, newmesh, newnodes);
                } else {
                    double error = options->error_control == BVP_DEFECT   ? dn
                                   : options->error_control == BVP_GLOBAL ? global
                                   : options->error_control == BVP_DEFECT_AND_GLOBAL
                                       ? (global_phase ? global : dn)
                                   : options->error_control == BVP_HYBRID ? dn + global
                                                                          : MAX(dn, derivative);
                    size_t intervals = nodes - 1;
                    int regressed = intervals < previous && error > 2 * last_error;
                    failed_coarsening |= regressed;
                    size_t minimum = regressed           ? (size_t)ceil(1.1 * previous)
                                     : failed_coarsening ? (size_t)ceil(1.1 * intervals)
                                                         : 0;
                    previous = intervals;
                    stagnant = error >= 0.9 * last_error ? stagnant + 1 : 0;
                    last_error = error;
                    if (options->trace)
                        fprintf(stderr,
                                "BVP C: %zu intervals, defect=%.3e, global=%.3e, derivative=%.3e\n",
                                intervals, dn, global, derivative);
                    if (error <= options->tolerance) {
                        c->stats.defect_norm = dn;
                        c->stats.derivative_error = derivative;
                        c->stats.global_error = global;
                        if (options->estimate_error) {
                            double dc = isnan(global)
                                            ? correction(c, mesh, nodes, y, options, &t, NULL, 0)
                                            : global;
                            c->stats.deferred_correction_error = c->stats.global_error = dc;
                            c->stats.higher_order_error =
                                correction(c, mesh, nodes, y, options, &t, NULL, 1);
                            c->stats.richardson_error =
                                richardson(c, mesh, nodes, y, ks, options, &t, 0, NULL);
                        }
                        if (options->estimate_condition) {
                            c->stats.condition_estimate = conditioning(c, mesh, nodes, y, &t);
                            c->stats.conditioning_error = c->stats.condition_estimate * dn;
                        }
                        double *bc = bvp_alloc(c, q, sizeof(double));
                        boundary(c, y, y + (nodes - 1) * q, bc, bc + c->left);
                        c->stats.boundary_residual = bvp_maximum(bc, q);
                        c->stats.scaled_boundary_residual =
                            boundary_norm(c, y, y + (nodes - 1) * q);
                        c->stats.intervals = nodes - 1;
                        c->stats.mesh_refinements = refinement;
                        double *values = bvp_alloc(c, nodes * c->n, sizeof(double)),
                               *parameters = bvp_alloc(c, p->parameter_count, sizeof(double));
                        for (size_t i = 0; i < nodes; i++)
                            memcpy(values + i * c->n, y + i * q, c->n * sizeof(double));
                        if (p->parameter_count)
                            memcpy(parameters, y + c->n, p->parameter_count * sizeof(double));
                        *out = (bvp_solution){nodes,
                                              c->n,
                                              p->parameter_count,
                                              t.s + t.extra,
                                              options->method,
                                              mesh,
                                              values,
                                              parameters,
                                              y,
                                              ks,
                                              c->stats};
                        bvp_detach(c, mesh);
                        bvp_detach(c, values);
                        bvp_detach(c, parameters);
                        bvp_detach(c, y);
                        bvp_detach(c, ks);
                        bvp_release(c, NULL);
                        free(c);
                        return BVP_OK;
                    }
                    if (refinement == options->max_refinements)
                        bvp_fail(c, BVP_REFINEMENT_LIMIT);
                    bvp_error_control mode = options->error_control == BVP_DEFECT_AND_GLOBAL
                                                 ? (global_phase ? BVP_GLOBAL : BVP_DEFECT)
                                                 : options->error_control;
                    minimum = mode == BVP_DEFECT
                                  ? 0
                                  : MAX(minimum, stagnant >= 2 ? (size_t)ceil(1.1 * intervals) : 0);
                    newmesh = select_mesh(c, mesh, nodes, defect, global_sub, mode, options,
                                          minimum, &newnodes);
                    newy =
                        interpolate_mesh(c, options->method, mesh, nodes, y, ks, newmesh, newnodes);
                }
            }
        } else if (e == BVP_NEWTON_FAILED || e == BVP_SINGULAR || e == BVP_NONFINITE) {
            if (refinement == options->max_refinements)
                bvp_fail(c, e);
            newmesh = half_mesh(c, mesh, nodes, options->max_intervals);
            newy = half_values(c, old, nodes);
            newnodes = 2 * nodes - 1;
        } else
            bvp_fail(c, e);
        /* Preserve only the next mesh and state, releasing every refinement scratch. */
        bvp_detach(c, newmesh);
        bvp_detach(c, newy);
        bvp_release(c, iteration_mark);
        free_ptr(c, mesh);
        free_ptr(c, y);
        adopt(c, newmesh);
        adopt(c, newy);
        mesh = newmesh;
        y = newy;
        nodes = newnodes;
        project(c, y);
    }
    bvp_fail(c, BVP_REFINEMENT_LIMIT);
    return BVP_REFINEMENT_LIMIT;
}
bvp_status bvp_linspace(double a, double b, size_t n, double *out) {
    if (!out || n < 2 || !isfinite(a) || !isfinite(b) || a >= b || !isfinite(b - a))
        return BVP_INVALID_INPUT;
    double delta = (b - a) / (n - 1);
    for (size_t i = 0; i < n; i++)
        out[i] = a + delta * i;
    return valid_mesh(out, n) ? BVP_OK : BVP_INVALID_INPUT;
}
void bvp_guess_free(bvp_initial_guess *g) {
    if (g) {
        free(g->mesh);
        free(g->values);
        free(g->parameters);
        memset(g, 0, sizeof(*g));
    }
}
bvp_status bvp_guess_from_values(size_t points, size_t n, size_t np, const double *mesh,
                                 const double *values, const double *parameters,
                                 bvp_initial_guess *out) {
    if (!out)
        return BVP_INVALID_INPUT;
    memset(out, 0, sizeof(*out));
    if (!valid_mesh(mesh, points) || !n || !values || (np && !parameters) ||
        n > SIZE_MAX / points / sizeof(double) || np > SIZE_MAX / sizeof(double))
        return BVP_INVALID_INPUT;
    for (size_t i = 0; i < points * n; i++)
        if (!isfinite(values[i]))
            return BVP_INVALID_INPUT;
    for (size_t i = 0; i < np; i++)
        if (!isfinite(parameters[i]))
            return BVP_INVALID_INPUT;
    out->mesh = malloc(points * sizeof(double));
    out->values = malloc(points * n * sizeof(double));
    out->parameters = malloc((np ? np : 1) * sizeof(double));
    if (!out->mesh || !out->values || !out->parameters) {
        bvp_guess_free(out);
        return BVP_NO_MEMORY;
    }
    out->points = points;
    out->dimension = n;
    out->parameter_count = np;
    memcpy(out->mesh, mesh, points * sizeof(double));
    memcpy(out->values, values, points * n * sizeof(double));
    if (np)
        memcpy(out->parameters, parameters, np * sizeof(double));
    return BVP_OK;
}
bvp_status bvp_guess_from_function(size_t points, size_t n, size_t np, const double *mesh,
                                   const double *parameters,
                                   void (*guess)(void *, double, double *), void *user,
                                   bvp_initial_guess *out) {
    if (!out || !guess || !valid_mesh(mesh, points) || !n ||
        n > SIZE_MAX / MAX(points, 10) / sizeof(double))
        return BVP_INVALID_INPUT;
    size_t actual = points == 2 ? 10 : points;
    double *x = malloc(actual * sizeof(double)), *y = calloc(actual * n, sizeof(double));
    if (!x || !y) {
        free(x);
        free(y);
        return BVP_NO_MEMORY;
    }
    if (points == 2) {
        bvp_status e = bvp_linspace(mesh[0], mesh[1], actual, x);
        if (e) {
            free(x);
            free(y);
            return e;
        }
    } else
        memcpy(x, mesh, actual * sizeof(double));
    for (size_t i = 0; i < actual; i++)
        guess(user, x[i], y + i * n);
    bvp_status e = bvp_guess_from_values(actual, n, np, x, y, parameters, out);
    free(x);
    free(y);
    return e;
}
void bvp_solution_free(bvp_solution *s) {
    if (s) {
        free(s->mesh);
        free(s->values);
        free(s->parameters);
        free(s->augmented);
        free(s->stages);
        memset(s, 0, sizeof(*s));
    }
}
static bvp_status evaluate(const bvp_solution *s, double x, int derivative, double *out) {
    if (!s || !out || !s->mesh || !s->augmented || !s->stages || s->points < 2)
        return BVP_INVALID_INPUT;
    size_t i, q = s->dimension + s->parameter_count;
    double tau;
    bvp_status e = find_interval(s->mesh, s->points, x, &i, &tau);
    if (e)
        return e;
    double w[8], wp[8];
    size_t sc = bvp_weights(s->method, tau, w, wp);
    if (sc != s->stage_count)
        return BVP_INVALID_INPUT;
    double h = s->mesh[i + 1] - s->mesh[i];
    const double *weights = derivative ? wp : w;
    for (size_t j = 0; j < s->dimension; j++) {
        double sum = 0;
        for (size_t k = 0; k < sc; k++)
            sum += weights[k] * s->stages[(i * sc + k) * q + j];
        out[j] = derivative ? sum : s->augmented[i * q + j] + h * sum;
        if (!isfinite(out[j]))
            return BVP_NONFINITE;
    }
    return BVP_OK;
}
bvp_status bvp_evaluate(const bvp_solution *s, double x, double *out) {
    return evaluate(s, x, 0, out);
}
bvp_status bvp_derivative(const bvp_solution *s, double x, double *out) {
    return evaluate(s, x, 1, out);
}
bvp_status bvp_solution_initial_guess(const bvp_solution *s, bvp_initial_guess *out) {
    if (!s)
        return BVP_INVALID_INPUT;
    return bvp_guess_from_values(s->points, s->dimension, s->parameter_count, s->mesh, s->values,
                                 s->parameters, out);
}
bvp_status bvp_extend(const bvp_solution *s, double left, double right, int order,
                      bvp_initial_guess *out) {
    if (!s || !out || order < 0 || order > 1 || !isfinite(left) || !isfinite(right) ||
        left >= right || !valid_mesh(s->mesh, s->points))
        return BVP_INVALID_INPUT;
    size_t n = s->dimension, points = 2;
    for (size_t i = 0; i < s->points; i++)
        if (s->mesh[i] > left && s->mesh[i] < right)
            points++;
    double *mesh = malloc(points * sizeof(double)), *values = calloc(points * n, sizeof(double)),
           *da = calloc(n, sizeof(double)), *db = calloc(n, sizeof(double));
    if (!mesh || !values || !da || !db) {
        free(mesh);
        free(values);
        free(da);
        free(db);
        return BVP_NO_MEMORY;
    }
    mesh[0] = left;
    size_t k = 1;
    for (size_t i = 0; i < s->points; i++)
        if (s->mesh[i] > left && s->mesh[i] < right)
            mesh[k++] = s->mesh[i];
    mesh[k] = right;
    bvp_status e = BVP_OK;
    double a = s->mesh[0], b = s->mesh[s->points - 1];
    if (order) {
        e = bvp_derivative(s, a, da);
        if (!e)
            e = bvp_derivative(s, b, db);
    }
    for (size_t i = 0; !e && i < points; i++) {
        double x = mesh[i];
        if (x < a)
            for (size_t j = 0; j < n; j++)
                values[i * n + j] = s->values[j] + (x - a) * da[j];
        else if (x > b)
            for (size_t j = 0; j < n; j++)
                values[i * n + j] = s->values[(s->points - 1) * n + j] + (x - b) * db[j];
        else
            e = bvp_evaluate(s, x, values + i * n);
    }
    if (!e)
        e = bvp_guess_from_values(points, n, s->parameter_count, mesh, values, s->parameters, out);
    free(mesh);
    free(values);
    free(da);
    free(db);
    return e;
}
bvp_status bvp_extend_with_values(const bvp_solution *s, double left, double right,
                                  const double *lv, const double *rv, bvp_initial_guess *out) {
    if (!lv || !rv)
        return BVP_INVALID_INPUT;
    bvp_status e = bvp_extend(s, left, right, 0, out);
    if (e)
        return e;
    if (left < s->mesh[0])
        memcpy(out->values, lv, s->dimension * sizeof(double));
    if (right > s->mesh[s->points - 1])
        memcpy(out->values + (out->points - 1) * s->dimension, rv, s->dimension * sizeof(double));
    for (size_t i = 0; i < out->points * out->dimension; i++)
        if (!isfinite(out->values[i])) {
            bvp_guess_free(out);
            return BVP_INVALID_INPUT;
        }
    return BVP_OK;
}
bvp_status bvp_save(const bvp_solution *s, const char *path) {
    if (!s || !path)
        return BVP_INVALID_INPUT;
    FILE *f = fopen(path, "w");
    if (!f)
        return BVP_IO_ERROR;
    size_t q = s->dimension + s->parameter_count;
    int ok = fprintf(f, "BVP_SOLVER_C 1\n%zu %zu %d %zu\n", s->dimension, s->parameter_count,
                     s->method, s->points) > 0;
    for (size_t i = 0; i < s->points && ok; i++) {
        ok = fprintf(f, "%.17e", s->mesh[i]) > 0;
        for (size_t j = 0; j < q && ok; j++)
            ok = fprintf(f, " %.17e", s->augmented[i * q + j]) > 0;
        ok = ok && fputc('\n', f) != EOF;
    }
    for (size_t i = 0; i < (s->points - 1) * s->stage_count && ok; i++) {
        for (size_t j = 0; j < q && ok; j++)
            ok = fprintf(f, "%.17e ", s->stages[i * q + j]) > 0;
        ok = ok && fputc('\n', f) != EOF;
    }
    if (fclose(f))
        ok = 0;
    return ok ? BVP_OK : BVP_IO_ERROR;
}
bvp_status bvp_load(const char *path, bvp_solution *s) {
    if (!path || !s)
        return BVP_INVALID_INPUT;
    memset(s, 0, sizeof(*s));
    FILE *f = fopen(path, "r");
    if (!f)
        return BVP_IO_ERROR;
    char magic[32];
    unsigned version;
    int method;
    size_t n, np, points;
    if (fscanf(f, "%31s %u %zu %zu %d %zu", magic, &version, &n, &np, &method, &points) != 6 ||
        strcmp(magic, "BVP_SOLVER_C") || version != 1 || !n || points < 2 ||
        (method != 2 && method != 4 && method != 6) || n > SIZE_MAX - np) {
        fclose(f);
        return BVP_IO_ERROR;
    }
    size_t q = n + np, sc = method == 2 ? 2 : method == 4 ? 4 : 8;
    if (q > SIZE_MAX / sizeof(double) / sc / points) {
        fclose(f);
        return BVP_IO_ERROR;
    }
    *s = (bvp_solution){points, n, np, sc, (bvp_method)method, NULL, NULL, NULL, NULL, NULL, {0}};
    s->mesh = malloc(points * sizeof(double));
    s->values = malloc(points * n * sizeof(double));
    s->parameters = malloc((np ? np : 1) * sizeof(double));
    s->augmented = malloc(points * q * sizeof(double));
    s->stages = malloc((points - 1) * sc * q * sizeof(double));
    if (!s->mesh || !s->values || !s->parameters || !s->augmented || !s->stages) {
        fclose(f);
        bvp_solution_free(s);
        return BVP_NO_MEMORY;
    }
    int ok = 1;
    for (size_t i = 0; i < points && ok; i++) {
        ok = fscanf(f, "%lf", s->mesh + i) == 1;
        for (size_t j = 0; j < q && ok; j++)
            ok = fscanf(f, "%lf", s->augmented + i * q + j) == 1 &&
                 isfinite(s->augmented[i * q + j]);
    }
    for (size_t i = 0; i < (points - 1) * sc * q && ok; i++)
        ok = fscanf(f, "%lf", s->stages + i) == 1 && isfinite(s->stages[i]);
    char extra;
    if (fscanf(f, " %c", &extra) == 1)
        ok = 0;
    fclose(f);
    if (!ok || !valid_mesh(s->mesh, points)) {
        bvp_solution_free(s);
        return BVP_IO_ERROR;
    }
    for (size_t i = 0; i < points; i++)
        memcpy(s->values + i * n, s->augmented + i * q, n * sizeof(double));
    if (np)
        memcpy(s->parameters, s->augmented + n, np * sizeof(double));
    init_statistics(&s->statistics);
    s->statistics.intervals = points - 1;
    return BVP_OK;
}
