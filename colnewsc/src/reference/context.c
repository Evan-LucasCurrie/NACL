/* Reference solver context, callback adapters and solve entry point. */
#include "reference_internal.h"
static void *context_alloc(cnref_context *c, size_t count, size_t size) {
    if (size && count > SIZE_MAX / size) {
        c->allocation_failed = 1;
        return NULL;
    }
    cnref_allocation *a = malloc(sizeof(*a));
    if (!a) {
        c->allocation_failed = 1;
        return NULL;
    }
    a->data = calloc(count ? count : 1, size ? size : 1);
    if (!a->data) {
        free(a);
        c->allocation_failed = 1;
        return NULL;
    }
    a->next = c->allocations;
    c->allocations = a;
    return a->data;
}
static void context_free(cnref_context *c) {
    if (!c)
        return;
    while (c->allocations) {
        cnref_allocation *a = c->allocations;
        c->allocations = a->next;
        free(a->data);
        free(a);
    }
    free(c->mesh_history);
    free(c);
}
static cnref_context *context_new(const colnewsc_reference_problem *p, int maximum) {
    cnref_context *c = calloc(1, sizeof(*c));
    if (!c)
        return NULL;
    int d = (int)p->base.dimension;
    size_t points = (size_t)maximum + 1;
    c->problem = p;
    c->iprint = 1;
    c->ncomp = c->mstar = d;
    c->mmax = 1;
    c->nmax = maximum;
    c->aright = 1;
    c->nonlin = 1;
    c->limit = 40;
    c->iguess = 1;
    c->ntol = d;
    c->maximum_corrections = 100;
    c->maximum_refinements = 20;
#define REAL(field, count) c->field = context_alloc(c, (size_t)(count), sizeof(double))
#define INTEGER(field, count) c->field = context_alloc(c, (size_t)(count), sizeof(int))
    REAL(rho, 7);
    REAL(coef, 49);
    INTEGER(m, 20);
    REAL(zeta, 40);
    REAL(wgtmsh, 40);
    REAL(wgterr, 40);
    REAL(tolin, 40);
    REAL(root, 40);
    INTEGER(jtol, 40);
    INTEGER(ltol, 40);
    REAL(b, 28);
    REAL(acol, 28 * 7);
    REAL(asave, 28 * 7);
    REAL(y_meshnew, points * d);
    REAL(y_meshold, points * d);
    REAL(zold, points * d);
    REAL(dmzold, (size_t)maximum * 7 * d);
    REAL(meshold, points);
    REAL(y_super, points * d);
    REAL(sk_discrete_super, (size_t)maximum * 10 * d);
    REAL(sk_interp_super, (size_t)maximum * 10 * d);
    REAL(c_star, 10);
    REAL(v_star, 10);
    REAL(x_star, 100);
#undef REAL
#undef INTEGER
    if (c->allocation_failed) {
        context_free(c);
        return NULL;
    }
    return c;
}
static int finite_vector(const double *v, int n) {
    for (int i = 0; i < n; i++)
        if (!isfinite(v[i]))
            return 0;
    return 1;
}
static void call_rhs(cnref_context *c, double x, R z, R f, int interpolation) {
    double zv[40], fv[20];
    for (int i = 0; i < c->ncomp; i++)
        fv[i] = NAN;
    for (int i = 0; i < c->mstar; i++)
        zv[i] = z[i];
    int valid = isfinite(x) && finite_vector(zv, c->mstar),
        preserve = interpolation && c->problem->preserve_non_finite_interpolation;
    if (valid || preserve) {
        if (c->problem->reference_rhs)
            c->problem->reference_rhs(c->problem->base.context, x, zv, fv);
        else
            c->problem->base.rhs(c->problem->base.context, x, zv, fv);
    }
    if (!isfinite(x) || !finite_vector(zv, c->mstar) || !finite_vector(fv, c->ncomp)) {
        if (preserve) {
            if (c->diagnostics.non_finite_interpolation_evaluations < SIZE_MAX)
                c->diagnostics.non_finite_interpolation_evaluations++;
            if (!c->diagnostics.has_first_non_finite_interpolation_point) {
                c->diagnostics.has_first_non_finite_interpolation_point = 1;
                c->diagnostics.first_non_finite_interpolation_point = x;
            }
        } else
            c->invalid_evaluation = 1;
    }
    for (int i = 0; i < c->ncomp; i++)
        f[i] = fv[i];
}
void cnref_fsub(cnref_context *c, double x, R z, R f) {
    call_rhs(c, x, z, f, 0);
}
void cnref_fsub_interpolant(cnref_context *c, double x, R z, R f) {
    call_rhs(c, x, z, f, 1);
}
double cnref_gsub(cnref_context *c, int condition, R z) {
    double zv[40];
    for (int i = 0; i < c->mstar; i++)
        zv[i] = z[i];
    if (!finite_vector(zv, c->mstar)) {
        c->invalid_evaluation = 1;
        return NAN;
    }
    double out =
        c->problem->reference_boundary
            ? c->problem->reference_boundary(c->problem->base.context, (size_t)condition - 1, zv)
            : c->problem->base.boundary(c->problem->base.context, (size_t)condition - 1, zv);
    if (!isfinite(out) || !finite_vector(zv, c->mstar))
        c->invalid_evaluation = 1;
    return out;
}
void cnref_guess(cnref_context *c, double x, R z, R dmz) {
    double zv[40], dv[20];
    for (int i = 0; i < c->mstar; i++)
        zv[i] = NAN;
    for (int i = 0; i < c->ncomp; i++)
        dv[i] = NAN;
    if (c->problem->reference_guess)
        c->problem->reference_guess(c->problem->base.context, x, zv, dv);
    else {
        if (c->problem->base.initial_guess)
            c->problem->base.initial_guess(c->problem->base.context, x, zv);
        else
            for (int i = 0; i < c->mstar; i++)
                zv[i] = 0;
        for (int i = 0; i < c->ncomp; i++)
            dv[i] = 0;
    }
    if (!finite_vector(zv, c->mstar) || !finite_vector(dv, c->ncomp))
        c->invalid_evaluation = 1;
    for (int i = 0; i < c->mstar; i++)
        z[i] = zv[i];
    for (int i = 0; i < c->ncomp; i++)
        dmz[i] = dv[i];
}
void cnref_dfdy(cnref_context *c, double x, R z, R df) {
    double fz[20] = {0}, fd[20] = {0}, sqrtu = sqrt(DBL_EPSILON);
    cnref_fsub(c, x, z, fz);
    for (int j = 1; j <= c->mstar; j++) {
        double save = GET(z, j), delta = save == 0 ? sqrtu : sqrtu * fabs(save);
        if (delta == 0 || !isfinite(delta) || !isfinite(save + delta) || save + delta == save) {
            c->invalid_evaluation = 1;
            return;
        }
        SET(z, j, save + delta);
        cnref_fsub(c, x, z, fd);
        for (int i = 1; i <= c->ncomp; i++) {
            double derivative = (GET(fd, i) - GET(fz, i)) / delta;
            if (!isfinite(derivative))
                c->invalid_evaluation = 1;
            SET2(df, i, j, c->ncomp, derivative);
        }
        SET(z, j, save);
    }
}
void cnref_dbcdy(cnref_context *c, int condition, R z, R dg) {
    double sqrtu = sqrt(DBL_EPSILON), gz = cnref_gsub(c, condition, z);
    for (int j = 1; j <= c->mstar; j++) {
        double save = GET(z, j), delta = save == 0 ? sqrtu : sqrtu * fabs(save);
        if (delta == 0 || !isfinite(delta) || !isfinite(save + delta) || save + delta == save) {
            c->invalid_evaluation = 1;
            return;
        }
        SET(z, j, save + delta);
        double gd = cnref_gsub(c, condition, z), derivative = (gd - gz) / delta;
        if (!isfinite(derivative))
            c->invalid_evaluation = 1;
        SET(dg, j, derivative);
        SET(z, j, save);
    }
}
void cnref_history(cnref_context *c, int n) {
    if (c->mesh_history_count == c->mesh_history_capacity) {
        size_t capacity = c->mesh_history_capacity ? c->mesh_history_capacity * 2 : 16;
        if (capacity < c->mesh_history_capacity || capacity > SIZE_MAX / sizeof(int)) {
            c->allocation_failed = c->invalid_evaluation = 1;
            return;
        }
        int *history = realloc(c->mesh_history, capacity * sizeof(int));
        if (!history) {
            c->allocation_failed = c->invalid_evaluation = 1;
            return;
        }
        c->mesh_history = history;
        c->mesh_history_capacity = capacity;
    }
    c->mesh_history[c->mesh_history_count++] = n;
}
int cnref_valid_mesh(R xi, int n) {
    if (n <= 0)
        return 0;
    for (int i = 0; i < n; i++) {
        double a = xi[i], b = xi[i + 1], mid = (a + b) * 0.5;
        if (!isfinite(a) || !isfinite(b) || !isfinite(mid) || !(a < mid && mid < b))
            return 0;
    }
    return 1;
}
double cnref_powi(double base, int exponent) {
    unsigned int e = exponent < 0 ? (unsigned int)(-(exponent + 1)) + 1u : (unsigned int)exponent;
    double result = 1;
    while (e) {
        if (e & 1u)
            result *= base;
        e >>= 1;
        if (e)
            base *= base;
    }
    return exponent < 0 ? 1 / result : result;
}
static int initial_valid(const double *mesh, size_t count, double a, double b, size_t max) {
    if (count < 2 || count - 1 > max || mesh[0] != a || mesh[count - 1] != b)
        return 0;
    for (size_t i = 0; i < count - 1; i++)
        if (!isfinite(mesh[i]) || mesh[i] >= mesh[i + 1])
            return 0;
    return 1;
}
colnewsc_status colnewsc_reference_solve(const colnewsc_reference_problem *p,
                                         const colnewsc_options *given,
                                         colnewsc_reference_solution **out) {
    if (!out)
        return COLNEWSC_INVALID_INPUT;
    *out = NULL;
    colnewsc_options defaults;
    if (!given) {
        defaults = colnewsc_default_options();
        given = &defaults;
    }
    const colnewsc_options *o = given;
    if (!p || !p->base.rhs || !p->base.boundary || p->base.dimension < 1 ||
        p->base.dimension > 20 || o->collocation_points < 1 || o->collocation_points > 7 ||
        o->max_subintervals > INT_MAX || !isfinite(o->tolerance) || o->tolerance <= 0 ||
        o->max_newton_iterations == 0)
        return COLNEWSC_INVALID_INPUT;
    int d = (int)p->base.dimension, k = (int)o->collocation_points, max = (int)o->max_subintervals;
    if (!o->initial_mesh &&
        (o->initial_subintervals == 0 || o->initial_subintervals > o->max_subintervals))
        return COLNEWSC_INVALID_INPUT;
    if (o->component_tolerances && o->component_tolerance_count != (size_t)d)
        return COLNEWSC_INVALID_INPUT;
    size_t fper = 4 + 3 * d + (5 + k * d) * ((k + 1) * d) + 4 * d * d, iper = 3 + (k + 1) * d;
    if ((size_t)max > INT_MAX / fper || (size_t)max > INT_MAX / iper)
        return COLNEWSC_INVALID_INPUT;
    int ndimf = (int)((size_t)max * fper), ndimi = (int)((size_t)max * iper);
    double left = p->base.a, right = p->base.b;
    if (!isfinite(left) || !isfinite(right) || left >= right || !isfinite(right - left) ||
        p->base.left_boundary_conditions > (size_t)d)
        return COLNEWSC_INVALID_INPUT;
    size_t initial_count = o->initial_mesh ? o->initial_mesh_count : o->initial_subintervals + 1;
    if (initial_count < 2 || initial_count - 1 > (size_t)max)
        return COLNEWSC_INVALID_INPUT;
    double *initial = malloc(initial_count * sizeof(double));
    if (!initial)
        return COLNEWSC_ALLOCATION_FAILED;
    if (o->initial_mesh)
        memcpy(initial, o->initial_mesh, initial_count * sizeof(double));
    else {
        double delta = (right - left) / o->initial_subintervals;
        for (size_t i = 0; i < initial_count; i++)
            initial[i] = left + delta * i;
        initial[initial_count - 1] = right;
    }
    if (!initial_valid(initial, initial_count, left, right, (size_t)max)) {
        free(initial);
        return COLNEWSC_INVALID_INPUT;
    }
    for (int i = 0; i < d; i++)
        if (o->component_tolerances &&
            (!isfinite(o->component_tolerances[i]) || o->component_tolerances[i] <= 0)) {
            free(initial);
            return COLNEWSC_INVALID_INPUT;
        }
    if (o->fixed_point_count > INT_MAX || (o->fixed_point_count && !o->fixed_points)) {
        free(initial);
        return COLNEWSC_INVALID_INPUT;
    }
    for (size_t i = 0; i < o->fixed_point_count; i++)
        if (!isfinite(o->fixed_points[i]) || o->fixed_points[i] <= left ||
            o->fixed_points[i] >= right || (i && o->fixed_points[i - 1] >= o->fixed_points[i])) {
            free(initial);
            return COLNEWSC_INVALID_INPUT;
        }
    cnref_context *c = context_new(p, max);
    if (!c) {
        free(initial);
        return COLNEWSC_ALLOCATION_FAILED;
    }
    c->maximum_corrections = o->max_newton_iterations;
    c->maximum_refinements = o->max_refinements;
    c->iprint = o->trace ? -1 : 1;
    c->aleft = left;
    c->aright = right;
    c->k = k;
    c->kd = k * d;
    c->n = (int)initial_count - 1;
    c->precis = 1;
    do {
        c->precis /= 2;
    } while (c->precis + 1 > 1);
    c->precis *= 100;
    for (int i = 1; i <= d; i++) {
        SET(c->m, i, 1);
        SET(c->zeta, i, i <= (int)p->base.left_boundary_conditions ? left : right);
        SET(c->ltol, i, i);
        SET(c->tolin, i, o->component_tolerances ? o->component_tolerances[i - 1] : o->tolerance);
    }
    c->mshlmt = 3;
    c->mshnum = 1;
    c->mshalt = 1;
    c->limit = 40;
    int nrec = d - (int)p->base.left_boundary_conditions;
    c->nmax = CN_MIN((ndimf - (nrec * (2 * d) + 5 * d + 3)) /
                         (4 + 3 * d + (c->kd + 5) * (c->kd + d) + (2 * d - nrec) * 2 * d),
                     (ndimi - d) / (3 + c->kd + d));
    colnewsc_status status = COLNEWSC_OK;
    if (c->nmax < c->n) {
        status = COLNEWSC_INVALID_INPUT;
        goto failed;
    }
    R fspace = context_alloc(c, (size_t)ndimf, sizeof(double));
    I ispace = context_alloc(c, (size_t)ndimi, sizeof(int));
    if (c->allocation_failed) {
        status = COLNEWSC_ALLOCATION_FAILED;
        goto failed;
    }
    for (size_t i = 0; i < initial_count; i++)
        fspace[i] = initial[i];
    int lxi = 1, lg = lxi + c->nmax + 1, lxiold = lg + 2 * d * (c->nmax * (2 * d - nrec) + nrec),
        lw = lxiold + c->nmax + 1, lv = lw + c->kd * c->kd * c->nmax, lz = lv + d * c->kd * c->nmax,
        ldmz = lz + d * (c->nmax + 1), ldelz = ldmz + c->kd * c->nmax,
        ldeldz = ldelz + d * (c->nmax + 1), ldqz = ldeldz + c->kd * c->nmax,
        ldqdmz = ldqz + d * (c->nmax + 1), lrhs = ldqdmz + c->kd * c->nmax,
        lvalst = lrhs + c->kd * c->nmax + d, lslope = lvalst + 4 * d * c->nmax,
        laccum = lslope + c->nmax, lscl = laccum + c->nmax + 1, ldscl = lscl + d * (c->nmax + 1),
        lpvtg = 1, lpvtw = lpvtg + d * (c->nmax + 1), linteg = lpvtw + c->kd * c->nmax;
    int nfx = (int)o->fixed_point_count;
    R fix = context_alloc(c, o->fixed_point_count, sizeof(double));
    if (c->allocation_failed) {
        status = COLNEWSC_ALLOCATION_FAILED;
        goto failed;
    }
    if (nfx)
        memcpy(fix, o->fixed_points, (size_t)nfx * sizeof(double));
    cnref_consts(c, k, c->rho, c->coef);
    if (!finite_vector(c->wgtmsh, c->ntol)) {
        status = COLNEWSC_NON_FINITE;
        goto failed;
    }
    double dummy[1] = {0};
    int mode = !o->initial_mesh && nfx > 0 ? 3 : 4;
    if (nfx >= c->nmax) {
        status = COLNEWSC_INVALID_INPUT;
        goto failed;
    }
    if (o->initial_mesh)
        for (int i = 0; i < nfx; i++) {
            int found = 0;
            for (size_t j = 0; j < initial_count; j++)
                if (initial[j] == fix[i])
                    found = 1;
            if (!found) {
                status = COLNEWSC_INVALID_INPUT;
                goto failed;
            }
        }
    cnref_newmsh(c, &mode, VIEW(fspace, lxi), VIEW(fspace, lxiold), dummy, dummy, dummy, dummy,
                 dummy, nfx, fix);
    if (!cnref_valid_mesh(VIEW(fspace, lxi), c->n)) {
        status = COLNEWSC_NON_FINITE;
        goto failed;
    }
    for (int i = 1; i <= c->n + 1; i++)
        SET(fspace, lxiold + i - 1, GET(fspace, lxi + i - 1));
    c->nold = c->n;
    int flag = -3;
    cnref_contrl(c, VIEW(fspace, lxi), VIEW(fspace, lxiold), VIEW(fspace, lz), VIEW(fspace, ldmz),
                 VIEW(fspace, lrhs), VIEW(fspace, ldelz), VIEW(fspace, ldeldz), VIEW(fspace, ldqz),
                 VIEW(fspace, ldqdmz), VIEW(fspace, lg), VIEW(fspace, lw), VIEW(fspace, lv),
                 VIEW(fspace, lvalst), VIEW(fspace, lslope), VIEW(fspace, lscl),
                 VIEW(fspace, ldscl), VIEW(fspace, laccum), VIEW(ispace, lpvtg),
                 VIEW(ispace, linteg), VIEW(ispace, lpvtw), nfx, fix, &flag);
    if (c->allocation_failed) {
        status = COLNEWSC_ALLOCATION_FAILED;
        goto failed;
    }
    if (flag != 1) {
        status = flag == -4   ? COLNEWSC_NON_FINITE
                 : flag == -5 ? COLNEWSC_REFINEMENT_LIMIT
                 : flag == 0  ? COLNEWSC_SINGULAR_JACOBIAN
                 : flag == -1 ? COLNEWSC_MESH_LIMIT
                 : flag == -2 ? COLNEWSC_NEWTON_NOT_CONVERGED
                              : COLNEWSC_INVALID_INPUT;
        goto failed;
    }
    SET(ispace, 1, c->n);
    SET(ispace, 2, k);
    SET(ispace, 3, d);
    SET(ispace, 4, d);
    SET(ispace, 5, 1);
    SET(ispace, 6, c->nz + c->ndmz + c->n + 2);
    SET(ispace, 7, GET(ispace, 6) + k * k - 1);
    for (int i = 1; i <= d; i++)
        SET(ispace, 7 + i, 1);
    for (int i = 1; i <= c->nz; i++)
        SET(fspace, c->n + 1 + i, GET(fspace, lz - 1 + i));
    int idmz = c->n + 1 + c->nz;
    for (int i = 1; i <= c->ndmz; i++)
        SET(fspace, idmz + i, GET(fspace, ldmz - 1 + i));
    int ic = idmz + c->ndmz;
    for (int i = 1; i <= k * k; i++)
        SET(fspace, ic + i, GET(c->coef, i));
    colnewsc_reference_solution *solution = calloc(1, sizeof(*solution));
    if (!solution) {
        status = COLNEWSC_ALLOCATION_FAILED;
        goto failed;
    }
    solution->context = c;
    solution->fspace = fspace;
    solution->ispace = ispace;
    solution->mesh_count = (size_t)c->n + 1;
    solution->mesh = context_alloc(c, solution->mesh_count, sizeof(double));
    solution->output = context_alloc(c, (size_t)d, sizeof(double));
    if (c->allocation_failed) {
        free(solution);
        status = COLNEWSC_ALLOCATION_FAILED;
        goto failed;
    }
    memcpy(solution->mesh, fspace, solution->mesh_count * sizeof(double));
    *out = solution;
    free(initial);
    return COLNEWSC_OK;
failed:
    free(initial);
    context_free(c);
    return status;
}
void colnewsc_reference_free(colnewsc_reference_solution *s) {
    if (s) {
        context_free(s->context);
        free(s);
    }
}
size_t colnewsc_reference_mesh_count(const colnewsc_reference_solution *s) {
    return s ? s->mesh_count : 0;
}
const double *colnewsc_reference_mesh(const colnewsc_reference_solution *s) {
    return s ? s->mesh : NULL;
}
const int *colnewsc_reference_mesh_history(const colnewsc_reference_solution *s, size_t *count) {
    if (count)
        *count = s ? s->context->mesh_history_count : 0;
    return s ? s->context->mesh_history : NULL;
}
const double *colnewsc_reference_representation(const colnewsc_reference_solution *s,
                                                size_t *count) {
    if (count)
        *count = s ? (size_t)GET(s->ispace, 7) : 0;
    return s ? s->fspace : NULL;
}
colnewsc_reference_diagnostics
colnewsc_reference_get_diagnostics(const colnewsc_reference_solution *s) {
    colnewsc_reference_diagnostics empty = {0, 0, 0};
    return s ? s->context->diagnostics : empty;
}
colnewsc_status colnewsc_reference_evaluate(colnewsc_reference_solution *s, double x, double *out,
                                            size_t count) {
    if (!s || !out || count != (size_t)GET(s->ispace, 4))
        return COLNEWSC_INVALID_INPUT;
    if (!isfinite(x) || x < s->mesh[0] || x > s->mesh[s->mesh_count - 1])
        return COLNEWSC_OUTSIDE_INTERVAL;
    cnref_context *c = s->context;
    c->invalid_evaluation = 0;
    cnref_appsln(c, x, s->output, s->fspace, s->ispace);
    if (c->invalid_evaluation || !finite_vector(s->output, GET(s->ispace, 4)))
        return COLNEWSC_NON_FINITE;
    memcpy(out, s->output, count * sizeof(double));
    return COLNEWSC_OK;
}
