#include "integrator.h"
#include "test_data/integrator_fixtures.h"
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void decay(void *ctx, double t, const double *y, const double *yp, double cj, double *r) {
    (void)t;
    (void)cj;
    int quadratic = ctx ? *(int *)ctx : 0;
    r[0] = yp[0] + (quadratic ? y[0] * y[0] : y[0]);
}
static int jac(void *ctx, double t, const double *y, const double *yp, double cj, BacBand *m,
               double *dense) {
    (void)t;
    (void)yp;
    int q = ctx ? *(int *)ctx : 0;
    double v = (q ? 2. * y[0] : 1.) + cj;
    if (m)
        bac_band_set(m, 0, 0, v);
    else
        dense[0] = v;
    return 1;
}
static uint64_t bits(double v) {
    uint64_t r;
    memcpy(&r, &v, sizeof r);
    return r;
}
static void fixture(const char *name, double h, int quadratic, int analytic) {
    FILE *f = bac_integrator_fixture(name);
    assert(f);
    TimeOptions o = bi_daskr_options(1e-6, 1e-6);
    o.initial_step = h;
    double y = 1., yp = -1.;
    BdfIntegrator b;
    assert(!bi_new(&b, 0., &y, &yp, 1, &o));
    DaeSystem s = {.n = 1, .ctx = &quadratic, .residual = decay, .jacobian = analytic ? jac : NULL};
    char line[2048];
    size_t checks = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n')
            continue;
        char *tok = strtok(line, " \n");
        assert(tok);
        size_t index = strtoull(tok, NULL, 10);
        tok = strtok(NULL, " \n");
        size_t used = strtoull(tok, NULL, 10);
        tok = strtok(NULL, " \n");
        size_t next = strtoull(tok, NULL, 10);
        StepOutcome step;
        if (bi_step(&b, &s, 1., 0, &step)) {
            fprintf(stderr, "%s step %zu: %s\n", name, index, b.error);
            abort();
        }
        assert(step.order == used && b.order == next);
        double values[19] = {b.t,    step.step_size, b.h, b.previous_cj, step.error_norm,
                             b.y[0], b.yp[0]};
        for (size_t j = 0; j < 6; j++)
            values[7 + j] = b.phi[j];
        for (size_t j = 0; j < 6; j++)
            values[13 + j] = b.psi[j];
        for (size_t j = 0; j < 19; j++) {
            tok = strtok(NULL, " \n");
            assert(tok);
            uint64_t expected = strtoull(tok, NULL, 16);
            if (bits(values[j]) != expected) {
                double v;
                memcpy(&v, &expected, sizeof v);
                fprintf(stderr,
                        "%s step %zu field %zu got %016" PRIX64 " %.17g expected %016" PRIX64
                        " %.17g\n",
                        name, index, j, bits(values[j]), values[j], expected, v);
                abort();
            }
            checks++;
        }
    }
    fclose(f);
    bi_free(&b);
    printf("%s: %zu exact binary64 fields\n", name, checks);
}
static void stiff(void *ctx, double t, const double *y, const double *yp, double cj, double *r) {
    (void)ctx;
    (void)cj;
    r[0] = yp[0] + 1000. * (y[0] - cos(t)) + sin(t);
}
static int stiff_jac(void *ctx, double t, const double *y, const double *yp, double cj, BacBand *m,
                     double *dense) {
    (void)ctx;
    (void)t;
    (void)y;
    (void)yp;
    if (m)
        bac_band_set(m, 0, 0, 1000. + cj);
    else
        dense[0] = 1000. + cj;
    return 1;
}
static void boundary(void *ctx, double t, const double *y, const double *yp, double cj, double *r) {
    (void)ctx;
    (void)cj;
    r[0] = yp[0] + y[0] - y[1];
    r[1] = y[1] - cos(t);
}
static void unit(void *ctx, double t, const double *y, const double *yp, double cj, double *r) {
    (void)ctx;
    (void)t;
    (void)y;
    (void)cj;
    r[0] = yp[0] - 1.;
}
static int roots(void *ctx, double t, const double *y, const double *yp, double *r, size_t n) {
    (void)ctx;
    (void)t;
    (void)yp;
    assert(n == 3);
    r[0] = y[0] - .3;
    r[1] = 2. * y[0] - .6;
    r[2] = .7 - y[0];
    return 0;
}
static int project(void *ctx, const double *v, size_t n, double *out, size_t m) {
    assert(n == m);
    double scale = ctx ? *(double *)ctx : 1.;
    for (size_t i = 0; i < n; i++)
        out[i] = v[i] * scale;
    return 0;
}
int main(void) {
    fixture("daskr_scalar.trace", 0., 0, 1);
    fixture("daskr_scalar_rejected.trace", .1, 0, 1);
    fixture("daskr_scalar_quadratic.trace", 0., 1, 1);
    fixture("daskr_scalar_numeric.trace", 0., 0, 0);
    fixture("daskr_scalar_quadratic_numeric.trace", 0., 1, 0);
    for (size_t order = 1; order <= 5; order++) {
        TimeOptions o = bi_options(1e-8, 1e-7);
        o.max_order = order;
        double y = 1., yp = 0.;
        BdfIntegrator b;
        assert(!bi_new(&b, 0., &y, &yp, 1, &o));
        DaeSystem s = {.n = 1, .residual = stiff, .jacobian = stiff_jac};
        StepOutcome step;
        while (b.t < 3.) {
            int rc = bi_step(&b, &s, 3., 1, &step);
            if (rc)
                fprintf(stderr, "%s\n", b.error);
            assert(!rc);
        }
        assert(fabs(b.y[0] - cos(3.)) < 3e-6);
        assert(b.statistics.maximum_order == order);
        double a, c;
        assert(!bi_interpolate(&b, (b.t + b.last_start) * .5, &a, &c));
        assert(fabs(a - cos((b.t + b.last_start) * .5)) < 4e-6);
        BdfIntegrator clone;
        assert(!bi_clone(&clone, &b));
        assert(clone.t == b.t && clone.statistics.accepted_steps == b.statistics.accepted_steps);
        bi_free(&clone);
        bi_free(&b);
    }
    {
        TimeOptions o = bi_options(1e-8, 1e-7);
        double y[] = {.5, 9.}, yp[] = {7., 0.};
        unsigned char d[] = {1, 0};
        BdfIntegrator b;
        assert(!bi_new(&b, .3, y, yp, 2, &o));
        DaeSystem s = {.n = 2, .residual = boundary};
        assert(!bi_initialize_consistent(&b, &s, d));
        assert(fabs(b.y[0] - .5) < 1e-14 && fabs(b.y[1] - cos(.3)) < 1e-10 &&
               fabs(b.yp[1] + sin(.3)) < 1e-6);
        StepOutcome step;
        while (b.t < 2.)
            assert(!bi_step(&b, &s, 2., 1, &step));
        double exact = .5 * (cos(2.) + sin(2.)) + (.5 - .5 * (cos(.3) + sin(.3))) * exp(.3 - 2.);
        assert(fabs(b.y[0] - exact) < 3e-6 && fabs(b.y[1] - cos(2.)) < 1e-10);
        bi_free(&b);
    }
    for (int method = 0; method < 2; method++) {
        TimeOptions o = method ? bi_daskr_options(1e-9, 1e-9) : bi_options(1e-9, 1e-9);
        o.initial_step = 1.;
        double y = 0., yp = 1.;
        BdfIntegrator b;
        assert(!bi_new(&b, 0., &y, &yp, 1, &o));
        DaeSystem s = {.n = 1, .residual = unit};
        StepOutcome step;
        assert(!bi_step(&b, &s, 1., 1, &step));
        int dirs[] = {BI_RISING, BI_ANY, BI_FALLING};
        RootEvent event;
        assert(bi_locate_roots(&b, roots, NULL, 3, dirs, 1., &event) == 1);
        assert(fabs(event.time - .3) < 1e-12 && event.count == 2 && event.indices[0] == 0 &&
               event.indices[1] == 1);
        assert(!bi_reset(&b, event.time, event.state, event.derivative, 1));
        bi_root_free(&event);
        assert(!bi_step(&b, &s, 1., 1, &step));
        int falling[] = {BI_FALLING, BI_FALLING, BI_FALLING};
        assert(bi_locate_roots(&b, roots, NULL, 3, falling, 1., &event) == 1);
        assert(fabs(event.time - .7) < 1e-12 && event.count == 1 && event.indices[0] == 2);
        bi_root_free(&event);
        bi_free(&b);
    }
    for (int method = 0; method < 2; method++) {
        TimeOptions o = method ? bi_daskr_options(1e-7, 1e-7) : bi_options(1e-7, 1e-7);
        double y = 1., yp = -1.;
        BdfIntegrator b, q;
        assert(!bi_new(&b, 0., &y, &yp, 1, &o));
        assert(!bi_new(&q, 0., &y, &yp, 1, &o));
        DaeSystem s = {.n = 1, .residual = decay, .jacobian = jac};
        StepOutcome step;
        for (size_t i = 0; i < 30; i++)
            assert(!bi_step(&b, &s, 1., 0, &step));
        double scale = 1.3;
        assert(bi_project_history(&q, &b, 1, &o, project, &scale) == 1);
        assert(fabs(q.y[0] - scale * b.y[0]) < 1e-14);
        assert(q.statistics.accepted_steps == 0);
        assert(!bi_step(&q, &s, 1., 0, &step));
        bi_free(&q);
        bi_free(&b);
    }
    {
        TimeOptions o = bi_options(1e-9, 1e-9);
        o.initial_step = .1;
        o.min_step = .1;
        double y = 0., yp = 1.;
        BdfIntegrator b;
        DaeSystem s = {.n = 1, .residual = unit};
        StepOutcome step;
        assert(!bi_new(&b, 0., &y, &yp, 1, &o));
        assert(bi_step(&b, &s, .05, 1, &step));
        assert(b.t == 0. && b.y[0] == 0.);
        bi_free(&b);
    }
    puts("BDF orders 1-5, stiff DAE, initialization, roots/restart, snapshots, warm transfer and "
         "invalid step checks passed");
    return 0;
}
