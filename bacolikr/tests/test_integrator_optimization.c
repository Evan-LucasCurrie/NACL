#include "integrator.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static void residual(void *ctx, double t, const double *y, const double *yp, double cj,
                     double *out) {
    (void)ctx;
    (void)t;
    (void)cj;
    out[0] = yp[0] + y[0] + .1 * y[1];
    out[1] = yp[1] + 2. * y[1];
}
static int jacobian(void *ctx, double t, const double *y, const double *yp, double cj,
                    BacBand *band, double *dense) {
    (void)ctx;
    (void)t;
    (void)y;
    (void)yp;
    if (band) {
        bac_band_set(band, 0, 0, cj + 1.);
        bac_band_set(band, 0, 1, .1);
        bac_band_set(band, 1, 0, 0.);
        bac_band_set(band, 1, 1, cj + 2.);
    } else {
        dense[0] = cj + 1.;
        dense[1] = .1;
        dense[2] = 0.;
        dense[3] = cj + 2.;
    }
    return 1;
}
static int project(void *ctx, const double *in, size_t n, double *out, size_t count) {
    int *calls = ctx;
    ++*calls;
    if (*calls < 0)
        return -1;
    assert(n == 2 && count == 3);
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[0] + 2. * in[1];
    return 0;
}
static void equal_state(const BdfIntegrator *a, const BdfIntegrator *b) {
    assert(a->n == b->n && a->t == b->t && a->order == b->order);
    assert(a->last_order == b->last_order && a->h == b->h);
    assert(!memcmp(a->y, b->y, a->n * sizeof(double)));
    assert(!memcmp(a->yp, b->yp, a->n * sizeof(double)));
    assert(!memcmp(a->phi, b->phi, 7 * a->n * sizeof(double)));
    assert(!memcmp(a->psi, b->psi, sizeof a->psi));
    assert(!memcmp(&a->statistics, &b->statistics, sizeof a->statistics));
}
int main(void) {
    for (int method = BI_VERIFIED; method <= BI_DASKR; ++method) {
        double atol[] = {1e-8, 2e-8}, rtol[] = {1e-7, 2e-7};
        TimeOptions options = method ? bi_daskr_options(1e-8, 1e-7) : bi_options(1e-8, 1e-7);
        options.atol = atol;
        options.rtol = rtol;
        options.atol_len = options.rtol_len = 2;
        double y[] = {1., .2}, yp[] = {-1.02, -.4};
        BdfIntegrator source, snapshot = {0};
        assert(!bi_new(&source, 0., y, yp, 2, &options));
        DaeSystem system = {.n = 2, .residual = residual, .jacobian = jacobian};
        StepOutcome step;
        assert(!bi_step(&source, &system, 1., 0, &step));
        assert(!bi_clone_from(&snapshot, &source));
        equal_state(&snapshot, &source);
        assert(snapshot.newton_cache && snapshot.newton_cache == source.newton_cache);
        assert(snapshot.matrix.values == source.matrix.values);
        assert(snapshot.factor.values == source.factor.values);
        assert(snapshot.y != source.y && !snapshot.step_work && !snapshot.newton_work);

        /* Reusing a checkpoint must not overwrite its independently allocated
         * state, tolerance or scratch-buffer storage. */
        assert(!bi_step(&snapshot, &system, 1., 0, &step));
        double *saved_y = snapshot.y, *saved_phi = snapshot.phi, *saved_atol = snapshot.atol_owned,
               *saved_rtol = snapshot.rtol_owned, *saved_step = snapshot.step_work,
               *saved_newton = snapshot.newton_work;
        for (size_t i = 0; i < 20; ++i) {
            assert(!bi_step(&source, &system, 1., 0, &step));
            assert(!bi_clone_from(&snapshot, &source));
            equal_state(&snapshot, &source);
            assert(snapshot.y == saved_y && snapshot.phi == saved_phi);
            assert(snapshot.atol_owned == saved_atol && snapshot.rtol_owned == saved_rtol);
            assert(snapshot.step_work == saved_step && snapshot.newton_work == saved_newton);
            assert(snapshot.newton_cache == source.newton_cache);
        }
        /* Invalidating/freeing either owner leaves the other's immutable
         * factors usable. A later refresh gets a different owner. */
        bi_invalidate_jacobian(&source);
        assert(snapshot.has_matrix && snapshot.newton_cache);
        bi_free(&source);
        assert(!bi_step(&snapshot, &system, 1., 0, &step));
        assert(snapshot.step_work == saved_step && snapshot.newton_work == saved_newton);

        BdfIntegrator clone = {0};
        assert(!bi_clone(&clone, &snapshot));
        double *old_y = snapshot.y;
        double old_time = snapshot.t;
        int calls = -10;
        TimeOptions projected_options =
            method ? bi_daskr_options(1e-8, 1e-7) : bi_options(1e-8, 1e-7);
        assert(bi_project_history(&snapshot, &clone, 3, &projected_options, project, &calls) == -1);
        assert(snapshot.y == old_y && snapshot.t == old_time && snapshot.n == 2);
        equal_state(&snapshot, &clone);
        calls = 0;
        assert(bi_project_history(&snapshot, &clone, 3, &projected_options, project, &calls) == 1);
        assert(calls > 0 && snapshot.n == 3 && !snapshot.has_matrix);
        assert(snapshot.step_work == NULL && snapshot.newton_work == NULL);
        const double *accepted = method ? clone.source_states : clone.y;
        assert(snapshot.y[0] == accepted[0] && snapshot.y[1] == accepted[1]);
        assert(snapshot.y[2] == accepted[0] + 2. * accepted[1]);
        assert(!bi_step(&clone, &system, 1., 0, &step));
        saved_step = clone.step_work;
        saved_newton = clone.newton_work;
        assert(!bi_reset(&clone, clone.t, clone.y, clone.yp, clone.n));
        assert(!clone.has_matrix && !clone.initialized);
        assert(clone.step_work == saved_step && clone.newton_work == saved_newton);
        assert(!bi_step(&clone, &system, 1., 0, &step));
        assert(clone.step_work == saved_step && clone.newton_work == saved_newton);
        bi_free(&snapshot);
        bi_free(&clone);
    }
    puts("BDF optimization: shared-cache lifetime, reusable snapshots/workspaces and transactional "
         "history changes passed");
    return 0;
}
