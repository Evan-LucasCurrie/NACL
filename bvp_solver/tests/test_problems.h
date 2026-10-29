#ifndef BVP_TEST_PROBLEMS_H
#define BVP_TEST_PROBLEMS_H
#include "bvp.h"
#include <math.h>
#include <string.h>
/* Analytic equations for solver regression tests. */
typedef struct {
    int kind, analytic;
} test_problem;
static void tp_rhs(void *v, double x, const double *y, const double *p, double *f) {
    test_problem *d = v;
    switch (d->kind) {
    case 0:
        f[0] = y[1];
        f[1] = y[0];
        break;
    case 1:
        f[0] = y[1];
        f[1] = (y[0] + y[0] * y[0] - exp(-20 * x)) / 0.01;
        break;
    case 2:
        f[0] = y[1];
        f[1] = -pow(y[0], 5);
        break;
    case 3:
        f[0] = p[0];
        break;
    case 4:
        f[0] = y[1];
        f[1] = -exp(y[0]);
        break;
    default:
        f[0] = 2;
    }
}
static void tp_boundary(void *v, const double *a, const double *b, const double *p, double *l,
                        double *r) {
    (void)p;
    test_problem *d = v;
    switch (d->kind) {
    case 0:
        l[0] = a[0] - 1;
        r[0] = b[0] - exp(1.0);
        break;
    case 1:
        l[0] = a[0] - 1;
        r[0] = b[0] - exp(-10);
        break;
    case 2:
        l[0] = a[1];
        r[0] = b[0] - sqrt(0.75);
        break;
    case 3:
        l[0] = a[0] - 1;
        r[0] = b[0] - 4;
        break;
    case 4:
        l[0] = a[0];
        r[0] = b[0];
        break;
    default:
        l[0] = a[0] - 1;
    }
}
static int tp_jac(void *v, double x, const double *y, const double *p, double *j) {
    (void)x;
    (void)p;
    test_problem *d = v;
    if (!d->analytic)
        return 0;
    if (d->kind == 3) {
        j[0] = 0;
        j[1] = 1;
        return 1;
    }
    if (d->kind == 5) {
        j[0] = 0;
        return 1;
    }
    j[0] = 0;
    j[1] = 1;
    j[3] = 0;
    j[2] = d->kind == 0   ? 1
           : d->kind == 1 ? (1 + 2 * y[0]) / 0.01
           : d->kind == 2 ? -5 * pow(y[0], 4)
                          : -exp(y[0]);
    return 1;
}
static int tp_bcjac(void *v, const double *a, const double *b, const double *p, double *l,
                    double *r) {
    (void)a;
    (void)b;
    (void)p;
    test_problem *d = v;
    if (!d->analytic)
        return 0;
    if (d->kind == 5) {
        l[0] = 1;
        return 1;
    }
    l[0] = d->kind == 2 ? 0 : 1;
    l[1] = d->kind == 2 ? 1 : 0;
    r[0] = 1;
    r[1] = 0;
    return 1;
}
static void tp_guess(void *v, double x, double *y) {
    test_problem *d = v;
    switch (d->kind) {
    case 0:
        y[0] = 1 + x * (exp(1.0) - 1);
        y[1] = 1;
        break;
    case 1:
        y[0] = 0.5;
        y[1] = 0;
        break;
    case 2:
        y[0] = sqrt(0.75);
        y[1] = 1e-4;
        break;
    case 3:
        y[0] = 1 + 3 * x;
        break;
    case 4:
        y[0] = y[1] = 0;
        break;
    default:
        y[0] = 1;
    }
}
static bvp_problem tp_make(test_problem *d) {
    static const double singular[] = {0, 0, 0, -2};
    bvp_problem p = {d->kind >= 3 && d->kind != 4 ? 1 : 2,
                     d->kind == 3 ? 1 : 0,
                     1,
                     0,
                     d,
                     tp_rhs,
                     tp_boundary,
                     tp_jac,
                     tp_bcjac,
                     d->kind == 2 ? singular : NULL};
    return p;
}
static double tp_exact(test_problem *d, double x) {
    if (d->kind == 0)
        return exp(x);
    if (d->kind == 1)
        return exp(-10 * x);
    if (d->kind == 2)
        return 1 / sqrt(1 + x * x / 3);
    if (d->kind == 3)
        return 1 + 3 * x;
    if (d->kind == 4) {
        double theta = 1.517164599050754;
        return -2 * log(cosh(theta * (x - 0.5) / 2) / cosh(theta / 4));
    }
    return 1 + 2 * x;
}
#endif
