#ifndef COLNEWSC_WORKLOADS_H
#define COLNEWSC_WORKLOADS_H
#include "colnewsc.h"
#include <math.h>
#include <stdint.h>
static void workload_rhs(void *ctx, double x, const double *y, double *f) {
    (void)x;
    int id = (int)(uintptr_t)ctx;
    if (id == 0)
        f[0] = y[0];
    else {
        f[0] = y[1];
        f[1] = id == 1 ? -y[0] : -exp(y[0]);
    }
}
static double workload_bc(void *ctx, size_t i, const double *y) {
    int id = (int)(uintptr_t)ctx;
    return y[0] - (id == 0 ? 1 : id == 1 && i ? sin(1) : 0);
}
static void workload_guess(void *ctx, double x, double *y) {
    int id = (int)(uintptr_t)ctx;
    if (id == 0)
        y[0] = 1 + x;
    else {
        y[0] = id == 1 ? x * sin(1) : 0;
        y[1] = id == 1 ? sin(1) : 0;
    }
}
static int workload_jac(void *ctx, double x, const double *y, double *j) {
    (void)x;
    int id = (int)(uintptr_t)ctx;
    if (!id)
        j[0] = 1;
    else {
        j[0] = 0;
        j[1] = 1;
        j[2] = id == 1 ? -1 : -exp(y[0]);
        j[3] = 0;
    }
    return 1;
}
static int workload_grad(void *ctx, size_t i, const double *y, double *j) {
    (void)i;
    (void)y;
    int id = (int)(uintptr_t)ctx;
    j[0] = 1;
    if (id)
        j[1] = 0;
    return 1;
}
static colnewsc_problem workload_problem(int id) {
    colnewsc_problem p = {0};
    p.dimension = id ? 2 : 1;
    p.left_boundary_conditions = 1;
    p.a = 0;
    p.b = 1;
    p.context = (void *)(uintptr_t)id;
    p.rhs = workload_rhs;
    p.boundary = workload_bc;
    p.initial_guess = workload_guess;
    p.jacobian = workload_jac;
    p.boundary_gradient = workload_grad;
    return p;
}
static double workload_exact(int id, double x) {
    if (id == 0)
        return exp(x);
    if (id == 1)
        return sin(x);
    double k = 0.7585822995253773;
    return -2 * log(cosh(k * (x - 0.5)) / cosh(k / 2));
}
#endif
