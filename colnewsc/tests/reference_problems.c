/* Callback variants for reference solver regression tests. */
#include "reference/reference_internal.h"
#include "problems.h"
static double power5(double v) {
    double v2 = v * v, v3 = v * v2;
    return v3 * v2;
}
static double power6(double v) {
    double v3 = (v * v) * v;
    return v3 * v3;
}
static void reference_guess(void *ctx, double x, double *z, double *dmz) {
    colnewsc_test_case id = (colnewsc_test_case)(uintptr_t)ctx;
    colnewsc_problem base = colnewsc_test_problem(id);
    base.initial_guess(ctx, x, z);
    memset(dmz, 0, base.dimension * sizeof(double));
    switch (id) {
    case CN_SIWA22_EX1:
        z[2] = 0.5 * x - (double)(float)0.025 * cnref_powi(x, 2);
        break;
    case CN_THVE16_EX1:
    case CN_KHWA14_EX4:
        z[2] = x + cnref_powi(x, 2) / 2;
        if (id == CN_THVE16_EX1)
            z[4] = z[2];
        break;
    case CN_SNK14_EX1:
    case CN_SNK14_EX2:
        z[2] = cnref_powi(x, 2) / 2 + x;
        if (id == CN_SNK14_EX2)
            z[3] = x * z[2];
        break;
    case CN_KHWA14_EX2:
        z[2] = x - cnref_powi(x, 2) / 2;
        break;
    default:
        break;
    }
}
static void reference_rhs(void *ctx, double x, const double *z, double *f) {
    colnewsc_test_case id = (colnewsc_test_case)(uintptr_t)ctx;
    colnewsc_problem base = colnewsc_test_problem(id);
    base.rhs(ctx, x, z, f);
    switch (id) {
    case CN_TP21:
        f[1] = (z[0] + cnref_powi(z[0], 2) - exp(-2 * x / sqrt(0.01))) / 0.01;
        break;
    case CN_SHOCK_WAVE: {
        double gamma = 1.4, cons1 = 1 / (0.05 * (1 + cnref_powi(x, 2))),
               cons2 = 1.0 / 2 + gamma / 2 - 2 * 0.05 * x;
        f[1] = (cons1 / z[0]) *
               (cons2 * z[0] * z[1] - z[1] / z[0] -
                (2 * x / (1 + cnref_powi(x, 2))) * (1 - 0.5 * (gamma - 1) * cnref_powi(z[0], 2)));
        break;
    }
    case CN_TMLPDE: {
        double k = 1.0 / (10 + 1), omega = 20;
        f[1] = z[0] * z[1] + (z[0] - 0) / k - cos(omega * x) -
               k * cnref_powi(omega, 2) * cos(omega * x) +
               cnref_powi(k, 2) * omega * cos(omega * x) * sin(omega * x);
        for (int i = 2; i <= 10; i++) {
            double ik = i * k;
            f[2 * i - 1] = z[2 * i - 2] * z[2 * i - 1] + (z[2 * i - 2] - z[2 * i - 4]) / k -
                           cos(omega * x) - ik * cnref_powi(omega, 2) * cos(omega * x) +
                           cnref_powi(ik, 2) * omega * cos(omega * x) * sin(omega * x);
        }
        break;
    }
    case CN_SIWA22_EX1:
        f[1] = -0.5 / x * z[1] + 12 * power6(x) * power5(z[0]) -
               2 * (3 + 0.5) * cnref_powi(x, 2) * cnref_powi(z[0], 3);
        break;
    case CN_KHWA14_EX4: {
        double c = 3 / (4 * (2 * sqrt(2) - 2));
        f[1] = c * z[3] * power5(z[0]);
        break;
    }
    case CN_SNK14_EX1:
        f[2] = exp(-x) * cnref_powi(z[0], 2);
        break;
    case CN_SNK14_EX2:
        f[1] = (-2 * cnref_powi(x, 2) - cnref_powi(x, 3) / 6 - 1 / cnref_powi(4 + x, 2)) +
               x * z[2] - z[3];
        break;
    default:
        break;
    }
}
static double reference_boundary(void *ctx, size_t i, const double *z) {
    colnewsc_test_case id = (colnewsc_test_case)(uintptr_t)ctx;
    switch (id) {
    case CN_TP21:
        return z[0] - (i == 0 ? 1 : exp(-1 / sqrt(0.01)));
    case CN_TP4:
        return z[0] - 1 - (i == 0 ? exp(-2) : exp(-2 * (1 + 0.01) / 0.01));
    case CN_TMLPDE: {
        double k = 1.0 / (10 + 1);
        size_t j = i % 10 + 1;
        return z[2 * (j - 1)] - (j * k) * (i < 10 ? 1 : cos(20));
    }
    default: {
        colnewsc_problem base = colnewsc_test_problem(id);
        return base.boundary(ctx, i, z);
    }
    }
}
colnewsc_reference_problem colnewsc_reference_test_problem(colnewsc_test_case id) {
    colnewsc_reference_problem p = {0};
    p.base = colnewsc_test_problem(id);
    p.reference_guess = reference_guess;
    p.reference_rhs = reference_rhs;
    p.reference_boundary = reference_boundary;
    p.preserve_non_finite_interpolation = id == CN_SIWA22_EX1;
    return p;
}
