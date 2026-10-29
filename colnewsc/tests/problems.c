/* Thirteen equations for solver regression tests. */
#include "internal.h"
#include "problems.h"
static const double pi = 3.14159265358979323846264338327950288;
static double sq(double x) {
    return x * x;
}
static double p3(double x) {
    return x * x * x;
}
static double p5(double x) {
    double y = x * x;
    return y * y * x;
}
static size_t dimension(colnewsc_test_case id) {
    switch (id) {
    case CN_SWIRLING_FLOW:
    case CN_MEASLES:
    case CN_THVE16_EX1:
        return 6;
    case CN_TMLPDE:
        return 20;
    case CN_SIWA22_EX1:
    case CN_SNK14_EX1:
        return 3;
    case CN_SNK14_EX2:
    case CN_KHWA14_EX2:
    case CN_KHWA14_EX4:
        return 4;
    default:
        return 2;
    }
}
static void rhs(void *ctx, double x, const double *y, double *f) {
    colnewsc_test_case id = (colnewsc_test_case)(uintptr_t)ctx;
    memset(f, 0, dimension(id) * sizeof(double));
    switch (id) {
    case CN_SWIRLING_FLOW:
        f[0] = y[1];
        f[1] = (y[0] * y[3] - y[2] * y[1]) / 0.005;
        f[2] = y[3];
        f[3] = y[4];
        f[4] = y[5];
        f[5] = (-y[2] * y[5] - y[0] * y[1]) / 0.005;
        break;
    case CN_TP20:
        f[0] = y[1];
        f[1] = (1 - y[1] * y[1]) / 0.05;
        break;
    case CN_TP21:
        f[0] = y[1];
        f[1] = (y[0] + y[0] * y[0] - exp(-20 * x)) / 0.01;
        break;
    case CN_TP4:
        f[0] = y[1];
        f[1] = (1.01 * y[0] - y[1]) / 0.01;
        break;
    case CN_MEASLES: {
        double bt = 1575 * (1 + cos(2 * pi * x));
        f[0] = 0.02 - bt * y[0] * y[2];
        f[1] = bt * y[0] * y[2] - y[1] / 0.0279;
        f[2] = y[1] / 0.0279 - y[2] / 0.01;
        break;
    }
    case CN_TMLPDE: {
        double step = 1.0 / 11, c = cos(20 * x), sn = sin(20 * x);
        for (size_t i = 0; i < 10; i++) {
            double a = (i + 1) * step;
            f[2 * i] = y[2 * i + 1];
            f[2 * i + 1] = y[2 * i] * y[2 * i + 1] + (y[2 * i] - (i ? y[2 * i - 2] : 0)) / step -
                           c - a * 400 * c + a * a * 20 * c * sn;
        }
        break;
    }
    case CN_SHOCK_WAVE: {
        double c1 = 1 / (0.05 * (1 + x * x)), c2 = 1.2 - 0.1 * x;
        f[0] = y[1];
        f[1] = c1 / y[0] *
               (c2 * y[0] * y[1] - y[1] / y[0] - 2 * x / (1 + x * x) * (1 - 0.2 * y[0] * y[0]));
        break;
    }
    case CN_THVE16_EX1: {
        double a = 46441.0 / 14400 - 537.0 / 50 * exp(-1) + 29851.0 / 2400 * exp(-2) +
                   63.0 / 50 * exp(-3) + 49.0 / 1600 * exp(-4);
        f[0] = y[1];
        f[1] = -a * (2 * exp(x) - 1 + 4 * x - x * x) * exp(-x) / (sq(0.5 * y[3] + 0.1 * y[5]) + 1);
        f[2] = y[0];
        f[4] = y[1] * y[1];
        break;
    }
    case CN_SIWA22_EX1:
        f[0] = y[1];
        f[1] = x == 0 ? 0 : -0.5 / x * y[1] + 12 * sq(p3(x)) * p5(y[0]) - 7 * x * x * p3(y[0]);
        f[2] = x / 2 * y[0];
        break;
    case CN_SNK14_EX1:
        f[0] = y[1];
        f[1] = 1 + y[2];
        f[2] = exp(-x) * y[0] * y[0];
        break;
    case CN_SNK14_EX2:
        f[0] = y[1];
        f[1] = -2 * x * x - p3(x) / 6 - 1 / sq(4 + x) + x * y[2] - y[3];
        f[2] = exp(y[0]);
        f[3] = x * exp(y[0]);
        break;
    case CN_KHWA14_EX2:
        f[0] = y[1];
        f[1] = -0.75 * cos(2 * pi * x / 3) / sq(y[3]);
        f[2] = y[0];
        break;
    case CN_KHWA14_EX4: {
        double c = 3 / (4 * (2 * sqrt(2) - 2));
        f[0] = y[1];
        f[1] = c * y[3] * p5(y[0]);
        f[2] = y[0];
        break;
    }
    default:
        break;
    }
}
static double boundary(void *ctx, size_t i, const double *y) {
    colnewsc_test_case id = (colnewsc_test_case)(uintptr_t)ctx;
    switch (id) {
    case CN_SWIRLING_FLOW:
        return i == 0 ? y[0] + 1 : i == 3 ? y[0] - 1 : i == 1 || i == 4 ? y[2] : y[3];
    case CN_TP20:
        return y[0] - 1 - 0.05 * log(cosh(i == 0 ? -0.745 / 0.05 : 0.255 / 0.05));
    case CN_TP21:
        return y[0] - (i == 0 ? 1 : exp(-10));
    case CN_TP4:
        return y[0] - 1 - (i == 0 ? exp(-2) : exp(-202));
    case CN_MEASLES:
        return y[i % 3] - y[i % 3 + 3];
    case CN_TMLPDE:
        return y[2 * (i % 10)] - (double)(i % 10 + 1) / 11 * (i < 10 ? 1 : cos(20));
    case CN_SHOCK_WAVE:
        return y[0] - (i == 0 ? (double)0.9129f : 0.375);
    case CN_THVE16_EX1:
        return i == 0 || i == 3 ? y[0]
               : i == 1         ? y[2]
               : i == 2         ? y[4]
               : i == 4         ? y[2] - y[3]
                                : y[4] - y[5];
    case CN_SIWA22_EX1:
        return i == 0   ? y[0] - 0.5
               : i == 1 ? y[2]
                        : y[0] - y[2] - (1 / sqrt(5) - 0.25 * asinh(0.5));
    case CN_SNK14_EX1:
        return i == 0 ? y[0] - 1 : i == 1 ? y[2] : y[0] - exp(1);
    case CN_SNK14_EX2:
        return i == 0 ? y[0] - log(4) : i == 1 ? y[2] : i == 2 ? y[3] : y[0] - log(5);
    case CN_KHWA14_EX2:
        return i == 0 ? y[0] - 1 : i == 1 ? y[2] : i == 2 ? y[0] + 0.5 : y[2] - y[3];
    case CN_KHWA14_EX4:
        return i == 0 ? y[0] - 1 : i == 1 ? y[2] : i == 2 ? y[0] - sqrt(2) / 2 : y[2] - y[3];
    default:
        return NAN;
    }
}
static void guess(void *ctx, double x, double *y) {
    colnewsc_test_case id = (colnewsc_test_case)(uintptr_t)ctx;
    size_t d = dimension(id);
    memset(y, 0, d * sizeof(double));
    switch (id) {
    case CN_SWIRLING_FLOW:
        y[0] = x;
        y[1] = 1;
        break;
    case CN_TP20:
    case CN_TP21:
    case CN_TP4:
    case CN_MEASLES:
    case CN_TMLPDE:
        for (size_t j = 0; j < d; j++)
            y[j] = 1;
        break;
    case CN_SHOCK_WAVE:
        y[0] = (double)0.9129f - (double)0.5379f * x;
        break;
    case CN_THVE16_EX1:
    case CN_KHWA14_EX4:
        y[0] = 1 + x;
        y[1] = 1;
        y[2] = x + x * x / 2;
        y[3] = 1.5;
        if (id == CN_THVE16_EX1) {
            y[4] = y[2];
            y[5] = 1.5;
        }
        break;
    case CN_SIWA22_EX1:
        y[0] = 0.5 - (double)0.05f * x;
        y[1] = -(double)0.05f;
        y[2] = 0.5 * x - (double)0.025f * x * x;
        break;
    case CN_SNK14_EX1:
    case CN_SNK14_EX2:
        y[0] = x + 1;
        y[1] = 1;
        y[2] = x * x / 2 + x;
        if (id == CN_SNK14_EX2)
            y[3] = x * y[2];
        break;
    case CN_KHWA14_EX2:
        y[0] = 1 - x;
        y[1] = -1;
        y[2] = x - x * x / 2;
        y[3] = 0.5;
        break;
    default:
        break;
    }
}
const char *colnewsc_test_name(colnewsc_test_case id) {
    static const char *const names[] = {
        "swf_iii",    "tp20",       "tp21",      "measles",   "tmlpde",     "swave",     "tp4",
        "ThVe16-Ex1", "SiWa22-Ex1", "SNK14-Ex1", "SNK14-Ex2", "KhWa14-Ex2", "KhWa14-Ex4"};
    return (unsigned)id < CN_TEST_CASE_COUNT ? names[id] : "unknown";
}
colnewsc_problem colnewsc_test_problem(colnewsc_test_case id) {
    colnewsc_problem p = {0};
    if ((unsigned)id >= CN_TEST_CASE_COUNT)
        return p;
    p.dimension = dimension(id);
    p.left_boundary_conditions =
        id == CN_SWIRLING_FLOW || id == CN_MEASLES || id == CN_THVE16_EX1 || id == CN_SNK14_EX2 ? 3
        : id == CN_TMLPDE                                                                       ? 10
        : id == CN_SIWA22_EX1 || id == CN_SNK14_EX1 || id == CN_KHWA14_EX2 || id == CN_KHWA14_EX4
            ? 2
            : 1;
    p.a = id == CN_TP4 ? -1 : 0;
    p.b = 1;
    p.context = (void *)(uintptr_t)id;
    p.rhs = rhs;
    p.boundary = boundary;
    p.initial_guess = guess;
    return p;
}
colnewsc_options colnewsc_test_options(colnewsc_test_case id) {
    colnewsc_options o = colnewsc_default_options();
    o.tolerance = id >= CN_THVE16_EX1 ? 1e-6 : 1e-8;
    o.collocation_points = id == CN_TP21                   ? 2
                           : id == CN_TP20 || id == CN_TP4 ? 7
                           : id >= CN_THVE16_EX1           ? 4
                                                           : 3;
    return o;
}
double colnewsc_test_exact(colnewsc_test_case id, double x) {
    switch (id) {
    case CN_TP4:
        return exp(x - 1) + exp(-101 * (x + 1));
    case CN_TP20:
        return 1 + 0.05 * log(cosh((x - 0.745) / 0.05));
    case CN_TP21:
        return exp(-10 * x);
    case CN_TMLPDE:
        return cos(20 * x) / 11;
    case CN_THVE16_EX1:
        return (1 - x * x) * (1 - exp(-x));
    case CN_SIWA22_EX1:
        return 1 / sqrt(4 + sq(sq(x)));
    case CN_SNK14_EX1:
        return exp(x);
    case CN_SNK14_EX2:
        return log(4 + x);
    case CN_KHWA14_EX2:
        return cos(2 * pi * x / 3);
    case CN_KHWA14_EX4:
        return 1 / sqrt(1 + x);
    default:
        return NAN;
    }
}
