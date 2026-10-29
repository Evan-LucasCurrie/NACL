/* Continuous superconvergent extensions. */
#include "internal.h"
static size_t total(size_t k) {
    return k == 2 ? 2 : k == 3 ? 6 : k == 4 ? 9 : 0;
}
static const double x3[] = {
    -0.02496, 0.02304, 0.1608365339512893, 0.05205333333333333, -0.2109698672846227, 0.0,
};
static const double x4[] = {
    0.01259503554386166,
    -0.004194459027329296,
    0.007325540972670704,
    0.02790115799284125,
    -0.001606014582884851,
    0.009913985417115149,
    0.03142445823600003,
    0.140128203528661,
    0.02161038635693079,
    0.1278455645496104,
    -0.02647240969774761,
    0.08052540747398258,
    -0.02486766857825578,
    -0.1547489916166944,
    -0.04775117444496425,
    -0.1748985477440576,
    -0.03317969806100672,
    -0.1516975152327369,
    0.0,
    0.08007336945700194,
    0.08007336945700194,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
};
static const double coefficients[9][6] = {
    {
        -30.59376195485321,
        190.3161507410791,
        -506.6326577159306,
        677.9074502047026,
        -450.2127204516637,
        118.2155391766658,
    },
    {
        18.26042862151987,
        -122.3369840744124,
        316.7368243825973,
        -401.428283538036,
        249.6918871183304,
        -60.92387250999916,
    },
    {
        66.44072302985358,
        -465.9005706172424,
        1290.052895837388,
        -1760.013679899507,
        1182.484831564531,
        -312.8902724924552,
    },
    {
        58.05733615326762,
        -415.3918067568109,
        1181.297756609702,
        -1652.285224203471,
        1134.24023288076,
        -305.5922221060169,
    },
    {
        -28.71080418445304,
        182.2499138934506,
        -432.9611914049302,
        497.0100823644704,
        -278.2997171236686,
        61.03778903256205,
    },
    {
        -50.78725499866815,
        341.5424634806026,
        -890.8894610421601,
        1143.788821738507,
        -725.9253473216229,
        182.44470556591,
    },
    {
        -32.66666666666667,
        238.1228070175439,
        -700.6140350877193,
        1032.438596491228,
        -752.1929824561404,
        214.9122807017544,
    },
    {
        -79.11786699358921,
        588.9108452680835,
        -1685.448258473616,
        2349.773462964052,
        -1603.256196611642,
        429.1380138467126,
    },
    {
        79.11786699358921,
        -537.512818952294,
        1428.458126894669,
        -1887.191226121947,
        1243.470012401116,
        -326.3419612151337,
    },
};
static void weights(size_t k, double t, double *w) {
    if (k == 3) {
        double root = sqrt(15.0), t2 = t * t, t3 = t2 * t, tm = t - 1;
        w[0] = t * tm * tm * (5 * t2 - 10 * t + 3) / 3;
        w[1] = t2 * tm * (30 * t2 - 35 * t + 9) / 4;
        w[2] = (-5 + root) *
               (60 * t3 - 15 * root * t2 - 195 * t2 + 32 * t * root + 220 * t - 18 * root - 90) *
               t2 / 36;
        w[3] = -4.0 / 9 * t2 * (-18 + 62 * t - 75 * t2 + 30 * t3);
        w[4] = -(5 + root) *
               (60 * t3 - 195 * t2 + 15 * root * t2 + 220 * t - 32 * t * root + 18 * root - 90) *
               t2 / 36;
        w[5] = 125.0 / 12 * t2 * (2 * t - 1) * tm * tm;
    } else {
        for (size_t i = 0; i < 9; i++) {
            double a = 0;
            for (size_t j = 6; j-- > 0;)
                a = a * t + coefficients[i][j];
            w[i] = t * t * a;
        }
        w[0] += t;
    }
}
colnewsc_status cn_prepare(const colnewsc_problem *p, colnewsc_solution *s) {
    size_t k = s->basis.k, count = total(k), n = s->dimension, cells = s->count - 1;
    if (!count)
        return COLNEWSC_OK;
    if (cells > SIZE_MAX / (count * n))
        return COLNEWSC_ALLOCATION_FAILED;
    free(s->extra);
    s->extra = cn_alloc(cells * count * n, sizeof(double));
    if (!s->extra)
        return COLNEWSC_ALLOCATION_FAILED;
    double *y = cn_alloc(n, sizeof(double));
    if (!y)
        return COLNEWSC_ALLOCATION_FAILED;
    static const double c3[] = {0.6}, c4[] = {0.6889822365046136, 0.2, 0.8};
    const double *c = k == 3 ? c3 : c4, *x = k == 3 ? x3 : x4;
    size_t clen = k == 3 ? 1 : 3;
    colnewsc_status e = COLNEWSC_OK;
    for (size_t i = 0; i < cells; i++) {
        double h = s->mesh[i + 1] - s->mesh[i];
        size_t offset = i * count * n;
        e = cn_checked_rhs(p, s->mesh[i], s->values + i * n, s->extra + offset);
        if (e)
            break;
        e = cn_checked_rhs(p, s->mesh[i + 1], s->values + (i + 1) * n, s->extra + offset + n);
        if (e)
            break;
        if (k == 2)
            continue;
        memcpy(s->extra + offset + 2 * n, s->stages + i * k * n, k * n * sizeof(double));
        for (size_t r = 0; r < clen; r++) {
            for (size_t q = 0; q < n; q++) {
                y[q] = (1 - c[r]) * s->values[i * n + q] + c[r] * s->values[(i + 1) * n + q];
                for (size_t j = 0; j < k + 2 + r; j++)
                    y[q] += h * x[j * clen + r] * s->extra[offset + j * n + q];
            }
            size_t j = k + 2 + r;
            e = cn_checked_rhs(p, s->mesh[i] + h * c[r], y, s->extra + offset + j * n);
            if (e)
                break;
        }
        if (e)
            break;
    }
    free(y);
    return e;
}
void cn_evaluate_interval(const colnewsc_solution *s, size_t i, double x, double *out) {
    size_t k = s->basis.k, n = s->dimension, count = total(k), offset = i * count * n;
    double h = s->mesh[i + 1] - s->mesh[i], t = (x - s->mesh[i]) / h;
    memcpy(out, s->values + i * n, n * sizeof(double));
    if (count && s->extra) {
        if (k == 2) {
            double v = t * t * (3 - 2 * t), wl = t * (1 - t) * (1 - t), wr = t * t * (t - 1);
            for (size_t q = 0; q < n; q++)
                out[q] = (1 - v) * s->values[i * n + q] + v * s->values[(i + 1) * n + q] +
                         h * (wl * s->extra[offset + q] + wr * s->extra[offset + n + q]);
        } else {
            double w[9] = {0};
            weights(k, t, w);
            if (n == 1) {
                double value = out[0];
                for (size_t j = 0; j < count; j++)
                    value += h * w[j] * s->extra[offset + j];
                out[0] = value;
                return;
            }
            if (n == 2) {
                double first = out[0], second = out[1];
                for (size_t j = 0; j < count; j++) {
                    first += h * w[j] * s->extra[offset + 2 * j];
                    second += h * w[j] * s->extra[offset + 2 * j + 1];
                }
                out[0] = first;
                out[1] = second;
                return;
            }
            for (size_t j = 0; j < count; j++)
                for (size_t q = 0; q < n; q++)
                    out[q] += h * w[j] * s->extra[offset + j * n + q];
        }
    } else {
        double w[7];
        cn_basis_weights(&s->basis, t, w);
        if (n == 1) {
            double value = out[0];
            for (size_t j = 0; j < k; j++)
                value += h * w[j] * s->stages[i * k + j];
            out[0] = value;
            return;
        }
        if (n == 2) {
            double first = out[0], second = out[1];
            for (size_t j = 0; j < k; j++) {
                first += h * w[j] * s->stages[(i * k + j) * 2];
                second += h * w[j] * s->stages[(i * k + j) * 2 + 1];
            }
            out[0] = first;
            out[1] = second;
            return;
        }
        for (size_t j = 0; j < k; j++)
            for (size_t q = 0; q < n; q++)
                out[q] += h * w[j] * s->stages[(i * k + j) * n + q];
    }
}
