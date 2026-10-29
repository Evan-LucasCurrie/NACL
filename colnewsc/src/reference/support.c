/* Native C versions of COLNEWSC reference BLAS/LINPACK and checked helpers.
   The original operation order, first-maximum pivot rule and streaming alias
   behavior are retained. No external BLAS or Fortran runtime is required. */
#include "reference_internal.h"
int cnref_all_finite(cnref_context *c, R a, int n) {
    for (int i = 0; i < n; i++)
        if (!cnref_mesh_finite(c, a[i]))
            return 0;
    return 1;
}
int cnref_all_pair_finite(cnref_context *c, R a, R b, int n) {
    for (int i = 0; i < n; i++)
        if (!cnref_mesh_finite(c, a[i]) || !cnref_mesh_finite(c, b[i]))
            return 0;
    return 1;
}
int cnref_check_mesh(cnref_context *c, R xi, int n) {
    for (int i = 0; i < n; i++)
        if (!cnref_mesh_finite(c, xi[i]) || !cnref_mesh_finite(c, xi[i + 1]) || xi[i] >= xi[i + 1])
            return 0;
    return 1;
}
int cnref_plain_finite(R a, int n) {
    for (int i = 0; i < n; i++)
        if (!isfinite(a[i]))
            return 0;
    return 1;
}
double cnref_mesh_density(cnref_context *c, double d1, double d2, int j, double oneovh, double z) {
    double difference = fabs(d2 - d1), weighted = difference * GET(c->wgtmsh, j),
           scaled = weighted * oneovh, denominator = 1 + fabs(z), quotient = scaled / denominator,
           density = pow(quotient, GET(c->root, j));
    if (!(cnref_mesh_finite(c, difference) && cnref_mesh_finite(c, weighted) &&
          cnref_mesh_finite(c, scaled) && cnref_mesh_finite(c, denominator) &&
          cnref_mesh_finite(c, quotient) && cnref_mesh_finite(c, density)))
        return NAN;
    return density;
}
static void axpy_cells(size_t n, double t, const double *x, double *y) {
    if (!n || t == 0)
        return;
    size_t m = n % 4;
    for (size_t i = 0; i < m; i++)
        y[i] += t * x[i];
    for (size_t i = m; i < n; i += 4) {
        y[i] += t * x[i];
        y[i + 1] += t * x[i + 1];
        y[i + 2] += t * x[i + 2];
        y[i + 3] += t * x[i + 3];
    }
}
static void scal_cells(size_t n, double t, double *x) {
    size_t m = n % 5;
    for (size_t i = 0; i < m; i++)
        x[i] = t * x[i];
    for (size_t i = m; i < n; i += 5) {
        x[i] = t * x[i];
        x[i + 1] = t * x[i + 1];
        x[i + 2] = t * x[i + 2];
        x[i + 3] = t * x[i + 3];
        x[i + 4] = t * x[i + 4];
    }
}
static double dot_cells(size_t n, const double *x, const double *y) {
    double result = 0;
    size_t m = n % 5;
    for (size_t i = 0; i < m; i++)
        result += x[i] * y[i];
    for (size_t i = m; i < n; i += 5)
        result = result + x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] +
                 x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    return result;
}
void cnref_dgefa(cnref_context *c, R a, int lda, int n, I ipvt, int *info) {
    (void)c;
    size_t rows = (size_t)lda, count = (size_t)n;
    *info = 0;
    for (size_t k = 0; k + 1 < count; k++) {
        size_t kp1 = k + 1, l = k;
        double maximum = fabs(a[k + k * rows]);
        for (size_t i = kp1; i < count; i++) {
            if (fabs(a[i + k * rows]) <= maximum)
                continue;
            l = i;
            maximum = fabs(a[i + k * rows]);
        }
        ipvt[k] = (int)l + 1;
        if (a[l + k * rows] == 0) {
            *info = (int)k + 1;
            continue;
        }
        if (l != k) {
            double t = a[l + k * rows];
            a[l + k * rows] = a[k + k * rows];
            a[k + k * rows] = t;
        }
        double t = -1 / a[k + k * rows];
        scal_cells(count - k - 1, t, a + kp1 + k * rows);
        for (size_t j = kp1; j < count; j++) {
            double t = a[l + j * rows];
            if (l != k) {
                a[l + j * rows] = a[k + j * rows];
                a[k + j * rows] = t;
            }
            axpy_cells(count - k - 1, t, a + kp1 + k * rows, a + kp1 + j * rows);
        }
    }
    ipvt[count - 1] = n;
    if (a[count - 1 + (count - 1) * rows] == 0)
        *info = n;
}
void cnref_dgesl(cnref_context *c, R a, int lda, int n, I ipvt, R b, int job) {
    (void)c;
    size_t rows = (size_t)lda, count = (size_t)n;
    if (job == 0) {
        for (size_t k = 0; k + 1 < count; k++) {
            size_t l = (size_t)ipvt[k] - 1;
            double t = b[l];
            if (l != k) {
                b[l] = b[k];
                b[k] = t;
            }
            axpy_cells(count - k - 1, t, a + k + 1 + k * rows, b + k + 1);
        }
        for (size_t kb = 0; kb < count; kb++) {
            size_t k = count - 1 - kb;
            b[k] /= a[k + k * rows];
            double t = -b[k];
            axpy_cells(k, t, a + k * rows, b);
        }
    } else {
        for (size_t k = 0; k < count; k++) {
            double t = dot_cells(k, a + k * rows, b);
            b[k] = (b[k] - t) / a[k + k * rows];
        }
        for (size_t kb = 0; kb + 1 < count; kb++) {
            size_t k = count - 2 - kb;
            b[k] += dot_cells(count - k - 1, a + k + 1 + k * rows, b + k + 1);
            size_t l = (size_t)ipvt[k] - 1;
            if (l != k) {
                double t = b[l];
                b[l] = b[k];
                b[k] = t;
            }
        }
    }
}
void cnref_daxpy(cnref_context *c, int n, double da, R dx, int incx, R dy, int incy) {
    (void)c;
    if (n <= 0 || da == 0)
        return;
    if (incx == 1 && incy == 1) {
        axpy_cells((size_t)n, da, dx, dy);
        return;
    }
    int ix = incx < 0 ? (-n + 1) * incx : 0, iy = incy < 0 ? (-n + 1) * incy : 0;
    for (int k = 0; k < n; k++) {
        dy[iy] += da * dx[ix];
        ix += incx;
        iy += incy;
    }
}
double cnref_ddot(cnref_context *c, int n, R dx, int incx, R dy, int incy) {
    (void)c;
    if (n <= 0)
        return 0;
    if (incx == 1 && incy == 1)
        return dot_cells((size_t)n, dx, dy);
    double result = 0;
    int ix = incx < 0 ? (-n + 1) * incx : 0, iy = incy < 0 ? (-n + 1) * incy : 0;
    for (int k = 0; k < n; k++) {
        result += dx[ix] * dy[iy];
        ix += incx;
        iy += incy;
    }
    return result;
}
void cnref_dscal(cnref_context *c, int n, double da, R dx, int incx) {
    (void)c;
    if (n <= 0)
        return;
    if (incx == 1) {
        scal_cells((size_t)n, da, dx);
        return;
    }
    int nincx = n * incx, i = 1;
    while (incx >= 0 ? i <= nincx : i >= nincx) {
        dx[i - 1] = da * dx[i - 1];
        i += incx;
    }
}
int cnref_idamax(cnref_context *c, int n, R dx, int incx) {
    (void)c;
    if (n < 1)
        return 0;
    int index = 1, ix = 1;
    double dmax = fabs(dx[0]);
    if (incx != 1)
        ix += incx;
    for (int i = 2; i <= n; i++) {
        int j = incx == 1 ? i : ix;
        if (fabs(dx[j - 1]) > dmax) {
            index = i;
            dmax = fabs(dx[j - 1]);
        }
        if (incx != 1)
            ix += incx;
    }
    return index;
}
void cnref_dcopy(cnref_context *c, int n, R dx, int incx, R dy, int incy) {
    (void)c;
    if (n <= 0)
        return;
    if (incx == 1 && incy == 1) {
        int m = n % 7;
        for (int i = 0; i < m; i++)
            dy[i] = dx[i];
        for (int i = m; i < n; i += 7)
            for (int j = 0; j < 7; j++)
                dy[i + j] = dx[i + j];
        return;
    }
    int ix = incx < 0 ? (-n + 1) * incx : 0, iy = incy < 0 ? (-n + 1) * incy : 0;
    for (int k = 0; k < n; k++) {
        dy[iy] = dx[ix];
        ix += incx;
        iy += incy;
    }
}
void cnref_interval(cnref_context *c, int nsub, R mesh, double t, int *i) {
    if ((t < mesh[0] || t > mesh[nsub]) && c->iprint < 1)
        fprintf(stderr, "ERROR-T IS OUTSIDE OF PROBLEM INTERVAL\n");
    if (t == mesh[nsub]) {
        *i = nsub;
        return;
    }
    if (t == mesh[0]) {
        *i = 1;
        return;
    }
    int inew = nsub / 2;
    double meshhalfnew = mesh[inew];
    if (t == meshhalfnew) {
        *i = inew + 1;
        return;
    }
    int iold = inew, ileft, iright;
    if (t < meshhalfnew) {
        while (t < meshhalfnew) {
            iold = inew;
            inew = iold / 2;
            meshhalfnew = mesh[inew];
        }
        iright = iold;
        ileft = inew;
    } else {
        while (t > meshhalfnew) {
            iold = inew;
            inew = (int)roundf((float)(iold + (nsub - iold) / 2) + 0.5f);
            meshhalfnew = mesh[inew];
        }
        ileft = iold;
        iright = inew;
    }
    *i = iright;
    if (t == meshhalfnew) {
        *i = inew + 1;
        return;
    }
    while (abs(iright - ileft) != 1) {
        inew = (ileft + iright) / 2;
        if (mesh[inew] < t)
            ileft = inew;
        else if (mesh[inew] > t)
            iright = inew;
        else {
            *i = inew + 1;
            return;
        }
    }
    *i = iright;
}
