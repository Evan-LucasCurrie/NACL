/* Reference basis, tableau and continuous interpolation routines. */
#include "reference_internal.h"
/* Keep the original streaming behavior when a private caller deliberately
 * aliases the accumulation target with one of its source views. */
static int disjoint(R a, int na, R b, int nb) {
    if (na <= 0 || nb <= 0)
        return 1;
    uintptr_t x = (uintptr_t)a, y = (uintptr_t)b;
    return x <= y ? y - x >= (size_t)na * sizeof(double) : x - y >= (size_t)nb * sizeof(double);
}
void cnref_appsln(cnref_context *c, double x, R z, R fspace, I ispace) {
    double a_storage[28] = {0};
    R a = a_storage;
    double dummy_storage[1] = {0};
    R dummy = dummy_storage;
    c->imaxe = 314;
    int s_super;
    int s_star_super;
    if (GET(ispace, 2) == 2) {
        s_super = 4;
        s_star_super = 4;
    } else if (GET(ispace, 2) == 3) {
        s_super = 5;
        s_star_super = 6;
    } else if (GET(ispace, 2) == 4) {
        s_super = 6;
        s_star_super = 9;
    } else if (1) {
        s_super = 0;
        s_star_super = 0;
    }
    int is6;
    is6 = GET(ispace, 6);
    int is5;
    is5 = (GET(ispace, 1) + 2);
    int is4;
    is4 = (is5 + (GET(ispace, 4) * ((GET(ispace, 1) + 1))));
    int i;
    i = 1;
    if ((GET(ispace, 2) >= 2) && (GET(ispace, 2) <= 4)) {
        cnref_setup_phase1(c, GET(ispace, 1), VIEW(fspace, 1), a, VIEW(fspace, is5),
                           VIEW(fspace, is4), GET(ispace, 2), GET(ispace, 3), GET(ispace, 5),
                           VIEW(ispace, 8), GET(ispace, 4), c->y_super, c->sk_discrete_super,
                           c->sk_interp_super);
        cnref_evaluate(c, GET(ispace, 4), GET(ispace, 1), VIEW(fspace, 1), c->y_super, x, s_super,
                       s_star_super, z, c->sk_discrete_super, c->sk_interp_super);
    } else {
        {
            cnref_approx(c, &(i), &(x), z, a, VIEW(fspace, is6), VIEW(fspace, 1), GET(ispace, 1),
                         VIEW(fspace, is5), VIEW(fspace, is4), GET(ispace, 2), GET(ispace, 3),
                         GET(ispace, 5), VIEW(ispace, 8), GET(ispace, 4), 2, dummy, 0);
        }
    }
}
void cnref_approx(cnref_context *c, int *i, double *x, R zval, R a, R coef, R xi, int n, R z, R dmz,
                  int k, int ncomp, int mmax, I m, int mstar, int mode, R dmval, int modm) {
    double dm_storage[7] = {0};
    R dm = dm_storage;
    double bm_storage[4] = {0};
    R bm = bm_storage;
    if (mode == 1) {
        (*x) = GET(xi, (*i));
        int iz;
        iz = ((((*i) - 1)) * mstar);
        for (int j = 1; j <= mstar; j += 1) {
            (iz += 1);
            SET(zval, j, GET(z, iz));
        }
        return;
    }
    if (mode == 2) {
        if (!((((*x) >= (GET(xi, 1) - c->precis)) && ((*x) <= (GET(xi, (n + 1)) + c->precis))))) {
            if (c->iprint < 1) {
                fprintf(stderr,
                        " ****** DOMAIN ERROR IN APPROX ****** X = %.17g, ALEFT = %.17g, ARIGHT = "
                        "%.17g\n",
                        *x, GET(xi, 1), GET(xi, n + 1));
            }
            if ((*x) < GET(xi, 1)) {
                (*x) = GET(xi, 1);
            }
            if ((*x) > GET(xi, (n + 1))) {
                (*x) = GET(xi, (n + 1));
            }
        }
        if (((*i) > n) || ((*i) < 1)) {
            (*i) = (((n + 1)) / 2);
        }
        int ileft;
        ileft = (*i);
        if ((*x) >= GET(xi, ileft)) {
            for (int l = ileft; l <= n; l += 1) {
                (*i) = l;
                if ((*x) < GET(xi, (l + 1))) {
                    break;
                }
            }
        } else {
            {
                int iright;
                iright = (ileft - 1);
                for (int l = 1; l <= iright; l += 1) {
                    (*i) = ((iright + 1) - l);
                    if ((*x) >= GET(xi, (*i))) {
                        break;
                    }
                }
            }
        }
    }
    if ((mode == 2) || (mode == 3)) {
        double s;
        s = ((((*x) - GET(xi, (*i)))) / ((GET(xi, ((*i) + 1)) - GET(xi, (*i)))));
        cnref_rkbas(c, s, coef, k, mmax, a, dm, modm);
    }
    SET(bm, 1, ((*x) - GET(xi, (*i))));
    for (int l = 2; l <= mmax; l += 1) {
        SET(bm, l, (GET(bm, 1) / ((double)(l))));
    }
    int ir;
    ir = 1;
    int iz;
    iz = (((((*i) - 1)) * mstar) + 1);
    int idmz;
    idmz = (((((*i) - 1)) * k) * ncomp);
    for (int jcomp = 1; jcomp <= ncomp; jcomp += 1) {
        int mj;
        mj = GET(m, jcomp);
        (ir += mj);
        (iz += mj);
        for (int l = 1; l <= mj; l += 1) {
            int ind;
            ind = (idmz + jcomp);
            double zsum;
            zsum = 0.0;
            for (int j = 1; j <= k; j += 1) {
                (zsum += (GET2(a, j, l, 7) * GET(dmz, ind)));
                (ind += ncomp);
            }
            for (int ll = 1; ll <= l; ll += 1) {
                int lb;
                lb = ((l + 1) - ll);
                zsum = ((zsum * GET(bm, lb)) + GET(z, (iz - ll)));
            }
            SET(zval, (ir - l), zsum);
        }
    }
    if (modm == 0) {
        return;
    }
    for (int jcomp = 1; jcomp <= ncomp; jcomp += 1) {
        SET(dmval, jcomp, 0.0);
    }
    (idmz += 1);
    for (int j = 1; j <= k; j += 1) {
        double fact;
        fact = GET(dm, j);
        for (int jcomp = 1; jcomp <= ncomp; jcomp += 1) {
            SET(dmval, jcomp, (GET(dmval, jcomp) + (fact * GET(dmz, idmz))));
            (idmz += 1);
        }
    }
}
void cnref_rkbas(cnref_context *c, double s, R coef, int k, int m, R rkb, R dm, int mode) {
    (void)c;
    if (k == 1) {
        SET2(rkb, 1, 1, 7, 1.0);
        SET(dm, 1, 1.0);
        return;
    }
    double t_storage[10] = {0};
    R t = t_storage;
    for (int i = 1; i <= ((k + m) - 1); i += 1) {
        SET(t, i, (s / ((double)(i))));
    }
    for (int l = 1; l <= m; l += 1) {
        int lb;
        lb = ((k + l) + 1);
        for (int i = 1; i <= k; i += 1) {
            double p;
            p = GET2(coef, 1, i, k);
            for (int j = 2; j <= k; j += 1) {
                p = ((p * GET(t, (lb - j))) + GET2(coef, j, i, k));
            }
            SET2(rkb, i, l, 7, p);
        }
    }
    if (mode == 0) {
        return;
    }
    for (int i = 1; i <= k; i += 1) {
        double p;
        p = GET2(coef, 1, i, k);
        for (int j = 2; j <= k; j += 1) {
            p = ((p * GET(t, ((k + 1) - j))) + GET2(coef, j, i, k));
        }
        SET(dm, i, p);
    }
}
void cnref_vmonde(cnref_context *c, R rho, R coef, int k) {
    (void)c;
    if (k == 1) {
        return;
    }
    int km1;
    km1 = (k - 1);
    for (int i = 1; i <= km1; i += 1) {
        int kmi;
        kmi = (k - i);
        for (int j = 1; j <= kmi; j += 1) {
            SET(coef, j,
                (((GET(coef, (j + 1)) - GET(coef, j))) / ((GET(rho, (j + i)) - GET(rho, j)))));
        }
    }
    int ifac;
    ifac = 1;
    for (int i = 1; i <= km1; i += 1) {
        int kmi;
        kmi = ((k + 1) - i);
        for (int j = 2; j <= kmi; j += 1) {
            SET(coef, j, (GET(coef, j) - (GET(rho, ((j + i) - 1)) * GET(coef, (j - 1)))));
        }
        SET(coef, kmi, (((double)(ifac)) * GET(coef, kmi)));
        (ifac *= i);
    }
    SET(coef, 1, (((double)(ifac)) * GET(coef, 1)));
}
void cnref_horder(cnref_context *c, int i, R uhigh, double hi, R dmz, int ncomp, int k) {
    double dn;
    dn = (1.0 / cnref_powi(hi, (k - 1)));
    for (int id = 1; id <= ncomp; id += 1) {
        SET(uhigh, id, 0.0);
    }
    int kin;
    kin = 1;
    int idmz;
    idmz = (((((i - 1)) * k) * ncomp) + 1);
    for (int _j = 1; _j <= k; _j += 1) {
        double fact;
        fact = (dn * GET(c->coef, kin));
        for (int id = 1; id <= ncomp; id += 1) {
            SET(uhigh, id, (GET(uhigh, id) + (fact * GET(dmz, idmz))));
            (idmz += 1);
        }
        (kin += k);
    }
}
void cnref_evaluate(cnref_context *c, int neqns, int nsub, R mesh, R y, double t, int s, int s_star,
                    R z, R k_discrete, R k_interp) {
    int i;
    i = 0;
    double weights_storage[10] = {0};
    R weights = weights_storage;
    cnref_interval(c, nsub, mesh, t, &(i));
    double hi;
    hi = (GET(mesh, (i + 1)) - GET(mesh, i));
    double tau;
    tau = (((t - GET(mesh, i))) / hi);
    cnref_col_weights(c, s_star, tau, weights);
    cnref_sum_stages(c, neqns, hi, VIEW(y, ((((i - 1)) * neqns) + 1)), s,
                     VIEW(k_discrete, (((((i - 1)) * s) * neqns) + 1)), s_star,
                     VIEW(k_interp, (((((i - 1)) * ((s_star - s))) * neqns) + 1)), weights, z);
}
void cnref_setup(cnref_context *c, int _node, int neqns, I _orders, double tim1, double hi, R yim1,
                 R yi, int s, R ki_discrete, int s_star, R ki_interp, R y) {
    (void)_node;
    (void)_orders;
    if (neqns <= 2 && disjoint(y, neqns, yim1, neqns) && disjoint(y, neqns, yi, neqns) &&
        disjoint(y, neqns, ki_discrete, s * neqns) &&
        disjoint(y, neqns, ki_interp, (s_star - s) * neqns)) {
        for (int r = 0; r < s_star - s; r++) {
            double values[2] = {0, 0};
            for (int j = 0; j < s; j++) {
                double a = c->x_star[j * (s_star - s) + r];
                if (a != 0)
                    for (int i = 0; i < neqns; i++)
                        values[i] += a * ki_discrete[j * neqns + i];
            }
            for (int j = 0; j < r; j++) {
                double a = c->x_star[(j + s) * (s_star - s) + r];
                if (a != 0)
                    for (int i = 0; i < neqns; i++)
                        values[i] += a * ki_interp[j * neqns + i];
            }
            for (int i = 0; i < neqns; i++)
                values[i] = hi * values[i];
            double a = 1.0 - c->v_star[r];
            if (a != 0)
                for (int i = 0; i < neqns; i++)
                    values[i] += a * yim1[i];
            a = c->v_star[r];
            if (a != 0)
                for (int i = 0; i < neqns; i++)
                    values[i] += a * yi[i];
            for (int i = 0; i < neqns; i++)
                y[i] = values[i];
            cnref_fsub_interpolant(c, tim1 + c->c_star[r] * hi, y, ki_interp + r * neqns);
        }
        return;
    }
    double zero_storage[1] = {0};
    R zero = zero_storage;
    for (int r = 1; r <= (s_star - s); r += 1) {
        cnref_dcopy(c, neqns, zero, 0, y, 1);
        for (int j = 1; j <= s; j += 1) {
            cnref_daxpy(c, neqns, GET(c->x_star, ((((j - 1)) * ((s_star - s))) + r)),
                        VIEW(ki_discrete, ((((j - 1)) * neqns) + 1)), 1, y, 1);
        }
        for (int j = 1; j <= (r - 1); j += 1) {
            cnref_daxpy(c, neqns, GET(c->x_star, (((((j + s) - 1)) * ((s_star - s))) + r)),
                        VIEW(ki_interp, ((((j - 1)) * neqns) + 1)), 1, y, 1);
        }
        cnref_dscal(c, neqns, hi, y, 1);
        cnref_daxpy(c, neqns, (1.0 - GET(c->v_star, r)), yim1, 1, y, 1);
        cnref_daxpy(c, neqns, GET(c->v_star, r), yi, 1, y, 1);
        cnref_fsub_interpolant(c, (tim1 + (GET(c->c_star, r) * hi)), y,
                               VIEW(ki_interp, ((((r - 1)) * neqns) + 1)));
    }
}
void cnref_sum_stages(cnref_context *c, int neqns, double hi, R yim1, int s, R ki_discrete,
                      int s_star, R ki_interp, R weights, R z) {
    if (neqns <= 2 && disjoint(z, neqns, yim1, neqns) &&
        disjoint(z, neqns, ki_discrete, s * neqns) &&
        disjoint(z, neqns, ki_interp, (s_star - s) * neqns) &&
        disjoint(z, neqns, weights, s_star)) {
        double values[2] = {0, 0};
        for (int j = 0; j < s; j++) {
            double a = weights[j];
            if (a != 0)
                for (int i = 0; i < neqns; i++)
                    values[i] += a * ki_discrete[j * neqns + i];
        }
        for (int j = 0; j < s_star - s; j++) {
            double a = weights[s + j];
            if (a != 0)
                for (int i = 0; i < neqns; i++)
                    values[i] += a * ki_interp[j * neqns + i];
        }
        for (int i = 0; i < neqns; i++) {
            values[i] = hi * values[i];
            z[i] = values[i] + 1.0 * yim1[i];
        }
        return;
    }
    cnref_dcopy(c, neqns, (double[1]){0}, 0, z, 1);
    for (int j = 1; j <= s; j += 1) {
        cnref_daxpy(c, neqns, GET(weights, j), VIEW(ki_discrete, ((((j - 1)) * neqns) + 1)), 1, z,
                    1);
    }
    for (int j = 1; j <= (s_star - s); j += 1) {
        cnref_daxpy(c, neqns, GET(weights, (s + j)), VIEW(ki_interp, ((((j - 1)) * neqns) + 1)), 1,
                    z, 1);
    }
    cnref_dscal(c, neqns, hi, z, 1);
    cnref_daxpy(c, neqns, 1.0, yim1, 1, z, 1);
}
void cnref_fsub_first(cnref_context *c, int node, int _nsol, I orders, R z, R f, R f_first) {
    (void)c;
    (void)_nsol;
    SET(f_first, GET(orders, 1), GET(f, 1));
    for (int i = 1; i <= (GET(orders, 1) - 1); i += 1) {
        SET(f_first, i, GET(z, (i + 1)));
    }
    int orders_sum;
    orders_sum = GET(orders, 1);
    for (int j = 2; j <= node; j += 1) {
        (orders_sum += GET(orders, j));
        SET(f_first, orders_sum, GET(f, j));
        for (int i = ((orders_sum - GET(orders, j)) + 1); i <= (orders_sum - 1); i += 1) {
            SET(f_first, i, GET(z, (i + 1)));
        }
    }
}
void cnref_setup_phase1(cnref_context *c, int nsub, R mesh, R _asave, R zval, R dmz, int k_points,
                        int ncomp, int mmax, I m, int mstar, R y_super, R sk_discrete_super,
                        R sk_interp_super) {
    (void)_asave;
    int s;
    s = (k_points + 2);
    int s_star;
    if (k_points == 2) {
        s_star = (k_points + 2);
    } else if (k_points == 3) {
        s_star = ((k_points + 2) + 1);
    } else if (k_points == 4) {
        s_star = ((k_points + 2) + 3);
    } else if (1) {
        s_star = 0;
    }
    cnref_tableau_star(c, k_points);
    double xc_storage[4] = {0};
    R xc = xc_storage;
    if (k_points == 2) {
        {
            SET(xc, 1, 0.2113248654051871);
            SET(xc, 2, 0.7886751345948129);
        }
    } else if (k_points == 3) {
        {
            SET(xc, 1, 0.1127016653792583);
            SET(xc, 2, 0.5);
            SET(xc, 3, 0.8872983346207417);
        }
    } else if (k_points == 4) {
        {
            SET(xc, 1, 0.6943184420297371e-1);
            SET(xc, 2, 0.3300094782075719);
            SET(xc, 3, 0.6699905217924281);
            SET(xc, 4, 0.9305681557970263);
        }
    } else if (1) {
        {
        }
    }
    double z_storage[160] = {0};
    R z = z_storage;
    double f_storage[50] = {0};
    R f = f_storage;
    double dummy_storage[1] = {0};
    R dummy = dummy_storage;
    double yy_storage[160] = {0};
    R yy = yy_storage;
    double a_storage[28] = {0};
    R a = a_storage;
    for (int i = 1; i <= (nsub + 1); i += 1) {
        double x;
        x = GET(mesh, i);
        for (int j = 1; j <= mstar; j += 1) {
            SET(z, j, GET(zval, ((((i - 1)) * mstar) + j)));
            SET(y_super, ((((i - 1)) * mstar) + j), GET(z, j));
        }
        cnref_fsub_interpolant(c, x, z, f);
        if (i != 1) {
            for (int k = 1; k <= mstar; k += 1) {
                SET(sk_discrete_super, ((((((i - 2)) * s) * mstar) + mstar) + k), GET(f, k));
            }
        }
        if (i != (nsub + 1)) {
            for (int k = 1; k <= mstar; k += 1) {
                SET(sk_discrete_super, (((((i - 1)) * s) * mstar) + k), GET(f, k));
            }
        }
    }
    for (int ii = 1; ii <= nsub; ii += 1) {
        int i;
        i = ii;
        double hi;
        hi = (GET(mesh, (i + 1)) - GET(mesh, i));
        for (int k = 1; k <= k_points; k += 1) {
            double x;
            x = (GET(mesh, i) + (hi * GET(xc, k)));
            cnref_approx(c, &(i), &(x), z, a, c->coef, mesh, nsub, zval, dmz, k_points, ncomp, mmax,
                         m, mstar, 3, dummy, 0);
            cnref_fsub_interpolant(c, x, z, f);
            for (int l = 1; l <= mstar; l += 1) {
                SET(sk_discrete_super,
                    ((((((i - 1)) * s) * mstar) + ((((2 + k) - 1)) * mstar)) + l), GET(f, l));
            }
        }
    }
    if (k_points != 2) {
        for (int i = 1; i <= nsub; i += 1) {
            cnref_setup(
                c, ncomp, mstar, m, GET(mesh, i), (GET(mesh, (i + 1)) - GET(mesh, i)),
                VIEW(zval, ((((i - 1)) * mstar) + 1)), VIEW(zval, ((i * mstar) + 1)), s,
                VIEW(sk_discrete_super, (((((i - 1)) * mstar) * ((k_points + 2))) + 1)), s_star,
                VIEW(sk_interp_super, (((((i - 1)) * mstar) * (((s_star - k_points) - 2))) + 1)),
                yy);
        }
    }
}
void cnref_tableau_star(cnref_context *c, int method) {
    if (method == 2) {
        c->s_star = 4;
    }
    if (method == 3) {
        c->s_star = 6;
        SET(c->c_star, 1, 0.6e0);
        SET(c->v_star, 1, 0.6e0);
        SET(c->x_star, 1, (-0.2496e-1));
        SET(c->x_star, 2, 0.2304e-1);
        SET(c->x_star, 3, 0.1608365339512893e0);
        SET(c->x_star, 4, 0.5205333333333333e-1);
        SET(c->x_star, 5, (-0.2109698672846227e0));
        SET(c->x_star, 6, 0.0e0);
    }
    if (method == 4) {
        c->s_star = 9;
        SET(c->c_star, 1, 0.6889822365046136e0);
        SET(c->c_star, 2, (1.0e0 / 5.0e0));
        SET(c->c_star, 3, (4.0e0 / 5.0e0));
        SET(c->v_star, 1, 0.6889822365046136e0);
        SET(c->v_star, 2, (1.0e0 / 5.0e0));
        SET(c->v_star, 3, (4.0e0 / 5.0e0));
        SET(c->x_star, 1, 0.1259503554386166e-1);
        SET(c->x_star, 2, (-0.4194459027329296e-2));
        SET(c->x_star, 3, 0.7325540972670704e-2);
        SET(c->x_star, 4, 0.2790115799284125e-1);
        SET(c->x_star, 5, (-0.1606014582884851e-2));
        SET(c->x_star, 6, 0.9913985417115149e-2);
        SET(c->x_star, 7, 0.3142445823600003e-1);
        SET(c->x_star, 8, 0.140128203528661e0);
        SET(c->x_star, 9, 0.2161038635693079e-1);
        SET(c->x_star, 10, 0.1278455645496104e0);
        SET(c->x_star, 11, (-0.2647240969774761e-1));
        SET(c->x_star, 12, 0.8052540747398258e-1);
        SET(c->x_star, 13, (-0.2486766857825578e-1));
        SET(c->x_star, 14, (-0.1547489916166944e0));
        SET(c->x_star, 15, (-0.4775117444496425e-1));
        SET(c->x_star, 16, (-0.1748985477440576e0));
        SET(c->x_star, 17, (-0.3317969806100672e-1));
        SET(c->x_star, 18, (-0.1516975152327369e0));
        SET(c->x_star, 19, 0.0e0);
        SET(c->x_star, 20, 0.8007336945700194e-1);
        SET(c->x_star, 21, 0.8007336945700194e-1);
        SET(c->x_star, 22, 0.0e0);
        SET(c->x_star, 23, 0.0e0);
        SET(c->x_star, 24, 0.0e0);
        SET(c->x_star, 25, 0.0e0);
        SET(c->x_star, 26, 0.0e0);
        SET(c->x_star, 27, 0.0e0);
    }
}
void cnref_col_weights(cnref_context *c, int s_star, double t, R w) {
    if (s_star == 4) {
        if (c->imaxe != 314) {
            if (t < 0.2e0) {
                SET(w, 1, 0.9435739743742127e-1);
                SET(w, 2, (-0.1764804844828675e-2));
                SET(w, 3, 0.7407407407407407e-1);
                SET(w, 4, 0.0e0);
            }
            if ((t < 0.5e0) && (t > 0.2e0)) {
                SET(w, 1, 0.7330644658652999e-1);
                SET(w, 2, 0.7676274875440806e-3);
                SET(w, 3, 0.2592592592592593e0);
                SET(w, 4, 0.0e0);
            }
            if ((t > 0.5e0) && (t < 0.8e0)) {
                SET(w, 1, (-0.1397593589591207e0));
                SET(w, 2, 0.6568528488504658e-1);
                SET(w, 3, 0.7407407407407407e0);
                SET(w, 4, 0.0e0);
            }
            if (t > 0.8e0) {
                SET(w, 1, (-0.2441436431433453e0));
                SET(w, 2, 0.1515510505507527e0);
                SET(w, 3, 0.9259259259259259e0);
                SET(w, 4, 0.0e0);
            }
        }
        if (c->imaxe == 314) {
            SET(w, 1,
                (((((3.0 + sqrt((3.0e0)))) *
                   ((((((4.0 * cnref_powi(t, 2)) - (9.0 * t)) + (t * sqrt((3.0e0)))) + 6.0) -
                     (2.0 * sqrt((3.0e0)))))) *
                  t) /
                 12.0));
            SET(w, 2,
                (((((3.0 - sqrt((3.0e0)))) * ((((4.0 * t) - 3.0) + sqrt((3.0e0))))) *
                  cnref_powi(t, 2)) /
                 12.0));
            SET(w, 3, (cnref_powi(t, 2) * ((3.0 - (2.0 * t)))));
            SET(w, 4, 0.0);
        }
    }
    if (s_star == 6) {
        if (c->imaxe != 314) {
            if (t < 0.2e0) {
                SET(w, 1, 0.5679869684499314e-1);
                SET(w, 2, (-0.2314814814814815e-1));
                SET(w, 3, 0.9491267653789125e-1);
                SET(w, 4, 0.1186556927297668e0);
                SET(w, 5, 0.5340693937431725e-1);
                SET(w, 6, (-0.1339591906721536e0));
            }
            if ((t < 0.5e0) && (t > 0.2e0)) {
                SET(w, 1, 0.1097393689986283e-1);
                SET(w, 2, (-0.1234567901234568e-1));
                SET(w, 3, 0.2468582536682419e0);
                SET(w, 4, 0.2249657064471879e0);
                SET(w, 5, 0.3434887939074302e-1);
                SET(w, 6, (-0.1714677640603567e0));
            }
            if ((t > 0.5e0) && (t < 0.8e0)) {
                SET(w, 1, (-0.3566529492455418e-1));
                SET(w, 2, 0.3703703703703704e-1);
                SET(w, 3, 0.3496835855257842e0);
                SET(w, 4, 0.2194787379972565e0);
                SET(w, 5, (-0.7533516302921355e-1));
                SET(w, 6, 0.1714677640603567e0);
            }
            if (t > 0.8e0) {
                SET(w, 1, (-0.1436042524005487e-1));
                SET(w, 2, (-0.1929012345679012e-1));
                SET(w, 3, 0.3073823127306085e0);
                SET(w, 4, 0.3257887517146776e0);
                SET(w, 5, 0.9985362691273853e-1);
                SET(w, 6, 0.1339591906721536e0);
            }
        }
        if (c->imaxe == 314) {
            SET(w, 1,
                (((t * cnref_powi(((t - 1.0)), 2)) *
                  ((((5.0 * cnref_powi(t, 2)) - (10.0 * t)) + 3.0))) /
                 3.0));
            SET(w, 2,
                (((cnref_powi(t, 2) * ((t - 1.0))) *
                  ((((30.0 * cnref_powi(t, 2)) - (35.0 * t)) + 9.0))) /
                 4.0));
            SET(w, 3,
                ((((((-5.0) + sqrt((15.0e0)))) *
                   ((((((((60.0 * cnref_powi(t, 3)) -
                          ((15.0 * sqrt((15.0e0))) * cnref_powi(t, 2))) -
                         (195.0 * cnref_powi(t, 2))) +
                        ((32.0 * t) * sqrt((15.0e0)))) +
                       (220.0 * t)) -
                      (18.0 * sqrt((15.0e0)))) -
                     90.0))) *
                  cnref_powi(t, 2)) /
                 36.0));
            SET(w, 4,
                ((((-4.0e0) / 9.0e0) * cnref_powi(t, 2)) *
                 (((((-18.0) + (62.0 * t)) - (75.0 * cnref_powi(t, 2))) +
                   (30.0 * cnref_powi(t, 3))))));
            SET(w, 5,
                ((((-((5.0 + sqrt((15.0e0))))) *
                   ((((((((60.0 * cnref_powi(t, 3)) - (195.0 * cnref_powi(t, 2))) +
                         ((15.0 * sqrt((15.0e0))) * cnref_powi(t, 2))) +
                        (220.0 * t)) -
                       ((32.0 * t) * sqrt((15.0e0)))) +
                      (18.0 * sqrt((15.0e0)))) -
                     90.0))) *
                  cnref_powi(t, 2)) /
                 36.0));
            SET(w, 6,
                ((((125.0e0 / 12.0e0) * cnref_powi(t, 2)) * (((2.0 * t) - 1.0))) *
                 cnref_powi(((t - 1.0)), 2)));
        }
    }
    if (s_star == 9) {
        if (c->imaxe != 314) {
            if (t < 0.2e0) {
                SET(w, 1, (-0.1150348945591317e0));
                SET(w, 2, 0.1387650209689777e0);
                SET(w, 3, 0.4819275638059541e0);
                SET(w, 4, 0.4118225626662212e0);
                SET(w, 5, (-0.2296787116499782e0));
                SET(w, 6, (-0.3847718606383842e0));
                SET(w, 7, (-0.2281661873591156e0));
                SET(w, 8, (-0.5024271855143026e0));
                SET(w, 9, 0.5942303589464258e0);
            }
            if ((t < 0.5e0) && (t > 0.2e0)) {
                SET(w, 1, (-0.4573080643257935e-1));
                SET(w, 2, 0.7095561972323692e-1);
                SET(w, 3, 0.2894461668828698e0);
                SET(w, 4, 0.266482498906871e0);
                SET(w, 5, (-0.9382952923818732e-1));
                SET(w, 6, (-0.1973529088423723e0));
                SET(w, 7, (-0.1446506068555011e0));
                SET(w, 8, (-0.1205466853346099e0));
                SET(w, 9, 0.3085595845236023e0);
            }
            if ((t > 0.5e0) && (t < 0.8e0)) {
                SET(w, 1, 0.4961740908461741e-1);
                SET(w, 2, (-0.7484222237526824e-1));
                SET(w, 3, 0.3788614407585321e-1);
                SET(w, 4, 0.1731351307224674e0);
                SET(w, 5, 0.3063570544713461e0);
                SET(w, 6, 0.2178754430209323e0);
                SET(w, 7, 0.1446506068554992e0);
                SET(w, 8, 0.1999299094367667e0);
                SET(w, 9, (-0.3879428086256833e0));
            }
            if (t > 0.8e0) {
                SET(w, 1, 0.5142214021528287e-1);
                SET(w, 2, (-0.7515226662512127e-1));
                SET(w, 3, 0.3281638326198989e-1);
                SET(w, 4, 0.1665107508495575e0);
                SET(w, 5, 0.3034905529965711e0);
                SET(w, 6, 0.21788275870745e0);
                SET(w, 7, 0.2281661873590951e0);
                SET(w, 8, 0.2078410023235988e0);
                SET(w, 9, (-0.2996441757555419e0));
            }
        }
        if (c->imaxe == 314) {
            double t2;
            t2 = (t * t);
            double t3;
            t3 = (t * t2);
            double t5;
            t5 = (t3 * t2);
            double t6;
            t6 = (t3 * t3);
            SET(w, 1,
                (((((((1.0e0 * t) - (0.3059376195485321e2 * cnref_powi(t, 2))) +
                     (0.1903161507410791e3 * cnref_powi(t, 3))) -
                    (0.5066326577159306e3 * cnref_powi(t, 4))) +
                   (0.6779074502047026e3 * t5)) -
                  (0.4502127204516637e3 * t6)) +
                 (0.1182155391766658e3 * cnref_powi(t, 7))));
            SET(w, 2,
                ((((((0.1826042862151987e2 * cnref_powi(t, 2)) -
                     (0.1223369840744124e3 * cnref_powi(t, 3))) +
                    (0.3167368243825973e3 * cnref_powi(t, 4))) -
                   (0.401428283538036e3 * t5)) +
                  (0.2496918871183304e3 * t6)) -
                 (0.6092387250999916e2 * cnref_powi(t, 7))));
            SET(w, 3,
                ((((((0.6644072302985358e2 * cnref_powi(t, 2)) -
                     (0.4659005706172424e3 * cnref_powi(t, 3))) +
                    (0.1290052895837388e4 * cnref_powi(t, 4))) -
                   (0.1760013679899507e4 * t5)) +
                  (0.1182484831564531e4 * t6)) -
                 (0.3128902724924552e3 * cnref_powi(t, 7))));
            SET(w, 4,
                ((((((0.5805733615326762e2 * cnref_powi(t, 2)) -
                     (0.4153918067568109e3 * cnref_powi(t, 3))) +
                    (0.1181297756609702e4 * cnref_powi(t, 4))) -
                   (0.1652285224203471e4 * t5)) +
                  (0.113424023288076e4 * t6)) -
                 (0.3055922221060169e3 * cnref_powi(t, 7))));
            SET(w, 5,
                (((((((-0.2871080418445304e2) * cnref_powi(t, 2)) +
                     (0.1822499138934506e3 * cnref_powi(t, 3))) -
                    (0.4329611914049302e3 * cnref_powi(t, 4))) +
                   (0.4970100823644704e3 * t5)) -
                  (0.2782997171236686e3 * t6)) +
                 (0.6103778903256205e2 * cnref_powi(t, 7))));
            SET(w, 6,
                (((((((-0.5078725499866815e2) * cnref_powi(t, 2)) +
                     (0.3415424634806026e3 * cnref_powi(t, 3))) -
                    (0.8908894610421601e3 * cnref_powi(t, 4))) +
                   (0.1143788821738507e4 * t5)) -
                  (0.7259253473216229e3 * t6)) +
                 (0.18244470556591e3 * cnref_powi(t, 7))));
            SET(w, 7,
                (((((((-0.3266666666666667e2) * cnref_powi(t, 2)) +
                     (0.2381228070175439e3 * cnref_powi(t, 3))) -
                    (0.7006140350877193e3 * cnref_powi(t, 4))) +
                   (0.1032438596491228e4 * t5)) -
                  (0.7521929824561404e3 * t6)) +
                 (0.2149122807017544e3 * cnref_powi(t, 7))));
            SET(w, 8,
                (((((((-0.7911786699358921e2) * cnref_powi(t, 2)) +
                     (0.5889108452680835e3 * cnref_powi(t, 3))) -
                    (0.1685448258473616e4 * cnref_powi(t, 4))) +
                   (0.2349773462964052e4 * t5)) -
                  (0.1603256196611642e4 * t6)) +
                 (0.4291380138467126e3 * cnref_powi(t, 7))));
            SET(w, 9,
                ((((((0.7911786699358921e2 * cnref_powi(t, 2)) -
                     (0.537512818952294e3 * cnref_powi(t, 3))) +
                    (0.1428458126894669e4 * cnref_powi(t, 4))) -
                   (0.1887191226121947e4 * t5)) +
                  (0.1243470012401116e4 * t6)) -
                 (0.3263419612151337e3 * cnref_powi(t, 7))));
        }
    }
}
