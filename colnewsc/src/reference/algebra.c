/* Block assembly, elimination, LINPACK and BLAS numerical kernels. */
#include "reference_internal.h"
void cnref_lsyslv(cnref_context *c, int *msing, R xi, R xiold, R z, R dmz, R delz, R deldmz, R g,
                  R w, R v, R rhs, R dmzo, I integs, I ipvtg, I ipvtw, double *rnorm, int mode) {
    double zval_storage[40] = {0};
    R zval = zval_storage;
    double f_storage[40] = {0};
    R f = f_storage;
    double dgz_storage[40] = {0};
    R dgz = dgz_storage;
    double dmval_storage[20] = {0};
    R dmval = dmval_storage;
    double df_storage[800] = {0};
    R df = df_storage;
    double at_storage[28] = {0};
    R at = at_storage;
    double dummy_storage[1] = {0};
    R dummy = dummy_storage;
    if (mode != 4) {
        if (mode == 0) {
            for (int i = 1; i <= c->mstar; i += 1) {
                SET(zval, i, 0.0);
            }
        }
        int idmz;
        idmz = 1;
        int idmzo;
        idmzo = 1;
        int irhs;
        irhs = 1;
        int ig;
        ig = 1;
        int iw;
        iw = 1;
        int iv;
        iv = 1;
        c->izeta = 1;
        int lside;
        lside = 0;
        int iold;
        iold = 1;
        int ncol;
        ncol = (2 * c->mstar);
        (*rnorm) = 0.0;
        if (mode <= 1) {
            for (int i = 1; i <= c->n; i += 1) {
                SET2(integs, 2, i, 3, ncol);
                if (i == c->n) {
                    SET2(integs, 3, c->n, 3, ncol);
                    lside = c->mstar;
                } else {
                    {
                        SET2(integs, 3, i, 3, c->mstar);
                        while ((lside != c->mstar)) {
                            if (GET(c->zeta, (lside + 1)) >= (GET(xi, i) + c->precis)) {
                                break;
                            }
                            (lside += 1);
                        }
                    }
                }
                SET2(integs, 1, i, 3, (c->mstar + lside));
            }
        }
        if (mode != 2) {
            for (int l = 1; l <= ((c->kd * c->kd) * c->n); l += 1) {
                SET(w, l, 0.0);
            }
        }
        for (int i = 1; i <= c->n; i += 1) {
            double xii;
            xii = GET(xi, i);
            double h;
            h = (GET(xi, (i + 1)) - GET(xi, i));
            int nrow;
            nrow = GET2(integs, 1, i, 3);
            while (((c->izeta <= c->mstar) && (GET(c->zeta, c->izeta) <= (xii + c->precis)))) {
                int skip_rhs;
                skip_rhs = 0;
                if (mode != 0) {
                    if (c->iguess == 1) {
                        cnref_guess(c, xii, zval, dmval);
                    } else {
                        if (mode == 1) {
                            cnref_approx(c, &(iold), &(xii), zval, at, c->coef, xiold, c->nold, z,
                                         dmz, c->k, c->ncomp, c->mmax, c->m, c->mstar, 2, dummy, 0);
                        } else {
                            {
                                int ii;
                                ii = i;
                                cnref_approx(c, &(ii), &(xii), zval, at, dummy, xi, c->n, z, dmz,
                                             c->k, c->ncomp, c->mmax, c->m, c->mstar, 1, dummy, 0);
                                skip_rhs = (mode == 3);
                            }
                        }
                    }
                }
                if (!skip_rhs) {
                    double gval;
                    gval = cnref_gsub(c, c->izeta, zval);
                    SET(rhs, (c->ndmz + c->izeta), (-gval));
                    ((*rnorm) += cnref_powi(gval, 2));
                }
                if (mode != 2) {
                    cnref_gderiv(c, VIEW(g, ig), nrow, c->izeta, zval, dgz, 1);
                }
                (c->izeta += 1);
            }
            for (int j = 1; j <= c->k; j += 1) {
                double hrho;
                hrho = (h * GET(c->rho, j));
                double xcol;
                xcol = (xii + hrho);
                int build;
                if (mode == 0) {
                    cnref_fsub(c, xcol, zval, VIEW(rhs, irhs));
                    (irhs += c->ncomp);
                    build = 1;
                } else {
                    if ((c->iguess == 1) || (mode == 1)) {
                        if (c->iguess == 1) {
                            cnref_guess(c, xcol, zval, VIEW(dmzo, irhs));
                        } else {
                            {
                                cnref_approx(c, &(iold), &(xcol), zval, at, c->coef, xiold, c->nold,
                                             z, dmz, c->k, c->ncomp, c->mmax, c->m, c->mstar, 2,
                                             VIEW(dmzo, irhs), 1);
                            }
                        }
                        cnref_fsub(c, xcol, zval, f);
                        for (int jj = 1; jj <= c->ncomp; jj += 1) {
                            double value;
                            value = (GET(dmzo, irhs) - GET(f, jj));
                            SET(rhs, irhs, (-value));
                            ((*rnorm) += cnref_powi(value, 2));
                            (irhs += 1);
                        }
                        build = 1;
                    } else {
                        {
                            int ii;
                            ii = i;
                            cnref_approx(c, &(ii), &(xcol), zval,
                                         VIEW(c->acol, (1 + (((j - 1)) * 28))), c->coef, xi, c->n,
                                         z, dmz, c->k, c->ncomp, c->mmax, c->m, c->mstar, 4, dummy,
                                         0);
                            build = (mode == 3);
                            if (!build) {
                                cnref_fsub(c, xcol, zval, f);
                                for (int jj = 1; jj <= c->ncomp; jj += 1) {
                                    double value;
                                    value = (GET(dmz, irhs) - GET(f, jj));
                                    SET(rhs, irhs, (-value));
                                    ((*rnorm) += cnref_powi(value, 2));
                                    (irhs += 1);
                                }
                            }
                        }
                    }
                }
                if (build) {
                    cnref_vwblok(c, xcol, hrho, j, VIEW(w, iw), VIEW(v, iv), VIEW(ipvtw, idmz),
                                 c->kd, zval, df, VIEW(c->acol, (1 + (((j - 1)) * 28))),
                                 VIEW(dmzo, idmzo), c->ncomp, msing);
                    if ((*msing) != 0) {
                        return;
                    }
                }
            }
            if (mode != 2) {
                cnref_gblock(c, h, VIEW(g, ig), nrow, c->izeta, VIEW(w, iw), VIEW(v, iv), c->kd,
                             dummy, VIEW(deldmz, idmz), VIEW(ipvtw, idmz), 1);
            }
            if (i == c->n) {
                c->idum = c->izeta;
                while ((c->izeta <= c->mstar)) {
                    int skip_rhs;
                    skip_rhs = 0;
                    if (mode != 0) {
                        if (c->iguess == 1) {
                            cnref_guess(c, c->aright, zval, dmval);
                        } else {
                            if (mode == 1) {
                                int ii;
                                ii = (c->nold + 1);
                                double x;
                                x = c->aright;
                                cnref_approx(c, &(ii), &(x), zval, at, c->coef, xiold, c->nold, z,
                                             dmz, c->k, c->ncomp, c->mmax, c->m, c->mstar, 1, dummy,
                                             0);
                                c->aright = x;
                            } else {
                                {
                                    int ii;
                                    ii = (c->n + 1);
                                    double x;
                                    x = c->aright;
                                    cnref_approx(c, &(ii), &(x), zval, at, c->coef, xi, c->n, z,
                                                 dmz, c->k, c->ncomp, c->mmax, c->m, c->mstar, 1,
                                                 dummy, 0);
                                    c->aright = x;
                                    skip_rhs = (mode == 3);
                                }
                            }
                        }
                    }
                    if (!skip_rhs) {
                        double gval;
                        gval = cnref_gsub(c, c->izeta, zval);
                        SET(rhs, (c->ndmz + c->izeta), (-gval));
                        ((*rnorm) += cnref_powi(gval, 2));
                    }
                    if (mode != 2) {
                        cnref_gderiv(c, VIEW(g, ig), nrow, (c->izeta + c->mstar), zval, dgz, 2);
                    }
                    (c->izeta += 1);
                }
            } else {
                {
                    (ig += (nrow * ncol));
                    (iv += (c->kd * c->mstar));
                    (iw += (c->kd * c->kd));
                    (idmz += c->kd);
                    if (mode == 1) {
                        (idmzo += c->kd);
                    }
                }
            }
        }
        if ((mode != 0) && (mode != 3)) {
            (*rnorm) = sqrt((((*rnorm) / ((double)((c->nz + c->ndmz))))));
            if (mode == 2) {
                return;
            }
        }
        cnref_fcblok(c, g, integs, c->n, ipvtg, df, msing);
        (*msing) = (-(*msing));
        if ((*msing) != 0) {
            return;
        }
    }
    for (int l = 1; l <= c->ndmz; l += 1) {
        SET(deldmz, l, GET(rhs, l));
    }
    int iz;
    iz = 1;
    int idmz;
    idmz = 1;
    int iw;
    iw = 1;
    int izet;
    izet = 1;
    for (int i = 1; i <= c->n; i += 1) {
        int nrow;
        nrow = GET2(integs, 1, i, 3);
        c->izeta = ((nrow + 1) - c->mstar);
        if (i == c->n) {
            c->izeta = c->idum;
        }
        while ((izet != c->izeta)) {
            SET(delz, ((iz - 1) + izet), GET(rhs, (c->ndmz + izet)));
            (izet += 1);
        }
        double h;
        h = (GET(xi, (i + 1)) - GET(xi, i));
        cnref_gblock(c, h, g, nrow, c->izeta, VIEW(w, iw), v, c->kd, VIEW(delz, iz),
                     VIEW(deldmz, idmz), VIEW(ipvtw, idmz), 2);
        (iz += c->mstar);
        (idmz += c->kd);
        (iw += (c->kd * c->kd));
        if (i == c->n) {
            while ((izet <= c->mstar)) {
                SET(delz, ((iz - 1) + izet), GET(rhs, (c->ndmz + izet)));
                (izet += 1);
            }
        }
    }
    cnref_sbblok(c, g, integs, c->n, ipvtg, delz);
    cnref_dmzsol(c, c->kd, c->mstar, c->n, v, delz, deldmz);
    if (mode != 1) {
        return;
    }
    for (int l = 1; l <= c->ndmz; l += 1) {
        SET(dmz, l, GET(dmzo, l));
    }
    iz = 1;
    idmz = 1;
    iw = 1;
    izet = 1;
    for (int i = 1; i <= c->n; i += 1) {
        int nrow;
        nrow = GET2(integs, 1, i, 3);
        c->izeta = ((nrow + 1) - c->mstar);
        if (i == c->n) {
            c->izeta = c->idum;
        }
        while ((izet != c->izeta)) {
            SET(z, ((iz - 1) + izet), GET(dgz, izet));
            (izet += 1);
        }
        double h;
        h = (GET(xi, (i + 1)) - GET(xi, i));
        cnref_gblock(c, h, g, nrow, c->izeta, VIEW(w, iw), df, c->kd, VIEW(z, iz), VIEW(dmz, idmz),
                     VIEW(ipvtw, idmz), 2);
        (iz += c->mstar);
        (idmz += c->kd);
        (iw += (c->kd * c->kd));
        if (i == c->n) {
            while ((izet <= c->mstar)) {
                SET(z, ((iz - 1) + izet), GET(dgz, izet));
                (izet += 1);
            }
        }
    }
    cnref_sbblok(c, g, integs, c->n, ipvtg, z);
    cnref_dmzsol(c, c->kd, c->mstar, c->n, v, z, dmz);
}
void cnref_gderiv(cnref_context *c, R gi, int nrow, int irow, R zval, R dgz, int mode) {
    double dg_storage[40] = {0};
    R dg = dg_storage;
    for (int j = 1; j <= c->mstar; j += 1) {
        SET(dg, j, 0.0);
    }
    cnref_dbcdy(c, c->izeta, zval, dg);
    if ((c->nonlin != 0) && (c->iter <= 0)) {
        double dot;
        dot = 0.0;
        for (int j = 1; j <= c->mstar; j += 1) {
            (dot += (GET(dg, j) * GET(zval, j)));
        }
        SET(dgz, c->izeta, dot);
    }
    if (mode == 2) {
        for (int j = 1; j <= c->mstar; j += 1) {
            SET2(gi, irow, j, nrow, 0.0);
            SET2(gi, irow, (c->mstar + j), nrow, GET(dg, j));
        }
    } else {
        {
            for (int j = 1; j <= c->mstar; j += 1) {
                SET2(gi, irow, j, nrow, GET(dg, j));
                SET2(gi, irow, (c->mstar + j), nrow, 0.0);
            }
        }
    }
}
void cnref_vwblok(cnref_context *c, double xcol, double hrho, int jj, R wi, R vi, I ipvtw, int kd,
                  R zval, R df, R acol, R dmzo, int ncomp, int *msing) {
    double ha[29] = {0};
    double basm[6] = {0};
    if (jj <= 1) {
        for (int id = 1; id <= kd; id += 1) {
            SET2(wi, id, id, kd, 1.0);
        }
    }
    double fact;
    fact = 1.0;
    for (int l = 1; l <= c->mmax; l += 1) {
        fact = ((fact * hrho) / ((double)(l)));
        basm[((size_t)(l))] = fact;
        for (int j = 1; j <= c->k; j += 1) {
            ha[((size_t)(((j + (((l - 1)) * 7)))))] = (fact * GET2(acol, j, l, 7));
        }
    }
    for (int jcol = 1; jcol <= c->mstar; jcol += 1) {
        for (int ir = 1; ir <= ncomp; ir += 1) {
            SET2(df, ir, jcol, ncomp, 0.0);
        }
    }
    cnref_dfdy(c, xcol, zval, df);
    int i0;
    i0 = (((jj - 1)) * ncomp);
    int i1;
    i1 = (i0 + 1);
    int i2;
    i2 = (i0 + ncomp);
    if ((c->nonlin != 0) && (c->iter <= 0)) {
        for (int j = 1; j <= c->mstar; j += 1) {
            fact = (-GET(zval, j));
            for (int id = 1; id <= ncomp; id += 1) {
                SET(dmzo, (i0 + id), (GET(dmzo, (i0 + id)) + (fact * GET2(df, id, j, ncomp))));
            }
        }
    }
    for (int j = 1; j <= c->mstar; j += 1) {
        for (int id = 1; id <= ncomp; id += 1) {
            SET2(vi, (i0 + id), j, kd, GET2(df, id, j, ncomp));
        }
    }
    int jn;
    jn = 1;
    for (int jcomp = 1; jcomp <= ncomp; jcomp += 1) {
        int mj;
        mj = GET(c->m, jcomp);
        (jn += mj);
        for (int l = 1; l <= mj; l += 1) {
            int jv;
            jv = (jn - l);
            int jw;
            jw = jcomp;
            for (int j = 1; j <= c->k; j += 1) {
                double ajl;
                ajl = (-ha[((size_t)(((j + (((l - 1)) * 7)))))]);
                for (int iw = i1; iw <= i2; iw += 1) {
                    SET2(wi, iw, jw, kd, (GET2(wi, iw, jw, kd) + (ajl * GET2(vi, iw, jv, kd))));
                }
                (jw += ncomp);
            }
            if (l != mj) {
                for (int ll = (l + 1); ll <= mj; ll += 1) {
                    int jdf;
                    jdf = (jn - ll);
                    double bl;
                    bl = basm[((size_t)(((ll - l))))];
                    for (int iw = i1; iw <= i2; iw += 1) {
                        SET2(vi, iw, jv, kd, (GET2(vi, iw, jv, kd) + (bl * GET2(vi, iw, jdf, kd))));
                    }
                }
            }
        }
    }
    if (jj < c->k) {
        return;
    }
    (*msing) = 0;
    cnref_dgefa(c, wi, kd, kd, ipvtw, msing);
    if ((*msing) != 0) {
        return;
    }
    for (int j = 1; j <= c->mstar; j += 1) {
        cnref_dgesl(c, wi, kd, kd, ipvtw, VIEW(vi, (1 + (((j - 1)) * kd))), 0);
    }
}
void cnref_gblock(cnref_context *c, double h, R gi, int nrow, int irow, R wi, R vi, int kd, R rhsz,
                  R rhsdmz, I ipvtw, int mode) {
    double hb[29] = {0};
    double basm[6] = {0};
    double fact;
    fact = 1.0;
    basm[1] = 1.0;
    for (int l = 1; l <= c->mmax; l += 1) {
        fact = ((fact * h) / ((double)(l)));
        basm[((size_t)(((l + 1))))] = fact;
        for (int j = 1; j <= c->k; j += 1) {
            hb[((size_t)(((j + (((l - 1)) * 7)))))] = (fact * GET2(c->b, j, l, 7));
        }
    }
    if (mode == 1) {
        for (int j = 1; j <= c->mstar; j += 1) {
            for (int ir = 1; ir <= c->mstar; ir += 1) {
                SET2(gi, ((irow - 1) + ir), j, nrow, 0.0);
                SET2(gi, ((irow - 1) + ir), (c->mstar + j), nrow, 0.0);
            }
            SET2(gi, ((irow - 1) + j), (c->mstar + j), nrow, 1.0);
        }
        int ir;
        ir = irow;
        for (int icomp = 1; icomp <= c->ncomp; icomp += 1) {
            int mj;
            mj = GET(c->m, icomp);
            (ir += mj);
            for (int l = 1; l <= mj; l += 1) {
                int id;
                id = (ir - l);
                for (int jcol = 1; jcol <= c->mstar; jcol += 1) {
                    int ind;
                    ind = icomp;
                    double rsum;
                    rsum = 0.0;
                    for (int j = 1; j <= c->k; j += 1) {
                        (rsum -=
                         (hb[((size_t)(((j + (((l - 1)) * 7)))))] * GET2(vi, ind, jcol, kd)));
                        (ind += c->ncomp);
                    }
                    SET2(gi, id, jcol, nrow, rsum);
                }
                int jd;
                jd = (id - irow);
                for (int ll = 1; ll <= l; ll += 1) {
                    SET2(gi, id, (jd + ll), nrow,
                         (GET2(gi, id, (jd + ll), nrow) - basm[((size_t)(ll))]));
                }
            }
        }
    } else {
        {
            cnref_dgesl(c, wi, kd, kd, ipvtw, rhsdmz, 0);
            int ir;
            ir = irow;
            for (int jcomp = 1; jcomp <= c->ncomp; jcomp += 1) {
                int mj;
                mj = GET(c->m, jcomp);
                (ir += mj);
                for (int l = 1; l <= mj; l += 1) {
                    int ind;
                    ind = jcomp;
                    double rsum;
                    rsum = 0.0;
                    for (int j = 1; j <= c->k; j += 1) {
                        (rsum += (hb[((size_t)(((j + (((l - 1)) * 7)))))] * GET(rhsdmz, ind)));
                        (ind += c->ncomp);
                    }
                    SET(rhsz, (ir - l), rsum);
                }
            }
        }
    }
}
void cnref_dmzsol(cnref_context *c, int kd, int mstar, int n, R v, R z, R dmz) {
    (void)c;
    int jz;
    jz = 1;
    for (int i = 1; i <= n; i += 1) {
        for (int _j = 1; _j <= mstar; _j += 1) {
            double fact;
            fact = GET(z, jz);
            for (int l = 1; l <= kd; l += 1) {
                SET2(dmz, l, i, kd, (GET2(dmz, l, i, kd) + (fact * GET2(v, l, jz, kd))));
            }
            (jz += 1);
        }
    }
}
void cnref_fcblok(cnref_context *c, R bloks, I integs, int nbloks, I ipivot, R scrtch, int *info) {
    (*info) = 0;
    int indexx;
    indexx = 1;
    int indexn;
    indexn = 1;
    int i;
    i = 1;
    for (;;) {
        int index;
        index = indexn;
        int nrow;
        nrow = GET2(integs, 1, i, 3);
        int ncol;
        ncol = GET2(integs, 2, i, 3);
        int last;
        last = GET2(integs, 3, i, 3);
        cnref_factrb(c, VIEW(bloks, index), VIEW(ipivot, indexx), scrtch, nrow, ncol, last, info);
        if ((*info) != 0) {
            ((*info) += (indexx - 1));
            return;
        }
        if (i == nbloks) {
            return;
        }
        (i += 1);
        indexn = ((nrow * ncol) + index);
        (indexx += last);
        cnref_shiftb(c, VIEW(bloks, index), nrow, ncol, last, VIEW(bloks, indexn),
                     GET2(integs, 1, i, 3), GET2(integs, 2, i, 3));
    }
}
void cnref_factrb(cnref_context *c, R w, I ipivot, R d, int nrow, int ncol, int last, int *info) {
    (void)c;
    for (int i = 1; i <= nrow; i += 1) {
        SET(d, i, 0.0);
    }
    for (int j = 1; j <= ncol; j += 1) {
        for (int i = 1; i <= nrow; i += 1) {
            SET(d, i, fmax(GET(d, i), fabs(GET2(w, i, j, nrow))));
        }
    }
    int k;
    k = 1;
    for (;;) {
        if (GET(d, k) == 0.0) {
            (*info) = k;
            return;
        }
        if (k == nrow) {
            if ((fabs(GET2(w, nrow, nrow, nrow)) + GET(d, nrow)) > GET(d, nrow)) {
                return;
            }
            (*info) = k;
            return;
        }
        int l;
        l = k;
        int kp1;
        kp1 = (k + 1);
        double colmax;
        colmax = (fabs(GET2(w, k, k, nrow)) / GET(d, k));
        for (int i = kp1; i <= nrow; i += 1) {
            if (fabs(GET2(w, i, k, nrow)) <= (colmax * GET(d, i))) {
                continue;
            }
            colmax = (fabs(GET2(w, i, k, nrow)) / GET(d, i));
            l = i;
        }
        SET(ipivot, k, l);
        double t;
        t = GET2(w, l, k, nrow);
        double s;
        s = GET(d, l);
        if (l != k) {
            SET2(w, l, k, nrow, GET2(w, k, k, nrow));
            SET2(w, k, k, nrow, t);
            SET(d, l, GET(d, k));
            SET(d, k, s);
        }
        if ((fabs(t) + GET(d, k)) <= GET(d, k)) {
            (*info) = k;
            return;
        }
        t = ((-1.0) / t);
        for (int i = kp1; i <= nrow; i += 1) {
            SET2(w, i, k, nrow, (GET2(w, i, k, nrow) * t));
        }
        for (int j = kp1; j <= ncol; j += 1) {
            t = GET2(w, l, j, nrow);
            if (l != k) {
                SET2(w, l, j, nrow, GET2(w, k, j, nrow));
                SET2(w, k, j, nrow, t);
            }
            if (t == 0.0) {
                continue;
            }
            for (int i = kp1; i <= nrow; i += 1) {
                SET2(w, i, j, nrow, (GET2(w, i, j, nrow) + (GET2(w, i, k, nrow) * t)));
            }
        }
        k = kp1;
        if (k > last) {
            return;
        }
    }
}
void cnref_shiftb(cnref_context *c, R ai, int nrowi, int ncoli, int last, R ai1, int nrowi1,
                  int ncoli1) {
    (void)c;
    int mmax;
    mmax = (nrowi - last);
    int jmax;
    jmax = (ncoli - last);
    if ((mmax < 1) || (jmax < 1)) {
        return;
    }
    for (int j = 1; j <= jmax; j += 1) {
        for (int m = 1; m <= mmax; m += 1) {
            SET2(ai1, m, j, nrowi1, GET2(ai, (last + m), (last + j), nrowi));
        }
    }
    if (jmax == ncoli1) {
        return;
    }
    for (int j = (jmax + 1); j <= ncoli1; j += 1) {
        for (int m = 1; m <= mmax; m += 1) {
            SET2(ai1, m, j, nrowi1, 0.0);
        }
    }
}
void cnref_sbblok(cnref_context *c, R bloks, I integs, int nbloks, I ipivot, R x) {
    int index;
    index = 1;
    int indexx;
    indexx = 1;
    for (int i = 1; i <= nbloks; i += 1) {
        int nrow;
        nrow = GET2(integs, 1, i, 3);
        int last;
        last = GET2(integs, 3, i, 3);
        cnref_subfor(c, VIEW(bloks, index), VIEW(ipivot, indexx), nrow, last, VIEW(x, indexx));
        (index += (nrow * GET2(integs, 2, i, 3)));
        (indexx += last);
    }
    int nbp1;
    nbp1 = (nbloks + 1);
    for (int j = 1; j <= nbloks; j += 1) {
        int i;
        i = (nbp1 - j);
        int nrow;
        nrow = GET2(integs, 1, i, 3);
        int ncol;
        ncol = GET2(integs, 2, i, 3);
        int last;
        last = GET2(integs, 3, i, 3);
        (index -= (nrow * ncol));
        (indexx -= last);
        cnref_subbak(c, VIEW(bloks, index), nrow, ncol, last, VIEW(x, indexx));
    }
}
void cnref_subfor(cnref_context *c, R w, I ipivot, int nrow, int last, R x) {
    (void)c;
    if (nrow == 1) {
        return;
    }
    int lstep;
    lstep = CN_MIN(((nrow - 1)), last);
    for (int k = 1; k <= lstep; k += 1) {
        int kp1;
        kp1 = (k + 1);
        int ip;
        ip = GET(ipivot, k);
        double t;
        t = GET(x, ip);
        SET(x, ip, GET(x, k));
        SET(x, k, t);
        if (t == 0.0) {
            continue;
        }
        for (int i = kp1; i <= nrow; i += 1) {
            SET(x, i, (GET(x, i) + (GET2(w, i, k, nrow) * t)));
        }
    }
}
void cnref_subbak(cnref_context *c, R w, int nrow, int ncol, int last, R x) {
    (void)c;
    int lp1;
    lp1 = (last + 1);
    if (lp1 <= ncol) {
        for (int j = lp1; j <= ncol; j += 1) {
            double t;
            t = (-GET(x, j));
            if (t == 0.0) {
                continue;
            }
            for (int i = 1; i <= last; i += 1) {
                SET(x, i, (GET(x, i) + (GET2(w, i, j, nrow) * t)));
            }
        }
    }
    if (last != 1) {
        int lm1;
        lm1 = (last - 1);
        for (int kb = 1; kb <= lm1; kb += 1) {
            int km1;
            km1 = (last - kb);
            int k;
            k = (km1 + 1);
            SET(x, k, (GET(x, k) / GET2(w, k, k, nrow)));
            double t;
            t = (-GET(x, k));
            if (t == 0.0) {
                continue;
            }
            for (int i = 1; i <= km1; i += 1) {
                SET(x, i, (GET(x, i) + (GET2(w, i, k, nrow) * t)));
            }
        }
    }
    SET(x, 1, (GET(x, 1) / GET2(w, 1, 1, nrow)));
}
