/* Reference adaptive mesh selection and error estimates. */
#include "reference_internal.h"
void cnref_newmsh(cnref_context *c, int *mode, R xi, R xiold, R z, R dmz, R valstr, R slope,
                  R accum, int nfxpnt, R fixpnt) {
    double d1_storage[40] = {0};
    R d1 = d1_storage;
    double d2_storage[40] = {0};
    R d2 = d2_storage;
    double dummy_storage[1] = {0};
    R dummy = dummy_storage;
    int s_super;
    int s_star_super;
    if (c->k == 2) {
        s_super = 4;
        s_star_super = 4;
    } else if (c->k == 3) {
        s_super = 5;
        s_star_super = 6;
    } else if (c->k == 4) {
        s_super = 6;
        s_star_super = 9;
    } else if (1) {
        s_super = 0;
        s_star_super = 0;
    }
    int nfxp1;
    nfxp1 = (nfxpnt + 1);
    int n2;
    n2 = 0;
    int pc;
    if ((*mode) == 1) {
        pc = 180;
    } else if ((*mode) == 2) {
        pc = 100;
    } else if ((*mode) == 3) {
        pc = 50;
    } else if ((*mode) == 4) {
        pc = 20;
    } else if (1) {
        pc = 10;
    }
    for (;;) {
        if (pc == 10) {
            {
                c->mshlmt = 1;
                pc = 20;
            }
        } else if (pc == 20) {
            {
                if ((c->iguess >= 2) && (c->iguess == 3)) {
                    c->n = (c->nold / 2);
                    int i;
                    i = 0;
                    for (int j = 1; j <= c->nold; j += 2) {
                        (i += 1);
                        SET(xi, i, GET(xiold, j));
                    }
                }
                int np1;
                np1 = (c->n + 1);
                SET(xi, 1, c->aleft);
                SET(xi, np1, c->aright);
                pc = 320;
            }
        } else if (pc == 50) {
            {
                if (c->n < nfxp1) {
                    c->n = nfxp1;
                }
                int np1;
                np1 = (c->n + 1);
                SET(xi, 1, c->aleft);
                int ileft;
                ileft = 1;
                double xleft;
                xleft = c->aleft;
                for (int j = 1; j <= nfxp1; j += 1) {
                    double xright;
                    xright = c->aright;
                    int iright;
                    iright = np1;
                    if (j != nfxp1) {
                        xright = GET(fixpnt, j);
                        int nmin;
                        nmin = ((int)(((((((xright - c->aleft)) / ((c->aright - c->aleft))) *
                                         ((double)(c->n))) +
                                        1.5))));
                        if (nmin > ((c->n - nfxpnt) + j)) {
                            nmin = ((c->n - nfxpnt) + j);
                        }
                        iright = CN_MAX(((ileft + 1)), nmin);
                    }
                    SET(xi, iright, xright);
                    int nregn;
                    nregn = ((iright - ileft) - 1);
                    if (nregn != 0) {
                        double dx;
                        dx = (((xright - xleft)) / ((double)((nregn + 1))));
                        for (int i = 1; i <= nregn; i += 1) {
                            SET(xi, (ileft + i), (xleft + (((double)(i)) * dx)));
                        }
                    }
                    ileft = iright;
                    xleft = xright;
                }
                pc = 320;
            }
        } else if (pc == 100) {
            {
                n2 = (2 * c->n);
                if (n2 <= c->nmax) {
                    pc = 120;
                } else {
                    if ((*mode) == 2) {
                        c->n = n2;
                        return;
                    } else {
                        {
                            c->n = (c->nmax / 2);
                            pc = 220;
                        }
                    }
                }
            }
        } else if (pc == 120) {
            {
                if (c->mshflg == 0) {
                    pc = 140;
                    continue;
                }
                if ((c->k >= 2 && c->k <= 4) && (c->switch_ == 1)) {
                    int kstore;
                    kstore = 1;
                    for (int i = 1; i <= c->nold; i += 1) {
                        double hd6;
                        hd6 = (((GET(xiold, (i + 1)) - GET(xiold, i))) / 6.0);
                        double x;
                        x = (GET(xiold, i) + hd6);
                        cnref_evaluate(c, c->mstar, c->nold, xiold, c->y_super, x, s_super,
                                       s_star_super, VIEW(valstr, kstore), c->sk_discrete_super,
                                       c->sk_interp_super);
                        (x += (4.0 * hd6));
                        (kstore += (3 * c->mstar));
                        cnref_evaluate(c, c->mstar, c->nold, xiold, c->y_super, x, s_super,
                                       s_star_super, VIEW(valstr, kstore), c->sk_discrete_super,
                                       c->sk_interp_super);
                        (kstore += c->mstar);
                    }
                } else {
                    if ((c->k >= 2 && c->k <= 4) && (c->switch_ == 0)) {
                        c->ykstore = 1;
                        for (int i = 1; i <= c->nold; i += 1) {
                            cnref_mesh_approx(c, i, GET(xiold, i), VIEW(c->y_meshold, c->ykstore),
                                              5, xiold, c->nold, z, dmz, dummy);
                            (c->ykstore += c->mstar);
                            if (i == c->nold) {
                                cnref_mesh_approx(c, i, GET(xiold, (i + 1)),
                                                  VIEW(c->y_meshold, c->ykstore), 7, xiold, c->nold,
                                                  z, dmz, dummy);
                            }
                        }
                    } else {
                        {
                            int kstore;
                            kstore = 1;
                            for (int i = 1; i <= c->nold; i += 1) {
                                double hd6;
                                hd6 = (((GET(xiold, (i + 1)) - GET(xiold, i))) / 6.0);
                                double x;
                                x = (GET(xiold, i) + hd6);
                                cnref_mesh_approx(c, i, x, VIEW(valstr, kstore), 1, xiold, c->nold,
                                                  z, dmz, dummy);
                                (x += (4.0 * hd6));
                                (kstore += (3 * c->mstar));
                                cnref_mesh_approx(c, i, x, VIEW(valstr, kstore), 4, xiold, c->nold,
                                                  z, dmz, dummy);
                                (kstore += c->mstar);
                            }
                        }
                    }
                }
                pc = 160;
            }
        } else if (pc == 140) {
            {
                if ((c->k >= 2 && c->k <= 4) && (c->switch_ == 1)) {
                    cnref_setup_phase1(c, c->nold, xiold, VIEW(c->asave, 57), z, dmz, c->k,
                                       c->ncomp, c->mmax, c->m, c->mstar, c->y_super,
                                       c->sk_discrete_super, c->sk_interp_super);
                    int kstore;
                    kstore = 1;
                    for (int i = 1; i <= c->n; i += 1) {
                        double x;
                        x = GET(xi, i);
                        double hd6;
                        hd6 = (((GET(xi, (i + 1)) - GET(xi, i))) / 6.0);
                        for (int j = 1; j <= 4; j += 1) {
                            (x += hd6);
                            if (j == 3) {
                                (x += hd6);
                            }
                            cnref_evaluate(c, c->mstar, c->nold, xiold, c->y_super, x, s_super,
                                           s_star_super, VIEW(valstr, kstore), c->sk_discrete_super,
                                           c->sk_interp_super);
                            (kstore += c->mstar);
                        }
                    }
                } else {
                    if ((c->k >= 2 && c->k <= 4) && (c->switch_ == 0)) {
                        c->ykstore = 1;
                        for (int i = 1; i <= c->n; i += 1) {
                            cnref_mesh_approx(c, i, GET(xi, i), VIEW(c->y_meshold, c->ykstore), 5,
                                              xiold, c->nold, z, dmz, dummy);
                            (c->ykstore += c->mstar);
                            if (i == c->n) {
                                cnref_mesh_approx(c, i, GET(xi, (i + 1)),
                                                  VIEW(c->y_meshold, c->ykstore), 7, xiold, c->nold,
                                                  z, dmz, dummy);
                            }
                        }
                    } else {
                        {
                            int kstore;
                            kstore = 1;
                            for (int i = 1; i <= c->n; i += 1) {
                                double x;
                                x = GET(xi, i);
                                double hd6;
                                hd6 = (((GET(xi, (i + 1)) - GET(xi, i))) / 6.0);
                                for (int j = 1; j <= 4; j += 1) {
                                    (x += hd6);
                                    if (j == 3) {
                                        (x += hd6);
                                    }
                                    cnref_mesh_approx(c, i, x, VIEW(valstr, kstore), j, xiold,
                                                      c->nold, z, dmz, dummy);
                                    (kstore += c->mstar);
                                }
                            }
                        }
                    }
                }
                pc = 160;
            }
        } else if (pc == 160) {
            {
                c->mshflg = 0;
                c->mshnum = 1;
                (*mode) = 2;
                int j;
                j = 2;
                for (int i = 1; i <= c->n; i += 1) {
                    SET(xi, j, (((GET(xiold, i) + GET(xiold, (i + 1)))) / 2.0));
                    SET(xi, (j + 1), GET(xiold, (i + 1)));
                    (j += 2);
                }
                c->n = n2;
                pc = 320;
            }
        } else if (pc == 180) {
            {
                if ((c->nold == 1) || (c->nold <= (2 * nfxpnt))) {
                    pc = 100;
                    continue;
                }
                double hiold;
                hiold = (GET(xiold, 2) - GET(xiold, 1));
                cnref_horder(c, 1, d1, hiold, dmz, c->ncomp, c->k);
                if (!cnref_all_finite(c, d1, c->ncomp)) {
                    return;
                }
                hiold = (GET(xiold, 3) - GET(xiold, 2));
                cnref_horder(c, 2, d2, hiold, dmz, c->ncomp, c->k);
                if (!cnref_all_finite(c, d2, c->ncomp)) {
                    return;
                }
                SET(accum, 1, 0.0);
                SET(slope, 1, 0.0);
                double oneovh;
                oneovh = (2.0 / ((GET(xiold, 3) - GET(xiold, 1))));
                if (!cnref_mesh_finite(c, oneovh)) {
                    return;
                }
                for (int j = 1; j <= c->ntol; j += 1) {
                    int jj;
                    jj = GET(c->jtol, j);
                    int jz;
                    jz = GET(c->ltol, j);
                    double density;
                    density =
                        cnref_mesh_density(c, GET(d1, jj), GET(d2, jj), j, oneovh, GET(z, jz));
                    if (!isfinite(density)) {
                        return;
                    }
                    SET(slope, 1, fmax(GET(slope, 1), density));
                }
                double slphmx;
                slphmx = (GET(slope, 1) * ((GET(xiold, 2) - GET(xiold, 1))));
                if (!cnref_mesh_finite(c, slphmx)) {
                    return;
                }
                SET(accum, 2, slphmx);
                int iflip;
                iflip = 1;
                for (int i = 2; i <= c->nold; i += 1) {
                    hiold = (GET(xiold, (i + 1)) - GET(xiold, i));
                    if (iflip == (-1)) {
                        cnref_horder(c, i, d1, hiold, dmz, c->ncomp, c->k);
                    }
                    if (iflip == 1) {
                        cnref_horder(c, i, d2, hiold, dmz, c->ncomp, c->k);
                    }
                    if (!cnref_all_pair_finite(c, d1, d2, c->ncomp)) {
                        return;
                    }
                    oneovh = (2.0 / ((GET(xiold, (i + 1)) - GET(xiold, (i - 1)))));
                    if (!cnref_mesh_finite(c, oneovh)) {
                        return;
                    }
                    SET(slope, i, 0.0);
                    for (int j = 1; j <= c->ntol; j += 1) {
                        int jj;
                        jj = GET(c->jtol, j);
                        int jz;
                        jz = (GET(c->ltol, j) + (((i - 1)) * c->mstar));
                        double density;
                        density =
                            cnref_mesh_density(c, GET(d1, jj), GET(d2, jj), j, oneovh, GET(z, jz));
                        if (!isfinite(density)) {
                            return;
                        }
                        SET(slope, i, fmax(GET(slope, i), density));
                    }
                    double temp;
                    temp = (GET(slope, i) * ((GET(xiold, (i + 1)) - GET(xiold, i))));
                    if (!cnref_mesh_finite(c, temp)) {
                        return;
                    }
                    slphmx = fmax(slphmx, temp);
                    SET(accum, (i + 1), (GET(accum, i) + temp));
                    if (!cnref_mesh_finite(c, GET(accum, (i + 1)))) {
                        return;
                    }
                    iflip = (-iflip);
                }
                double avrg;
                avrg = (GET(accum, (c->nold + 1)) / ((double)(c->nold)));
                double degequ;
                degequ = (avrg / fmax(slphmx, c->precis));
                double prediction;
                prediction = (GET(accum, (c->nold + 1)) + 1.0);
                if (!(cnref_mesh_finite(c, avrg) && cnref_mesh_finite(c, degequ) &&
                      cnref_mesh_finite(c, prediction))) {
                    return;
                }
                int naccum;
                naccum = ((int)(prediction));
                if ((avrg < c->precis) || (degequ >= 0.5)) {
                    pc = 100;
                    continue;
                }
                int nmx;
                nmx = (CN_MAX(((c->nold + 1)), naccum) / 2);
                int nmax2;
                nmax2 = (c->nmax / 2);
                c->n = CN_MIN(CN_MIN(nmax2, c->nold), nmx);
                pc = 220;
            }
        } else if (pc == 220) {
            {
                int noldp1;
                noldp1 = (c->nold + 1);
                if (c->n < nfxp1) {
                    c->n = nfxp1;
                }
                (c->mshnum += 1);
                if (c->n < c->nold) {
                    c->mshnum = c->mshlmt;
                }
                if (c->n > (c->nold / 2)) {
                    c->mshalt = 1;
                }
                if (c->n == (c->nold / 2)) {
                    (c->mshalt += 1);
                }
                c->mshflg = 0;
                int inside;
                inside = 1;
                double accl;
                accl = 0.0;
                int lold;
                lold = 2;
                SET(xi, 1, c->aleft);
                SET(xi, (c->n + 1), c->aright);
                for (int i = 1; i <= nfxp1; i += 1) {
                    int lnew;
                    lnew = noldp1;
                    double accr;
                    int nregn;
                    if (i != nfxp1) {
                        for (int j = lold; j <= noldp1; j += 1) {
                            lnew = j;
                            if (GET(fixpnt, i) <= GET(xiold, j)) {
                                break;
                            }
                        }
                        double shift;
                        shift = (((GET(fixpnt, i) - GET(xiold, lnew))) * GET(slope, (lnew - 1)));
                        accr = (GET(accum, lnew) + shift);
                        double proposed;
                        proposed =
                            (((((accr - accl)) / GET(accum, noldp1)) * ((double)(c->n))) - 0.5);
                        if (!(cnref_mesh_finite(c, shift) && cnref_mesh_finite(c, accr) &&
                              cnref_mesh_finite(c, proposed))) {
                            return;
                        }
                        int proposed_i;
                        proposed_i = ((int)(proposed));
                        nregn = CN_MIN(proposed_i, (((c->n - inside) - nfxp1) + i));
                        SET(xi, ((inside + nregn) + 1), GET(fixpnt, i));
                    } else {
                        {
                            accr = GET(accum, noldp1);
                            nregn = (c->n - inside);
                        }
                    }
                    if (nregn != 0) {
                        double temp;
                        temp = accl;
                        double tsum;
                        tsum = (((accr - accl)) / ((double)((nregn + 1))));
                        if (!cnref_mesh_finite(c, tsum)) {
                            return;
                        }
                        for (int _j = 1; _j <= nregn; _j += 1) {
                            (inside += 1);
                            (temp += tsum);
                            if (!cnref_mesh_finite(c, temp)) {
                                return;
                            }
                            int lcarry;
                            lcarry = lold;
                            for (int l = lold; l <= lnew; l += 1) {
                                lcarry = l;
                                if (temp <= GET(accum, l)) {
                                    break;
                                }
                            }
                            lold = lcarry;
                            SET(xi, inside,
                                (GET(xiold, (lold - 1)) +
                                 (((temp - GET(accum, (lold - 1)))) / GET(slope, (lold - 1)))));
                        }
                    }
                    (inside += 1);
                    accl = accr;
                    lold = lnew;
                }
                (*mode) = 1;
                pc = 320;
            }
        } else if (pc == 320) {
            {
                if (!cnref_check_mesh(c, xi, c->n)) {
                    (c->invalid_evaluation = 1);
                    return;
                }
                c->nz = (c->mstar * ((c->n + 1)));
                c->ndmz = (c->kd * c->n);
                return;
            }
        } else if (1) {
            abort();
        }
    }
}
int cnref_mesh_finite(cnref_context *c, double value) {
    if (isfinite(value)) {
        return 1;
    } else {
        {
            (c->invalid_evaluation = 1);
            return 0;
        }
    }
}
int cnref_estimate_finite(cnref_context *c, double value, int preserve_scsa) {
    return ((isfinite(value) || preserve_scsa) || cnref_mesh_finite(c, value));
}
void cnref_consts(cnref_context *c, int k, R rho, R coef) {
    double cnsts1[28] = {0.25,      0.625e-1,  7.2169e-2, 1.8342e-2, 1.9065e-2, 5.8190e-2,
                         5.4658e-3, 5.3370e-3, 1.8890e-2, 2.7792e-2, 1.6095e-3, 1.4964e-3,
                         7.5938e-3, 5.7573e-3, 1.8342e-2, 4.673e-3,  4.150e-4,  1.919e-3,
                         1.468e-3,  6.371e-3,  4.610e-3,  1.342e-4,  1.138e-4,  4.889e-4,
                         4.177e-4,  1.374e-3,  1.654e-3,  2.863e-3};
    double cnsts2[28] = {1.25e-1,   2.604e-3,  8.019e-3,  2.170e-5,  7.453e-5,  5.208e-4,
                         9.689e-8,  3.689e-7,  3.100e-6,  2.451e-5,  2.691e-10, 1.120e-9,
                         1.076e-8,  9.405e-8,  1.033e-6,  5.097e-13, 2.290e-12, 2.446e-11,
                         2.331e-10, 2.936e-9,  3.593e-8,  7.001e-16, 3.363e-15, 3.921e-14,
                         4.028e-13, 5.646e-12, 7.531e-11, 1.129e-9};
    if (k == 2) {
        cnsts1[2] = cnsts1[1];
    }
    if (k == 3) {
        cnsts1[5] = cnsts1[3];
    }
    if (k == 4) {
        cnsts1[9] = cnsts1[6];
    }
    int koff;
    koff = ((k * ((k + 1))) / 2);
    int iz;
    iz = 1;
    for (int j = 1; j <= c->ncomp; j += 1) {
        int mj;
        mj = GET(c->m, j);
        for (int l = 1; l <= mj; l += 1) {
            SET(c->wgterr, iz, cnsts1[((size_t)(((((koff - mj) + l) - 1))))]);
            (iz += 1);
        }
    }
    int jcomp;
    jcomp = 1;
    int mtot;
    mtot = GET(c->m, 1);
    for (int i = 1; i <= c->ntol; i += 1) {
        int ltoli;
        ltoli = GET(c->ltol, i);
        while ((ltoli > mtot)) {
            (jcomp += 1);
            (mtot += GET(c->m, jcomp));
        }
        SET(c->jtol, i, jcomp);
        SET(c->wgtmsh, i,
            ((1e1 * cnsts2[((size_t)(((((koff + ltoli) - mtot) - 1))))]) / GET(c->tolin, i)));
        SET(c->root, i,
            ((k >= 2 && k <= 4) ? (1.0 / ((double)((2 * k))))
                                : (1.0 / ((double)((((k + mtot) - ltoli) + 1))))));
    }
    if (k == 1) {
        {
            SET(rho, 1, 0.0);
        }
    } else if (k == 2) {
        {
            SET(rho, 2, 0.57735026918962576451);
            SET(rho, 1, (-GET(rho, 2)));
        }
    } else if (k == 3) {
        {
            SET(rho, 3, 0.77459666924148337704);
            SET(rho, 2, 0.0);
            SET(rho, 1, (-GET(rho, 3)));
        }
    } else if (k == 4) {
        {
            SET(rho, 4, 0.86113631159405257523);
            SET(rho, 3, 0.33998104358485626480);
            SET(rho, 2, (-GET(rho, 3)));
            SET(rho, 1, (-GET(rho, 4)));
        }
    } else if (k == 5) {
        {
            SET(rho, 5, 0.90617984593866399280);
            SET(rho, 4, 0.53846931010568309104);
            SET(rho, 3, 0.0);
            SET(rho, 2, (-GET(rho, 4)));
            SET(rho, 1, (-GET(rho, 5)));
        }
    } else if (k == 6) {
        {
            SET(rho, 6, 0.93246951420315202781);
            SET(rho, 5, 0.66120938646626451366);
            SET(rho, 4, 0.23861918608319690863);
            SET(rho, 3, (-GET(rho, 4)));
            SET(rho, 2, (-GET(rho, 5)));
            SET(rho, 1, (-GET(rho, 6)));
        }
    } else if (k == 7) {
        {
            SET(rho, 7, 0.949107991234275852452);
            SET(rho, 6, 0.74153118559939443986);
            SET(rho, 5, 0.40584515137739716690);
            SET(rho, 4, 0.0);
            SET(rho, 3, (-GET(rho, 5)));
            SET(rho, 2, (-GET(rho, 6)));
            SET(rho, 1, (-GET(rho, 7)));
        }
    } else if (1) {
        abort();
    }
    for (int j = 1; j <= k; j += 1) {
        SET(rho, j, (0.5 * ((1.0 + GET(rho, j)))));
    }
    for (int j = 1; j <= k; j += 1) {
        for (int i = 1; i <= k; i += 1) {
            SET2(coef, i, j, k, 0.0);
        }
        SET2(coef, j, j, k, 1.0);
        cnref_vmonde(c, rho, VIEW(coef, (1 + (((j - 1)) * k))), k);
    }
    double dummy_storage[1] = {0};
    R dummy = dummy_storage;
    cnref_rkbas(c, 1.0, coef, k, c->mmax, c->b, dummy, 0);
    for (int i = 1; i <= k; i += 1) {
        cnref_rkbas(c, GET(rho, i), coef, k, c->mmax, VIEW(c->acol, (1 + (((i - 1)) * 28))), dummy,
                    0);
    }
    cnref_rkbas(c, (1.0 / 6.0), coef, k, c->mmax, VIEW(c->asave, 1), dummy, 0);
    cnref_rkbas(c, (1.0 / 3.0), coef, k, c->mmax, VIEW(c->asave, 29), dummy, 0);
    cnref_rkbas(c, (2.0 / 3.0), coef, k, c->mmax, VIEW(c->asave, 57), dummy, 0);
    cnref_rkbas(c, (5.0 / 6.0), coef, k, c->mmax, VIEW(c->asave, 85), dummy, 0);
    cnref_rkbas(c, 0.0, coef, k, c->mmax, VIEW(c->asave, 113), dummy, 0);
    cnref_rkbas(c, (1.0 / 2.0), coef, k, c->mmax, VIEW(c->asave, 141), dummy, 0);
    cnref_rkbas(c, 1.0, coef, k, c->mmax, VIEW(c->asave, 169), dummy, 0);
}
void cnref_errchk(cnref_context *c, R xi, R z, R dmz, R valstr, int *ifin) {
    double err_storage[40] = {0};
    R err = err_storage;
    int inherited_err[40] = {0};
    double errest_storage[40] = {0};
    R errest = errest_storage;
    double err1_storage[40] = {0};
    R err1 = err1_storage;
    double err3_storage[40] = {0};
    R err3 = err3_storage;
    double dummy_storage[1] = {0};
    R dummy = dummy_storage;
    int scsa;
    scsa = (c->k >= 2 && c->k <= 4);
    int s_super;
    int s_star_super;
    if (c->k == 2) {
        s_super = 4;
        s_star_super = 4;
    } else if (c->k == 3) {
        s_super = 5;
        s_star_super = 6;
    } else if (c->k == 4) {
        s_super = 6;
        s_star_super = 9;
    } else if (1) {
        s_super = 0;
        s_star_super = 0;
    }
    if (scsa && (c->switch_ == 0)) {
        (*ifin) = 1;
        c->mshflg = 1;
        for (int j = 1; j <= c->mstar; j += 1) {
            SET(errest, j, 0.0);
        }
        c->yknew = 1;
        for (int i = 1; i <= c->n; i += 2) {
            cnref_mesh_approx(c, i, GET(xi, i), VIEW(c->y_meshnew, c->yknew), 5, xi, c->n, z, dmz,
                              dummy);
            for (int l = 1; l <= c->mstar; l += 1) {
                double difference;
                difference = fabs(((GET(c->y_meshnew, c->yknew) - GET(c->y_meshold, c->yknew))));
                double estimate;
                estimate = (difference / cnref_powi(2.0, (2 * c->k)));
                if ((!cnref_mesh_finite(c, difference)) || (!cnref_mesh_finite(c, estimate))) {
                    (*ifin) = 0;
                    return;
                }
                SET(err1, l, estimate);
                (c->yknew += 1);
            }
            if (i == (c->n - 1)) {
                cnref_mesh_approx(c, (i + 1), GET(xi, (i + 2)), VIEW(c->y_meshnew, c->yknew), 7, xi,
                                  c->n, z, dmz, dummy);
                for (int l = 1; l <= c->mstar; l += 1) {
                    double difference;
                    difference =
                        fabs(((GET(c->y_meshnew, c->yknew) - GET(c->y_meshold, c->yknew))));
                    double estimate;
                    estimate = (difference / cnref_powi(2.0, (2 * c->k)));
                    if ((!cnref_mesh_finite(c, difference)) || (!cnref_mesh_finite(c, estimate))) {
                        (*ifin) = 0;
                        return;
                    }
                    SET(err3, l, estimate);
                    (c->yknew += 1);
                }
            }
            if (i != (c->n - 1)) {
                for (int l = 1; l <= c->mstar; l += 1) {
                    SET(errest, l, fmax(GET(errest, l), GET(err1, l)));
                }
            } else {
                {
                    for (int l = 1; l <= c->mstar; l += 1) {
                        SET(errest, l, fmax(fmax(GET(errest, l), GET(err1, l)), GET(err3, l)));
                    }
                }
            }
            if ((*ifin) != 0) {
                for (int j = 1; j <= c->ntol; j += 1) {
                    int ltolj;
                    ltolj = GET(c->ltol, j);
                    int ltjz;
                    ltjz = (ltolj + (((i - 1)) * c->mstar));
                    double threshold;
                    threshold = (GET(c->tolin, j) * ((fabs(GET(z, ltjz)) + 1.0)));
                    if (!cnref_mesh_finite(c, threshold)) {
                        (*ifin) = 0;
                        return;
                    }
                    if (GET(errest, ltolj) > threshold) {
                        (*ifin) = 0;
                    }
                }
            }
        }
        if ((*ifin) != 1) {
            return;
        }
        c->switch_ = 1;
        (*ifin) = 0;
        c->mshflg = 0;
        cnref_setup_phase1(c, c->nold, c->meshold, VIEW(c->asave, 57), c->zold, c->dmzold, c->k,
                           c->ncomp, c->mmax, c->m, c->mstar, c->y_super, c->sk_discrete_super,
                           c->sk_interp_super);
        int kstore;
        kstore = 1;
        for (int i = 1; i <= c->nold; i += 1) {
            double x;
            x = GET(c->meshold, i);
            double hd6;
            hd6 = (((GET(c->meshold, (i + 1)) - GET(c->meshold, i))) / 6.0);
            for (int j = 1; j <= 4; j += 1) {
                (x += hd6);
                if (j == 3) {
                    (x += hd6);
                }
                cnref_evaluate(c, c->mstar, c->nold, c->meshold, c->y_super, x, s_super,
                               s_star_super, VIEW(valstr, kstore), c->sk_discrete_super,
                               c->sk_interp_super);
                (kstore += c->mstar);
            }
        }
    }
    (*ifin) = 1;
    c->mshflg = 1;
    for (int j = 1; j <= c->mstar; j += 1) {
        SET(errest, j, 0.0);
    }
    if (scsa) {
        cnref_setup_phase1(c, c->n, xi, VIEW(c->asave, 57), z, dmz, c->k, c->ncomp, c->mmax, c->m,
                           c->mstar, c->y_super, c->sk_discrete_super, c->sk_interp_super);
    }
    int preserve_scsa;
    preserve_scsa = ((scsa && c->problem->preserve_non_finite_interpolation) &&
                     (c->diagnostics.non_finite_interpolation_evaluations > 0));
    for (int iback = 1; iback <= c->n; iback += 1) {
        int i;
        i = ((c->n + 1) - iback);
        int knew;
        knew = (((((4 * ((i - 1))) + 2)) * c->mstar) + 1);
        int kstore;
        kstore = (((((2 * ((i - 1))) + 1)) * c->mstar) + 1);
        double x;
        x = (GET(xi, i) + ((((GET(xi, (i + 1)) - GET(xi, i))) * 2.0) / 3.0));
        if (scsa) {
            cnref_evaluate(c, c->mstar, c->n, xi, c->y_super, x, s_super, s_star_super,
                           VIEW(valstr, knew), c->sk_discrete_super, c->sk_interp_super);
        } else {
            {
                cnref_mesh_approx(c, i, x, VIEW(valstr, knew), 3, xi, c->n, z, dmz, dummy);
            }
        }
        for (int l = 1; l <= c->mstar; l += 1) {
            double newer;
            newer = GET(valstr, knew);
            double older;
            older = GET(valstr, kstore);
            double difference;
            difference = fabs(((newer - older)));
            double estimate;
            estimate = (GET(c->wgterr, l) * difference);
            int inherited_difference;
            inherited_difference = (preserve_scsa && (((!isfinite(newer)) || (!isfinite(older)))));
            int inherited_estimate;
            inherited_estimate = (inherited_difference && (!isfinite(difference)));
            if ((!cnref_estimate_finite(c, difference, inherited_difference)) ||
                (!cnref_estimate_finite(c, estimate, inherited_estimate))) {
                (*ifin) = 0;
                return;
            }
            SET(err, l, estimate);
            inherited_err[((size_t)(((l - 1))))] = (inherited_estimate && (!isfinite(estimate)));
            (knew += 1);
            (kstore += 1);
        }
        knew = (((((4 * ((i - 1))) + 1)) * c->mstar) + 1);
        kstore = (((2 * ((i - 1))) * c->mstar) + 1);
        x = (GET(xi, i) + (((GET(xi, (i + 1)) - GET(xi, i))) / 3.0));
        if (scsa) {
            cnref_evaluate(c, c->mstar, c->n, xi, c->y_super, x, s_super, s_star_super,
                           VIEW(valstr, knew), c->sk_discrete_super, c->sk_interp_super);
        } else {
            {
                cnref_mesh_approx(c, i, x, VIEW(valstr, knew), 2, xi, c->n, z, dmz, dummy);
            }
        }
        for (int l = 1; l <= c->mstar; l += 1) {
            double newer;
            newer = GET(valstr, knew);
            double older;
            older = GET(valstr, kstore);
            double difference;
            difference = fabs(((newer - older)));
            double weighted;
            weighted = (GET(c->wgterr, l) * difference);
            double estimate;
            estimate = (GET(err, l) + weighted);
            int inherited_difference;
            inherited_difference = (preserve_scsa && (((!isfinite(newer)) || (!isfinite(older)))));
            int inherited_weighted;
            inherited_weighted = (inherited_difference && (!isfinite(difference)));
            int inherited_sum;
            inherited_sum = (inherited_err[((size_t)(((l - 1))))] ||
                             ((inherited_weighted && (!isfinite(weighted)))));
            if (((!cnref_estimate_finite(c, difference, inherited_difference)) ||
                 (!cnref_estimate_finite(c, weighted, inherited_weighted))) ||
                (!cnref_estimate_finite(c, estimate, inherited_sum))) {
                (*ifin) = 0;
                return;
            }
            SET(err, l, estimate);
            (knew += 1);
            (kstore += 1);
        }
        for (int l = 1; l <= c->mstar; l += 1) {
            SET(errest, l, fmax(GET(errest, l), GET(err, l)));
        }
        if ((*ifin) != 0) {
            for (int j = 1; j <= c->ntol; j += 1) {
                int ltolj;
                ltolj = GET(c->ltol, j);
                int ltjz;
                ltjz = (ltolj + (((i - 1)) * c->mstar));
                double threshold;
                threshold = (GET(c->tolin, j) * ((fabs(GET(z, ltjz)) + 1.0)));
                if (!cnref_mesh_finite(c, threshold)) {
                    (*ifin) = 0;
                    return;
                }
                if (GET(err, ltolj) > threshold) {
                    (*ifin) = 0;
                }
            }
        }
    }
}
void cnref_mesh_approx(cnref_context *c, int i, double x, R zval, int column, R xi, int n, R z,
                       R dmz, R dummy) {
    cnref_approx(c, &(i), &(x), zval, VIEW(c->asave, (1 + (((column - 1)) * 28))), dummy, xi, n, z,
                 dmz, c->k, c->ncomp, c->mmax, c->m, c->mstar, 4, dummy, 0);
}
