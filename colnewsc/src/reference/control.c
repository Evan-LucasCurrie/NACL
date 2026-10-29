/* Reference collocation iteration and error controller. */
#include "reference_internal.h"
void cnref_contrl(cnref_context *c, R xi, R xiold, R z, R dmz, R rhs, R delz, R deldmz, R dqz,
                  R dqdmz, R g, R w, R v, R valstr, R slope, R scale, R dscale, R accum, I ipvtg,
                  I integs, I ipvtw, int nfxpnt, R fixpnt, int *iflag) {
    double dummy_storage[1] = {0};
    R dummy = dummy_storage;
    double relmin;
    relmin = 1e-3;
    double rstart;
    rstart = 1e-2;
    int lmtfrz;
    lmtfrz = 4;
    double check;
    check = 0.0;
    for (int i = 1; i <= c->ntol; i += 1) {
        check = fmax(GET(c->tolin, i), check);
    }
    int imesh;
    imesh = 1;
    int iconv;
    if (c->nonlin == 0) {
        iconv = 1;
    } else {
        iconv = 0;
    }
    int noconv;
    noconv = 0;
    int msing;
    msing = 0;
    int ifreez;
    ifreez = 0;
    int ifrz;
    ifrz = 0;
    int ipred;
    ipred = 0;
    int ifin;
    ifin = 0;
    double relax;
    relax = 1.0;
    double rlxold;
    rlxold = 1.0;
    double rnold;
    rnold = 0.0;
    double rnorm;
    rnorm = 0.0;
    double anfix;
    anfix = 0.0;
    double anorm;
    anorm = 0.0;
    int pc;
    pc = 20;
    size_t corrections;
    corrections = 0;
    size_t refinements;
    refinements = 0;
    for (;;) {
        if (c->invalid_evaluation || !(isfinite(rnold) && isfinite(rnorm) && isfinite(relax) &&
                                       isfinite(rlxold) && isfinite(anfix) && isfinite(anorm))) {
            (*iflag) = (-4);
            return;
        }
        if (pc == 20) {
            {
                corrections = 0;
                cnref_history(c, c->n);
                if (c->iprint < 1) {
                    fprintf(stderr, "reference mesh: %d intervals\n", c->n);
                }
                c->iter = 0;
                if (c->nonlin > 0) {
                    pc = 50;
                } else {
                    {
                        cnref_lsyslv(c, &(msing), xi, xiold, dummy, dummy, z, dmz, g, w, v, rhs,
                                     dummy, integs, ipvtg, ipvtw, &(rnorm), 0);
                        if (msing == 0) {
                            pc = 400;
                        } else {
                            pc = 30;
                        }
                    }
                }
            }
        } else if (pc == 30) {
            {
                if (msing < 0) {
                    (*iflag) = 0;
                    return;
                }
                pc = 460;
            }
        } else if (pc == 50) {
            {
                if ((c->icare == 1) || (c->icare == (-1))) {
                    relax = rstart;
                } else {
                    relax = 1.0;
                }
                if (iconv == 0) {
                    pc = 160;
                } else {
                    {
                        ifreez = 0;
                        cnref_lsyslv(c, &(msing), xi, xiold, z, dmz, delz, deldmz, g, w, v, rhs,
                                     dqdmz, integs, ipvtg, ipvtw, &(rnold), 1);
                        pc = 70;
                    }
                }
            }
        } else if (pc == 60) {
            {
                rnold = rnorm;
                cnref_lsyslv(c, &(msing), xi, xiold, z, dmz, delz, deldmz, g, w, v, rhs, dummy,
                             integs, ipvtg, ipvtw, &(rnorm), (3 + ifreez));
                pc = 70;
            }
        } else if (pc == 70) {
            {
                if (msing != 0) {
                    pc = 30;
                    continue;
                }
                if (corrections >= c->maximum_corrections) {
                    (*iflag) = (-2);
                    return;
                }
                (corrections += 1);
                if (ifreez != 1) {
                    (c->iter += 1);
                    ifrz = 0;
                }
                for (int i = 1; i <= c->nz; i += 1) {
                    SET(z, i, (GET(z, i) + GET(delz, i)));
                }
                for (int i = 1; i <= c->ndmz; i += 1) {
                    SET(dmz, i, (GET(dmz, i) + GET(deldmz, i)));
                }
                cnref_lsyslv(c, &(msing), xi, xiold, z, dmz, delz, deldmz, g, w, v, rhs, dummy,
                             integs, ipvtg, ipvtw, &(rnorm), 2);
                if (rnorm < c->precis) {
                    pc = 390;
                } else {
                    if (rnorm > rnold) {
                        pc = 130;
                    } else {
                        if (ifreez == 1) {
                            pc = 110;
                        } else {
                            ifreez = 1;
                            pc = 60;
                        }
                    }
                }
            }
        } else if (pc == 110) {
            {
                (ifrz += 1);
                if ((ifrz >= lmtfrz) || (rnold < (4.0 * rnorm))) {
                    ifreez = 0;
                }
                int large;
                large = 0;
                for (int it = 1; it <= c->ntol; it += 1) {
                    int iz;
                    iz = GET(c->ltol, it);
                    while ((iz <= c->nz)) {
                        if (fabs(GET(delz, iz)) >
                            (GET(c->tolin, it) * ((fabs(GET(z, iz)) + 1.0)))) {
                            large = 1;
                            break;
                        }
                        (iz += c->mstar);
                    }
                    if (large) {
                        break;
                    }
                }
                if (large) {
                    pc = 60;
                } else {
                    pc = 400;
                }
            }
        } else if (pc == 130) {
            {
                iconv = 0;
                relax = rstart;
                for (int i = 1; i <= c->nz; i += 1) {
                    SET(z, i, (GET(z, i) - GET(delz, i)));
                }
                for (int i = 1; i <= c->ndmz; i += 1) {
                    SET(dmz, i, (GET(dmz, i) - GET(deldmz, i)));
                }
                for (int i = 1; i <= (c->n + 1); i += 1) {
                    SET(xiold, i, GET(xi, i));
                }
                c->nold = c->n;
                c->iter = 0;
                pc = 160;
            }
        } else if (pc == 160) {
            {
                cnref_lsyslv(c, &(msing), xi, xiold, z, dmz, delz, deldmz, g, w, v, rhs, dqdmz,
                             integs, ipvtg, ipvtw, &(rnold), 1);
                if (msing != 0) {
                    pc = 30;
                    continue;
                }
                if (c->iguess == 1) {
                    c->iguess = 0;
                }
                cnref_skale(c, c->n, c->mstar, c->kd, z, xi, scale, dscale);
                pc = 220;
            }
        } else if (pc == 170) {
            {
                rnold = rnorm;
                if (c->iter >= c->limit) {
                    pc = 450;
                    continue;
                }
                cnref_skale(c, c->n, c->mstar, c->kd, z, xi, scale, dscale);
                double anscl;
                anscl = 0.0;
                for (int i = 1; i <= c->nz; i += 1) {
                    (anscl += cnref_powi(((GET(delz, i) * GET(scale, i))), 2));
                }
                for (int i = 1; i <= c->ndmz; i += 1) {
                    (anscl += cnref_powi(((GET(deldmz, i) * GET(dscale, i))), 2));
                }
                anscl = sqrt(((anscl / ((double)((c->nz + c->ndmz))))));
                if (!isfinite(anscl)) {
                    (*iflag) = (-4);
                    return;
                }
                cnref_lsyslv(c, &(msing), xi, xiold, z, dmz, delz, deldmz, g, w, v, rhs, dummy,
                             integs, ipvtg, ipvtw, &(rnorm), 3);
                if (msing != 0) {
                    pc = 30;
                    continue;
                }
                double andif;
                andif = 0.0;
                for (int i = 1; i <= c->nz; i += 1) {
                    (andif += cnref_powi(((((GET(dqz, i) - GET(delz, i))) * GET(scale, i))), 2));
                }
                for (int i = 1; i <= c->ndmz; i += 1) {
                    (andif +=
                     cnref_powi(((((GET(dqdmz, i) - GET(deldmz, i))) * GET(dscale, i))), 2));
                }
                andif = sqrt((((andif / ((double)((c->nz + c->ndmz)))) + c->precis)));
                if (!isfinite(andif)) {
                    (*iflag) = (-4);
                    return;
                }
                relax = ((relax * anscl) / andif);
                if (relax > 1.0) {
                    relax = 1.0;
                }
                pc = 220;
            }
        } else if (pc == 220) {
            {
                if (corrections >= c->maximum_corrections) {
                    (*iflag) = (-2);
                    return;
                }
                (corrections += 1);
                rlxold = relax;
                ipred = 1;
                (c->iter += 1);
                for (int i = 1; i <= c->nz; i += 1) {
                    SET(z, i, (GET(z, i) + (relax * GET(delz, i))));
                }
                for (int i = 1; i <= c->ndmz; i += 1) {
                    SET(dmz, i, (GET(dmz, i) + (relax * GET(deldmz, i))));
                }
                pc = 250;
            }
        } else if (pc == 250) {
            {
                cnref_lsyslv(c, &(msing), xi, xiold, z, dmz, dqz, dqdmz, g, w, v, rhs, dummy,
                             integs, ipvtg, ipvtw, &(rnorm), 2);
                cnref_lsyslv(c, &(msing), xi, xiold, z, dmz, dqz, dqdmz, g, w, v, rhs, dummy,
                             integs, ipvtg, ipvtw, &(rnorm), 4);
                anorm = 0.0;
                anfix = 0.0;
                for (int i = 1; i <= c->nz; i += 1) {
                    (anorm += cnref_powi(((GET(delz, i) * GET(scale, i))), 2));
                    (anfix += cnref_powi(((GET(dqz, i) * GET(scale, i))), 2));
                }
                for (int i = 1; i <= c->ndmz; i += 1) {
                    (anorm += cnref_powi(((GET(deldmz, i) * GET(dscale, i))), 2));
                    (anfix += cnref_powi(((GET(dqdmz, i) * GET(dscale, i))), 2));
                }
                anorm = sqrt(((anorm / ((double)((c->nz + c->ndmz))))));
                anfix = sqrt(((anfix / ((double)((c->nz + c->ndmz))))));
                if ((anfix < c->precis) || (rnorm < c->precis)) {
                    pc = 390;
                } else {
                    if (anfix > anorm) {
                        pc = 300;
                    } else {
                        if (anfix <= check) {
                            pc = 350;
                        } else {
                            if (ipred != 1) {
                                pc = 170;
                            } else {
                                pc = 300;
                            }
                        }
                    }
                }
            }
        } else if (pc == 300) {
            {
                if (c->iter >= c->limit) {
                    pc = 450;
                    continue;
                }
                ipred = 0;
                double arg;
                arg = (((((anfix / anorm) - 1.0)) / relax) + 1.0);
                if (!isfinite(arg)) {
                    (*iflag) = (-4);
                    return;
                }
                if (arg < 0.0) {
                    pc = 170;
                    continue;
                }
                if (arg <= ((0.25 * relax) + (0.125 * cnref_powi(relax, 2)))) {
                    pc = 310;
                } else {
                    {
                        double factor;
                        factor = ((-1.0) + sqrt(((1.0 + (8.0 * arg)))));
                        if (!isfinite(factor)) {
                            (*iflag) = (-4);
                            return;
                        }
                        if (fabs(((factor - 1.0))) < (0.1 * factor)) {
                            pc = 170;
                            continue;
                        }
                        if (factor < 0.5) {
                            factor = 0.5;
                        }
                        (relax /= factor);
                        pc = 320;
                    }
                }
            }
        } else if (pc == 310) {
            {
                if (relax >= 0.9) {
                    pc = 170;
                } else {
                    {
                        relax = 1.0;
                        pc = 320;
                    }
                }
            }
        } else if (pc == 320) {
            {
                if (relax < relmin) {
                    pc = 450;
                    continue;
                }
                if (corrections >= c->maximum_corrections) {
                    (*iflag) = (-2);
                    return;
                }
                (corrections += 1);
                double fact;
                fact = (relax - rlxold);
                for (int i = 1; i <= c->nz; i += 1) {
                    SET(z, i, (GET(z, i) + (fact * GET(delz, i))));
                }
                for (int i = 1; i <= c->ndmz; i += 1) {
                    SET(dmz, i, (GET(dmz, i) + (fact * GET(deldmz, i))));
                }
                rlxold = relax;
                pc = 250;
            }
        } else if (pc == 350) {
            {
                int large;
                large = 0;
                for (int it = 1; it <= c->ntol; it += 1) {
                    int iz;
                    iz = GET(c->ltol, it);
                    while ((iz <= c->nz)) {
                        if (fabs(GET(dqz, iz)) > (GET(c->tolin, it) * ((fabs(GET(z, iz)) + 1.0)))) {
                            large = 1;
                            break;
                        }
                        (iz += c->mstar);
                    }
                    if (large) {
                        break;
                    }
                }
                if (large) {
                    pc = 170;
                    continue;
                }
                for (int i = 1; i <= c->nz; i += 1) {
                    SET(z, i, (GET(z, i) + GET(dqz, i)));
                }
                for (int i = 1; i <= c->ndmz; i += 1) {
                    SET(dmz, i, (GET(dmz, i) + GET(dqdmz, i)));
                }
                pc = 390;
            }
        } else if (pc == 390) {
            {
                iconv = 1;
                if (c->icare == (-1)) {
                    c->icare = 0;
                }
                pc = 400;
            }
        } else if (pc == 400) {
            {
                if (!cnref_plain_finite(z, c->nz) || !cnref_plain_finite(dmz, c->ndmz)) {
                    (*iflag) = (-4);
                    return;
                }
                if (imesh == 2) {
                    cnref_errchk(c, xi, z, dmz, valstr, &(ifin));
                }
                if (c->invalid_evaluation) {
                    (*iflag) = (-4);
                    return;
                }
                if ((imesh == 1) || (((ifin == 0) && (c->icare != 2)))) {
                    pc = 460;
                } else {
                    {
                        (*iflag) = 1;
                        return;
                    }
                }
            }
        } else if (pc == 450) {
            {
                (*iflag) = (-2);
                (noconv += 1);
                if ((c->icare == 2) && (noconv > 1)) {
                    return;
                }
                if (c->icare == 0) {
                    c->icare = (-1);
                }
                pc = 460;
            }
        } else if (pc == 460) {
            {
                if (refinements >= c->maximum_refinements) {
                    (*iflag) = (-5);
                    return;
                }
                (refinements += 1);
                for (int i = 1; i <= (c->n + 1); i += 1) {
                    SET(xiold, i, GET(xi, i));
                }
                c->nold = c->n;
                if (((iconv == 0) || (c->mshnum >= c->mshlmt)) || (c->mshalt >= c->mshlmt)) {
                    imesh = 2;
                } else {
                    imesh = 1;
                }
                if ((c->mshalt >= c->mshlmt) && (c->mshnum < c->mshlmt)) {
                    c->mshalt = 1;
                }
                if ((c->k >= 2 && c->k <= 4) && (c->switch_ == 1)) {
                    imesh = 2;
                }
                if (c->k >= 2 && c->k <= 4) {
                    for (int i = 1; i <= (c->nold + 1); i += 1) {
                        SET(c->meshold, i, GET(xiold, i));
                    }
                    for (int i = 1; i <= (((c->nold + 1)) * c->mstar); i += 1) {
                        SET(c->zold, i, GET(z, i));
                    }
                    for (int i = 1; i <= c->ndmz; i += 1) {
                        SET(c->dmzold, i, GET(dmz, i));
                    }
                }
                cnref_newmsh(c, &(imesh), xi, xiold, z, dmz, valstr, slope, accum, nfxpnt, fixpnt);
                if (c->n > c->nmax) {
                    (c->n /= 2);
                    (*iflag) = (-1);
                    return;
                }
                if (!cnref_valid_mesh(xi, c->n)) {
                    (*iflag) = (-4);
                    return;
                }
                if (iconv == 0) {
                    imesh = 1;
                }
                if (c->icare == 1) {
                    iconv = 0;
                }
                pc = 20;
            }
        } else if (1) {
            abort();
        }
    }
}
void cnref_skale(cnref_context *c, int n, int mstar, int kd, R z, R xi, R scale, R dscale) {
    double basm_storage[5] = {0};
    R basm = basm_storage;
    SET(basm, 1, 1.0);
    for (int j = 1; j <= n; j += 1) {
        int iz;
        iz = 1;
        double h;
        h = (GET(xi, (j + 1)) - GET(xi, j));
        for (int l = 1; l <= c->mmax; l += 1) {
            SET(basm, (l + 1), ((GET(basm, l) * h) / ((double)(l))));
        }
        for (int icomp = 1; icomp <= c->ncomp; icomp += 1) {
            double scal;
            scal =
                ((((fabs(GET2(z, iz, j, mstar)) + fabs(GET2(z, iz, (j + 1), mstar)))) * 0.5) + 1.0);
            int mj;
            mj = GET(c->m, icomp);
            if ((!isfinite(scal)) || (scal <= 0.0)) {
                (c->invalid_evaluation = 1);
            }
            for (int l = 1; l <= mj; l += 1) {
                SET2(scale, iz, j, mstar, (GET(basm, l) / scal));
                (iz += 1);
            }
            scal = (GET(basm, (mj + 1)) / scal);
            if ((!isfinite(scal)) || (scal <= 0.0)) {
                (c->invalid_evaluation = 1);
            }
            int idmz;
            idmz = icomp;
            while ((idmz <= kd)) {
                SET2(dscale, idmz, j, kd, scal);
                (idmz += c->ncomp);
            }
        }
    }
    for (int iz = 1; iz <= mstar; iz += 1) {
        SET2(scale, iz, (n + 1), mstar, GET2(scale, iz, n, mstar));
    }
}
