void cnref_lsyslv(cnref_context *c, int *msing, R xi, R xiold, R z, R dmz, R delz, R deldmz, R g,
                  R w, R v, R rhs, R dmzo, I integs, I ipvtg, I ipvtw, double *rnorm, int mode);
void cnref_gderiv(cnref_context *c, R gi, int nrow, int irow, R zval, R dgz, int mode);
void cnref_vwblok(cnref_context *c, double xcol, double hrho, int jj, R wi, R vi, I ipvtw, int kd,
                  R zval, R df, R acol, R dmzo, int ncomp, int *msing);
void cnref_gblock(cnref_context *c, double h, R gi, int nrow, int irow, R wi, R vi, int kd, R rhsz,
                  R rhsdmz, I ipvtw, int mode);
void cnref_dmzsol(cnref_context *c, int kd, int mstar, int n, R v, R z, R dmz);
void cnref_fcblok(cnref_context *c, R bloks, I integs, int nbloks, I ipivot, R scrtch, int *info);
void cnref_factrb(cnref_context *c, R w, I ipivot, R d, int nrow, int ncol, int last, int *info);
void cnref_shiftb(cnref_context *c, R ai, int nrowi, int ncoli, int last, R ai1, int nrowi1,
                  int ncoli1);
void cnref_sbblok(cnref_context *c, R bloks, I integs, int nbloks, I ipivot, R x);
void cnref_subfor(cnref_context *c, R w, I ipivot, int nrow, int last, R x);
void cnref_subbak(cnref_context *c, R w, int nrow, int ncol, int last, R x);
void cnref_appsln(cnref_context *c, double x, R z, R fspace, I ispace);
void cnref_approx(cnref_context *c, int *i, double *x, R zval, R a, R coef, R xi, int n, R z, R dmz,
                  int k, int ncomp, int mmax, I m, int mstar, int mode, R dmval, int modm);
void cnref_rkbas(cnref_context *c, double s, R coef, int k, int m, R rkb, R dm, int mode);
void cnref_vmonde(cnref_context *c, R rho, R coef, int k);
void cnref_horder(cnref_context *c, int i, R uhigh, double hi, R dmz, int ncomp, int k);
void cnref_evaluate(cnref_context *c, int neqns, int nsub, R mesh, R y, double t, int s, int s_star,
                    R z, R k_discrete, R k_interp);
void cnref_setup(cnref_context *c, int _node, int neqns, I _orders, double tim1, double hi, R yim1,
                 R yi, int s, R ki_discrete, int s_star, R ki_interp, R y);
void cnref_sum_stages(cnref_context *c, int neqns, double hi, R yim1, int s, R ki_discrete,
                      int s_star, R ki_interp, R weights, R z);
void cnref_fsub_first(cnref_context *c, int node, int _nsol, I orders, R z, R f, R f_first);
void cnref_setup_phase1(cnref_context *c, int nsub, R mesh, R _asave, R zval, R dmz, int k_points,
                        int ncomp, int mmax, I m, int mstar, R y_super, R sk_discrete_super,
                        R sk_interp_super);
void cnref_tableau_star(cnref_context *c, int method);
void cnref_col_weights(cnref_context *c, int s_star, double t, R w);
void cnref_newmsh(cnref_context *c, int *mode, R xi, R xiold, R z, R dmz, R valstr, R slope,
                  R accum, int nfxpnt, R fixpnt);
int cnref_mesh_finite(cnref_context *c, double value);
int cnref_estimate_finite(cnref_context *c, double value, int preserve_scsa);
void cnref_consts(cnref_context *c, int k, R rho, R coef);
void cnref_errchk(cnref_context *c, R xi, R z, R dmz, R valstr, int *ifin);
void cnref_mesh_approx(cnref_context *c, int i, double x, R zval, int column, R xi, int n, R z,
                       R dmz, R dummy);
void cnref_contrl(cnref_context *c, R xi, R xiold, R z, R dmz, R rhs, R delz, R deldmz, R dqz,
                  R dqdmz, R g, R w, R v, R valstr, R slope, R scale, R dscale, R accum, I ipvtg,
                  I integs, I ipvtw, int nfxpnt, R fixpnt, int *iflag);
void cnref_skale(cnref_context *c, int n, int mstar, int kd, R z, R xi, R scale, R dscale);
