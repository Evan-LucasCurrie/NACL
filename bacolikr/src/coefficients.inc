/* LOWCOEFF and SCICOEFF tables, preserving source arithmetic. */
static void low_coeff(double x,size_t kcol,double h[2],double hd[2],double g[10]) {
memset(h,0,2*sizeof(double));memset(hd,0,2*sizeof(double));memset(g,0,10*sizeof(double));
    double phi;
    double eta1 = (x - 1.0e0) * (x - 1.0e0);
    double eta2 = x * x;
    double eta = eta1 * eta2;
    if (kcol == 2 ){
        h[0] = (1e0 + 2e0 * x) * eta1;
        h[1] = (1e0 - 2e0 * (x - 1.0)) * eta2;
        hd[0] = x * eta1;
        hd[1] = (x - 1e0) * eta2;
    } else if (kcol == 3 ){
        phi = x - 0.5e0;
        h[0] = -2.0e0 * (1.0e0 + 4.0e0 * x) * eta1 * phi;
        h[1] = 2.0e0 * (5.0e0 - 4.0e0 * x) * eta2 * phi;
        hd[0] = -2.0e0 * x * eta1 * phi;
        hd[1] = 2.0e0 * (x - 1.0e0) * eta2 * phi;
        g[0] = 16.0e0 * eta;
    } else if (kcol == 4 ){
        phi = (x - 0.3110177634953864e0) * (x - 0.6889822365046136e0);
        h[0] = (0.4666666666666667e1 + 0.3111111111111111e2 * x) * eta1 * phi;
        h[1] = (0.3577777777777778e2 - 0.3111111111111111e2 * x) * eta2 * phi;
        hd[0] = 0.4666666666666667e1 * x * eta1 * phi;
        hd[1] = (0.4666666666666667e1 * x - 0.4666666666666667e1) * eta2 * phi;
        g[0] = -0.5761858410762886e2 * eta;
        g[1] = -g[0] * (x - 0.3110177634953864e0);
        g[0] = g[0] * (x - 0.6889822365046136e0);
    } else if (kcol == 5 ){
        phi = (x - 0.2113248654051871e0) * (x - 0.5e0) * (x - 0.7886751345948129e0);
        h[0] = (-0.12e2 - 0.12e3 * x) * eta1 * phi;
        h[1] = (0.132e3 - 0.12e3 * x) * eta2 * phi;
        hd[0] = -0.12e2 * x * eta1 * phi;
        hd[1] = (0.12e2 * x - 0.12e2) * eta2 * phi;
        g[0] = 0.216e3 * eta * phi;
        g[1] = -192.0e0 * eta * phi / (x - 0.5e0);
        g[2] = g[0] / (x - 0.7886751345948129e0);
        g[0] = g[0] / (x - 0.2113248654051871e0);
    } else if (kcol == 6 ){
        phi = (x - 0.1526267046965671e0)
            * (x - 0.3747185964571342e0)
            * (x - 0.6252814035428658e0)
            * (x - 0.8473732953034329e0);
        h[0] = (0.33e2 + 0.462e3 * x) * eta1 * phi;
        h[1] = (0.495e3 - 0.462e3 * x) * eta2 * phi;
        hd[0] = 0.33e2 * x * eta1 * phi;
        hd[1] = (-0.33e2 + 0.33e2 * x) * eta2 * phi;
        g[0] = -0.8197591840300222e3 * eta * phi;
        g[1] = 0.6925405260332655e3 * eta * phi;
        g[2] = -g[1] / (x - 0.6252814035428658e0);
        g[3] = -g[0] / (x - 0.8473732953034329e0);
        g[0] = g[0] / (x - 0.1526267046965671e0);
        g[1] = g[1] / (x - 0.3747185964571342e0);
    } else if (kcol == 7 ){
        phi = (x - 0.1152723378341063e0)
            * (x - 0.2895425974880943e0)
            * (x - 0.5e0)
            * (x - 0.7104574025119057e0)
            * (x - 0.8847276621658937e0);
        h[0] = (-0.9533333333333333e2 - 0.1779555555555556e4 * x) * eta1 * phi;
        h[1] = (0.1874888888888889e4 - 0.1779555555555556e4 * x) * eta2 * phi;
        hd[0] = -0.9533333333333333e2 * x * eta1 * phi;
        hd[1] = (-0.9533333333333333e2 + 0.9533333333333333e2 * x) * eta2 * phi;
        g[0] = 0.3131255164938139e4 * eta * phi;
        g[1] = -0.257196627604925e4 * eta * phi;
        g[2] = 0.2440533333333333e4 * eta * phi / (x - 0.5e0);
        g[3] = g[1] / (x - 0.7104574025119057e0);
        g[4] = g[0] / (x - 0.8847276621658937e0);
        g[0] = g[0] / (x - 0.1152723378341063e0);
        g[1] = g[1] / (x - 0.2895425974880943e0);
    } else if (kcol == 8 ){
        phi = (x - 0.9007700226825652e-1)
            * (x - 0.2296976813063206e0)
            * (x - 0.405661288754607e0)
            * (x - 0.594338711245393e0)
            * (x - 0.7703023186936794e0)
            * (x - 0.9099229977317435e0);
        h[0] = (0.286e3 + 0.6864e4 * x) * eta1 * phi;
        h[1] = (0.715e4 - 0.6864e4 * x) * eta2 * phi;
        hd[0] = 0.286e3 * x * eta1 * phi;
        hd[1] = (-0.286e3 + 0.286e3 * x) * eta2 * phi;
        g[0] = -0.1201314415336355e5 * eta * phi;
        g[1] = 0.9696023778080566e4 * eta * phi;
        g[2] = -0.8929458911121293e4 * eta * phi;
        g[3] = -g[2] / (x - 0.594338711245393e0);
        g[4] = -g[1] / (x - 0.7703023186936794e0);
        g[5] = -g[0] / (x - 0.9099229977317435e0);
        g[0] = g[0] / (x - 0.9007700226825652e-1);
        g[1] = g[1] / (x - 0.2296976813063206e0);
        g[2] = g[2] / (x - 0.405661288754607e0);
    } else if (kcol == 9 ){
        phi = (x - 0.7229898685756272e-1)
            * (x - 0.1863109301186906e0)
            * (x - 0.3341852231986051e0)
            * (x - 0.5e0)
            * (x - 0.6658147768013949e0)
            * (x - 0.8136890698813094e0)
            * (x - 0.9277010131424373e0);
        h[0] = (-0.884e3 - 0.2652e5 * x) * eta1 * phi;
        h[1] = (0.27404e5 - 0.2652e5 * x) * eta2 * phi;
        hd[0] = -0.884e3 * x * eta1 * phi;
        hd[1] = (-0.884e3 + 0.884e3 * x) * eta2 * phi;
        g[0] = 0.4624522420172525e5 * eta * phi;
        g[1] = -0.3688889968687672e5 * eta * phi;
        g[2] = 0.3332824691372289e5 * eta * phi;
        g[3] = -0.3232914285714286e5 * eta * phi / (x - 0.5e0);
        g[4] = g[2] / (x - 0.6658147768013949e0);
        g[5] = g[1] / (x - 0.8136890698813094e0);
        g[6] = g[0] / (x - 0.9277010131424373e0);
        g[0] = g[0] / (x - 0.7229898685756272e-1);
        g[1] = g[1] / (x - 0.1863109301186906e0);
        g[2] = g[2] / (x - 0.3341852231986051e0);
    } else {
        h[0] = 0.0e0;
        h[1] = 0.0e0;
        hd[0] = 0.0e0;
        hd[1] = 0.0e0;
        for (size_t i=1;i<=7;i++) {
            g[i - 1] = 0.0;
        }
    }
}
static void sci_coeff(size_t kcol,const double *mesh,size_t nint,const double *points,double (*h)[2],double (*hd)[2],double (*g)[10]) {
size_t quad=kcol+3;
    double eta[13][3]={{0}};
    double phi;
    double a;
    double aa;
    double b;
    double bb;
    double hx;
    double gamma1;
    double gamma2;
    for (size_t j=1;j<=kcol + 3;j++) {
        eta[j - 1][0] = (points[j - 1] - 1.0e0) * (points[j - 1] - 1.0e0);
        eta[j - 1][1] = points[j - 1] * points[j - 1];
        eta[j - 1][2] = eta[j - 1][0] * eta[j - 1][1];
    }
    if (kcol == 3 ){
        for (size_t i=1;i<=nint;i++) {
            hx = mesh[i + 1 - 1] - mesh[i - 1];
            if (i != 1 && i != nint ){
                a = -0.5e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                b = 1.0e0 + 0.5e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -1.0e0 / a - 4.0e0 - 1.0e0 / b;
                gamma2 = 1.0e0 / (1.0e0 - a) + 4.0e0 + 1.0e0 / (1.0e0 - b);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -2.0e0 / a / b;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 2.0e0 / (1.0e0 - a) / (1.0e0 - b);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] =
                    1.0e0 / ((a - 0.5e0) * (a - b) * a * a * (a - 1.0e0) * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 16.0e0 / (0.5e0 - a) / (0.5e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    1.0e0 / ((b - a) * (b - 0.5e0) * b * b * (b - 1.0e0) * (b - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - a) * (points[j - 1] - 0.5e0) * (points[j - 1] - b);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                }
            } else if (i == 1 ){
                b = 1.0e0 + 0.5e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                bb = 1.0e0 + (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -4.0e0 - 1.0e0 / b - 1.0e0 / bb;
                gamma2 = 4.0e0 + 1.0e0 / (1.0e0 - b) + 1.0e0 / (1.0e0 - bb);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -2.0e0 / b / bb;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 2.0e0 / (1.0e0 - b) / (1.0e0 - bb);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 16.0e0 / (0.5e0 - bb) / (0.5e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    1.0e0 / ((b - bb) * (b - 0.5e0) * b * b * (b - 1.0e0) * (b - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    1.0e0 / ((bb - b) * (bb - 0.5e0) * bb * bb * (bb - 1.0e0) * (bb - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - 0.5e0) * (points[j - 1] - b) * (points[j - 1] - bb);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - bb);
                }
            } else {
                aa = -(mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                a = -0.5e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                gamma1 = -1.0e0 / aa - 1.0e0 / a - 4.0e0;
                gamma2 = 1.0e0 / (1.0e0 - aa) + 1.0e0 / (1.0e0 - a) + 4.0e0;
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -2.0e0 / aa / a;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 2.0e0 / (1.0e0 - aa) / (1.0e0 - a);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] =
                    1.0e0 / ((aa - 0.5e0) * (aa - a) * aa * aa * (aa - 1.0e0) * (aa - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    1.0e0 / ((a - 0.5e0) * (a - aa) * a * a * (a - 1.0e0) * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = 16.0e0 / (0.5e0 - a) / (0.5e0 - aa);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - aa) * (points[j - 1] - a) * (points[j - 1] - 0.5e0);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - aa);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                }
            }
        }
    } else if (kcol == 4 ){
        for (size_t i=1;i<=nint;i++) {
            hx = mesh[i + 1 - 1] - mesh[i - 1];
            if (i != 1 && i != nint ){
                a = -0.3110177634953864e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                b = 1.0e0 + 0.3110177634953864e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -1.0e0 / a - 0.6666666666666667e1 - 1.0e0 / b;
                gamma2 = 1.0e0 / (1.0e0 - a) + 0.6666666666666667e1 + 1.0e0 / (1.0e0 - b);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.4666666666666667e1 / a / b;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.4666666666666667e1 / (1.0e0 - a) / (1.0e0 - b);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / (a - 0.3110177634953864e0)
                    / (a - 0.6889822365046136e0)
                    / (a - b)
                    / a
                    / a
                    / (a - 1.0e0)
                    / (a - 1.0e0);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    -0.5761858410762886e2 / (0.3110177634953864e0 - a) / (0.3110177634953864e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.5761858410762886e2 / (0.6889822365046136e0 - a) / (0.6889822365046136e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] = 1.0e0
                    / (b - a)
                    / (b - 0.3110177634953864e0)
                    / (b - 0.6889822365046136e0)
                    / b
                    / b
                    / (b - 1.0e0)
                    / (b - 1.0e0);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - a)
                        * (points[j - 1] - 0.3110177634953864e0)
                        * (points[j - 1] - 0.6889822365046136e0)
                        * (points[j - 1] - b);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3110177634953864e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6889822365046136e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                }
            } else if (i == 1 ){
                b = 1.0e0 + 0.3110177634953864e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                bb = 1.0e0 + 0.6889822365046136e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -0.6666666666666667e1 - 1.0e0 / b - 1.0e0 / bb;
                gamma2 = 0.6666666666666667e1 + 1.0e0 / (1.0e0 - b) + 1.0e0 / (1.0e0 - bb);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.4666666666666667e1 / b / bb;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.4666666666666667e1 / (1.0e0 - b) / (1.0e0 - bb);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.5761858410762886e2
                    / (0.3110177634953864e0 - b)
                    / (0.3110177634953864e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.5761858410762886e2 / (0.6889822365046136e0 - b) / (0.6889822365046136e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = 1.0e0
                    / (b - bb)
                    / (b - 0.3110177634953864e0)
                    / (b - 0.6889822365046136e0)
                    / b
                    / b
                    / (b - 1.0e0)
                    / (b - 1.0e0);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] = 1.0e0
                    / (bb - b)
                    / (bb - 0.3110177634953864e0)
                    / (bb - 0.6889822365046136e0)
                    / bb
                    / bb
                    / (bb - 1.0e0)
                    / (bb - 1.0e0);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - 0.3110177634953864e0)
                        * (points[j - 1] - 0.6889822365046136e0)
                        * (points[j - 1] - b)
                        * (points[j - 1] - bb);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3110177634953864e0);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6889822365046136e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - bb);
                }
            } else {
                aa = -0.6889822365046136e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                a = -0.3110177634953864e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                gamma1 = -1.0e0 / aa - 1.0e0 / a - 0.6666666666666667e1;
                gamma2 = 1.0e0 / (1.0e0 - aa) + 1.0e0 / (1.0e0 - a) + 0.6666666666666667e1;
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.4666666666666667e1 / aa / a;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.4666666666666667e1 / (1.0e0 - aa) / (1.0e0 - a);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / (aa - a)
                    / (aa - 0.3110177634953864e0)
                    / (aa - 0.6889822365046136e0)
                    / aa
                    / aa
                    / (aa - 1.0e0)
                    / (aa - 1.0e0);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 1.0e0
                    / (a - aa)
                    / (a - 0.3110177634953864e0)
                    / (a - 0.6889822365046136e0)
                    / a
                    / a
                    / (a - 1.0e0)
                    / (a - 1.0e0);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = -0.5761858410762886e2
                    / (0.3110177634953864e0 - a)
                    / (0.3110177634953864e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.5761858410762886e2 / (0.6889822365046136e0 - a) / (0.6889822365046136e0 - aa);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - aa)
                        * (points[j - 1] - a)
                        * (points[j - 1] - 0.3110177634953864e0)
                        * (points[j - 1] - 0.6889822365046136e0);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - aa);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3110177634953864e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6889822365046136e0);
                }
            }
        }
    } else if (kcol == 5 ){
        for (size_t i=1;i<=nint;i++) {
            hx = mesh[i + 1 - 1] - mesh[i - 1];
            if (i != 1 && i != nint ){
                a = -0.2113248654051871e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                b = 1.0e0 + 0.2113248654051871e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -1.0e0 / a - 0.1e2 - 1.0e0 / b;
                gamma2 = 1.0e0 / (1.0e0 - a) + 0.1e2 + 1.0e0 / (1.0e0 - b);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.12e2 / a / b;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.12e2 / (1.0e0 - a) / (1.0e0 - b);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((a - 0.2113248654051871e0)
                        * (a - 0.5e0)
                        * (a - 0.7886751345948129e0)
                        * (a - b)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.216e3 / (0.2113248654051871e0 - a) / (0.2113248654051871e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = -192.0e0 / (0.5e0 - a) / (0.5e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.216e3 / (0.7886751345948129e0 - a) / (0.7886751345948129e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] = 1.0e0
                    / ((b - a)
                        * (b - 0.2113248654051871e0)
                        * (b - 0.5e0)
                        * (b - 0.7886751345948129e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - a)
                        * (points[j - 1] - 0.2113248654051871e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.7886751345948129e0)
                        * (points[j - 1] - b);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2113248654051871e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7886751345948129e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                }
            } else if (i == 1 ){
                b = 1.0e0 + 0.2113248654051871e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                bb = 1.0e0 + 0.5e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -0.1e2 - 1.0e0 / b - 1.0e0 / bb;
                gamma2 = 0.1e2 + 1.0e0 / (1.0e0 - b) + 1.0e0 / (1.0e0 - bb);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.12e2 / b / bb;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.12e2 / (1.0e0 - b) / (1.0e0 - bb);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] =
                    0.216e3 / (0.2113248654051871e0 - b) / (0.2113248654051871e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = -192.0e0 / (0.5e0 - b) / (0.5e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.216e3 / (0.7886751345948129e0 - b) / (0.7886751345948129e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] = 1.0e0
                    / ((b - 0.2113248654051871e0)
                        * (b - 0.5e0)
                        * (b - 0.7886751345948129e0)
                        * (b - bb)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][4] = 1.0e0
                    / ((bb - 0.2113248654051871e0)
                        * (bb - 0.5e0)
                        * (bb - 0.7886751345948129e0)
                        * (bb - b)
                        * bb
                        * bb
                        * (bb - 1.0e0)
                        * (bb - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - 0.2113248654051871e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.7886751345948129e0)
                        * (points[j - 1] - b)
                        * (points[j - 1] - bb);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2113248654051871e0);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7886751345948129e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - bb);
                }
            } else {
                aa = -0.5e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                a = -0.2113248654051871e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                gamma1 = -1.0e0 / aa - 1.0e0 / a - 0.1e2;
                gamma2 = 1.0e0 / (1.0e0 - aa) + 1.0e0 / (1.0e0 - a) + 0.1e2;
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.12e2 / aa / a;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.12e2 / (1.0e0 - aa) / (1.0e0 - a);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((aa - 0.2113248654051871e0)
                        * (aa - 0.5e0)
                        * (aa - 0.7886751345948129e0)
                        * (aa - a)
                        * aa
                        * aa
                        * (aa - 1.0e0)
                        * (aa - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 1.0e0
                    / ((a - 0.2113248654051871e0)
                        * (a - 0.5e0)
                        * (a - 0.7886751345948129e0)
                        * (a - aa)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.216e3 / (0.2113248654051871e0 - a) / (0.2113248654051871e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] = -192.0e0 / (0.5e0 - a) / (0.5e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    0.216e3 / (0.7886751345948129e0 - a) / (0.7886751345948129e0 - aa);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - aa)
                        * (points[j - 1] - a)
                        * (points[j - 1] - 0.2113248654051871e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.7886751345948129e0);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - aa);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2113248654051871e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7886751345948129e0);
                }
            }
        }
    } else if (kcol == 6 ){
        for (size_t i=1;i<=nint;i++) {
            hx = mesh[i + 1 - 1] - mesh[i - 1];
            if (i != 1 && i != nint ){
                a = -0.1526267046965671e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                b = 1.0e0 + 0.1526267046965671e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -1.0e0 / a - 0.14e2 - 1.0e0 / b;
                gamma2 = 1.0e0 / (1.0e0 - a) + 0.14e2 + 1.0e0 / (1.0e0 - b);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.33e2 / a / b;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.33e2 / (1.0e0 - a) / (1.0e0 - b);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((a - 0.1526267046965671e0)
                        * (a - 0.3747185964571342e0)
                        * (a - 0.6252814035428658e0)
                        * (a - 0.8473732953034329e0)
                        * (a - b)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    -0.8197591840300222e3 / (0.1526267046965671e0 - a) / (0.1526267046965671e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.6925405260332655e3 / (0.3747185964571342e0 - a) / (0.3747185964571342e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    -0.6925405260332655e3 / (0.6252814035428658e0 - a) / (0.6252814035428658e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    0.8197591840300222e3 / (0.8473732953034329e0 - a) / (0.8473732953034329e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] = 1.0e0
                    / ((b - a)
                        * (b - 0.1526267046965671e0)
                        * (b - 0.3747185964571342e0)
                        * (b - 0.6252814035428658e0)
                        * (b - 0.8473732953034329e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - a)
                        * (points[j - 1] - 0.1526267046965671e0)
                        * (points[j - 1] - 0.3747185964571342e0)
                        * (points[j - 1] - 0.6252814035428658e0)
                        * (points[j - 1] - 0.8473732953034329e0)
                        * (points[j - 1] - b);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1526267046965671e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3747185964571342e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6252814035428658e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8473732953034329e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                }
            } else if (i == 1 ){
                b = 1.0e0 + 0.1526267046965671e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                bb = 1.0e0 + 0.3747185964571342e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -0.14e2 - 1.0e0 / b - 1.0e0 / bb;
                gamma2 = 0.14e2 + 1.0e0 / (1.0e0 - b) + 1.0e0 / (1.0e0 - bb);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.33e2 / b / bb;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.33e2 / (1.0e0 - b) / (1.0e0 - bb);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.8197591840300222e3
                    / (0.1526267046965671e0 - b)
                    / (0.1526267046965671e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.6925405260332655e3 / (0.3747185964571342e0 - b) / (0.3747185964571342e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = -0.6925405260332655e3
                    / (0.6252814035428658e0 - b)
                    / (0.6252814035428658e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.8197591840300222e3 / (0.8473732953034329e0 - b) / (0.8473732953034329e0 - bb);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] = 1.0e0
                    / ((b - bb)
                        * (b - 0.1526267046965671e0)
                        * (b - 0.3747185964571342e0)
                        * (b - 0.6252814035428658e0)
                        * (b - 0.8473732953034329e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][5] = 1.0e0
                    / ((bb - b)
                        * (bb - 0.1526267046965671e0)
                        * (bb - 0.3747185964571342e0)
                        * (bb - 0.6252814035428658e0)
                        * (bb - 0.8473732953034329e0)
                        * bb
                        * bb
                        * (bb - 1.0e0)
                        * (bb - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - 0.1526267046965671e0)
                        * (points[j - 1] - 0.3747185964571342e0)
                        * (points[j - 1] - 0.6252814035428658e0)
                        * (points[j - 1] - 0.8473732953034329e0)
                        * (points[j - 1] - b)
                        * (points[j - 1] - bb);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1526267046965671e0);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3747185964571342e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6252814035428658e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8473732953034329e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - bb);
                }
            } else {
                aa = -0.3747185964571342e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                a = -0.1526267046965671e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                gamma1 = -1.0e0 / aa - 1.0e0 / a - 0.14e2;
                gamma2 = 1.0e0 / (1.0e0 - aa) + 1.0e0 / (1.0e0 - a) + 0.14e2;
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.33e2 / aa / a;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.33e2 / (1.0e0 - aa) / (1.0e0 - a);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((aa - 0.1526267046965671e0)
                        * (aa - 0.3747185964571342e0)
                        * (aa - 0.6252814035428658e0)
                        * (aa - 0.8473732953034329e0)
                        * (aa - a)
                        * aa
                        * aa
                        * (aa - 1.0e0)
                        * (aa - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 1.0e0
                    / ((a - 0.1526267046965671e0)
                        * (a - 0.3747185964571342e0)
                        * (a - 0.6252814035428658e0)
                        * (a - 0.8473732953034329e0)
                        * (a - aa)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = -0.8197591840300222e3
                    / (0.1526267046965671e0 - a)
                    / (0.1526267046965671e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.6925405260332655e3 / (0.3747185964571342e0 - a) / (0.3747185964571342e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] = -0.6925405260332655e3
                    / (0.6252814035428658e0 - a)
                    / (0.6252814035428658e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    0.8197591840300222e3 / (0.8473732953034329e0 - a) / (0.8473732953034329e0 - aa);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - aa)
                        * (points[j - 1] - a)
                        * (points[j - 1] - 0.1526267046965671e0)
                        * (points[j - 1] - 0.3747185964571342e0)
                        * (points[j - 1] - 0.6252814035428658e0)
                        * (points[j - 1] - 0.8473732953034329e0);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - aa);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1526267046965671e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3747185964571342e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6252814035428658e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8473732953034329e0);
                }
            }
        }
    } else if (kcol == 7 ){
        for (size_t i=1;i<=nint;i++) {
            hx = mesh[i + 1 - 1] - mesh[i - 1];
            if (i != 1 && i != nint ){
                a = -0.1152723378341063e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                b = 1.0e0 + 0.1152723378341063e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -1.0e0 / a - 0.1866666666666667e2 - 1.0e0 / b;
                gamma2 = 1.0e0 / (1.0e0 - a) + 0.1866666666666667e2 + 1.0e0 / (1.0e0 - b);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.9533333333333333e2 / a / b;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.9533333333333333e2 / (1.0e0 - a) / (1.0e0 - b);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((a - 0.1152723378341063e0)
                        * (a - 0.2895425974880943e0)
                        * (a - 0.5e0)
                        * (a - 0.7104574025119057e0)
                        * (a - 0.8847276621658937e0)
                        * (a - b)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.3131255164938139e4 / (0.1152723378341063e0 - a) / (0.1152723378341063e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    -0.257196627604925e4 / (0.2895425974880943e0 - a) / (0.2895425974880943e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.2440533333333333e4 / (0.5e0 - a) / (0.5e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    -0.257196627604925e4 / (0.7104574025119057e0 - a) / (0.7104574025119057e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    0.3131255164938139e4 / (0.8847276621658937e0 - a) / (0.8847276621658937e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] = 1.0e0
                    / ((b - a)
                        * (b - 0.1152723378341063e0)
                        * (b - 0.2895425974880943e0)
                        * (b - 0.5e0)
                        * (b - 0.7104574025119057e0)
                        * (b - 0.8847276621658937e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - a)
                        * (points[j - 1] - 0.1152723378341063e0)
                        * (points[j - 1] - 0.2895425974880943e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.7104574025119057e0)
                        * (points[j - 1] - 0.8847276621658937e0)
                        * (points[j - 1] - b);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1152723378341063e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2895425974880943e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7104574025119057e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8847276621658937e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                }
            } else if (i == 1 ){
                b = 1.0e0 + 0.1152723378341063e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                bb = 1.0e0 + 0.2895425974880943e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -0.1866666666666667e2 - 1.0e0 / b - 1.0e0 / bb;
                gamma2 = 0.1866666666666667e2 + 1.0e0 / (1.0e0 - b) + 1.0e0 / (1.0e0 - bb);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.9533333333333333e2 / b / bb;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.9533333333333333e2 / (1.0e0 - b) / (1.0e0 - bb);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] =
                    0.3131255164938139e4 / (0.1152723378341063e0 - bb) / (0.1152723378341063e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    -0.257196627604925e4 / (0.2895425974880943e0 - bb) / (0.2895425974880943e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.2440533333333333e4 / (0.5e0 - bb) / (0.5e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    -0.257196627604925e4 / (0.7104574025119057e0 - bb) / (0.7104574025119057e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    0.3131255164938139e4 / (0.8847276621658937e0 - bb) / (0.8847276621658937e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] = 1.0e0
                    / ((b - bb)
                        * (b - 0.1152723378341063e0)
                        * (b - 0.2895425974880943e0)
                        * (b - 0.5e0)
                        * (b - 0.7104574025119057e0)
                        * (b - 0.8847276621658937e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][6] = 1.0e0
                    / ((bb - b)
                        * (bb - 0.1152723378341063e0)
                        * (bb - 0.2895425974880943e0)
                        * (bb - 0.5e0)
                        * (bb - 0.7104574025119057e0)
                        * (bb - 0.8847276621658937e0)
                        * bb
                        * bb
                        * (bb - 1.0e0)
                        * (bb - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - 0.1152723378341063e0)
                        * (points[j - 1] - 0.2895425974880943e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.7104574025119057e0)
                        * (points[j - 1] - 0.8847276621658937e0)
                        * (points[j - 1] - b)
                        * (points[j - 1] - bb);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1152723378341063e0);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2895425974880943e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7104574025119057e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8847276621658937e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - bb);
                }
            } else {
                aa = -0.2895425974880943e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                a = -0.1152723378341063e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                gamma1 = -1.0e0 / aa - 1.0e0 / a - 0.1866666666666667e2;
                gamma2 = 1.0e0 / (1.0e0 - aa) + 1.0e0 / (1.0e0 - a) + 0.1866666666666667e2;
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.9533333333333333e2 / aa / a;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.9533333333333333e2 / (1.0e0 - aa) / (1.0e0 - a);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((aa - 0.1152723378341063e0)
                        * (aa - 0.2895425974880943e0)
                        * (aa - 0.5e0)
                        * (aa - 0.7104574025119057e0)
                        * (aa - 0.8847276621658937e0)
                        * (aa - a)
                        * aa
                        * aa
                        * (aa - 1.0e0)
                        * (aa - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 1.0e0
                    / ((a - 0.1152723378341063e0)
                        * (a - 0.2895425974880943e0)
                        * (a - 0.5e0)
                        * (a - 0.7104574025119057e0)
                        * (a - 0.8847276621658937e0)
                        * (a - aa)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.3131255164938139e4 / (0.1152723378341063e0 - a) / (0.1152723378341063e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    -0.257196627604925e4 / (0.2895425974880943e0 - a) / (0.2895425974880943e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    0.2440533333333333e4 / (0.5e0 - a) / (0.5e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    -0.257196627604925e4 / (0.7104574025119057e0 - a) / (0.7104574025119057e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] =
                    0.3131255164938139e4 / (0.8847276621658937e0 - a) / (0.8847276621658937e0 - aa);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - aa)
                        * (points[j - 1] - a)
                        * (points[j - 1] - 0.1152723378341063e0)
                        * (points[j - 1] - 0.2895425974880943e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.7104574025119057e0)
                        * (points[j - 1] - 0.8847276621658937e0);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - aa);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1152723378341063e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2895425974880943e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7104574025119057e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8847276621658937e0);
                }
            }
        }
    } else if (kcol == 8 ){
        for (size_t i=1;i<=nint;i++) {
            hx = mesh[i + 1 - 1] - mesh[i - 1];
            if (i != 1 && i != nint ){
                a = -0.9007700226825652e-1 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                b = 1.0e0 + 0.9007700226825652e-1 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -1.0e0 / a - 0.24e2 - 1.0e0 / b;
                gamma2 = 1.0e0 / (1.0e0 - a) + 0.24e2 + 1.0e0 / (1.0e0 - b);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.286e3 / a / b;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.286e3 / (1.0e0 - a) / (1.0e0 - b);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((a - 0.9007700226825652e-1)
                        * (a - 0.2296976813063206e0)
                        * (a - 0.405661288754607e0)
                        * (a - 0.594338711245393e0)
                        * (a - 0.7703023186936794e0)
                        * (a - 0.9099229977317435e0)
                        * (a - b)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = -0.1201314415336355e5
                    / (0.9007700226825652e-1 - a)
                    / (0.9007700226825652e-1 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.9696023778080566e4 / (0.2296976813063206e0 - a) / (0.2296976813063206e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    -0.8929458911121293e4 / (0.405661288754607e0 - a) / (0.405661288754607e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    0.8929458911121293e4 / (0.594338711245393e0 - a) / (0.594338711245393e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    -0.9696023778080565e4 / (0.7703023186936794e0 - a) / (0.7703023186936794e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] =
                    0.1201314415336355e5 / (0.9099229977317435e0 - a) / (0.9099229977317435e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][7] = 1.0e0
                    / ((b - a)
                        * (b - 0.9007700226825652e-1)
                        * (b - 0.2296976813063206e0)
                        * (b - 0.405661288754607e0)
                        * (b - 0.594338711245393e0)
                        * (b - 0.7703023186936794e0)
                        * (b - 0.9099229977317435e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - a)
                        * (points[j - 1] - 0.9007700226825652e-1)
                        * (points[j - 1] - 0.2296976813063206e0)
                        * (points[j - 1] - 0.405661288754607e0)
                        * (points[j - 1] - 0.594338711245393e0)
                        * (points[j - 1] - 0.7703023186936794e0)
                        * (points[j - 1] - 0.9099229977317435e0)
                        * (points[j - 1] - b);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9007700226825652e-1);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2296976813063206e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.405661288754607e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.594338711245393e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7703023186936794e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9099229977317435e0);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                }
            } else if (i == 1 ){
                b = 1.0e0 + 0.9007700226825652e-1 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                bb = 1.0e0 + 0.2296976813063206e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -0.24e2 - 1.0e0 / b - 1.0e0 / bb;
                gamma2 = 0.24e2 + 1.0e0 / (1.0e0 - b) + 1.0e0 / (1.0e0 - bb);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.286e3 / b / bb;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.286e3 / (1.0e0 - b) / (1.0e0 - bb);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.1201314415336355e5
                    / (0.9007700226825652e-1 - bb)
                    / (0.9007700226825652e-1 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.9696023778080566e4 / (0.2296976813063206e0 - bb) / (0.2296976813063206e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    -0.8929458911121293e4 / (0.405661288754607e0 - bb) / (0.405661288754607e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.8929458911121293e4 / (0.594338711245393e0 - bb) / (0.594338711245393e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] = -0.9696023778080565e4
                    / (0.7703023186936794e0 - bb)
                    / (0.7703023186936794e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    0.1201314415336355e5 / (0.9099229977317435e0 - bb) / (0.9099229977317435e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] = 1.0e0
                    / ((b - bb)
                        * (b - 0.9007700226825652e-1)
                        * (b - 0.2296976813063206e0)
                        * (b - 0.405661288754607e0)
                        * (b - 0.594338711245393e0)
                        * (b - 0.7703023186936794e0)
                        * (b - 0.9099229977317435e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][7] = 1.0e0
                    / ((bb - b)
                        * (bb - 0.9007700226825652e-1)
                        * (bb - 0.2296976813063206e0)
                        * (bb - 0.405661288754607e0)
                        * (bb - 0.594338711245393e0)
                        * (bb - 0.7703023186936794e0)
                        * (bb - 0.9099229977317435e0)
                        * bb
                        * bb
                        * (bb - 1.0e0)
                        * (bb - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - 0.9007700226825652e-1)
                        * (points[j - 1] - 0.2296976813063206e0)
                        * (points[j - 1] - 0.405661288754607e0)
                        * (points[j - 1] - 0.594338711245393e0)
                        * (points[j - 1] - 0.7703023186936794e0)
                        * (points[j - 1] - 0.9099229977317435e0)
                        * (points[j - 1] - b)
                        * (points[j - 1] - bb);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9007700226825652e-1);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2296976813063206e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.405661288754607e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.594338711245393e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7703023186936794e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9099229977317435e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - bb);
                }
            } else {
                aa = -0.2296976813063206e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                a = -0.9007700226825652e-1 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                gamma1 = -1.0e0 / aa - 1.0e0 / a - 0.24e2;
                gamma2 = 1.0e0 / (1.0e0 - aa) + 1.0e0 / (1.0e0 - a) + 0.24e2;
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.286e3 / aa / a;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.286e3 / (1.0e0 - aa) / (1.0e0 - a);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((aa - 0.9007700226825652e-1)
                        * (aa - 0.2296976813063206e0)
                        * (aa - 0.405661288754607e0)
                        * (aa - 0.594338711245393e0)
                        * (aa - 0.7703023186936794e0)
                        * (aa - 0.9099229977317435e0)
                        * (aa - a)
                        * aa
                        * aa
                        * (aa - 1.0e0)
                        * (aa - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 1.0e0
                    / ((a - 0.9007700226825652e-1)
                        * (a - 0.2296976813063206e0)
                        * (a - 0.405661288754607e0)
                        * (a - 0.594338711245393e0)
                        * (a - 0.7703023186936794e0)
                        * (a - 0.9099229977317435e0)
                        * (a - aa)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = -0.1201314415336355e5
                    / (0.9007700226825652e-1 - a)
                    / (0.9007700226825652e-1 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.9696023778080566e4 / (0.2296976813063206e0 - a) / (0.2296976813063206e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    -0.8929458911121293e4 / (0.405661288754607e0 - a) / (0.405661288754607e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    0.8929458911121293e4 / (0.594338711245393e0 - a) / (0.594338711245393e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] = -0.9696023778080565e4
                    / (0.7703023186936794e0 - a)
                    / (0.7703023186936794e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][7] =
                    0.1201314415336355e5 / (0.9099229977317435e0 - a) / (0.9099229977317435e0 - aa);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - aa)
                        * (points[j - 1] - a)
                        * (points[j - 1] - 0.9007700226825652e-1)
                        * (points[j - 1] - 0.2296976813063206e0)
                        * (points[j - 1] - 0.405661288754607e0)
                        * (points[j - 1] - 0.594338711245393e0)
                        * (points[j - 1] - 0.7703023186936794e0)
                        * (points[j - 1] - 0.9099229977317435e0);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - aa);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9007700226825652e-1);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2296976813063206e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.405661288754607e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.594338711245393e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7703023186936794e0);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9099229977317435e0);
                }
            }
        }
    } else if (kcol == 9 ){
        for (size_t i=1;i<=nint;i++) {
            hx = mesh[i + 1 - 1] - mesh[i - 1];
            if (i != 1 && i != nint ){
                a = -0.7229898685756272e-1 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                b = 1.0e0 + 0.7229898685756272e-1 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -1.0e0 / a - 0.3e2 - 1.0e0 / b;
                gamma2 = 1.0e0 / (1.0e0 - a) + 0.3e2 + 1.0e0 / (1.0e0 - b);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.884e3 / a / b;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.884e3 / (1.0e0 - a) / (1.0e0 - b);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((a - 0.7229898685756272e-1)
                        * (a - 0.1863109301186906e0)
                        * (a - 0.3341852231986051e0)
                        * (a - 0.5e0)
                        * (a - 0.6658147768013949e0)
                        * (a - 0.8136890698813094e0)
                        * (a - 0.9277010131424373e0)
                        * (a - b)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.4624522420172525e5
                    / (0.7229898685756272e-1 - a)
                    / (0.7229898685756272e-1 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    -0.3688889968687672e5 / (0.1863109301186906e0 - a) / (0.1863109301186906e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.3332824691372289e5 / (0.3341852231986051e0 - a) / (0.3341852231986051e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    -0.3232914285714286e5 / (0.5e0 - a) / (0.5e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    0.3332824691372289e5 / (0.6658147768013949e0 - a) / (0.6658147768013949e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] =
                    -0.3688889968687672e5 / (0.8136890698813094e0 - a) / (0.8136890698813094e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][7] =
                    0.4624522420172525e5 / (0.9277010131424373e0 - a) / (0.9277010131424373e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][8] = 1.0e0
                    / ((b - a)
                        * (b - 0.7229898685756272e-1)
                        * (b - 0.1863109301186906e0)
                        * (b - 0.3341852231986051e0)
                        * (b - 0.5e0)
                        * (b - 0.6658147768013949e0)
                        * (b - 0.8136890698813094e0)
                        * (b - 0.9277010131424373e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - a)
                        * (points[j - 1] - 0.7229898685756272e-1)
                        * (points[j - 1] - 0.1863109301186906e0)
                        * (points[j - 1] - 0.3341852231986051e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.6658147768013949e0)
                        * (points[j - 1] - 0.8136890698813094e0)
                        * (points[j - 1] - 0.9277010131424373e0)
                        * (points[j - 1] - b);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7229898685756272e-1);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1863109301186906e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3341852231986051e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6658147768013949e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8136890698813094e0);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9277010131424373e0);
                    g[(i - 1) * quad + (j - 1)][8] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][8] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                }
            } else if (i == 1 ){
                b = 1.0e0 + 0.7229898685756272e-1 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                bb = 1.0e0 + 0.1863109301186906e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -0.3e2 - 1.0e0 / b - 1.0e0 / bb;
                gamma2 = 0.3e2 + 1.0e0 / (1.0e0 - b) + 1.0e0 / (1.0e0 - bb);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.884e3 / b / bb;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.884e3 / (1.0e0 - b) / (1.0e0 - bb);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.4624522420172525e5
                    / (0.7229898685756272e-1 - bb)
                    / (0.7229898685756272e-1 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = -0.3688889968687672e5
                    / (0.1863109301186906e0 - bb)
                    / (0.1863109301186906e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.3332824691372289e5 / (0.3341852231986051e0 - bb) / (0.3341852231986051e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    -0.3232914285714286e5 / (0.5e0 - bb) / (0.5e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    0.3332824691372289e5 / (0.6658147768013949e0 - bb) / (0.6658147768013949e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] = -0.3688889968687672e5
                    / (0.8136890698813094e0 - bb)
                    / (0.8136890698813094e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] =
                    0.4624522420172525e5 / (0.9277010131424373e0 - bb) / (0.9277010131424373e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][7] = 1.0e0
                    / ((b - bb)
                        * (b - 0.7229898685756272e-1)
                        * (b - 0.1863109301186906e0)
                        * (b - 0.3341852231986051e0)
                        * (b - 0.5e0)
                        * (b - 0.6658147768013949e0)
                        * (b - 0.8136890698813094e0)
                        * (b - 0.9277010131424373e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][8] = 1.0e0
                    / ((bb - b)
                        * (bb - 0.7229898685756272e-1)
                        * (bb - 0.1863109301186906e0)
                        * (bb - 0.3341852231986051e0)
                        * (bb - 0.5e0)
                        * (bb - 0.6658147768013949e0)
                        * (bb - 0.8136890698813094e0)
                        * (bb - 0.9277010131424373e0)
                        * bb
                        * bb
                        * (bb - 1.0e0)
                        * (bb - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - 0.7229898685756272e-1)
                        * (points[j - 1] - 0.1863109301186906e0)
                        * (points[j - 1] - 0.3341852231986051e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.6658147768013949e0)
                        * (points[j - 1] - 0.8136890698813094e0)
                        * (points[j - 1] - 0.9277010131424373e0)
                        * (points[j - 1] - b)
                        * (points[j - 1] - bb);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7229898685756272e-1);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1863109301186906e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3341852231986051e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6658147768013949e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8136890698813094e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9277010131424373e0);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                    g[(i - 1) * quad + (j - 1)][8] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][8] * eta[j - 1][2] * phi
                            / (points[j - 1] - bb);
                }
            } else {
                aa = -0.1863109301186906e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                a = -0.7229898685756272e-1 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                gamma1 = -1.0e0 / aa - 1.0e0 / a - 0.3e2;
                gamma2 = 1.0e0 / (1.0e0 - aa) + 1.0e0 / (1.0e0 - a) + 0.3e2;
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.884e3 / aa / a;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] = 0.884e3 / (1.0e0 - aa) / (1.0e0 - a);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((aa - 0.7229898685756272e-1)
                        * (aa - 0.1863109301186906e0)
                        * (aa - 0.3341852231986051e0)
                        * (aa - 0.5e0)
                        * (aa - 0.6658147768013949e0)
                        * (aa - 0.8136890698813094e0)
                        * (aa - 0.9277010131424373e0)
                        * (aa - a)
                        * aa
                        * aa
                        * (aa - 1.0e0)
                        * (aa - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 1.0e0
                    / ((a - 0.7229898685756272e-1)
                        * (a - 0.1863109301186906e0)
                        * (a - 0.3341852231986051e0)
                        * (a - 0.5e0)
                        * (a - 0.6658147768013949e0)
                        * (a - 0.8136890698813094e0)
                        * (a - 0.9277010131424373e0)
                        * (a - aa)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = 0.4624522420172525e5
                    / (0.7229898685756272e-1 - a)
                    / (0.7229898685756272e-1 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] = -0.3688889968687672e5
                    / (0.1863109301186906e0 - a)
                    / (0.1863109301186906e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    0.3332824691372289e5 / (0.3341852231986051e0 - a) / (0.3341852231986051e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    -0.3232914285714286e5 / (0.5e0 - a) / (0.5e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] =
                    0.3332824691372289e5 / (0.6658147768013949e0 - a) / (0.6658147768013949e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][7] = -0.3688889968687672e5
                    / (0.8136890698813094e0 - a)
                    / (0.8136890698813094e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][8] =
                    0.4624522420172525e5 / (0.9277010131424373e0 - a) / (0.9277010131424373e0 - aa);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - aa)
                        * (points[j - 1] - a)
                        * (points[j - 1] - 0.7229898685756272e-1)
                        * (points[j - 1] - 0.1863109301186906e0)
                        * (points[j - 1] - 0.3341852231986051e0)
                        * (points[j - 1] - 0.5e0)
                        * (points[j - 1] - 0.6658147768013949e0)
                        * (points[j - 1] - 0.8136890698813094e0)
                        * (points[j - 1] - 0.9277010131424373e0);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - aa);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7229898685756272e-1);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1863109301186906e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.3341852231986051e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.6658147768013949e0);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8136890698813094e0);
                    g[(i - 1) * quad + (j - 1)][8] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][8] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.9277010131424373e0);
                }
            }
        }
    } else if (kcol == 10 ){
        for (size_t i=1;i<=nint;i++) {
            hx = mesh[i + 1 - 1] - mesh[i - 1];
            if (i != 1 && i != nint ){
                a = -0.5929571219129399e-1 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                b = 1.0e0 + 0.5929571219129399e-1 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -1.0e0 / a - 0.3666666666666667e2 - 1.0e0 / b;
                gamma2 = 1.0e0 / (1.0e0 - a) + 0.3666666666666667e2 + 1.0e0 / (1.0e0 - b);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.2799333333333333e4 / a / b;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.2799333333333333e4 / (1.0e0 - a) / (1.0e0 - b);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((a - 0.5929571219129399e-1)
                        * (a - 0.1539696908715823e0)
                        * (a - 0.2792835119457421e0)
                        * (a - 0.4241841678533667e0)
                        * (a - 0.5758158321466333e0)
                        * (a - 0.7207164880542579e0)
                        * (a - 0.8460303091284177e0)
                        * (a - 0.940704287808706e0)
                        * (a - b)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = -0.1785201442945021e6
                    / (0.5929571219129399e-1 - a)
                    / (0.5929571219129399e-1 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] =
                    0.1412225097382635e6 / (0.1539696908715823e0 - a) / (0.1539696908715823e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    -0.1259240998337037e6 / (0.2792835119457421e0 - a) / (0.2792835119457421e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] =
                    0.1197516586185479e6 / (0.4241841678533667e0 - a) / (0.4241841678533667e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    -0.1197516586185479e6 / (0.5758158321466333e0 - a) / (0.5758158321466333e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] =
                    0.1259240998337037e6 / (0.7207164880542579e0 - a) / (0.7207164880542579e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][7] =
                    -0.1412225097382635e6 / (0.8460303091284177e0 - a) / (0.8460303091284177e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][8] =
                    0.1785201442945021e6 / (0.940704287808706e0 - a) / (0.940704287808706e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][9] = 1.0e0
                    / ((b - a)
                        * (b - 0.5929571219129399e-1)
                        * (b - 0.1539696908715823e0)
                        * (b - 0.2792835119457421e0)
                        * (b - 0.4241841678533667e0)
                        * (b - 0.5758158321466333e0)
                        * (b - 0.7207164880542579e0)
                        * (b - 0.8460303091284177e0)
                        * (b - 0.940704287808706e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - a)
                        * (points[j - 1] - 0.5929571219129399e-1)
                        * (points[j - 1] - 0.1539696908715823e0)
                        * (points[j - 1] - 0.2792835119457421e0)
                        * (points[j - 1] - 0.4241841678533667e0)
                        * (points[j - 1] - 0.5758158321466333e0)
                        * (points[j - 1] - 0.7207164880542579e0)
                        * (points[j - 1] - 0.8460303091284177e0)
                        * (points[j - 1] - 0.940704287808706e0)
                        * (points[j - 1] - b);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5929571219129399e-1);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1539696908715823e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2792835119457421e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.4241841678533667e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5758158321466333e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7207164880542579e0);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8460303091284177e0);
                    g[(i - 1) * quad + (j - 1)][8] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][8] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.940704287808706e0);
                    g[(i - 1) * quad + (j - 1)][9] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][9] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                }
            } else if (i == 1 ){
                b = 1.0e0 + 0.5929571219129399e-1 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                bb = 1.0e0 + 0.1539696908715823e0 * (mesh[i + 2 - 1] - mesh[i + 1 - 1]) / hx;
                gamma1 = -0.3666666666666667e2 - 1.0e0 / b - 1.0e0 / bb;
                gamma2 = 0.3666666666666667e2 + 1.0e0 / (1.0e0 - b) + 1.0e0 / (1.0e0 - bb);
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.2799333333333333e4 / b / bb;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.2799333333333333e4 / (1.0e0 - b) / (1.0e0 - bb);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = -0.1785201442945021e6
                    / (0.5929571219129399e-1 - bb)
                    / (0.5929571219129399e-1 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.1412225097382635e6 / (0.1539696908715823e0 - bb) / (0.1539696908715823e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = -0.1259240998337037e6
                    / (0.2792835119457421e0 - bb)
                    / (0.2792835119457421e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.1197516586185479e6 / (0.4241841678533667e0 - bb) / (0.4241841678533667e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] = -0.1197516586185479e6
                    / (0.5758158321466333e0 - bb)
                    / (0.5758158321466333e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    0.1259240998337037e6 / (0.7207164880542579e0 - bb) / (0.7207164880542579e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] = -0.1412225097382635e6
                    / (0.8460303091284177e0 - bb)
                    / (0.8460303091284177e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][7] =
                    0.1785201442945021e6 / (0.940704287808706e0 - bb) / (0.940704287808706e0 - b);
                g[(i - 1) * quad + (kcol + 3 - 1)][8] = 1.0e0
                    / ((b - bb)
                        * (b - 0.5929571219129399e-1)
                        * (b - 0.1539696908715823e0)
                        * (b - 0.2792835119457421e0)
                        * (b - 0.4241841678533667e0)
                        * (b - 0.5758158321466333e0)
                        * (b - 0.7207164880542579e0)
                        * (b - 0.8460303091284177e0)
                        * (b - 0.940704287808706e0)
                        * b
                        * b
                        * (b - 1.0e0)
                        * (b - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][9] = 1.0e0
                    / ((bb - b)
                        * (bb - 0.5929571219129399e-1)
                        * (bb - 0.1539696908715823e0)
                        * (bb - 0.2792835119457421e0)
                        * (bb - 0.4241841678533667e0)
                        * (bb - 0.5758158321466333e0)
                        * (bb - 0.7207164880542579e0)
                        * (bb - 0.8460303091284177e0)
                        * (bb - 0.940704287808706e0)
                        * bb
                        * bb
                        * (bb - 1.0e0)
                        * (bb - 1.0e0));
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - 0.5929571219129399e-1)
                        * (points[j - 1] - 0.1539696908715823e0)
                        * (points[j - 1] - 0.2792835119457421e0)
                        * (points[j - 1] - 0.4241841678533667e0)
                        * (points[j - 1] - 0.5758158321466333e0)
                        * (points[j - 1] - 0.7207164880542579e0)
                        * (points[j - 1] - 0.8460303091284177e0)
                        * (points[j - 1] - 0.940704287808706e0)
                        * (points[j - 1] - b)
                        * (points[j - 1] - bb);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5929571219129399e-1);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1539696908715823e0);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2792835119457421e0);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.4241841678533667e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5758158321466333e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7207164880542579e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8460303091284177e0);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.940704287808706e0);
                    g[(i - 1) * quad + (j - 1)][8] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][8] * eta[j - 1][2] * phi
                            / (points[j - 1] - b);
                    g[(i - 1) * quad + (j - 1)][9] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][9] * eta[j - 1][2] * phi
                            / (points[j - 1] - bb);
                }
            } else {
                aa = -0.1539696908715823e0 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                a = -0.5929571219129399e-1 * (mesh[i - 1] - mesh[i - 1 - 1]) / hx;
                gamma1 = -1.0e0 / aa - 1.0e0 / a - 0.3666666666666667e2;
                gamma2 = 1.0e0 / (1.0e0 - aa) + 1.0e0 / (1.0e0 - a) + 0.3666666666666667e2;
                h[(i - 1) * quad + (kcol + 3 - 1)][0] = 0.2799333333333333e4 / aa / a;
                h[(i - 1) * quad + (kcol + 3 - 1)][1] =
                    0.2799333333333333e4 / (1.0e0 - aa) / (1.0e0 - a);
                hd[(i - 1) * quad + (kcol + 3 - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0] * hx;
                hd[(i - 1) * quad + (kcol + 3 - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1] * hx;
                g[(i - 1) * quad + (kcol + 3 - 1)][0] = 1.0e0
                    / ((aa - 0.5929571219129399e-1)
                        * (aa - 0.1539696908715823e0)
                        * (aa - 0.2792835119457421e0)
                        * (aa - 0.4241841678533667e0)
                        * (aa - 0.5758158321466333e0)
                        * (aa - 0.7207164880542579e0)
                        * (aa - 0.8460303091284177e0)
                        * (aa - 0.940704287808706e0)
                        * (aa - a)
                        * aa
                        * aa
                        * (aa - 1.0e0)
                        * (aa - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][1] = 1.0e0
                    / ((a - 0.5929571219129399e-1)
                        * (a - 0.1539696908715823e0)
                        * (a - 0.2792835119457421e0)
                        * (a - 0.4241841678533667e0)
                        * (a - 0.5758158321466333e0)
                        * (a - 0.7207164880542579e0)
                        * (a - 0.8460303091284177e0)
                        * (a - 0.940704287808706e0)
                        * (a - aa)
                        * a
                        * a
                        * (a - 1.0e0)
                        * (a - 1.0e0));
                g[(i - 1) * quad + (kcol + 3 - 1)][2] = -0.1785201442945021e6
                    / (0.5929571219129399e-1 - a)
                    / (0.5929571219129399e-1 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][3] =
                    0.1412225097382635e6 / (0.1539696908715823e0 - a) / (0.1539696908715823e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][4] = -0.1259240998337037e6
                    / (0.2792835119457421e0 - a)
                    / (0.2792835119457421e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][5] =
                    0.1197516586185479e6 / (0.4241841678533667e0 - a) / (0.4241841678533667e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][6] = -0.1197516586185479e6
                    / (0.5758158321466333e0 - a)
                    / (0.5758158321466333e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][7] =
                    0.1259240998337037e6 / (0.7207164880542579e0 - a) / (0.7207164880542579e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][8] = -0.1412225097382635e6
                    / (0.8460303091284177e0 - a)
                    / (0.8460303091284177e0 - aa);
                g[(i - 1) * quad + (kcol + 3 - 1)][9] =
                    0.1785201442945021e6 / (0.940704287808706e0 - a) / (0.940704287808706e0 - aa);
                for (size_t j=1;j<=kcol + 3;j++) {
                    phi = (points[j - 1] - aa)
                        * (points[j - 1] - a)
                        * (points[j - 1] - 0.5929571219129399e-1)
                        * (points[j - 1] - 0.1539696908715823e0)
                        * (points[j - 1] - 0.2792835119457421e0)
                        * (points[j - 1] - 0.4241841678533667e0)
                        * (points[j - 1] - 0.5758158321466333e0)
                        * (points[j - 1] - 0.7207164880542579e0)
                        * (points[j - 1] - 0.8460303091284177e0)
                        * (points[j - 1] - 0.940704287808706e0);
                    h[(i - 1) * quad + (j - 1)][0] = h[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * (1.0e0 - points[j - 1] * gamma1)
                        * eta[j - 1][0]
                        * phi;
                    h[(i - 1) * quad + (j - 1)][1] = h[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (1.0e0 + gamma2 - points[j - 1] * gamma2)
                        * eta[j - 1][1]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][0] = hd[(i - 1) * quad + (kcol + 3 - 1)][0]
                        * points[j - 1]
                        * eta[j - 1][0]
                        * phi;
                    hd[(i - 1) * quad + (j - 1)][1] = hd[(i - 1) * quad + (kcol + 3 - 1)][1]
                        * (points[j - 1] - 1.0e0)
                        * eta[j - 1][1]
                        * phi;
                    g[(i - 1) * quad + (j - 1)][0] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][0] * eta[j - 1][2] * phi
                            / (points[j - 1] - aa);
                    g[(i - 1) * quad + (j - 1)][1] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][1] * eta[j - 1][2] * phi
                            / (points[j - 1] - a);
                    g[(i - 1) * quad + (j - 1)][2] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][2] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5929571219129399e-1);
                    g[(i - 1) * quad + (j - 1)][3] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][3] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.1539696908715823e0);
                    g[(i - 1) * quad + (j - 1)][4] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][4] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.2792835119457421e0);
                    g[(i - 1) * quad + (j - 1)][5] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][5] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.4241841678533667e0);
                    g[(i - 1) * quad + (j - 1)][6] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][6] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.5758158321466333e0);
                    g[(i - 1) * quad + (j - 1)][7] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][7] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.7207164880542579e0);
                    g[(i - 1) * quad + (j - 1)][8] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][8] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.8460303091284177e0);
                    g[(i - 1) * quad + (j - 1)][9] =
                        g[(i - 1) * quad + (kcol + 3 - 1)][9] * eta[j - 1][2] * phi
                            / (points[j - 1] - 0.940704287808706e0);
                }
            }
        }
    }
}
