/* Integrated Lagrange basis on Gauss-Legendre nodes, orders 1 through 7. */
#include "internal.h"
void cn_basis_weights(const cn_basis *b, double t, double *out) {
    for (size_t j = 0; j < b->k; j++) {
        double v = 0;
        for (size_t l = b->k + 1; l-- > 0;)
            v = v * t + b->integrated[j][l];
        out[j] = v;
    }
}
void cn_basis_new(cn_basis *b, size_t k) {
    static const double roots[7][7] = {
        {0},
        {-0.5773502691896258, 0.5773502691896258},
        {-0.7745966692414834, 0, 0.7745966692414834},
        {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526},
        {-0.906179845938664, -0.5384693101056831, 0, 0.5384693101056831, 0.906179845938664},
        {-0.932469514203152, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969,
         0.6612093864662645, 0.932469514203152},
        {-0.9491079912342759, -0.7415311855993945, -0.4058451513773972, 0, 0.4058451513773972,
         0.7415311855993945, 0.9491079912342759}};
    memset(b, 0, sizeof *b);
    b->k = k;
    for (size_t j = 0; j < k; j++)
        b->nodes[j] = 0.5 * (1 + roots[k - 1][j]);
    for (size_t j = 0; j < k; j++) {
        double p[8] = {1};
        size_t count = 1;
        for (size_t i = 0; i < k; i++)
            if (i != j) {
                double q[8] = {0}, den = b->nodes[j] - b->nodes[i];
                for (size_t l = 0; l < count; l++) {
                    q[l] -= p[l] * b->nodes[i] / den;
                    q[l + 1] += p[l] / den;
                }
                count++;
                memcpy(p, q, sizeof p);
            }
        for (size_t l = 0; l < count; l++)
            b->integrated[j][l + 1] = p[l] / (double)(l + 1);
    }
    cn_basis_weights(b, 1, b->b);
    for (size_t j = 0; j < k; j++)
        cn_basis_weights(b, b->nodes[j], b->a + j * k);
}
