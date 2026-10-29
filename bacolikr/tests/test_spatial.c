/* Direct inclusion exposes private arithmetic kernels to fixture validation. */
#include "spatial.c"
#include "test_data/spatial_fixture.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
static uint64_t bits(double x) {
    uint64_t u;
    memcpy(&u, &x, sizeof u);
    return u;
}
int main(void) {
    FILE *f = bac_spatial_fixture();
    assert(f);
    char line[512];
    double mesh[] = {-1.3, -.91, -.1, .11, 1.2},
           points[] = {-1.3, -1.2, -.91, -.3, -.1, .11, .7, 1.2};
    size_t seen = 0;
    while (fgets(line, sizeof line, f)) {
        char record[512];
        memcpy(record, line, sizeof line);
        char *token = strtok(line, " \t\r\n"), kind[32];
        assert(token);
        snprintf(kind, sizeof kind, "%s", token);
        size_t i[8] = {0}, count = 0;
        char *fields[10];
        while ((token = strtok(NULL, " \t\r\n")))
            fields[count++] = token;
        assert(count > 0);
        for (size_t j = 0; j < count - 1; j++)
            i[j] = (size_t)strtoul(fields[j], NULL, 10);
        uint64_t expected = strtoull(fields[count - 1], NULL, 16);
        double actual = NAN;
        if (!strcmp(kind, "gauss")) {
            double d[128], w[128];
            assert(!bac_gauss(i[0], i[1] == 4, d, w));
            actual = i[3] ? w[i[2]] : d[i[2]];
        } else if (!strcmp(kind, "basis")) {
            BacGrid g;
            assert(!bac_grid_new(&g, mesh, 4, i[0]));
            BacBasis b;
            assert(!bac_grid_basis(&g, points[i[1]], &b));
            double v[13][13];
            source_basis(&g, points[i[1]], b.start, i[2] + 1, v);
            actual = v[i[2]][i[3]];
            bac_grid_free(&g);
        } else if (!strcmp(kind, "low")) {
            double d[13], h[2], hd[2], gg[10];
            assert(!bac_gauss(i[0] + 2, 0, d, NULL));
            low_coeff(d[i[1]], i[0] - 1, h, hd, gg);
            actual = i[2] == 0 ? h[i[3]] : i[2] == 1 ? hd[i[3]] : gg[i[3]];
        } else if (!strcmp(kind, "sci")) {
            size_t q = i[0] + 3;
            double d[13], h[52][2] = {{0}}, hd[52][2] = {{0}}, gg[52][10] = {{0}};
            assert(!bac_gauss(q, 0, d, NULL));
            sci_coeff(i[0], mesh, 4, d, h, hd, gg);
            size_t j = i[1] * q + i[2];
            actual = i[3] == 0 ? h[j][i[4]] : i[3] == 1 ? hd[j][i[4]] : gg[j][i[4]];
        } else if (!strcmp(kind, "remesh")) {
            BacGrid g;
            assert(!bac_grid_new(&g, mesh, 4, i[0]));
            double ratios[] = {.04, .31, 1.13, 2.91}, ind[] = {.371, .233, .151, .847};
            BacSpatialError e = {0};
            e.intervals = 4;
            e.ratio = ratios[i[2]];
            e.indicators = ind;
            e.estimator = i[1] ? BAC_SCI : BAC_LOI;
            double *newmesh;
            size_t ni;
            assert(!bac_grid_refined_mesh(&g, &e, 100, 0, 0, &newmesh, &ni));
            assert(ni == i[3]);
            actual = newmesh[i[4]];
            free(newmesh);
            bac_grid_free(&g);
        } else
            assert(0);
        if (bits(actual) != expected) {
            fprintf(stderr,
                    "spatial mismatch: %sactual=%.17g bits=%016" PRIx64 " expected=%016" PRIx64
                    "\n",
                    record, actual, bits(actual), expected);
            return 1;
        }
        seen++;
    }
    fclose(f);
    assert(seen == 6404);
    printf("spatial: %zu exact binary64 fields (Gauss, DBSPVD, LOI, SCI, remesh)\n", seen);
    return 0;
}
