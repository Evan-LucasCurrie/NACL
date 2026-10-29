#ifndef COLNEWSC_TEST_PROBLEMS_H
#define COLNEWSC_TEST_PROBLEMS_H
#include "colnewsc_reference.h"
typedef enum {
    CN_SWIRLING_FLOW,
    CN_TP20,
    CN_TP21,
    CN_MEASLES,
    CN_TMLPDE,
    CN_SHOCK_WAVE,
    CN_TP4,
    CN_THVE16_EX1,
    CN_SIWA22_EX1,
    CN_SNK14_EX1,
    CN_SNK14_EX2,
    CN_KHWA14_EX2,
    CN_KHWA14_EX4,
    CN_TEST_CASE_COUNT
} colnewsc_test_case;
const char *colnewsc_test_name(colnewsc_test_case);
colnewsc_problem colnewsc_test_problem(colnewsc_test_case);
colnewsc_options colnewsc_test_options(colnewsc_test_case);
/* Returns NAN when no analytic first component is available. */
double colnewsc_test_exact(colnewsc_test_case, double);
colnewsc_reference_problem colnewsc_reference_test_problem(colnewsc_test_case);
#endif
