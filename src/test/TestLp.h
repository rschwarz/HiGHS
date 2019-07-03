
#ifndef TEST_TESTLP_H_
#define TEST_TESTLP_H_

#include "lp_data/HighsLp.h"

HighsLp generateTestLpMini();
double getOptimalObjectiveMini();
std::vector<double> getOptimalSolutionMini();

#endif