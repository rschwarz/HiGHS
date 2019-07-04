#ifndef PRESOVLE_FINDFEASIBILITY_H_
#define PRESOVLE_FINDFEASIBILITY_H_

#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"

double getQuadraticObjective(const std::vector<double> cost,
                             const std::vector<double>& x,
                             std::vector<double>& r, const double mu,
                             const std::vector<double> lambda);

enum class MinimizationType {
  kComponentWise,
  kComponentWisePenalty,
  kComponentWiseAdmm,
  kComponentWiseBreakpoints,
  kExact,
  kExactPenalty,
  kExactAdmm
};

HighsStatus runFeasibility(const HighsLp& lp, HighsSolution& solution,
                           const MinimizationType type, const double initial_weight);

#endif