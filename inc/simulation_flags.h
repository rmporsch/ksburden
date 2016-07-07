#ifndef simulation_flags_H
#define simulation_flags_H

#include <gflags/gflags.h>

DECLARE_string(genotypes);
DECLARE_string(variant);
DECLARE_string(gene);
DECLARE_int32(threads);
DECLARE_int32(verbose);
DECLARE_int32(subjects);
DECLARE_int32(powerIter);
DECLARE_int32(iter);
DECLARE_bool(trackPerformance);
DECLARE_bool(storeSim);
DECLARE_bool(pvalues);
DECLARE_string(output);
DECLARE_string(outPath);
DECLARE_bool(ks);
DECLARE_bool(burden);
DECLARE_bool(CMC);
DECLARE_double(maxEffect);
DECLARE_double(effectSteps);
DECLARE_double(percentageSteps);
DECLARE_double(lifetimerisk);
DECLARE_double(maxPercentage);
DECLARE_int32(causalVar);
DECLARE_double(probSteps);
DECLARE_double(maxProb);
void sanity_check(int argc, char* argv[]);

#endif
