#ifndef simulation_flags_H
#define simulation_flags_H

#include <gflags/gflags.h>

DECLARE_string(genotypes);
DECLARE_string(variant);
DECLARE_string(gene);
DECLARE_string(tests);
DECLARE_string(simmat);


DECLARE_int32(threads);
DECLARE_int32(verbose);
DECLARE_int32(numcases);
DECLARE_int32(numcontrols);
DECLARE_int32(powerIter);
DECLARE_int32(iter);

DECLARE_bool(trackPerformance);
DECLARE_bool(storeSim);
DECLARE_bool(pvalues);
DECLARE_string(out);
DECLARE_string(path);

DECLARE_double(maxEffect);
DECLARE_double(minEffect);
DECLARE_double(effectSteps);
DECLARE_double(percentageSteps);
DECLARE_double(lifetimerisk);
DECLARE_double(maxPercentage);
DECLARE_int32(causalVar);
DECLARE_double(probSteps);
DECLARE_double(maxProb);
void sanity_check(int argc, char* argv[]);

#endif
