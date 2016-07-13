#include <simulation_flags.h>
#include <gflags/gflags.h>
#include <iostream>
#include <stdexcept>

//required
DEFINE_string(genotypes, "/path/to/genotypes", 
		"the genotypes in vcf formate.");
DEFINE_string(variant, "/path/to/varinants", 
		"This is either a variant file containing the position of each variant in the first column and the genes in the second. Otherwise this is a standardized genotype matrix");
DEFINE_string(gene, "somegene", "gene name");
DEFINE_string(tests, "tests", "a string of tests seperated by comma");

//optional
DEFINE_int32(threads, 1, "number of threads");
DEFINE_int32(powerIter, 100, "number of iterations to compute the empirical power");
DEFINE_int32(verbose, 1, "verbose level");
DEFINE_int32(iter, 1000, "number of iterations for tests");
DEFINE_int32(numcases, 500, "number of iterations for tests");
DEFINE_int32(numcontrols, 500, "number of iterations for tests");

DEFINE_bool(pvalues, true, "save p-values for each gene");
DEFINE_bool(storeSim, false, "save simulation IDs for each iteration");
DEFINE_string(out, "power.txt", "file postfix name of the outputfile");

DEFINE_double(maxEffect, 0.0100, "maximal effect size");
DEFINE_double(minEffect, 0.0010, "effect size steps");
DEFINE_double(effectSteps, 0.0010, "effect size steps");

DEFINE_double(percentageSteps, 0.00, "percentage steps in which the program increases the causal cluster");
DEFINE_double(maxPercentage, 1.00, "maximal percentage steps");
DEFINE_double(probSteps, 0.00, "probability steps in which the program increases");
DEFINE_double(maxProb, 1.00, "maximal probability");
DEFINE_double(lifetimerisk, 0.10, "life time risk in percentage");
DEFINE_int32(causalVar, 0, "Number of causal mutations, currently mutual exclusive with 'precentageSteps'");
//
// general stuff
void sanity_check(int argc, char *argv[]) {

  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (::google::GetCommandLineFlagInfoOrDie("genotypes").is_default) {
    throw std::runtime_error("a file path for genotypes needs to be specficed");
  }
  if (::google::GetCommandLineFlagInfoOrDie("variant").is_default) {
    throw std::runtime_error("a file path for variants needs to be specficed");
  }
  if (::google::GetCommandLineFlagInfoOrDie("percentageSteps").is_default &
      ::google::GetCommandLineFlagInfoOrDie("causalVar").is_default &
      ::google::GetCommandLineFlagInfoOrDie("probSteps").is_default) {
    throw std::runtime_error(
        "at least one causal model needs to be specificed");
  }
  if (!(::google::GetCommandLineFlagInfoOrDie("percentageSteps").is_default) &
      !(::google::GetCommandLineFlagInfoOrDie("causalVar").is_default) &
      !(::google::GetCommandLineFlagInfoOrDie("probSteps").is_default)) {
    throw std::runtime_error("only select one causal model");
  }
}

