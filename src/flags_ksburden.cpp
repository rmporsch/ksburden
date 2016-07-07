#include <gflags/gflags.h>
#include <iostream>
#include <stdexcept>
#include <flags_ksburden.h>

DEFINE_string(vcf, "/path/to/bcfFile", 
		"The genotypes. This can be either a vcf file or a non-standardized genotype matrix If its a vcf file then one has also to provide the variant file");
DEFINE_string(variant, "/path/to/varinants", 
		"This is either a variant file containing the position of each variant in the first column and the genes in the second. Otherwise this is a standardized genotype matrix");
DEFINE_string(ped, "/path/to/varinants", 
		"This is either a variant file containing the position of each variant in the first column and the genes in the second. Otherwise this is a standardized genotype matrix");

DEFINE_string(out, "/path/to/outputfile", 
		"path of output file");

DEFINE_int32(threads, 1, "number of threads");
DEFINE_int32(iter, 1000, "number of iterations for tests");
DEFINE_int32(verbose, 1, "verbose level");

DEFINE_bool(ks, true, "perform KS test");
DEFINE_bool(burden, true, "perform burden test");
DEFINE_bool(CMC, true, "perform CMC test");

void sanity_check(int argc, char *argv[]) {

  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (::google::GetCommandLineFlagInfoOrDie("vcf").is_default) {
    throw std::runtime_error("a file path for genotypes needs to be specficed");
  }
  if (::google::GetCommandLineFlagInfoOrDie("variant").is_default) {
    throw std::runtime_error("a file path for variants needs to be specficed");
  }
  if (::google::GetCommandLineFlagInfoOrDie("ped").is_default) {
    throw std::runtime_error("a file path for the pedigree needs to be specficed");
  }
}

