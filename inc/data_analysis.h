#ifndef data_analysis_H
#define data_analysis_H

#include <armadillo>
#include <load_vcf.h>
#include <vector>
#include <string>

class analysis : public LoadVCF {

public:
  arma::uvec id_genotype_include;
  arma::Col<int> phenotype;
  std::vector <std::string> exclude;
  std::vector<std::vector<std::string>> ped;
  int num_subjects = ped.size();
  std::vector<std::string> name_subjects;

  void pedigree(std::string ped_file, char sep, std::vector<std::string> vcfID);
  std::vector<std::string> sampleNames(VcfHeader& header);

  analysis(std::string ped_file, std::string vcfFile, std::string variant_file)
      : LoadVCF(vcf_file, variant_file) {
    name_subjects = sampleNames(header);
    pedigree(ped_file, '\t', name_subjects);
  };
  analysis(){};
};
#endif

