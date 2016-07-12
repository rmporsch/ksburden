#include <gtest/gtest.h>
#include <data_analysis.h>
#include <load_vcf.h>

class test_analysis : public testing::Test {
protected:
  std::string ped_file;
  std::string vcf_file;
  std::string variant_file;
  std::vector<std::string> name_subjects;
  test_analysis() {
    ped_file = "./src/tests/test_data/pedigree.ped";
    vcf_file = "./src/tests/test_data/small_AMPD2.recode.vcf.gz";
    variant_file = "./src/tests/test_data/small_AMPD2.txt";
  };
};

TEST_F(test_analysis, subjects) {
  LoadVCF gene_loc(vcf_file, variant_file);
  analysis ana;

  name_subjects = ana.sampleNames(gene_loc.header);
  int num_subjects =  name_subjects.size();
  EXPECT_EQ(915, num_subjects);
}

TEST_F(test_analysis, pedigree) {
  LoadVCF gene_loc(vcf_file, variant_file);
  analysis ana;

  name_subjects = ana.sampleNames(gene_loc.header);
  ana.pedigree(ped_file, '\t', name_subjects);
  EXPECT_EQ(915, ana.phenotype.size());
  arma::uvec cases = arma::find(ana.phenotype == 2);
  int num_cases = cases.size();
  double per_cases = (double)num_cases / (double)ana.phenotype.size();
  EXPECT_NEAR(0.5, per_cases, 0.1);
}
