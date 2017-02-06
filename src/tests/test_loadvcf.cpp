#include <gtest/gtest.h>
#include <load_vcf.h>

class test_loadvcf : public testing::Test
{
protected:
  std::string vcf_file;
  std::string variant_file;
  LoadVCF instance;
  std::string gene_name;
  test_loadvcf() {
    vcf_file = "./src/tests/test_data/small_AMPD2.recode.vcf.gz";
    variant_file = "./src/tests/test_data/small_AMPD2.txt";
    gene_name = "AMPD2";
  };
};


TEST_F(test_loadvcf, variant_location) {
  std::vector<std::vector<std::string>> variant_loc;
  variant_loc = instance.variant_location(variant_file, '\t');
  int size_variant_loc = variant_loc.size();
  EXPECT_EQ(3, size_variant_loc);
}

TEST_F(test_loadvcf, num_genes) {
  std::vector<std::vector<std::string>> variant_loc;
  variant_loc = instance.variant_location(variant_file, '\t');
  std::vector<std::string> genes;
  genes = instance.get_all_genes(variant_loc);
  int num_genes = genes.size();
  EXPECT_EQ(1, num_genes);
  EXPECT_EQ(gene_name, genes[0]);
}

TEST_F(test_loadvcf, get_gene_loc) {
  LoadVCF gene_loc(vcf_file, variant_file);
  std::vector<std::vector<std::string>> current_gene_loc;
  current_gene_loc = gene_loc.get_gene_loc(gene_name);
  std::string chrom = "chr1";
  std::string pos = "110170502";
  EXPECT_EQ(chrom, current_gene_loc[0][0]);
  EXPECT_EQ(pos, current_gene_loc[0][1]);
  EXPECT_EQ(3, current_gene_loc.size());
}

TEST_F(test_loadvcf, get_matrix) {
  LoadVCF gene_mat(vcf_file, variant_file);
  gene_mat.load_gene(gene_name);
  EXPECT_EQ(915, gene_mat.genotype_matrix.n_rows);
  EXPECT_EQ(3, gene_mat.genotype_matrix.n_cols);
  EXPECT_EQ(4, arma::accu(gene_mat.genotype_matrix));
}
