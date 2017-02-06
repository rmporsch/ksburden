#include <gtest/gtest.h>
#include <load_plink.h>
#include <armadillo>

class test_loadplink : public testing::Test
{
protected:
  std::string plink_file;
  std::string variant_file;
  LoadPlink instance;
  std::string gene_name;
  test_loadplink() { plink_file = "./src/tests/test_data/dummy_gene";
    variant_file = "./src/tests/test_data/dummy_gene.variant.list";
    gene_name = "BMS1P17(0)|BMS1P18(0)";
  };
};

TEST_F(test_loadplink, variant_location) {
  std::vector<std::vector<std::string>> variant_loc;
  variant_loc = instance.variant_location(variant_file, '\t');
  int size_variant_loc = variant_loc.size();
  EXPECT_EQ(9, size_variant_loc);
}

TEST_F(test_loadplink, num_genes) {
  std::vector<std::vector<std::string>> variant_loc;
  variant_loc = instance.variant_location(variant_file, '\t');
  std::vector<std::string> genes;
  genes = instance.get_all_genes(variant_loc, 1);
  int num_genes = genes.size();
  EXPECT_EQ(1, num_genes);
  EXPECT_EQ(gene_name, genes[0]);
}

TEST_F(test_loadplink, openPlinkBinaryFile) {
  std::string bedfile = plink_file + ".bed";
  std::ifstream bfile;
  bool opened = instance.openPlinkBinaryFile(bedfile, bfile);
  EXPECT_TRUE(opened);
}

TEST_F(test_loadplink, get_gene_loc) {
  std::string bedfile = plink_file + ".bed";
  std::string bimfile = plink_file + ".bim";
  std::string famfile = plink_file + ".fam";
  std::vector<std::vector<std::string>> current_gene_loc;
  LoadPlink gene_loc(famfile, bimfile, bedfile, variant_file);
  current_gene_loc = gene_loc.get_gene_loc(gene_name); 
  EXPECT_EQ("rs575268414", current_gene_loc[1][0]);
}

TEST_F(test_loadplink, getcol_skip) {
  std::string bedfile = plink_file + ".bed";
  std::string bimfile = plink_file + ".bim";
  std::string famfile = plink_file + ".fam";
  std::vector<std::vector<std::string>> current_gene_loc;
  LoadPlink gene_loc(famfile, bimfile, bedfile, variant_file);
  current_gene_loc = gene_loc.get_gene_loc(gene_name); 
  arma::Col<int> variants_select = gene_loc.get_col_skip(current_gene_loc);
  EXPECT_EQ(9, arma::accu(variants_select));
}

TEST_F(test_loadplink, get_genotype_matrix) {
  std::string bedfile = plink_file + ".bed";
  std::string bimfile = plink_file + ".bim";
  std::string famfile = plink_file + ".fam";
  std::vector<std::vector<std::string>> current_gene_loc;
  LoadPlink gene_loc(famfile, bimfile, bedfile, variant_file);
  current_gene_loc = gene_loc.get_gene_loc(gene_name); 
  arma::Col<int> variants_select = gene_loc.get_col_skip(current_gene_loc);
  arma::Col<int> row_skip(gene_loc.N, arma::fill::ones);
  arma::mat genotypeMatrix = gene_loc.get_genotype_matrix(bedfile, variants_select, 
      row_skip);
  EXPECT_EQ(9, genotypeMatrix.n_cols);
  EXPECT_EQ(gene_loc.N, genotypeMatrix.n_rows);
  std::cout << arma::accu(genotypeMatrix) << std::endl;
}
