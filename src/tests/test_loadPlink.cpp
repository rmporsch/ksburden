#include <gtest/gtest.h>
#include <load_plink.h>

class test_loadplink : public testing::Test
{
protected:
  std::string plink_file;
  std::string variant_file;
  LoadPlink instance;
  std::string gene_name;
  test_loadplink() {
    plink_file = "./src/tests/test_data/dummy_gene";
    variant_file = "./src/tests/test_data/dummy_gene.variant.list";
    gene_name = "BMS1P17(0)|BMS1P18(0)";
  };
};

TEST_F(test_loadplink, variant_location) {
  std::vector<std::vector<std::string>> variant_loc;
  variant_loc = instance.variant_location(variant_file, '\t');
  int size_variant_loc = variant_loc.size();
  EXPECT_EQ(3, size_variant_loc);
}

TEST_F(test_loadplink, num_genes) {
  std::vector<std::vector<std::string>> variant_loc;
  variant_loc = instance.variant_location(variant_file, '\t');
  std::vector<std::string> genes;
  genes = instance.get_all_genes(variant_loc);
  int num_genes = genes.size();
  EXPECT_EQ(1, num_genes);
  EXPECT_EQ(gene_name, genes[0]);
}
