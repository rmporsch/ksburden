#include <gtest/gtest.h>
#include <models.h>

class models_test : public testing::Test {
protected:
  int num_subjects = 10;
  int num_variants = 2;
  arma::mat somemat;
  arma::uvec test_subjects;
  arma::Col<int> phenotype;
  arma::Mat<int> test_matrix;
  models instance;

models_test() {
  test_matrix.set_size(num_subjects, num_variants);
  test_matrix.zeros();
  test_matrix(0, 0) = 1;
  test_matrix(7, 1) = 1;

  test_subjects = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  phenotype.set_size(10);
  phenotype.ones();
  arma::Col<int> cases(5);
  cases.zeros();
  cases = cases - 1;
  phenotype(arma::span(5, 9)) = cases;

}
};

TEST_F(models_test, ecdf) {

  arma::vec test_out = instance.edf(test_matrix, test_subjects);
  EXPECT_EQ(1.5, sum(test_out));
  EXPECT_EQ(2, test_out.size());
  EXPECT_EQ(0, sum(phenotype));
}

TEST_F(models_test, ksburden) {

  EXPECT_EQ(1.0, instance.ksburden(test_matrix, phenotype));
}

TEST_F(models_test, fisher) {
  arma::vec pvalues = {0.01, 0.05};
  EXPECT_FLOAT_EQ(0.004300451, instance.fisher(pvalues));
}

TEST_F(models_test, large_fisher) {
  arma::mat lpvalues(2,2);
  lpvalues = {{0.01, 0.05}, {0.02, 0.01}};
  arma::vec out = instance.large_fisher(lpvalues);
  EXPECT_FLOAT_EQ(0.004300451, out(0));
  EXPECT_FLOAT_EQ(0.001903439, out(1));
}

TEST_F(models_test, cmc) {
  EXPECT_EQ(10, test_matrix.n_rows);
  EXPECT_EQ(2, test_matrix.n_cols);
  EXPECT_EQ(0, sum(phenotype));
  EXPECT_EQ(0.0, instance.cmc(test_matrix, phenotype));
}

TEST_F(models_test, burden) {
  EXPECT_EQ(0.0, instance.burden(test_matrix, phenotype));
}

TEST_F(models_test, permutation) {
  int num_tests = 3;
  EXPECT_GE(1, instance.permutation(instance.model_array[0], 1, test_matrix,
                                     phenotype, 1));
  EXPECT_GE(1, instance.permutation(instance.model_array[1], 1, test_matrix,
                                     phenotype, 1));
  EXPECT_GE(1, instance.permutation(instance.model_array[2], 1, test_matrix,
                                      phenotype, 1));
}

TEST_F(models_test, actualdata) {
  std::string ret_genotype_file = "./src/tests/test_data/saved_KIAA0556.csv";
  arma::Mat<int> ret_mat;

  ret_mat.load(ret_genotype_file);

  int num_cases = 458;
  int num_controls = 915 - num_cases;
  arma::Col<int> phenotype(915);
  phenotype.ones();
  arma::Col<int> case_vec(num_cases);
  case_vec.zeros();
  case_vec = case_vec - 1;
  phenotype.subvec(num_controls, (915 -1)) = case_vec;

  double test_stat_ks = instance.ksburden(ret_mat, phenotype);
  double test_stat_burden = instance.burden(ret_mat, phenotype);
  double test_stat_cmc = instance.cmc(ret_mat, phenotype);
  EXPECT_NEAR(0.2441, test_stat_ks, 0.0001);
  EXPECT_FLOAT_EQ(196, test_stat_burden);
  EXPECT_NEAR(0.000557, test_stat_cmc, 0.0001);

  EXPECT_NEAR(0.133, instance.permutation(instance.model_array[0], 1000,
                                          ret_mat, phenotype, 1001),
              0.03);
  EXPECT_NEAR(0.118, instance.permutation(instance.model_array[1], 1000,
                                          ret_mat, phenotype, 1001),
              0.03);
  EXPECT_NEAR(0.184, instance.permutation(instance.model_array[2], 1000,
                                          ret_mat, phenotype, 1001),
              0.03);
}
