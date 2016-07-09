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
  test_matrix(2, 1) = 1;

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
  EXPECT_GE(1, instance.permutation(instance.model_array[1], 1, test_matrix,
                                     phenotype, 1));
  EXPECT_GE(1, instance.permutation(instance.model_array[1], 1, test_matrix,
                                     phenotype, 1));
  EXPECT_GE(1, instance.permutation(instance.model_array[2], 1, test_matrix,
                                      phenotype, 1));
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
