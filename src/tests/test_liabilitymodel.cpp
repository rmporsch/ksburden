#include <gtest/gtest.h>
#include <liability_model.h>

class test_liabilitymodel: public testing::Test
{
protected:
  LiabilityModel liabmodel;
  arma::Mat<int> test_mat;
  test_liabilitymodel() {
	  test_mat.set_size(100, 10);
	  test_mat.zeros();
	  test_mat(5, 0) = 1;
	  test_mat(10, 1) = 1;
	  test_mat(15, 2) = 1;
	  test_mat(20, 3) = 1;
	  test_mat(25, 4) = 1;
	  test_mat(30, 5) = 1;
	  test_mat(35, 6) = 1;
	  test_mat(50, 7) = 1;
	  test_mat(60, 8) = 1;
	  test_mat(70, 9) = 1;
  }
};

TEST_F(test_liabilitymodel, normal_random) {
  arma::vec normalr = liabmodel.normal_random(10000, 0, 1);
  double mean = arma::mean(normalr);
  double var = arma::var(normalr);
  EXPECT_NEAR(0, mean, 0.05);
  EXPECT_NEAR(1, var, 0.05);
}

TEST_F(test_liabilitymodel, uniform_random) {
  arma::vec unifr = liabmodel.uniform_random(10);
  int num_n = unifr.size();
  EXPECT_EQ(10, num_n);
}

TEST_F(test_liabilitymodel, generate_causal_variants) {
  liabmodel.num_variants = 10;
  liabmodel.size_cluster = 1;
  arma::Col<int> causal_var = liabmodel.generate_causal_variants(true);
  EXPECT_EQ(10, causal_var.size());
  EXPECT_EQ(1, sum(causal_var));
}


TEST_F(test_liabilitymodel, effect_generation) {
  liabmodel.num_variants = 10;
  liabmodel.size_cluster = 1;
  liabmodel.wished_effect = 0.1;
  liabmodel.genotype_matrix = test_mat;
  arma::Col<int> causal_var = liabmodel.generate_causal_variants(true);
  arma::vec effect = liabmodel.effect_generation(causal_var);
  EXPECT_NEAR(0.1, pow(sum(effect),2), 0.001);
  EXPECT_EQ(10, effect.size());
}

TEST_F(test_liabilitymodel, standardize_matrix) {
  liabmodel.num_variants = 10;
  liabmodel.size_cluster = 1;
  liabmodel.genotype_matrix = test_mat;
  liabmodel.standardize_matrix();

  EXPECT_NEAR(0, mean(liabmodel.genotype_matrix_standarized.col(0)), 0.0001);
  EXPECT_NEAR(1, var(liabmodel.genotype_matrix_standarized.col(0)), 0.0001);
  EXPECT_NEAR(0, mean(liabmodel.genotype_matrix_standarized.col(1)), 0.0001);
  EXPECT_NEAR(1, var(liabmodel.genotype_matrix_standarized.col(1)), 0.0001);
}

TEST_F(test_liabilitymodel, simulate_data) {
  liabmodel.num_variants = 10;
  liabmodel.size_cluster = 1;
  liabmodel.wished_effect = 0.1;
  liabmodel.num_cases = 10;
  liabmodel.life_time_risk = 0.1;
  liabmodel.num_controls = 10;
  liabmodel.genotype_matrix = test_mat;
  liabmodel.standardize_matrix();
  arma::Col<int> causal_var = liabmodel.generate_causal_variants(true);
  arma::uvec sim_id = liabmodel.simulate_data(causal_var);
  double mean_risk = mean(liabmodel.liability_dist);
  double var_risk = var(liabmodel.liability_dist);
  EXPECT_EQ(20, sim_id.size());
  EXPECT_NEAR(0, mean_risk, 0.05);
  EXPECT_NEAR(1, var_risk, 0.05);
}
