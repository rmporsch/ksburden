#include <gtest/gtest.h>
#include <liability_model.h>
#include <models.h>

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
  arma::vec normalr = liabmodel.normal_random(10000, 0, sqrt(1));
  double mean = arma::mean(normalr);
  double var = arma::var(normalr);
  EXPECT_NEAR(0, mean, 0.02);
  EXPECT_NEAR(1, var, 0.02);
}

TEST_F(test_liabilitymodel, uniform_random) {
  arma::vec unifr = liabmodel.uniform_random(10);
  int num_n = unifr.size();
  EXPECT_EQ(10, num_n);
}

TEST_F(test_liabilitymodel, generate_causal_variants) {
  liabmodel.num_variants = 100;
  liabmodel.size_cluster = 5;
  liabmodel.num_cluster = 5;
  arma::Col<int> causal_var = liabmodel.generate_causal_variants(false);
  EXPECT_EQ(100, causal_var.size());
  EXPECT_EQ(25, sum(causal_var));
  causal_var.print(); 
  //EXPECT_EQ(1, causal_var(0));
}


TEST_F(test_liabilitymodel, effect_generation) {
  liabmodel.num_variants = 10;
  liabmodel.size_cluster = 1;
  liabmodel.num_cases = 1000;
  liabmodel.num_controls = 1000;
  liabmodel.num_subjects = 2000;

  liabmodel.genotype_matrix = test_mat;
  liabmodel.wished_effect = 0.1;
  arma::Col<int> causal_var = liabmodel.generate_causal_variants(true);
  arma::vec effect = liabmodel.effect_generation(causal_var);
  liabmodel.standardize_matrix();
  arma::vec risk = liabmodel.genotype_matrix_standarized * effect;
  EXPECT_NEAR(0.0, sum(risk), 0.001); 

  EXPECT_NEAR(0.1, pow(sum(effect),2), 0.001);
  EXPECT_EQ(10, effect.size());
  liabmodel.wished_effect = 0.0;
  effect = liabmodel.effect_generation(causal_var);
  EXPECT_NEAR(0.0, pow(sum(effect),2), 0.001);

  liabmodel.wished_effect = 0.001;
  effect = liabmodel.effect_generation(causal_var);
  effect.print();
  EXPECT_NEAR(0.001, pow(sum(effect),2), 0.0000001);
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
  liabmodel.wished_effect = 0.0;
  liabmodel.num_cases = 1000;
  liabmodel.num_controls = 1000;
  liabmodel.num_subjects = 2000;
  liabmodel.life_time_risk = 0.1;
  liabmodel.genotype_matrix = test_mat;
  liabmodel.standardize_matrix();
  int num_sub = liabmodel.num_subjects;
  arma::Col<int> causal_var = liabmodel.generate_causal_variants(true);
  arma::uvec sim_id = liabmodel.simulate_data(causal_var);
  double mean_risk = mean(liabmodel.liability_dist);
  double var_risk = var(liabmodel.liability_dist);
  EXPECT_EQ(num_sub, sim_id.size());
  EXPECT_EQ(100, liabmodel.liability_dist.size());

  EXPECT_NEAR(0, mean_risk, 0.10);
  EXPECT_NEAR(1, var_risk, 0.10);
}

//TEST_F(test_liabilitymodel, realdata) {
//  std::string ret_genotype_file = "./src/tests/test_data/saved_KIAA0556.csv";
//  arma::Mat<int> ret_mat;
//
//  ret_mat.load(ret_genotype_file);
//
//
//  liabmodel.wished_effect = 0.0;
//  liabmodel.num_cases = 1000;
//  liabmodel.num_controls = 1000;
//  liabmodel.num_subjects = 2000;
//  liabmodel.life_time_risk = 0.1;
//  liabmodel.genotype_matrix = ret_mat;
//  liabmodel.standardize_matrix();
//
//  arma::Col<int> phenotpe(liabmodel.num_subjects);
//  phenotpe.ones();
//  arma::Col<int> case_vec(liabmodel.num_cases);
//  case_vec.zeros();
//  case_vec = case_vec - 1;
//  phenotpe.subvec(liabmodel.num_controls, (liabmodel.num_subjects -1)) = case_vec;
//
//  models instance;
//  arma::mat pvalues_output(100,3);
//  pvalues_output.zeros();
//
//  auto causal_variants = liabmodel.generate_causal_variants(true);
//  // power calculations
//  int i, m;
//  for (i = 0; i < 100; ++i) {
//
//    arma::uvec simulated_pheno_id = liabmodel.simulate_data(causal_variants);
//    arma::Mat<int> temp_genotypes = liabmodel.genotype_matrix.rows(simulated_pheno_id);
//
//    // run models
//    for (m = 0; m < 3; ++m) {
//      pvalues_output(i, m) = instance.permutation(
//          instance.model_array[m], 1000, temp_genotypes, phenotpe, 1001);
//    }
//  }
//}
