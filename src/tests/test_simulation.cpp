#include <gtest/gtest.h>
#include <simulation.h>

class test_simulation : public testing::Test
{
public:
	Simulation sim;
        arma::Mat<int> test_mat;
        test_simulation() {

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

	  sim.genotype_matrix = test_mat;
        };
};

TEST_F(test_simulation, select_test)
{
	std::string test_string = "ksburden,burden";
	sim.select_test(test_string);
	int num_tests = sim.id_perform_models.size();
	EXPECT_EQ(2, num_tests);
	EXPECT_EQ(0, sim.id_perform_models[0]);
	EXPECT_EQ(1, sim.id_perform_models[1]);

	test_string = "ksburden";
	sim.select_test(test_string);
	num_tests = sim.id_perform_models.size();
	
	EXPECT_EQ(1, num_tests);
	EXPECT_EQ(0, sim.id_perform_models[0]);

	test_string = "ksburden,something";
	sim.select_test(test_string);
	num_tests = sim.id_perform_models.size();
	
	EXPECT_EQ(1, num_tests);
	EXPECT_EQ(0, sim.id_perform_models[0]);
}

TEST_F(test_simulation, num_causal_var) {
	sim.fixed_causal_var = 1;
	sim.num_variants = 10;
	bool executed = sim.num_causal_var();
	EXPECT_EQ(1, sim.size_cluster); 
	EXPECT_TRUE(executed);
	sim.fixed_causal_var = 10;
	sim.num_variants = 1;
	executed = sim.num_causal_var();
	EXPECT_FALSE(executed);
}

TEST_F(test_simulation, power_calculation) {
  sim.num_cases = 50;
  sim.num_controls = 50;
  sim.num_subjects = 100;
  sim.num_variants = 10;
  sim.size_cluster = 1;
  sim.wished_effect = 0.0;

  std::string test_string = "ksburden,burden,cmc";
  sim.select_test(test_string);

  sim.life_time_risk = 0.1;
  arma::Col<int> phenotype;
  phenotype.set_size(100);
  phenotype.ones();
  arma::Col<int> case_vec(50, arma::fill::zeros);
  case_vec = case_vec - 1;
  phenotype.subvec(50, (100 - 1)) = case_vec;
  sim.phenotype = phenotype;
  arma::Col<int> causal_variants = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  sim.standardize_matrix();
  
  sim.power_calculation(100, causal_variants);

  EXPECT_EQ(100, sim.pvalues_output.n_rows);
  EXPECT_EQ(4, sim.pvalues_output.n_cols);
  EXPECT_EQ(4, sim.power.size());
  EXPECT_NEAR(0.05, sim.power(0), 0.04);
  EXPECT_NEAR(0.05, sim.power(1), 0.04);
  EXPECT_NEAR(0.05, sim.power(2), 0.04);
}

TEST_F(test_simulation, writehead) {
	std::string test_string = "ksburden,burden";
	sim.select_test(test_string);
	FILE * pFile = std::fopen("sometest.txt", "w");
	sim.writehead(pFile);
	std::cout << sim.id_perform_models.size() << std::endl;
	std::cout << sim.output_file_body << std::endl;
	std::cout << sim.output_file_header << std::endl;
}
