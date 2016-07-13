#include <liability_model.h>
#include <models.h>
#include <simulation.h>
#include <string>
#include <easylogging++.h>
#include <armadillo>
#include <models.h>

void Simulation::select_test(std::string tests_input) {
	std::istringstream ss(tests_input);
	std::string token;
	std::vector<std::string> perfrom_tests;

	while (std::getline(ss, token, ',')) {
		perfrom_tests.push_back(token);
	}

        for (int i = 0; i < available_tests.size(); ++i) {

          if (std::find(perfrom_tests.begin(), perfrom_tests.end(), available_tests[i]) !=
              perfrom_tests.end()) {
		  id_perform_models.push_back(i);
          }
        }
}

void Simulation::power_calculation(int power_iter, arma::Col<int> causal_variants) {
// power calculations
int i;
std::vector<int>::iterator m;
arma::mat pvalues_output(power_iter, id_perform_models.size(), arma::fill::zeros);
#pragma omp parallel for private(m)
  for (i = 0; i < power_iter; ++i) {

    arma::uvec simulated_pheno_id = simulate_data(causal_variants);
    arma::Mat<int> temp_genotypes =
        genotype_matrix.rows(simulated_pheno_id);

    // run models
    for (m=id_perform_models.begin(); m<id_perform_models.end(); m++) {
      pvalues_output(i, *m) =
          permutation(model_array[*m], test_iteration,
                               temp_genotypes, phenotype, max_test_iteration);
    }
  }

  // calculate power
  arma::mat::col_iterator a = pvalues_output.begin_col(0);
  arma::mat::col_iterator b =
      pvalues_output.end_col((pvalues_output.n_cols - 1));

  int t_test = 0;
  for(t_test = 0; t_test < id_perform_models.size(); ++t_test){
    arma::uvec temp = arma::find(pvalues_output.col(t_test) <= 0.05);
    power(t_test) = temp.size() / (double)power_iter;
    ++t_test;
  }
}

bool Simulation::num_causal_var() {

  if (fixed_causal_var > 0) {
    size_cluster = fixed_causal_var;
  }

  if (fixed_causal_var == 0) {
    int temp_size_cluster =
        std::floor((num_variants * current_percentage) + 0.0001);

    while (temp_size_cluster < size_cluster) {
      current_percentage += steps_percentage;
      int temp_size_cluster =
          std::floor((num_variants * current_percentage) + 0.0001);
    }
    size_cluster = temp_size_cluster;
    if (size_cluster >= num_variants) {
      return false;
    } else {
      return true;
    }
  }
}
