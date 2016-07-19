#include <liability_model.h>
#include <models.h>
#include <simulation.h>
#include <string>
#include <easylogging++.h>
#include <armadillo>
#include <models.h>
#include <algorithm>

void Simulation::select_test(std::string tests_input) {
	std::istringstream ss(tests_input);
	std::string token;
	perform_tests.resize(0);

	while (std::getline(ss, token, ',')) {
		perform_tests.push_back(token);
	}
	id_perform_models.resize(0);

        for (int i = 0; i < perform_tests.size(); ++i) {

          if (std::find(available_tests.begin(), available_tests.end(), perform_tests[i]) !=
              available_tests.end()) {
		  id_perform_models.push_back(i);
          }
        }
}

void Simulation::power_calculation(int power_iter, arma::Col<int> causal_variants) {
// power calculations
int i;
std::vector<int>::iterator m;
pvalues_output.set_size(power_iter, id_perform_models.size());
pvalues_output.ones();
omp_set_dynamic(threads);
#pragma omp parallel for private(m)
for (i = 0; i < power_iter; ++i) {

    arma::uvec simulated_pheno_id = simulate_data(causal_variants);
    arma::Mat<int> temp_genotypes =
        genotype_matrix.rows(simulated_pheno_id);
    // run models
    for (m=id_perform_models.begin(); m<id_perform_models.end(); ++m) {
      pvalues_output(i, *m) =
          permutation(model_array[*m], test_iteration,
                               temp_genotypes, phenotype, max_test_iteration);
    }
  }

  // calculate power
  arma::mat::col_iterator a = pvalues_output.begin_col(0);
  arma::mat::col_iterator b =
      pvalues_output.end_col((pvalues_output.n_cols - 1));

  power.set_size(id_perform_models.size());
  for(int t_test = 0; t_test < id_perform_models.size(); ++t_test){
    arma::uvec temp = arma::find(pvalues_output.col(t_test) <= 0.05);
    power(t_test) = temp.size() / (double)power_iter;
  }
}

bool Simulation::num_causal_var() {

  if (fixed_causal_var > 0) {
    size_cluster = fixed_causal_var;
  }

  if (fixed_causal_var == 0) {
    int temp_size_cluster =
        std::floor((num_variants * current_percentage) + 0.0001);

    while (temp_size_cluster <= size_cluster) {
      current_percentage += steps_percentage;
      temp_size_cluster =
          std::floor((num_variants * current_percentage) + 0.0001);
    }
    size_cluster = temp_size_cluster;
    std::cout << size_cluster << std::endl;

    if (size_cluster >= num_variants) {
      return false;
    } else {
      return true;
    }
  }
}

void Simulation::writehead(FILE *pFile) {
    for (int i = 0; i < perform_tests.size(); ++i) {
      output_file_body = "%f\t" + output_file_body;
      output_file_header = "%s\t" + output_file_header;
    }
    switch (perform_tests.size()) {
      case 1
          : fprintf(pFile, output_file_header.c_str(), perform_tests[0].c_str(),
                    "percentage", "num_car", "effect_size", "simID");
      break;
    case 2:
      fprintf(pFile, output_file_header.c_str(), perform_tests[0].c_str(),
              perform_tests[1].c_str(), "percentage", "num_car", "effect_size",
              "simID");
      break;
    case 3:
      fprintf(pFile, output_file_header.c_str(), perform_tests[0].c_str(),
              perform_tests[1].c_str(), perform_tests[2].c_str(), "percentage",
              "num_car", "effect_size", "simID");
      break;
    }
}

void Simulation::writeoutput(FILE *pFile, ...) {
  if (!wrote_head) {
    writehead(pFile);
    wrote_head = true;
  } else {
  va_list argptr;
  va_start(argptr, output_file_body.c_str());
  vfprintf(pFile, output_file_body.c_str(), argptr);
  va_end(argptr);
  }
}
