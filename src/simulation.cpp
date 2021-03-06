#include <liability_model.h>
#include <models.h>
#include <simulation.h>
#include <string>
#include <easylogging++.h>
#include <armadillo>
#include <math.h>
#include <models.h>
#include <algorithm>

/*! \brief helper function fo select tests
 * \param test_input a string of tests seperate by comma
 */
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
  if (perform_tests.size()==3) {
    perform_tests.push_back("cksburden");
  }
}

/*! computes power for a number of models
 *
 * Calculates the empirical power for each model. Results are saved into the power vector
 *
 * \param power_iter number of iterations to compute power
 * \param causal_variants a vector of causal variants
 *
 */
void Simulation::power_calculation(int power_iter, arma::Col<int> causal_variants) {
  // power calculations
  int i;
  std::vector<int>::iterator m;
  pvalues_output.set_size(power_iter, id_perform_models.size());
  pvalues_output.ones();
  omp_set_num_threads(threads);
  saveSim.set_size(power_iter, num_subjects);
  VLOG(9) << "sum of causal variants is " << arma::sum(causal_variants);
  VLOG(9) << causal_variants;

#pragma omp parallel for private(m)
  for (i = 0; i < power_iter; ++i) {

    // VLOG(9) << "simulating iteration " << i;
    arma::uvec simulated_pheno_id = simulate_data(causal_variants);
    arma::Mat<int> temp_genotypes =
      genotype_matrix.rows(simulated_pheno_id);
    saveSim.row(i) = arma::conv_to<arma::Row<int>>::from(simulated_pheno_id);
    // run models
    for (m=id_perform_models.begin(); m<id_perform_models.end(); ++m) {
      pvalues_output(i, *m) =
        permutation(model_array[*m], test_iteration,
            temp_genotypes, phenotype, max_test_iteration);
    }
  }

  if (id_perform_models.size() >= 2) {
    pvalues_output.resize(power_iter, (id_perform_models.size() + 1));
    arma::mat cksburden = pvalues_output.cols(0, 1);
    pvalues_output.col(id_perform_models.size()) =
      large_fisher(cksburden);
  }


  power.set_size(pvalues_output.n_cols);
  for(int t_test = 0; t_test < pvalues_output.n_cols; ++t_test){
    arma::uvec temp = arma::find(pvalues_output.col(t_test) <= pvalue_threshold);
    power(t_test) = temp.size() / (double)power_iter;
    VLOG(9) << "Power for test " << (t_test + 1) << " is " << power(t_test);
  }
}


/*! \brief define number of causal varinats
 *
 * This is a helper function which can take the current and desired 
 * percentage of a gene covered by causal mutations and defines the raw count of
 * causal mutations
 */
bool Simulation::num_causal_var() {

  if (fixed_causal_var > 0) {
    size_cluster = fixed_causal_var;
  }
  if (num_variants < 1) {
    throw std::runtime_error("no variant present");}


  if (fixed_causal_var == 0 & steps_percentage>0) {
    int temp_size_cluster =
      std::floor((num_variants * current_percentage) + 0.0001);

    VLOG(9) << "computing new cluster size";
    VLOG(9) << "cluster size is set to " << size_cluster;
    VLOG(9) << "number of variants is " << num_variants;
    while (temp_size_cluster <= size_cluster) {
      current_percentage += steps_percentage;
      VLOG(9) << "percentage is now (temp) " << current_percentage;
      temp_size_cluster =
        std::floor((num_variants * current_percentage) + 0.0001);
      VLOG(9) << "temp cluster size is " << temp_size_cluster;
      if (current_percentage > 1) {
        VLOG(9) << "Current percentage is too large with " << current_percentage;
        throw std::runtime_error("percentage is higher than 1");
      }
    }
    size_cluster = temp_size_cluster;
    VLOG(9) << "cluster size is now " << size_cluster;

    if (size_cluster >= num_variants) {
      return false;
    } else {
      return true;
    }
  } else {
    size_cluster = std::floor((num_variants * current_percentage) + 0.0001);
    return true;
  }
}

/*! \brief writes the header of the outputfile
 *
 * \param pFile a stream
 */
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
    case 4:
    fprintf(pFile, output_file_header.c_str(), perform_tests[0].c_str(),
        perform_tests[1].c_str(), perform_tests[2].c_str(), perform_tests[3].c_str(), "percentage",
        "num_car", "effect_size", "simID");
    break;
  }
}

/*! \brief writes the output file
 *
 * \param pFile a stream
 * \param ... additional things to print into the outputfile
 *
 * This is a warpper around writehead. It takes a file name and additional paramters 
 * to print.  Important is that this function needs to be called first without any additional 
 * arguments in order to write the header.
 */
void Simulation::writeoutput(FILE *pFile, ...) {
  if (!wrote_head) {
    writehead(pFile);
    wrote_head = true;
  } else {
    std::cout << "writing now" << std::endl;
    va_list argptr;
    va_start(argptr, output_file_body.c_str());
    vfprintf(pFile, output_file_body.c_str(), argptr);
    va_end(argptr);
  }
}
