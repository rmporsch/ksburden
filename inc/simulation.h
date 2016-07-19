#ifndef SIMULATION_H
#define SIMULATION_H

#include <liability_model.h>
#include <models.h>
#include <armadillo>

class Simulation : public LiabilityModel, public models {
private:
public:
  int test_iteration = 1000;
  int max_test_iteration = 100000;
  arma::Col<int> phenotype;
  std::vector<int> id_perform_models;
  std::vector<std::string> perform_tests; // vector containing the test names
  arma::vec power;

  int fixed_causal_var;
  double current_percentage;
  double max_percentage;
  double min_percentage;
  double steps_percentage;
  arma::mat pvalues_output;

  int input_num_cases;
  int input_num_controls;
  int threads = 1;

  Simulation(std::string vcf_file, std::string variant_file,
             int input_num_cases, int input_num_controls)
      : LiabilityModel(vcf_file, variant_file) {

    num_cases = input_num_cases;
    num_controls = input_num_controls;
    num_subjects = num_cases + num_controls;
    phenotype.set_size(input_num_cases + input_num_controls);
    phenotype.ones();

    arma::Col<int> case_vec(num_cases);
    case_vec.zeros();
    case_vec = case_vec - 1;
    phenotype.subvec(num_controls, (num_subjects - 1)) = case_vec;
  };
  Simulation() : LiabilityModel(){};

  void select_test(std::string tests_input);
  void power_calculation(int power_iter, arma::Col<int> causal_variants);
  bool num_causal_var();
  void writeoutput(FILE *pFile, ...);
  void writehead(FILE *pFile);
  bool wrote_head = false;
  std::string output_file_body = "%f\t%i\t%f\t%i\n";
  std::string output_file_header = "%s\t%s\t%s\t%s\n";

private:
};
#endif /* SIMULATION_H */
