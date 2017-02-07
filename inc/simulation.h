#ifndef SIMULATION_H
#define SIMULATION_H

#include <liability_model.h>
#include <models.h>
#include <armadillo>

/*! \brief simulates various scenarios
 *
 * Simulation class to simulate a number of different situations.
 * The class mainly contains helper functions to run the simulations 
 * wile the LiabilityModel does most of the heavy lifting.
 */
class Simulation : public LiabilityModel, public models {
public:
  int test_iteration = 1000; /*!< number of MC iterations for each model*/
  int max_test_iteration = 100000; /*!< max. number of iteration for each model*/
  arma::Col<int> phenotype; /*!< integer vector for phenotype */
  std::vector<int> id_perform_models; /*!< integer vector of models to perform */
  std::vector<std::string> perform_tests;  /*!< vector with the names of models to compute */
  arma::vec power; /*!< output vector containing the stat. power of a simulation*/
  arma::Mat<int> saveSim; /*! matrix which holds all simulations */

  int fixed_causal_var = 0; /*!< fixed number of causal variants */
  double current_percentage; /*!< current percentage of the gene coverd by causal mutations*/
  double max_percentage; /*!< max percentage of the gene to cover in causal mutations */
  double min_percentage; /*!< min percentage of the gene to cover in causal mutations */
  double steps_percentage; /*!< steps of percentages */
  arma::mat pvalues_output; /*!< matrix of p-values */

  int input_num_cases; /*!< induced number of cases*/
  int input_num_controls; /*!< induced number of controls */
  int threads = 1; /*!< number of threads to use */

  /*! Simulation contructor
   *
   * \param vcf_file path to vcf file
   * \param variant_file path to variant file
   * \param input_num_cases number of cases to simulate
   * \param input_num_controls number of controls to simulate
   */
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
   Simulation(
       std::string fam_file,
       std::string bim_file,
       std::string bam_file,
       std::string variant_file,
       int input_num_cases,
       int input_num_controls)
      : LiabilityModel(fam_file, bim_file, bam_file, variant_file) {

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

  Simulation(int input_num_cases, int input_num_controls)
      : LiabilityModel() {
    num_cases = input_num_cases;
    num_controls = input_num_controls;
    num_subjects = num_cases + num_controls;
    phenotype.set_size(input_num_cases + input_num_controls);
    phenotype.ones();

    arma::Col<int> case_vec(num_cases);
    case_vec.zeros();
    case_vec = case_vec - 1;
    phenotype.subvec(num_controls, (num_subjects - 1)) = case_vec;
  }
  void select_test(std::string tests_input);
  void power_calculation(int power_iter, arma::Col<int> causal_variants);
  bool num_causal_var();
  void writeoutput(FILE *pFile, ...);
  void writehead(FILE *pFile);
  bool wrote_head = false; /*!< alrady wrote the head of the file? */
  std::string output_file_body = "%f\t%i\t%f\t%i\n"; /*!< basic output body format */
  std::string output_file_header = "%s\t%s\t%s\t%s\n"; /*!< basic output header format */

private:
};
#endif /* SIMULATION_H */
