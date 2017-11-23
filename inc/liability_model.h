#ifndef liability_model_H
#define liability_model_H

#include <load_vcf.h>
#include <load_plink.h>
#include <armadillo>
#include <random>
#include <easylogging++.h>
#include <boost/math/distributions/normal.hpp>

/*! \brief Generates a liability model given paramters 
 *
 * This class does the heavy lifting of the simulation by
 * generating disease affected and unaffected subjects based on a given 
 * genotype matrix. It is fairly simple writen. In addtion there are also
 * some helper function included. The main function is simulate_data which does 
 * most of the work.
 */
class LiabilityModel : public LoadVCF, public LoadPlink
{
private:
	

public:
  arma::mat genotype_matrix_standarized; /*!< standardized genotype matrix */
	int num_cluster = 1; /*!< number of clusters */
	int size_cluster = 1; /*!< size of eachj cluster*/
	double causal_probability = 1.0; /*!< mutation probability in cluster */
  int real_num_causal; /*!< number of real causal mutations used*/
	int num_cases; /*!< number of cases to simulate */
	int num_controls; /*!< number of controls to simulate */
	int num_subjects = num_cases + num_controls; /*!< number of subjects to simulate */
	arma::vec liability_dist; /*!< vector of the liabilities of each subject */
	double wished_effect; /*!< desired effect size */
	double life_time_risk; /*!< life time risk in percentage */

	arma::vec normal_random(int n, double mean, double stdev);
	arma::vec uniform_random(int n);
  arma::Col<int> generate_causal_variants(bool EmptyStart);
  arma::vec effect_generation(arma::Col<int> causal);
  arma::uvec simulate_data(arma::Col<int> causal);
  arma::colvec standardize(arma::Col<int> variant);
  void standardize_matrix();

  LiabilityModel(std::string vcf_file, std::string variant_file)
    : LoadVCF(vcf_file, variant_file) {};
  LiabilityModel(
       std::string fam,
       std::string bim,
       std::string bam, std::string variant_file)
    : LoadPlink( fam, bim, bam, variant_file) {};
  LiabilityModel() : LoadVCF(){};
};

#endif /* LIABILITY_MODEL_H */
