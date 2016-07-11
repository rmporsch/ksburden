#ifndef liability_model_H
#define liability_model_H

#include <load_vcf.h>
#include <armadillo>
#include <random>
#include <easylogging++.h>
#include <boost/math/distributions/normal.hpp>

class LiabilityModel : public LoadVCF
{
private:
	

public:
	int num_cluster = 1;
	int size_cluster = 1;
	double causal_probability = 1.0;
	int num_cases;
	int num_controls;
	int num_subjects = num_cases + num_controls;
	arma::vec liability_dist;
	double wished_effect;
	double life_time_risk;
	arma::vec normal_random(int n, double mean, double stdev);
	arma::vec uniform_random(int n);
	arma::Col<int> generate_causal_variants(bool EmptyStart);
	arma::vec effect_generation(arma::Col<int> causal);
        arma::uvec simulate_data(arma::Col<int> causal);
        arma::colvec standardize(arma::Col<int> variant);
        void standardize_matrix();

        LiabilityModel(std::string vcf_file, std::string variant_file)
		: LoadVCF(vcf_file, variant_file) {};
	LiabilityModel() : LoadVCF(){};
};

#endif /* LIABILITY_MODEL_H */
