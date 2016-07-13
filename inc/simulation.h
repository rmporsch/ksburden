#ifndef SIMULATION_H
#define SIMULATION_H

#include <liability_model.h>
#include <models.h>
#include <armadillo>

class Simulation : public LiabilityModel, public models
{
private:
	

public:
	int test_iteration = 1000;
	int max_test_iteration = 100000;
	arma::Col<int> phenotype;
	std::vector<int> id_perform_models;
	arma::vec power;

	int fixed_causal_var;
	double current_percentage;
	double max_percentage;
	double min_percentage;
	double steps_percentage;


        Simulation(std::string vcf_file, std::string variant_file)
            : LiabilityModel(vcf_file, variant_file){};
        Simulation() : LiabilityModel(){};

        void select_test(std::string tests_input);
	void power_calculation(int power_iter, arma::Col<int> causal_variants);
        bool num_causal_var();

};
#endif /* SIMULATION_H */
