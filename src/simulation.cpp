#include <liability_model.h>
#include <simulation_flags.h>
#include <string>
#include <easylogging++.h>
#include <armadillo>
#include <models.h>

INITIALIZE_EASYLOGGINGPP

int main(int argc, char *argv[]) {

  std::string vcf_file = FLAGS_genotypes;
  std::string variant_file = FLAGS_variant;
  LiabilityModel sim(vcf_file, variant_file);
  auto gene_loc = sim.get_gene_loc(FLAGS_gene);
  sim.get_gene_matrix(gene_loc);
  sim.standardize_matrix();

  sim.life_time_risk = fLD::FLAGS_lifetimerisk;
  sim.num_controls = 500;
  sim.num_cases = 500;


  int num_models = 3; int m;
  models instance;

  arma::mat pvalues_output(fLI::FLAGS_powerIter, num_models);
  pvalues_output.zeros();

  auto causal_variants = sim.generate_causal_variants(true);
  int i = 0;
  arma::Col<int> phenotpe(sim.num_subjects);
  phenotpe.zeros();
  phenotpe.subvec(0, (sim.num_cases - 1)) = 1;

  // power calculations
  for (i = 0; i < fLI::FLAGS_powerIter; ++i) {

    arma::uvec simulated_pheno_id = sim.simulate_data(causal_variants);
    arma::Mat<int> temp_genotypes = sim.genotype_matrix.rows(simulated_pheno_id);

    // run models
    for (m = 0; m < num_models; ++m) {
      pvalues_output(i, m) = instance.permutation(
          instance.model_array[m], fLI::FLAGS_iter, temp_genotypes, phenotpe, 100000);
    }
  }

  return 0;
}
