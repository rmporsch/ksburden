#include <liability_model.h>
#include <simulation_flags.h>
#include <string>
#include <easylogging++.h>
#include <armadillo>
#include <models.h>

INITIALIZE_EASYLOGGINGPP

int main(int argc, char *argv[]) {

  gflags::ParseCommandLineFlags(&argc, &argv, true);
  el::Loggers::setVerboseLevel(FLAGS_verbose);
  START_EASYLOGGINGPP(argc, argv);

  LOG(INFO) << "Start simulation";

  std::string vcf_file = FLAGS_genotypes;
  std::string variant_file = FLAGS_variant;
  LiabilityModel sim(vcf_file, variant_file);
  sim.life_time_risk = fLD::FLAGS_lifetimerisk;
  sim.num_controls = 500;
  sim.num_cases = 500;
  sim.num_subjects =  1000;
  sim.life_time_risk = 0.1;
  sim.wished_effect = 0.0;
  sim.size_cluster = 1;

  auto gene_loc = sim.get_gene_loc(FLAGS_gene);
  sim.get_gene_matrix(gene_loc);
  sim.standardize_matrix();

  int num_models = 3; int m;
  models instance;

  // output storage
  arma::mat pvalues_output(fLI::FLAGS_powerIter, num_models);
  pvalues_output.zeros();

  auto causal_variants = sim.generate_causal_variants(true);
  int i = 0;

  // dummy phenotype
  arma::Col<int> phenotpe(sim.num_subjects);
  phenotpe.ones();
  arma::Col<int> case_vec(sim.num_cases);
  case_vec.zeros();
  case_vec = case_vec - 1;
  phenotpe.subvec(sim.num_controls, (sim.num_subjects -1)) = case_vec;

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

  pvalues_output.save("sometest.csv", arma::csv_ascii);
  LOG(INFO) << "Simulation finished";

  return 0;
}
