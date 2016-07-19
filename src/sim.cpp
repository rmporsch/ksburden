#include <simulation.h>
#include <simulation_flags.h>
#include <string>
#include <easylogging++.h>

INITIALIZE_EASYLOGGINGPP

int main(int argc, char *argv[]) {

  gflags::ParseCommandLineFlags(&argc, &argv, true);
  el::Loggers::setVerboseLevel(FLAGS_verbose);
  START_EASYLOGGINGPP(argc, argv);
  LOG(INFO) << "Start simulation";

  std::string vcf_file = FLAGS_genotypes;
  std::string variant_file = FLAGS_variant;

  Simulation sim(vcf_file, variant_file, FLAGS_numcases, FLAGS_numcontrols);
  sim.life_time_risk = fLD::FLAGS_lifetimerisk;
  sim.num_controls = FLAGS_numcontrols;
  sim.num_cases = FLAGS_numcontrols;
  sim.num_subjects =  FLAGS_numcontrols + FLAGS_numcontrols;
  sim.wished_effect = FLAGS_minEffect;
  sim.size_cluster = 1;
  sim.min_percentage = 0.1;
  sim.current_percentage = 0.1;
  sim.steps_percentage = FLAGS_percentageSteps;
  sim.max_percentage = 1.0;

  auto gene_loc = sim.get_gene_loc(FLAGS_gene);
  sim.get_gene_matrix(gene_loc);
  sim.standardize_matrix();

  sim.select_test(FLAGS_tests);
  arma::Col<int> causal_variants = sim.generate_causal_variants(true);
  int i = 0;

  std::string power_output = FLAGS_path + "power_" + FLAGS_out;
  FILE *pFile = std::fopen(power_output.c_str(), "w");
  sim.writeoutput(pFile);

  sim.threads = fLI::FLAGS_threads;

  // cluster size loop
  bool cover_full_gene = true;
  int sim_id = 0;
  while (cover_full_gene) {
    cover_full_gene = sim.num_causal_var();
    causal_variants = sim.generate_causal_variants(true);
    std::cout << sim.current_percentage << std::endl;
    sim.wished_effect = FLAGS_minEffect;
    while (sim.wished_effect <= FLAGS_maxEffect) {
      sim.power_calculation(FLAGS_powerIter, causal_variants);
      sim.wished_effect += FLAGS_effectSteps;
      std::string output = FLAGS_path + std::to_string(sim_id) + FLAGS_out;
      sim.pvalues_output.save(output, arma::csv_ascii);
      sim_id += 1;

      sim.writeoutput(pFile, sim.power[0], sim.power[1], sim.power[2],
                      sim.current_percentage, sim.size_cluster,
                      sim.wished_effect, sim_id);
    }
  }
  std::fclose(pFile);

  LOG(INFO) << "Simulation finished";

  return 0;
}
