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
  std::string simmat = FLAGS_simmat;

  Simulation* sim = NULL;
  if (FLAGS_simmat == "simmat") {
    sim = new Simulation(FLAGS_numcases, FLAGS_numcontrols);
    auto gene_loc = sim->get_gene_loc(FLAGS_gene);
    sim->get_gene_matrix(gene_loc);

  } else if(FLAGS_simmat == "plink") {

    std::string bed = vcf_file+".bed";
    std::string bim = vcf_file+".bim";
    std::string fam = vcf_file+".fam";
    LOG(INFO) << "using a plink file";

    LOG(INFO) << "Input files are:\n"
      << bed << '\n'
      << bim << '\n'
      << fam << std::endl;

    sim = new Simulation(fam, bim, bed, variant_file,
        FLAGS_numcases, FLAGS_numcontrols);
    auto gene_loc = sim->get_gene_loc(FLAGS_gene);
    arma::Col<int> variants_select = sim->get_col_skip(gene_loc);
    arma::Col<int> row_select(sim->N, arma::fill::ones);
    sim->get_genotype_matrix(bed, variants_select, row_select);

  } else if (FLAGS_simmat == "vcf") {

    sim = new Simulation(vcf_file, variant_file, FLAGS_numcases, FLAGS_numcontrols);
    auto gene_loc = sim->get_gene_loc(FLAGS_gene);
    sim->get_gene_matrix(gene_loc);

  }

  if (sim->genotype_matrix.n_cols < 1) {
    throw std::runtime_error("genotype matrix has not columns"); }
  if (sim->genotype_matrix.n_rows < 1) {
    throw std::runtime_error("genotype matrix has not rows"); }
  if (arma::accu(sim->genotype_matrix) < 1) {
    throw std::runtime_error("genotype matrix has no variants"); }

  sim->standardize_matrix();
  sim->fixed_causal_var = 0;
  sim->max_test_iteration = FLAGS_iter;
  sim->life_time_risk = fLD::FLAGS_lifetimerisk;
  sim->num_controls = FLAGS_numcontrols;
  sim->num_cases = FLAGS_numcontrols;
  sim->num_subjects =  FLAGS_numcontrols + FLAGS_numcontrols;
  sim->wished_effect = FLAGS_minEffect;
  sim->num_cluster = FLAGS_numCluster;
  sim->min_percentage = 0;
  sim->current_percentage = FLAGS_percentageSteps;
  sim->steps_percentage = FLAGS_percentageSteps;
  sim->max_percentage = FLAGS_maxPercentage;
  sim->threads = fLI::FLAGS_threads;


  sim->select_test(FLAGS_tests);
  std::string power_output = FLAGS_path + "power_" + FLAGS_out;
  FILE *pFile = std::fopen(power_output.c_str(), "w");

  VLOG(9) << "Effect size: "<< sim->wished_effect << " to " << FLAGS_maxEffect;
  VLOG(9) << "Cluster size: " << sim->min_percentage << " to "
          << sim->max_percentage << " by " << sim->steps_percentage;
  VLOG(9) << "Number of cases: " << sim->num_cases;
  VLOG(9) << "Number of controls: " << sim->num_controls;
  VLOG(9) << "Lifetimerisk: " << sim->life_time_risk;
  VLOG(9) << "Max. iteration: " << sim->max_test_iteration;
  VLOG(9) << "min. iteration: " << sim->test_iteration;
  VLOG(9) << "runnig with " << sim->threads << " threads";
  VLOG(9) << "Loaded Matrix with " << sim->genotype_matrix.n_cols
          << " variants and " << sim->genotype_matrix.n_rows << " subjects";
  VLOG(9) << "Using  " << sim->current_percentage << " percentage of variants ";

  sim->writeoutput(pFile);
  sim->genotype_matrix.save("genotypes_" + FLAGS_gene +".txt", arma::csv_ascii);


  // cluster size loop
  VLOG(9) << "Finished setting up the simulation and will start computing now";
  bool cover_full_gene = true; int sim_id = 0;
  while(sim->current_percentage <= sim->max_percentage) {

    VLOG(1) << "Current Cluster size: " << sim->current_percentage;
    sim->num_causal_var();
    arma::Col<int> causal_variants(sim->num_variants, arma::fill::zeros);
    VLOG(1) << "Double check percentage as: " << sim->current_percentage;
    if (sim->current_percentage >= 0.99) {
      VLOG(9) << "setting causal variants to max";
      causal_variants = sim->generate_causal_variants(true);
    } else {
      VLOG(9) << "setting variants to cluster size";
      causal_variants = sim->generate_causal_variants(false);
    }
    sim->wished_effect = FLAGS_minEffect;
    while (sim->wished_effect <= FLAGS_maxEffect) {

      VLOG(1) << "Current Efffect size: " << sim->wished_effect;
      VLOG(1) << "Current Cluster size: " << sim->current_percentage;
      sim->power_calculation(FLAGS_powerIter, causal_variants);
      std::string output = FLAGS_path + std::to_string(sim_id) + FLAGS_out;
      sim->pvalues_output.save(output, arma::csv_ascii);
      std::string output_sim =
        FLAGS_path + std::to_string(sim_id) + "_sim_" + FLAGS_out;
      sim->saveSim.save(output_sim, arma::csv_ascii);
      std::cout << sim_id << std::endl;

      sim->writeoutput(pFile, sim->power[0], sim->power[1], sim->power[2],
          sim->power[3], sim->current_percentage,
          sim->real_num_causal,
          sim->wished_effect, sim_id);
      sim->wished_effect += FLAGS_effectSteps;
      sim_id += 1;
      VLOG(9) << sim_id;
    }
//    sim->current_percentage += sim->steps_percentage;
  }
  std::fclose(pFile);

  LOG(INFO) << "Simulation finished";

  return 0;
}
