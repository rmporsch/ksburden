#include <models.h>
#include <flags_ksburden.h>
#include <data_analysis.h>
#include <armadillo>
#include <easylogging++.h>

INITIALIZE_EASYLOGGINGPP

int main(int argc, char *argv[])
{
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  el::Loggers::setVerboseLevel(FLAGS_verbose);
  START_EASYLOGGINGPP(argc, argv);

  std::string ped_file= fLS::FLAGS_ped;
  std::string vcf_file = fLS::FLAGS_vcf;
  std::string variant_file = fLS::FLAGS_variant;
  std::string fileName = fLS::FLAGS_out;
  VLOG(9) << "Got input from: \n ped file: " << ped_file
          << "\n variant file: " << variant_file << "\n vcf file: " << vcf_file
          << "\n output is writen to: " << fileName;

  FILE *fout = fopen(fileName.c_str(), "w");
  fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\n","Gene", "Ks","Burden","CMC","ksburden", "NumVar");

  analysis data(ped_file, vcf_file, variant_file);
  models instance;
  int num_models = 3; int m;

  VLOG(9) << "starting compuation";
  for (auto g = data.genes.begin(); g != data.genes.end(); ++g) {
    VLOG(9) << *g;
    std::vector<std::vector<std::string>> gene_loc = data.get_gene_loc(*g);
    data.get_gene_matrix(gene_loc);
    int num_var_gene = data.genotype_matrix.n_cols;
    //data.genotype_matrix = data.genotype_matrix.rows(data.id_genotype_include);
    // run models
    std::vector<double> output;
    output.resize(num_models);
    for (m = 0; m < num_models; ++m) {
      output[m] = instance.permutation(
          instance.model_array[m], FLAGS_iter,
          data.genotype_matrix.rows(data.id_genotype_include), data.phenotype,
          100000);
    }

    arma::vec fisher_comb = {output[0], output[1]};
    double p_fisher = instance.fisher(fisher_comb);



    fprintf(fout, "%s\t%f\t%f\t%f\t%f\t%i\n",
		    g->c_str(),
		    output[0],
		    output[1],
		    output[2],
		    p_fisher,
		    num_var_gene);
  }
  return 0;
}
