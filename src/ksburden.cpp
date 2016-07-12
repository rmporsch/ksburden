#include <models.h>
#include <flags_ksburden.h>
#include <data_analysis.h>
#include <armadillo>
#include <easylogging++.h>

INITIALIZE_EASYLOGGINGPP

int main(int argc, char *argv[])
{
  START_EASYLOGGINGPP(argc, argv);
  sanity_check(argc, argv);

  std::string ped_file= fLS::FLAGS_ped;
  std::string vcf_file = fLS::FLAGS_vcf;
  std::string variant_file = fLS::FLAGS_variant;
  std::string fileName = fLS::FLAGS_out;

  FILE *fout = fopen(fileName.c_str(), "w");
  fprintf(fout, "%s\t%s\t%s\t%s\t%s\n","Gene", "Ks","Burden","CMC", "NumVar");

  analysis data(ped_file, vcf_file, variant_file);
  models instance;
  int num_models = 3; int m;

  for (auto g = data.genes.begin(); g != data.genes.end(); ++g) {
    std::vector<std::vector<std::string>> gene_loc = data.get_gene_loc(*g);
    data.get_gene_matrix(gene_loc);
    int num_var_gene = data.genotype_matrix.n_cols;
    // run models
    std::vector<double> output;
    output.resize(num_models);
    for (m = 0; m < num_models; ++m) {
      output[m] = instance.permutation(
          instance.model_array[m], FLAGS_iter, data.genotype_matrix, data.phenotype, 100000);
    }
    fprintf(fout, "%s\t%f\t%f\t%f\t%i\n",
		    g->c_str(),
		    output[0],
		    output[1],
		    output[2],
		    num_var_gene);
  }
  return 0;
}
