#ifndef loadVCF_H
#define loadVCF_H
 
#include <armadillo>
#include <random>
#include <string>
#include <vector>
#include <VcfFileReader.h>

using namespace std;

class LoadVCF {
	public:
		int num_variants;
		arma::Mat<int> genotype_matrix;
		arma::mat genotype_matrix_standarized;
                string vcf_file;
                string variant_file;
                VcfFileReader reader;
                VcfHeader header;
                std::random_device rd;
                vector<string> genes;
                vector<vector<string>> gene_loc;

                LoadVCF(string vcf_file, string variant_file) {
                  gene_loc = variant_location(variant_file);
                  genes = get_all_genes(gene_loc);
                  reader.open(vcf_file.c_str(), header);
                  reader.readVcfIndex();
                };

                vector<vector<string>>
                variant_location(string VarFile, char sep = '\t');
                vector<string>
                get_all_genes(vector<vector<string>> &data);
                vector<vector<string>>
                get_gene_loc(string gene);

                int get_gene_matrix(vector<vector<string>> geneLoc);
		void load_gene(string gene_name);
                int generate_genotype_matrix(int subjects, int variants);
		arma::Col<int> bernulli(double p, int n);
};

#endif
