#ifndef loadVCF_H
#define loadVCF_H
 
#include <armadillo>
#include <random>
#include <string>
#include <vector>
#include <VcfFileReader.h>
#include <easylogging++.h>

using namespace std;

/*! \brief Loads the genotypes from a vcf and variant file 
 *
 * This is a simple class which loads a genotype matrix from a vcf file.
 * The class requires a variant file as well which is a tab delim file
 * with chromosom,location, and gene name as columns. No headers are
 * allowed.
 */
class LoadVCF {
	public:
		int num_variants; /*!< number of variants in gene*/
		arma::Mat<int> genotype_matrix; /*!< genotype matrix */
		arma::mat genotype_matrix_standarized; /*!< standardized genotype matrix */
                string vcf_file; /*!< path of vcf file */
                string variant_file; /*!< path of variant file */
                VcfFileReader reader; /*!< reader of vcf file*/
                VcfHeader header; /*!< header of vcf file */
                std::random_device rd; /*!< random device */
                vector<string> genes; /*!< vector of genes in variant file*/
                vector<vector<string>> gene_loc; /*!< database of variants in gene*/

		/*! The LoadVCF constructor
		 * \param vcf_file a string with the vcf file path
		 * \param variant_file a string with the varaint file path
		 * */
                LoadVCF(string vcf_file, string variant_file) {
                  gene_loc = variant_location(variant_file);
                  genes = get_all_genes(gene_loc);
                  reader.open(vcf_file.c_str(), header);
                  reader.readVcfIndex();
                };
		/*! A simple LoadVCF construct
		*/
		LoadVCF(){};

                vector<vector<string>>
                variant_location(string VarFile, char sep = '\t');
                vector<string>
                get_all_genes(vector<vector<string>> &data);
                vector<vector<string>> get_gene_loc(string gene);

                int get_gene_matrix(vector<vector<string>> geneLoc);
		void load_gene(string gene_name);
                int generate_genotype_matrix(int subjects, int variants);
		arma::Col<int> bernulli(double p, int n);
};

#endif
