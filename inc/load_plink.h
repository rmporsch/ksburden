#ifndef loadPlink_H
#define loadPlink_H
 
#include <armadillo>
#include <random>
#include <string>
#include <vector>
#include <load_vcf.h>
#include <easylogging++.h>

class LoadPlink {
  public:
		int num_variants; /*!< number of variants in gene*/
		arma::Mat<int> genotype_matrix; /*!< genotype matrix */
    string vcf_file; /*!< path of vcf file */
    string variant_file; /*!< path of variant file */
    std::random_device rd; /*!< random device */
    vector<string> genes; /*!< vector of genes in variant file*/
    vector<vector<string>> gene_loc; /*!< database of variants in gene*/

    LoadPlink(string variant_file) {
      VLOG(9) << "loading variant file";
      gene_loc = variant_location(variant_file);
      genes = get_all_genes(gene_loc);
    };
    /*! A simple LoadVCF construct
    */
    LoadPlink(){};

    vector<vector<string>>
      variant_location(string VarFile, char sep = '\t');
    vector<string>
      get_all_genes(vector<vector<string>> &data);

    vector<vector<string>> get_gene_loc(string gene);
    bool openPlinkBinaryFile(const std::string s, std::ifstream &BIT);
    arma::mat read_genotype_matrix(const std::string fileName,
        int N, int P, const arma::mat input, 
        arma::Col<int> col_skip_pos, arma::Col<int> col_skip, 
        arma::Col<int> keepbytes, arma::Col<int> keepoffset);
    int countlines(const char* fileName);
};

#endif

