#ifndef loadPlink_H
#define loadPlink_H

#include <armadillo>
#include <random>
#include <string>
#include <vector>
#include <easylogging++.h>
#include <load_variantFile.h>

using namespace std;

class LoadPlink: virtual public VariantFile {
  public:
    int P; /*!< number of variants in gene*/
    int N;
    string plink_file; /*!< path of vcf file */
    vector<string> genes; /*!< vector of genes in variant file*/
    vector<vector<string>> gene_loc; /*!< database of variants in gene*/
    vector<vector<string>> bim_file; /*!< bim file content*/
    vector<vector<string>> fam_file; /*!< bim file content*/

    LoadPlink(std::string fam,
        std::string bim,
        std::string bam,
        std::string variant_file): VariantFile(variant_file) {
      LOG(INFO) << "loading bim and fam file";
      bim_file = variant_location(bim, '\t');
      fam_file = variant_location(fam, ' ');
      P = bim_file.size();
      N = fam_file.size();
      LOG(INFO) << "read " << P << " variants in bim file";
      LOG(INFO) << "read " << N << " subjects in fam file";
    };
    /*! A simple LoadPlink construct
    */
    LoadPlink(): VariantFile() {};

    bool openPlinkBinaryFile(const std::string s, std::ifstream &BIT);

    arma::Col<int> get_col_skip(vector<vector<string>> gene_loc);

    int get_genotype_matrix(const std::string fileName,
        arma::Col<int> col_skip,
        arma::Col<int> row_skip);

    int countlines(const char* fileName);
};

#endif
