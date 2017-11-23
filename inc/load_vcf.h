#ifndef loadVCF_H
#define loadVCF_H
 
#include <armadillo>
#include <random>
#include <string>
#include <vector>
#include <VcfFileReader.h>
#include <easylogging++.h>
#include <load_variantFile.h>

using namespace std;

/*! \brief Loads the genotypes from a vcf and variant file 
 *
 * This is a simple class which loads a genotype matrix from a vcf file.
 * The class requires a variant file as well which is a tab delim file
 * with chromosom,location, and gene name as columns. No headers are
 * allowed.
 */
class LoadVCF: virtual public VariantFile {
	public:
    string vcf_file; /*!< path of vcf file */
    string variant_file; /*!< path of variant file */
    random_device rd; /*!< random device */
    VcfFileReader reader; /*!< reader of vcf file*/
    VcfHeader header; /*!< header of vcf file */

    /*! The LoadVCF constructor
     * \param vcf_file a string with the vcf file path
     * \param variant_file a string with the varaint file path
     * */
    //LoadVCF(string vcf_file, string variant_file): VariantFile(variant_file) {
    LoadVCF(string vcf_file, string variant_file) {
      VLOG(9) << "loading vcf file";
      reader.open(vcf_file.c_str(), header);
      reader.readVcfIndex();
    };
    /*! A simple LoadVCF construct
    */
    LoadVCF(): VariantFile() {};

    int get_gene_matrix(vector<vector<string>> geneLoc);
    void load_gene(string gene_name);
    int generate_genotype_matrix(int subjects, int variants);
    arma::Col<int> bernulli(double p, int n);
};

#endif
