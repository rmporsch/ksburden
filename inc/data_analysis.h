#ifndef data_analysis_H
#define data_analysis_H

#include <armadillo>
#include <load_vcf.h>
#include <vector>
#include <string>

/*! \brief performs basic data analysis
 *
 * This class was designed to perfom the standard task of
 * analysing a sample of data given a vcf file, variant file and
 * pedigree file containing the case-control status
 */
class analysis : public LoadVCF {

public:
  arma::uvec id_genotype_include; /*!< IDs of the genotype matrix to include in the analysis */
  arma::Col<int> phenotype; /*!< imported phenotype from the ped file */
  std::vector <std::string> exclude; /*!< subject names to exclude */
  std::vector<std::vector<std::string>> ped; /*!< contant of the pedigree file*/
  std::string ped_file;
  std::string vcfFile;
  std::string variant_file;
  int num_subjects; /*!< number of subjects based on ped file */
  std::vector<std::string> name_subjects; /*!< name of subjects from ped file */

  void pedigree(std::string ped_file, char sep, std::vector<std::string> vcfID);
  std::vector<std::string> sampleNames(VcfHeader& header);

  /*! \brief Constructor of analysis
   *
   * \param ped_file path to the ped file
   * \param vcfFile path to the vcf file
   * \param variant_file path to the variant file
   */
  analysis(std::string ped_file, std::string vcfFile, std::string variant_file)
      : LoadVCF(vcfFile, variant_file) {
    VLOG(9) << "load subjects";
    name_subjects = sampleNames(header);
    pedigree(ped_file, '\t', name_subjects);
    num_subjects = ped.size();
  };
  /*! \brief Simple constructor of analysis
   */
  analysis(){};
};
#endif

