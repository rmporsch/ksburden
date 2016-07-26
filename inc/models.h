#ifndef models_H
#define models_H

#include <armadillo>
#include <easylogging++.h>

/*! \brief Models
 *
 * Containing all the different models including a permutation function.
 * Each model can be called by its pointer
 */
class models {
private:
public:
  std::vector<std::string> available_tests; /*!< a vector of strings with available models */

  /*! specifying the typdef for each model
   *
   * This is used as a paramter in the permutation function 
   */
  typedef double (models::*model_members)(const arma::Mat<int> &genotypes,
                                          const arma::Col<int> &phenotpe);

  arma::vec edf(const arma::Mat<int> &genotypes, const arma::uvec &subjects);

  double ksburden(const arma::Mat<int> &genotypes,
                  const arma::Col<int> &phenotype);

  double burden(const arma::Mat<int> &genotypes,
                const arma::Col<int> &phenotype);

  double cmc(const arma::Mat<int> &genotypes, const arma::Col<int> &phenotype);

  double permutation(model_members model, int iteration,
                     const arma::Mat<int> &genotypes, arma::Col<int> phenotype,
                     int max_iteration);
  double fisher(arma::vec pvalues);

  arma::vec large_fisher(arma::mat & pvalues);

  model_members model_array[3] = {NULL}; /*!< array of pointers to each model function */
  /*! Model Constructir
   */
  models() {

    model_array[0] = &models::ksburden;
    model_array[1] = &models::burden;
    model_array[2] = &models::cmc;

    available_tests.resize(3);
    available_tests[0] = "ksburden";
    available_tests[1] = "burden";
    available_tests[2] = "cmc";
  };
};

#endif /* MODELS_H */
