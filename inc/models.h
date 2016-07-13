#ifndef models_H
#define models_H

#include <armadillo>
#include <easylogging++.h>

class models {
private:
public:
  std::vector<std::string> available_tests;
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

  model_members model_array[3] = {NULL};
  models() {
    model_array[0] = &models::ksburden;
    model_array[1] = &models::burden;
    model_array[2] = &models::cmc;

    available_tests[0] = "ksburden";
    available_tests[1] = "burden";
    available_tests[2] = "cmc";
  };
};

#endif /* MODELS_H */
