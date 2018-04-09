#include <models.h>
#include <boost/math/distributions/chi_squared.hpp>

using namespace arma;

/*! \brief Function to compute cumulative distribution function on genotypes
 *
 * \param genotypes a genotype matrix
 * \param subjects vector of row IDs of genotype matrix
 * \returns a vector of the ecdf
 */
vec models::edf(const Mat<int> &genotypes, const uvec &subjects) {
  Row<int> variantSum(genotypes.n_cols, fill::zeros);
  variantSum = sum(genotypes.rows(subjects), 0);

  Row<int> ecdf = cumsum(variantSum);
  int AllVariants = sum(variantSum);
  vec out = conv_to<vec>::from(ecdf) / as_scalar(AllVariants);
  return out;
}

/*! \brief ks-test
 *
 * \param genotypes a genotype matrix
 * \param phenotype the phenotype
 * \returns the ks test statistics
 */
double models::ksburden(const Mat<int> &genotypes, const Col<int> &phenotype)
{
  double test_statistic;
  uvec cases_ID = find(phenotype == 1);
  uvec control_ID = find(phenotype == -1);
  vec cases = edf(genotypes, cases_ID);
  vec controls = edf(genotypes, control_ID);

  vec diff_ksburden = cases - controls;
  test_statistic = max(abs(diff_ksburden));

  return test_statistic;
}

/*! \brief burden-test
 *
 * \param genotypes a genotype matrix
 * \param phenotype the phenotype
 * \returns the burden test statistics
 */
double models::burden(const Mat<int> &genotypes, const Col<int> &phenotype)
{
  double test_statistic = pow(sum(phenotype.t() * genotypes),2);
  return test_statistic;
}

/*! \brief cmc-test
 *
 * \param genotypes a genotype matrix
 * \param phenotype the phenotype
 * \returns the cmc test statistics
 */
double models::cmc(const Mat<int> &genotypes, const Col<int> &phenotype)
{
  double test_statistic;

  uvec cases_ID = find(phenotype > 0);
  uvec control_ID = find(phenotype < 0);
  int number_cases = cases_ID.size();
  int number_controls = control_ID.size();

  Col<int> test_vec = arma::sign(phenotype % sum(genotypes, 1));
  uvec larger_ID = find(test_vec > 0);
  uvec smaller_ID = find(test_vec < 0);

  test_statistic =
      ((double)larger_ID.size() / number_cases) - ((double)smaller_ID.size() / number_controls);
  test_statistic = pow(test_statistic, 2);
  return test_statistic;
}

/*! computes p-values given a model
 *
 * \param model a pointer to the model function
 * \param iteration number of iteration
 * \param genotypes a genotype matrix
 * \param phenotype the phenotype
 * \param max_iteration the amximal number of iterations
 *
 * \return a p-value
 */
double models::permutation(model_members model, int iteration,
                           const Mat<int> &genotypes, Col<int> phenotype,
                           int max_iteration) {
  double p_value;
  int larger = 0; int i=0;
  Col<int> temp_phenotype(phenotype.size());

  double test_statistic = (this->*model)(genotypes, phenotype);

  while (i <= iteration || larger < 1) {
    temp_phenotype = shuffle(phenotype);
    double iteration_test_statistic = (this->*model)(genotypes, temp_phenotype);
    if (iteration_test_statistic >= test_statistic)
      larger += 1;
    i += 1;
    if ( i >= max_iteration) break;
  }

  p_value = (double)(larger+1) / (i+1);
  return p_value;
}

/*! \brief Fisher's p-value combination
 *
 * This is a function which compues Fisher's p-value combination
 * on the assumption the tests are independent.
 *
 * \param pvalues a vector of p-values
 * \returns a combined p-value
 */
double models::fisher(arma::vec pvalues)
{
  arma::vec container(pvalues.size(), arma::fill::zeros);
  int i;
  for (i = 0; i < pvalues.size(); ++i) {
    container(i) = log(pvalues(i));
  }
  double chi = -2 * arma::sum(container);
  boost::math::chi_squared_distribution<double> chi_dist(2 * pvalues.size());
  double p = 1 - boost::math::cdf(chi_dist, chi);
  return p;
}

/*! \brief computes fisher's method on a matrix of pvalues
 *
 * \param a aramdillo matrix with pvalues
 * \return a vector of pvalues
 */
arma::vec models::large_fisher(arma::mat & pvalues)
{
  arma::vec coutput;
  coutput.set_size(pvalues.n_rows);
  int i;
  for (i = 0; i < pvalues.n_rows; ++i) {
	  coutput.at(i) = fisher(pvalues.row(i).t());
  }
  return coutput;
}


/*! \brief kssum-test
 *
 * \param genotypes a genotype matrix
 * \param phenotype the phenotype
 * \returns the kssum test statistics
 */
double models::kssum(const Mat<int> &genotypes, const Col<int> &phenotype)
{
  double test_statistic;
  uvec cases_ID = find(phenotype == 1);
  uvec control_ID = find(phenotype == -1);
  vec cases = edf(genotypes, cases_ID);
  vec controls = edf(genotypes, control_ID);

  vec diff_ksburden = cases - controls;
  test_statistic = sum(abs(diff_ksburden));

  return test_statistic;
}
