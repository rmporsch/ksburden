#include <models.h>

using namespace arma;

/**
 * Function to compute cumulative distribution function on genotypes
 */
vec models::edf(const Mat<int> &genotypes, const uvec &subjects) {
  Row<int> variantSum(genotypes.n_cols, fill::zeros);
  variantSum = sum(genotypes.rows(subjects), 0);

  Row<int> ecdf = cumsum(variantSum);
  int AllVariants = sum(variantSum);
  vec out = conv_to<vec>::from(ecdf) / as_scalar(AllVariants);
  return out;
}

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

double models::burden(const Mat<int> &genotypes, const Col<int> &phenotype)
{
  double test_statistic = pow(sum(phenotype.t() * genotypes),2);
  return test_statistic;
}

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

