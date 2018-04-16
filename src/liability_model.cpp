#include <armadillo>
#include <random>
#include <liability_model.h>
#include <load_vcf.h>
#include <easylogging++.h>
#include <boost/math/distributions/normal.hpp>
#include <stdlib.h>
#include <stdio.h>


using namespace arma;

/*! \brief normal distribution
 *
 * \param n size of the distribition
 * \param mean the mean
 * \param stdev the standard deviation
 */
vec LiabilityModel::normal_random(int n, double mean, double stdev) {
  vec dat(n);
  normal_distribution<double> distribution(mean, stdev);
  for (int i = 0; i < n; ++i) {
    dat(i) = distribution(rd);
  }
  return dat;
}

/*! \brief unifrom distribution
 *
 * \param n size of the distribition
 */
vec LiabilityModel::uniform_random(int n) {
  vec dat(n);
  uniform_real_distribution<double> distribution(1, 6);
  for (int i = 0; i < n; ++i) {
    dat(i) = distribution(rd);
  }
  return dat;
}

/*! \brief declares the causal varinats
 *
 * \param EmptyStart should there be no causal variants at the beginning of the gene
 * \return a vector of integers with either 1 or 0
 */
Col<int> LiabilityModel::generate_causal_variants(bool EmptyStart) {

  Col<int> causal(num_variants, fill::zeros);
  int sides = (size_cluster / 2) / num_cluster;
  std::cout << "cluster size: " << size_cluster << std::endl;
  int gamp_min = 1;
  int pos;
  std::uniform_int_distribution<> dis(0+(sides+gamp_min),
      (num_variants-(1+sides+gamp_min)));

  Col<int> c_vec(2*sides+1, fill::ones);
  real_num_causal = 2*sides+1;

  if (EmptyStart) {
    pos = 0;
    causal.ones();
  } else {
    // for each causal cluster
    for (int i = 0; i != num_cluster; ++i) {
      int counter = 0;
      int ff;
      do {
        // choose random position
        pos = dis(rd);
        // check if there is any overlap
        ff = arma::sum(causal.subvec(pos-(sides+gamp_min), pos+(sides+gamp_min)));
        counter += 1;
      } while (ff > 0 & counter <= 100);
      // asign causal cluster
      if (counter <=100) {
        causal.subvec(pos-(sides), pos+(sides)) = c_vec;
        std::cout << "current number of causal mutations" << std::endl;
        std::cout << arma::sum(causal) << std::endl;
      } else {
        printf ("Error: could not find suitable location for causal cluster");
        exit (EXIT_FAILURE);
      }
    }
  }

  VLOG(9) << "simulated " << sum(causal) << " causal mutations"
          << " with probability " << causal_probability;
  return causal;
}

/*! \brief makes effect size
 *
 * \param causal vector of causal mutations
 * \return effect
 *
 * This function standardizes the effect size of all causal variants 
 * to be unform across the causal mutations. Hence the indicated effect size from wished_effect
 * indicates the $R^2$ of the intire gene.
 */
vec LiabilityModel::effect_generation(Col<int> causal) {
  vec effect(num_variants);
  effect.ones();
  effect = causal % effect;
  effect = sqrt(pow(effect, 2) / as_scalar(sum(pow(effect, 2)) / wished_effect));
  return effect;
}

/*! \brief generates cases and controls
 *
 * \param causal a vector of causal mutations
 * \returns a vector of IDs from the original genotype matrix
 *
 * This is the work horse of the class. It generates a liability distribution
 * based on the effect size. Case-control status is then decided on the life_time_risk 
 * This process is repeated until enough cases and controls have been collected.
 */
uvec LiabilityModel::simulate_data(Col<int> causal) {
  // generate effect size
  vec effectsize = effect_generation(causal);
  vec risk =  genotype_matrix_standarized * effectsize;

  int counter_num_cases = 0;
  int counter_num_controls = 0;
  uvec id_num_cases(num_cases);
  uvec id_num_controls(num_controls);

  boost::math::normal dist(0.0, 1.0);
  double q = quantile(dist, 1-life_time_risk);

  while (counter_num_controls < num_controls || counter_num_cases < num_cases) {

    liability_dist =
        risk + normal_random(risk.size(), 0, sqrt(1 - wished_effect));

    mat::iterator start = liability_dist.begin();
    mat::iterator end = liability_dist.end();

    for (vec::iterator i = start; i != end; ++i) {
      if (*i >= q && counter_num_cases < num_cases) {
        id_num_cases(counter_num_cases) = std::distance(start, i);
        ++counter_num_cases;
      }
      if (*i < q && counter_num_controls < num_controls) {
        id_num_controls(counter_num_controls) = std::distance(start, i);
        ++counter_num_controls;
      }
    }
  }
  uvec out(num_cases + num_controls);
  out.subvec(0, num_cases-1) = id_num_cases;
  out.subvec(num_cases, num_cases+num_controls-1) = id_num_controls;
  return out;
}

/*! \brief simple standardization function
 * 
 * \param variant a vector of binary varinats
 * \returns a standardized vector of variants
 *
 * If the input vector contains only 0 the output vector will as well
 */
colvec LiabilityModel::standardize(Col<int> variant) {
  colvec genomean(variant.n_elem);
  genomean.fill(mean(conv_to<Col<double>>::from(variant)));

  double genosd = stddev(conv_to<Col<double>>::from(variant));
  if (genosd == 0) {
    colvec out(variant.n_elem, fill::zeros);
    return out;
  } else {
    colvec out = (variant - genomean) / as_scalar(genosd);
    return out;
  }
}

/*! \brief standardization of the whole matrix
 */
void LiabilityModel::standardize_matrix()
{
  genotype_matrix_standarized.set_size(genotype_matrix.n_rows,
                                       genotype_matrix.n_cols);
  for (int i = 0; i < genotype_matrix.n_cols; ++i) {
    genotype_matrix_standarized.col(i) = standardize(genotype_matrix.col(i));
  }
}
