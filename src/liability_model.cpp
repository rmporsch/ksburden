#include <armadillo>
#include <random>
#include <liability_model.h>
#include <load_vcf.h>
#include <easylogging++.h>
#include <boost/math/distributions/normal.hpp>

using namespace arma;

vec LiabilityModel::normal_random(int n, double mean, double stdev) {
  vec dat(n);
  normal_distribution<double> distribution(mean, stdev);
  for (int i = 0; i < n; ++i) {
    dat(i) = distribution(rd);
  }
  return dat;
}

vec LiabilityModel::uniform_random(int n) {
  vec dat(n);
  uniform_real_distribution<double> distribution(1, 6);
  for (int i = 0; i < n; ++i) {
    dat(i) = distribution(rd);
  }
  return dat;
}

Col<int> LiabilityModel::generate_causal_variants(bool EmptyStart) {

  Col<int> causal(num_variants, fill::zeros);
  Col<int> cluster(size_cluster);

  int start;
  int emptySpace = num_variants - size_cluster * num_cluster;
  int sizeEmpty = size_cluster-1;

  if (EmptyStart) {
    start = sizeEmpty;
  } else {
    start = 0;
  }

  int i = 0;
  while (i < num_cluster) {
    causal.subvec(start, start + size_cluster - 1) =
        bernulli(causal_probability, size_cluster);

    i++;
    start = start + size_cluster;
    start = start + sizeEmpty;
  }
  VLOG(9) << "simulated " << sum(causal) << " causal mutations"
          << " with probability " << causal_probability;
  return causal;
}

vec LiabilityModel::effect_generation(Col<int> causal) {
  vec effect(num_variants);
  effect.ones();
  effect = causal % effect;
  effect = sqrt(pow(effect, 2) / as_scalar(sum(pow(effect, 2)) / wished_effect));
  //std::cout << effect.subvec(0, 20) << std::endl;
  return effect;
}

uvec LiabilityModel::simulate_data(Col<int> causal) {
  int num_variants = causal.size();
  double q;
  // generate effect size
  vec effectsize = effect_generation(causal);
  vec risk =  genotype_matrix_standarized * effectsize;

  int counter_num_cases = 0;
  int counter_num_controls = 0;
  uvec id_num_cases(num_cases);
  uvec id_num_controls(num_controls);

  boost::math::normal dist(0.0, 1.0);
  q = quantile(dist, 1-life_time_risk);

  while (counter_num_controls < num_controls || counter_num_cases < num_cases) {

    liability_dist =
        risk + normal_random(risk.size(), 0, sqrt(1 - wished_effect));

    mat::iterator start = liability_dist.begin();
    mat::iterator end = liability_dist.end();
    for (vec::iterator i = start; i != end; ++i) {
      if (liability_dist(*i) >= q && counter_num_cases < num_cases) {
        id_num_cases(counter_num_cases) = std::distance(start, i);
        ++counter_num_cases;
      }
      if (liability_dist(*i) < q && counter_num_controls < num_controls) {
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

void LiabilityModel::standardize_matrix()
{
  genotype_matrix_standarized.set_size(genotype_matrix.n_rows,
                                       genotype_matrix.n_cols);
  for (int i = 0; i < genotype_matrix.n_cols; ++i) {
    genotype_matrix_standarized.col(i) = standardize(genotype_matrix.col(i));
  }
}
