#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <easylogging++.h>
#include <vector>
#include <armadillo>
#include <VcfFileReader.h>
#include <load_vcf.h>
#include <load_variantFile.h>

using namespace std;
using namespace arma;


/*! \brief getGene
 *
 * extracts the genotype matrix of a given gene
 *
 * \param geneLoc the result from getGeneLoc()
 *
 */
int LoadVCF::get_gene_matrix(vector<vector<string>> geneLoc) {
  VcfRecord record;
  int num_samples = header.getNumSamples();
  int i,s;
  int invalidGenotypeCount = 0;
  int missingGT = 0;
  int u = 0;
  num_variants = geneLoc.size();
  if (num_variants < 3)
    std::cout << "WARNING: Gene has less than 3 variants" << std::endl;
  genotype_matrix.set_size(num_samples, num_variants);
  genotype_matrix.zeros();

  for(auto i = geneLoc.begin(); i!=geneLoc.end(); ++i) {
    // the last positon is the gene name
    string chrom = i->at(0);
    string pos = i->at(1);
    reader.set1BasedReadSection(chrom.c_str(), (stoi(pos)), (stoi(pos)+1));
    if(reader.isEOF()) throw runtime_error("end of file reached");

    while(reader.readRecord(record)) {
      for (s = 0; s < num_samples; ++s) 
      {
        auto R =  record.getGT(s, 0);
        auto A =  record.getGT(s, 1);
        if (R == VcfGenotypeSample::INVALID_GT or
            A == VcfGenotypeSample::INVALID_GT) {
          LOG(WARNING) << "error in sample ID " << s
            << " ignored genotype at position "
            << chrom << ":" << pos;
          invalidGenotypeCount +=1;
          continue;
        }
        if (R == VcfGenotypeSample::MISSING_GT or
            A == VcfGenotypeSample::MISSING_GT) {
          missingGT +=1;
          continue;
        } else {
          genotype_matrix(s, u) = (R + A);
        }
      }
    }
    if (missingGT/arma::as_scalar(num_samples)  > 0.1) {
      LOG(WARNING) << "missingness of " << chrom << ":" << pos
        << " is higher that 10%";
    }
    u += 1;
  }
  if (invalidGenotypeCount > 0) {
    LOG(WARNING) << "ignored " << invalidGenotypeCount
      << " genotypes due to invalid coding";
  }

  return 0;
}

/*! \brief loads the genotype matrix
 *
 * \param gene_name a string with a gene name
 *
 * its a simple warpper function for get_gene_loc and get_gene_matrix
 */
void LoadVCF::load_gene(string gene_name)
{
	auto gene_loc = get_gene_loc(gene_name);
	get_gene_matrix(gene_loc);
}


/*! \brief generates a genotype matrix
 *
 * \param subjects number of subjects to simulate
 * \param variants number of variants
 */
int LoadVCF::generate_genotype_matrix(int subjects, int variants) {
  int i;
  double variantMAF;
  num_variants = variants;

  genotype_matrix.set_size(subjects, variants);
  arma::vec MAF = {0.0001, 0.001, 0.01};
  arma::vec FreqMAFtype = {0.9, 0.09, 0.01};
  string MatrixFileName = "genotype_matrixtypeMatrix" + std::to_string(subjects) +
                          std::to_string(variants) + ".mat";
  string MatrixFileNameStand = "genotype_matrixtypeMatrixSTAND" +
                               std::to_string(subjects) +
                               std::to_string(variants) + ".mat";

  for (i = 0; i < variants; ++i) {
    double val = (double)rand() / RAND_MAX;
    if (val < 0.01) {
      variantMAF = 0.01;
    } else if (val < 0.09) {
      variantMAF = 0.001;
    } else if (val < 0.9) {
      variantMAF = 0.0001;
    }

    genotype_matrix.col(i) = conv_to<Col<int>>::from(bernulli(variantMAF, subjects));
  }
  genotype_matrix.save(MatrixFileName);
  return 0;
}

/*! \brief generates a bernulli distribution
 *
 * \param p probability
 * \param n size of the distribution
 */
Col<int> LoadVCF::bernulli(double p, int n) {
  Col<int> dat(n);
  bernoulli_distribution distribution(p);
  for (int i = 0; i < n; ++i) {
    dat(i) = distribution(rd);
  }
  return dat;
}
