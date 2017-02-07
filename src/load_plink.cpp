#include <stdio.h>
#include <string>
#include <armadillo>
#include <easylogging++.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <bitset>
#include <load_plink.h>
#include <load_vcf.h>

using namespace arma;

bool LoadPlink::openPlinkBinaryFile(const std::string s, std::ifstream &BIT) {
  BIT.open(s.c_str(), std::ios::in | std::ios::binary);
  if (!BIT.is_open()) {
    std::cerr << "Cannot open the bed file";
  }

  // 2) else check for 0.99 SNP/Ind coding
  // 3) else print warning that file is too old
  char ch[1];
  BIT.read(ch, 1);
  std::bitset<8> b;
  b = ch[0];
  bool bfile_SNP_major = false;
  bool v1_bfile = true;
  // If v1.00 file format
  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
  // std::std::cerr << "check magic number" << std::endl;
  if ((b[2] && b[3] && b[5] && b[6]) && !(b[0] || b[1] || b[4] || b[7])) {
    // Next number
    BIT.read(ch, 1);
    b = ch[0];
    if ((b[0] && b[1] && b[3] && b[4]) && !(b[2] || b[5] || b[6] || b[7])) {
      // Read SNP/Ind major coding
      BIT.read(ch, 1);
      b = ch[0];
      if (b[0])
        bfile_SNP_major = true;
      else
        bfile_SNP_major = false;

      // if (bfile_SNP_major) std::std::cerr << "Detected that binary PED file is
      // v1.00 SNP-major mode" << std::endl;
      // else std::std::cerr << "Detected that binary PED file is v1.00
      // individual-major mode" << std::endl;

    } else
      v1_bfile = false;

  } else
    v1_bfile = false;
  // Reset file if < v1
  if (!v1_bfile) {
    std::cerr << "Warning, old BED file <v1.00 : will try to recover..."
      << std::endl;
   std::cerr << "  but you should --make-bed from PED )" << std::endl;
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
    BIT.read(ch, 1);
    b = ch[0];
  }
  // If 0.99 file format
  if ((!v1_bfile) && (b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7])) {
    std::cerr << std::endl
      << " *** Possible problem: guessing that BED is < v0.99      *** "
      << std::endl;
    std::cerr << " *** High chance of data corruption, spurious results    *** "
      << std::endl;
    std::cerr
      << " *** Unless you are _sure_ this really is an old BED file *** "
      << std::endl;
    std::cerr << " *** you should recreate PED -> BED                      *** "
      << std::endl
      << std::endl;
    bfile_SNP_major = false;
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
  } else if (!v1_bfile) {
    if (b[0])
      bfile_SNP_major = true;
    else
      bfile_SNP_major = false;
    std::cerr << "Binary PED file is v0.99" << std::endl;
    if (bfile_SNP_major)
      std::cerr << "Detected that binary PED file is in SNP-major mode"
        << std::endl;
    else
      std::cerr << "Detected that binary PED file is in individual-major mode"
        << std::endl;
  }
  return bfile_SNP_major;
}

int LoadPlink::countlines(const char* fileName) {
  // Stolen from http://stackoverflow.com/questions/3482064/counting-the-number-of-lines-in-a-text-file
  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(fileName);

  while (std::getline(myfile, line))
    ++number_of_lines;
  return number_of_lines;
}

arma::mat LoadPlink::get_genotype_matrix(const std::string fileName,
    arma::Col<int> col_skip,
    arma::Col<int> row_skip) {

  std::ifstream bedFile;
  bool snpMajor = openPlinkBinaryFile(fileName, bedFile);

  if (!snpMajor)
    throw std::runtime_error("We currently have no plans of implementing the "
        "individual-major mode. Please use the snp-major "
        "format");

  int i = 0;
  int ii = 0;
  int num_snps = arma::accu(col_skip);
  const bool colskip = (num_snps != P);

  int Nbytes = ceil(N / 4.0);
  int n = arma::accu(row_skip);
  const bool selectrow = (n != N);

  int j, jj;

  arma::mat genotypes = arma::mat(n, num_snps, arma::fill::zeros);
  std::bitset<8> b; // Initiate the bit array
  char ch[Nbytes];

  i=0;
  int curr_snp = 0;
  while (i < P && i < num_snps) {
    if (col_skip[i] == 1) {
      bedFile.seekg(i * Nbytes, bedFile.beg);
      bedFile.read(ch, Nbytes); // Read the information
      if (!bedFile)
        throw std::runtime_error(
            "Problem with the BED file...has the FAM/BIM file been changed?");

      j = 0; 
      for (jj = 0; jj < Nbytes; jj++) {
        b = ch[jj];
        int c = 0;
        while (c < 7 && j < N && j < n) { //from the original PLINK: 7 because of 8 bits
          if (row_skip[j] == 1) {
            int first = b[c++];
            int second = b[c++];
            if (first == 0) {
              genotypes(j,curr_snp) = (2 - second);
            } // since genotypeMatrix is filled with 1
            if (first == 1 && second == 0) {
              genotypes(j, curr_snp) = 0;
            } // in case of missing values
          }
          j++;
        }
      }
      curr_snp++;
    }
    i++;
  }
  return genotypes;
}

/*! \brief generates vector of variants to include 
 *
 * This functions takes the gene location with rsids and
 * looks those up in the bed file. It then returns a vector of
 * size p (size of all snps) with 1 for variant to include
 * and 0 for variant not to include
 *
 * \param gene_loc the gene location vector
 *
 * \return Col<int> with 1 and 0
 */
arma::Col<int> LoadPlink::get_col_skip(vector<vector<string>> gene_loc)
{
  // bim file
  arma::Col<int> col_skip(bim_file.size(), arma::fill::zeros);

  int i = 0;
  for (auto g = gene_loc.begin(); g != gene_loc.end(); g++ ) {
    while (i < P) {
      if(find(bim_file[i].begin(), bim_file[i].end(), g->at(0)) !=bim_file[i].end() )
      {
        col_skip(i) = 1;
        break;
      }
      i++;
    }
  }
  if(arma::accu(col_skip) == 0)
    throw std::runtime_error("Gene or positions not in file");

  return col_skip;
}
