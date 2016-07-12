#include <armadillo>
#include <easylogging++.h>
#include <load_vcf.h>
#include <data_analysis.h>
#include <VcfFileReader.h>
#include <string>
#include <vector>
#include <iostream>

using namespace std;
using namespace arma;


void analysis::pedigree(string ped_file, char sep, vector<string> vcfID) {

  ped = variant_location(ped_file, sep); // can use the same function
  vector<int> vcfOrder;
  vector<int> pheno;
  vector<string>::iterator it;
  int num_subjects_ped = ped.size();

  LOG(INFO) << "Number of samples in ped file: " << num_subjects_ped << "\n";
  vector<string> excludeSubjects;
  int i = 0;
  for (auto row = ped.begin(); row != ped.end(); ++row) {
    if (row->at(5) == "-9") {
      excludeSubjects.push_back(row->at(2));
      continue;
    } else {
      it = find(vcfID.begin(), vcfID.end(), row->at(1));
      if (it != vcfID.end()) {
	int pos = distance(vcfID.begin(), it);
        vcfOrder.push_back(pos);
        pheno.push_back(stoi(row->at(5)));
      } else {
        excludeSubjects.push_back(row->at(2));
        LOG(INFO) << "subject " << row->at(1) << " not found in vcf file";
      }
    }
  }
  LOG(INFO) << "Excluding " << excludeSubjects.size()
            << " subjects due to missing phenotype";

  exclude = excludeSubjects;
  id_genotype_include = arma::conv_to< Col<uword> >::from(vcfOrder);
  phenotype.set_size(pheno.size());
  phenotype = conv_to<Col<int> >::from(pheno);
  LOG(INFO) << "Proceeding with " << vcfOrder.size() << " subjects";
}

vector<string> analysis::sampleNames(VcfHeader& header) {
  vector<string> sNames;
  int nSamples = header.getNumSamples();
  for (int i = 0; i < nSamples; ++i) {
    sNames.push_back(header.getSampleName(i));
  }
  ofstream file("sampleNames.txt");
  for(unsigned int i=0; i<sNames.size(); i++)
	  file << sNames[i] << endl;
  file.close();
  return sNames;
}
