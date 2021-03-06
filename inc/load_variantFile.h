#ifndef loadVariantFile_H
#define loadVariantFile_H

#include <string>
#include <vector>
#include <easylogging++.h>
#include <armadillo>

using namespace std;

class VariantFile
{
  public:
    vector<vector<string>> gene_loc; /*!< database of variants in gene*/
    vector<string> genes; /*!< vector of genes in variant file*/
    arma::Mat<int> genotype_matrix; /*!< genotype matrix */
		int num_variants; /*!< number of variants in gene*/

    vector <string> get_all_genes(vector <vector <string> >& data, int r = 2);
    vector <vector <string> > variant_location(string VarFile, char sep);
    vector <vector <string> > get_gene_loc(string gene);

    VariantFile(string variant_file) {
      VLOG(9) << "loading variant file";
      gene_loc = variant_location(variant_file, '\t');
      int gene_col = (gene_loc[0].size() -1);
      genes = get_all_genes(gene_loc, gene_col);
      std::cout << "Found: " << genes.size() << " genes" << std::endl;
    };
    VariantFile() {
      VLOG(9) << "loading variant file but empty";
    };
};

#endif /* loadVariantFile_H */
