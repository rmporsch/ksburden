#ifndef loadVariantFile_H
#define loadVariantFile_H

#include <string>
#include <vector>

using namespace std;

class VariantFile
{
public:
    vector<vector<string>> gene_loc; /*!< database of variants in gene*/
    vector<string> genes; /*!< vector of genes in variant file*/

  vector <string> get_all_genes(vector <vector <string> >& data);
  vector <vector <string> > variant_location(string VarFile, char sep);
  vector <vector <string> > get_gene_loc(string gene);
  VariantFile(string variant_file) {
      VLOG(9) << "loading variant file";
      gene_loc = variant_location(variant_file, '\t');
      genes = get_all_genes(gene_loc);
  };
  VariantFile() {};
};

#endif /* loadVariantFile_H */
