#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <easylogging++.h>
#include <load_variantFile.h>

using namespace std;

/*! \brief GetAllGenes
 * extracts all genes or pathways mentioned in data
 *
 * \param data the result from variantLocation()
 * \returns a vector of gene names
 */
vector <string> VariantFile::get_all_genes(vector <vector <string> >& data, int r)
{
	vector <string> genes;
	for(auto row = data.begin(); row!=data.end(); ++row) 
	{
		if(find(genes.begin(), genes.end(), row->at(r)) !=genes.end() )
		{
			continue;
		} else {
			genes.push_back(row->at(r));
		}
	}
	return genes;
}

/*! \brief Function to read the varinat location file
 *
 * \param VarFile File of containing the variants
 * \param sep seperator of VarFile default is '\t'
 * \returns a vector of vector<string> containing the content of the file
 *
 * Imput file requires that the first two columns contain the chr and pos.
 * The third column contains the gene name or pathway
 */
vector <vector <string> >VariantFile::variant_location(std::string VarFile, char sep)
{
	vector <vector <string> > data;
	ifstream infile(VarFile);

	while(infile) {
		string s;
		if (!getline(infile, s)) break;

		istringstream ss(s);
		vector <string> record;

		while (ss) {
			string s;
			if (!getline(ss, s, sep)) break;
			record.push_back(s);
		}
		data.push_back(record);
	}
	if (!infile.eof()) {
		cerr <<"Foey!\n";
	}

	return data;
}

/*! \brief getGeneLoc
 *
 * extracts location of a given gene
 *
 * \param gene the gene name
 *
 */
vector <vector <string> >VariantFile::get_gene_loc(string gene)
{
	vector <vector <string> > positionGene;
	for(auto i = gene_loc.begin(); i!=gene_loc.end(); ++i) {
		if(find(i->begin(), i->end(), gene) !=i->end() )
		{
			positionGene.push_back(*i);
		}
	}
	if(positionGene.size() == 0) throw std::runtime_error("Gene not in file");
	return positionGene;
}
