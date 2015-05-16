#include <string>
#include <set>
#include <vector>
#include <map>

// OrderedMap is template, cannot forward delcaration
#include "OrderedMap.h"

void loadGeneFile(const std::string& geneFile,
                  const std::vector<std::string>& geneName,
                  std::vector<std::string>* geneRange);

void loadGeneFile(const std::string& geneFile, const std::vector<std::string>& geneName, OrderedMap< std::string, std::string>* geneRangePtr);
void loadGeneFile(const std::string& geneFile, const std::set<std::string>& geneName, OrderedMap< std::string, std::string>* geneRangePtr);
void loadGeneFile(const std::string& geneFile, const std::string& geneName, OrderedMap< std::string, std::string>* geneRangePtr);

#if 0

void loadGeneFile(const std::string& geneFile, const std::set<std::string>& geneName, std::map< std::string, std::string>* geneRangePtr);
void loadGeneFile(const std::string& geneFile, const std::string& geneName, std::map< std::string, std::string>* geneRangePtr);
void loadGeneFile(const std::string& geneFile, const std::vector<std::string>& geneName, std::map< std::string, std::string>* geneRangePtr);

#endif
