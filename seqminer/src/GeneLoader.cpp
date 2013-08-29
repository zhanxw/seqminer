#include "GeneLoader.h"

#include <string>
#include <map>
#include <vector>
#include <set>

#include "R.h"

#include "IO.h"
#include "Utils.h"

/**
 * @param geneRangePtr: store results here, key is gene name, value is comma-separated range
 */
void loadGeneFile(const std::string& geneFile, const std::set<std::string>& geneName, std::map< std::string, std::string>* geneRangePtr){
  std::map< std::string, std::string>& geneRange = *geneRangePtr;
  geneRange.clear();

  // load gene ranges
  if (!geneName.empty() ){
    if (geneFile.size() == 0) {
      //REprintf("Have to provide --geneFile to extract by gene.\n");
      //REprintf("Critical error happening!"); //abort();
      error("Have to provide --geneFile to extract by gene.\n");
    }
    LineReader lr(geneFile);
    std::vector< std::string > fd;
    while (lr.readLineBySep(&fd, "\t ")) {
      if (!geneName.count(fd[0]) ) continue;
      fd[2] = chopChr(fd[2]); // chop "chr1" to "1"
      if (geneRange.find(fd[0])  == geneRange.end()) {
        geneRange[fd[0]] = fd[2] + ":" + fd[4] + "-" + fd[5];
      } else {
        geneRange[fd[0]] += "," + fd[2] + ":" + fd[4] + "-" + fd[5];
      }
    };
  }
}
void loadGeneFile(const std::string& geneFile, const std::string& geneName, std::map< std::string, std::string>* geneRangePtr){
  if (geneName.empty()) return;
  std::set<std::string> g;
  g.insert(geneName);
  loadGeneFile(geneFile, g, geneRangePtr);
}
void loadGeneFile(const std::string& geneFile, const std::vector<std::string>& geneName, std::map< std::string, std::string>* geneRangePtr){
  if (geneName.empty()) return;
  std::set<std::string> g;
  for (size_t i = 0; i < geneName.size(); ++i) {
    if (geneName[i].empty()) continue;
    g.insert(geneName[i]);
  }
  loadGeneFile(geneFile, g, geneRangePtr);
}

//////////////////////////////////////////////////
/**
 * @param geneRangePtr: store results here, key is gene name, value is comma-separated range
 */
void loadGeneFile(const std::string& geneFile, const std::set<std::string>& geneName, OrderedMap< std::string, std::string>* geneRangePtr){
  OrderedMap< std::string, std::string>& geneRange = *geneRangePtr;
  geneRange.clear();

  // load gene ranges
  if (!geneName.empty() ){
    if (geneFile.size() == 0) {
      //REprintf("Have to provide --geneFile to extract by gene.\n");
      //REprintf("Critical error happening!"); //abort();
      error("Have to provide --geneFile to extract by gene.\n");
    }
    LineReader lr(geneFile);
    std::vector< std::string > fd;
    while (lr.readLineBySep(&fd, "\t ")) {
      if (!geneName.count(fd[0]) ) continue;
      fd[2] = chopChr(fd[2]); // chop "chr1" to "1"
      if (!geneRange.find(fd[0])) {
        geneRange[fd[0]] = fd[2] + ":" + fd[4] + "-" + fd[5];
      } else {
        geneRange[fd[0]] += "," + fd[2] + ":" + fd[4] + "-" + fd[5];
      }
    };
  }
}
void loadGeneFile(const std::string& geneFile, const std::string& geneName, OrderedMap< std::string, std::string>* geneRangePtr){
  if (geneName.empty()) return;
  std::set<std::string> g;
  g.insert(geneName);
  loadGeneFile(geneFile, g, geneRangePtr);
}
void loadGeneFile(const std::string& geneFile, const std::vector<std::string>& geneName, OrderedMap< std::string, std::string>* geneRangePtr){
  if (geneName.empty()) return;
  std::set<std::string> g;
  for (size_t i = 0; i < geneName.size(); ++i) {
    if (geneName[i].empty()) continue;
    g.insert(geneName[i]);
  }
  loadGeneFile(geneFile, g, geneRangePtr);
}
