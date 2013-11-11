#include "GeneLoader.h"

#include <string>
#include <map>
#include <vector>
#include <set>

#include "R.h"

#include "IO.h"
#include "Utils.h"
#include "OrderedMap.h"
#include "RangeList.h"

/**
 * Load @param geneFile (assuming refFlat format)
 * for those exists in @param geneName (allow duplicated gene names)
 * store its range (after merging multiple transcripts) to @param geneRange
 * e.g.
 * geneName = {A, B, C, A} and
 * range for A is rangeA1, rangeA2
 * range for B is rangeB
 * range for C is rangeC
 * geneRange = {rangeForA1:rangeA2, rangeForB, rangeForC, rangeForA1:rangeA2}
 */
void loadGeneFile(const std::string& geneFile,
                  const std::vector<std::string>& geneName,
                  std::vector<std::string>* geneRange) {
  if (!geneRange) return;
  
  int n = geneName.size();
  std::map<std::string, std::vector<int> > data;
  std::vector<RangeList> ranges(n);

  // store gene <-> idx
  for (int i = 0; i < n; ++i) {
    data[geneName[i]].push_back(i);
  }

  // store ranges for each gene
  LineReader lr(geneFile);
  std::vector< std::string > fd;
  while (lr.readLineBySep(&fd, "\t ")) {
    if (!data.count(fd[0]) ) continue;
    fd[2] = chopChr(fd[2]); // chop "chr1" to "1"

    const int firstIndex = data[fd[0]][0];
    RangeList& rl = ranges[firstIndex];

    int beg = atoi(fd[4]);
    int end = atoi(fd[5]);
    if (beg <= 0 || end <= 0) {
      continue;
    }
    rl.addRange(fd[2], beg, end);
  }

  // store results
  std::vector<std::string>& output = *geneRange;
  output.resize(geneName.size());
  std::map<std::string, std::vector<int> >::iterator it;
  for(it = data.begin(); it != data.end(); ++it){
    const std::string& name = it->first;
    const std::vector<int>& indice = it->second;
    int n = indice.size();
    ranges[indice[0]].sort();
    const std::string& r = ranges[indice[0]].toString();
    for (int i = 0; i < n; ++i) {
      (output)[indice[i]] = r;  
    }
  }

  // debug outputs
  assert(geneName.size() == output.size());
  for (size_t i = 0; i < geneName.size(); ++i) {
    REprintf("range of [ %s ] is [ %s ]\n", geneName[i].c_str(), output[i].c_str());
  }
}

static void setToVector(const std::set<std::string>& input,
                        std::vector<std::string>* output) {
  int n = input.size();
  std::set<std::string>::const_iterator iter;
  for (iter = input.begin(); iter != input.end(); ++iter){
    output->push_back(*iter);
  }
}

void loadGeneFile(const std::string& geneFile,
                  const std::vector<std::string>& geneName,
                  OrderedMap< std::string, std::string>* geneRangePtr){
  std::vector<std::string> v;
  loadGeneFile(geneFile, geneName, &v);
  int n = v.size();
  for (int i = 0; i < n; ++i) {
    (*geneRangePtr)[geneName[i]] = v[i];
  }
}

void loadGeneFile(const std::string& geneFile,
                  const std::set<std::string>& geneName,
                  OrderedMap< std::string, std::string>* geneRangePtr){
  
  std::vector<std::string> v;
  setToVector(geneName, &v);
  loadGeneFile(geneFile, v, geneRangePtr);
}

void loadGeneFile(const std::string& geneFile,
                  const std::string& geneName,
                  OrderedMap< std::string, std::string>* geneRangePtr){
  
  std::vector<std::string> v;
  v.push_back(geneName);
  loadGeneFile(geneFile, v, geneRangePtr);
}

#if 0
/**
 * @param geneRangePtr: store results here, key is gene name, value is comma-separated range
 */
void loadGeneFile(const std::string& geneFile,
                  const std::set<std::string>& geneName,
                  OrderedMap< std::string, std::string>* geneRangePtr){
  std::vector<std::string> v;
  std::vector<std::string> v;
  setToVector(geneName, &v);
  
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

/**
 * 
 */
void loadGeneFile(const std::string& geneFile, const std::vector<std::string>& geneName, std::map< std::string, std::string>* geneRangePtr){
  if (geneName.empty()) return;
  std::set<std::string> g;
  for (size_t i = 0; i < geneName.size(); ++i) {
    if (geneName[i].empty()) continue;
    g.insert(geneName[i]);
  }
  loadGeneFile(geneFile, g, geneRangePtr);
}

#endif
