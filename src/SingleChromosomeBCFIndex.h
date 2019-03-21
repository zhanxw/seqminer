#ifndef _SINGLECHROMOSOMEBCFINDEX_H_
#define _SINGLECHROMOSOMEBCFINDEX_H_

#include <string>
#include <vector>
#include "bgzf.h"

class SingleChromosomeBCFIndex {
 public:
  SingleChromosomeBCFIndex(const std::string& bcfFile,
                           const std::string& indexFile);
  virtual ~SingleChromosomeBCFIndex();

  void close();
  void closeIndex();

  // @return 0 for success
  int createIndex();

  // open index
  int openIndex();

  // @return 1 if found, 0 not found, -1 if error
  int query(int chromPos, int64_t* voffset);

  // @return number of offsets found, -1 if error
  int query(int chromPosBeg, int chromPosEnd, int64_t* voffset);

  // @return length of @param line if success, or -1 if error
  int readLine(int64_t pos,
               uint32_t* l_shared,
               uint32_t* l_indiv,
               std::vector<char>* line);

  int nextLine(uint32_t* l_shared,
               uint32_t* l_indiv,
               std::vector<char>* line);

 private:
  std::string bcfFile_;  // must be bgzFile
  std::string indexFile_;
  void* data_;  // store indices
  BGZF* fBcfFile_;
};

#endif /* _SINGLECHROMOSOMEBCFINDEX_H_ */
