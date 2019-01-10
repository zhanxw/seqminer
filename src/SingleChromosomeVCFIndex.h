#ifndef _SINGLECHROMOSOMEVCFINDEX_H_
#define _SINGLECHROMOSOMEVCFINDEX_H_

#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#include <R.h>
#include <Rinternals.h>

#include <string>
#include "bgzf.h"

class SingleChromosomeVCFIndex {
 public:
  SingleChromosomeVCFIndex(const std::string& vcfFile,
                           const std::string& indexFile);
  virtual ~SingleChromosomeVCFIndex();

  void close();
  // @return 0 for success
  int createIndex();

  // @return 1 if found, 0 not found, -1 if error
  int query(int chromPos, int64_t* voffset);

  // @return number of offsets found, -1 if error
  int query(int chromPosBeg,
            int chromPosEnd,
            int64_t* voffset);
  
  // @return length of @param line if success, or -1 if error
  int readLine(int64_t pos, std::string* line);

  int nextLine(std::string* line);
  
 private:
  std::string vcfFile_;  // must be bgzFile
  std::string indexFile_;
  kstring_t* str;
  BGZF* fVcfFile_;
};

#endif /* _SINGLECHROMOSOMEVCFINDEX_H_ */
