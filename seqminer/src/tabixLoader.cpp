#include "tabixLoader.h"

#include "R_CPP_interface.h"
#include "TabixReader.h"

#include <string>
#include <vector>

SEXP impl_readTabixByRange(SEXP arg_tabixFile, SEXP arg_range) {
  std::vector<std::string> FLAG_tabixFile;
  std::vector<std::string> FLAG_range;
  extractStringArray(arg_tabixFile, &FLAG_tabixFile);
  extractStringArray(arg_range, &FLAG_range);

  if (FLAG_tabixFile.size() != 1) {
    Rprintf("Read the first tabix file: %s\n", FLAG_tabixFile[0].c_str() );
  }
  TabixReader tr(FLAG_tabixFile[0]);
  for (size_t i = 0; i < FLAG_range.size(); ++i) {
    tr.addRange(FLAG_range[i]);  
  }

  std::string line;
  std::vector<std::string> res;
  while( tr.readLine(&line)) {
    res.push_back(line);
  }

  SEXP ret = R_NilValue;
  storeResult(res, &ret);
  UNPROTECT(1);
  return ret;
}
