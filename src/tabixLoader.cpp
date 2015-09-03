#include "tabixLoader.h"

#include "R_CPP_interface.h"
#include "TabixReader.h"

#include <string>
#include <vector>

SEXP impl_readTabixByRange(SEXP arg_tabixFile, SEXP arg_range) {
  SEXP ret = R_NilValue;

  std::vector<std::string> FLAG_tabixFile;
  std::vector<std::string> FLAG_range;
  extractStringArray(arg_tabixFile, &FLAG_tabixFile);
  extractStringArray(arg_range, &FLAG_range);

  if (FLAG_tabixFile.size() != 1) {
    Rprintf("Read the first tabix file: %s\n", FLAG_tabixFile[0].c_str() );
  }
  TabixReader tr(FLAG_tabixFile[0]);
  if (!tr.good()) {
    REprintf("Cannot open specified tabix file: %s\n", FLAG_tabixFile[0].c_str());
    return ret;
  }

  for (size_t i = 0; i < FLAG_range.size(); ++i) {
    int ret = tr.addRange(FLAG_range[i]);
    // REprintf("add range %s , ret = %d\n", FLAG_range[i].c_str(), ret);
  }


  std::string line;
  std::vector<std::string> res;
  while( tr.readLine(&line)) {
    res.push_back(line);
  }

  storeResult(res, &ret);
  UNPROTECT(1);
  return ret;
}


SEXP impl_readTabixSkippedLine(SEXP arg_tabixFile) {
  SEXP ret = R_NilValue;
  
  std::vector<std::string> FLAG_tabixFile;
  extractStringArray(arg_tabixFile, &FLAG_tabixFile);
  TabixReader tr(FLAG_tabixFile[0]);
  if (!tr.good()) {
    REprintf("Cannot open specified tabix file: %s\n", FLAG_tabixFile[0].c_str());
    return ret;
  }

  std::vector<std::string> headers;
  stringTokenize(stringStrip(tr.getSkippedLine()), "\n", &headers);
  
  storeResult(headers, &ret);
  UNPROTECT(1);
  return ret;
}

SEXP impl_readTabixHeader(SEXP arg_tabixFile) {
  SEXP ret = R_NilValue;

  std::vector<std::string> FLAG_tabixFile;
  extractStringArray(arg_tabixFile, &FLAG_tabixFile);
  TabixReader tr(FLAG_tabixFile[0]);
  if (!tr.good()) {
    REprintf("Cannot open specified tabix file: %s\n", FLAG_tabixFile[0].c_str());
    return ret;
  }

  std::vector<std::string> headers;
  stringTokenize(stringStrip(tr.getHeader()), "\n", &headers);
  
  storeResult(headers, &ret);
  UNPROTECT(1);
  return ret;
}

SEXP impl_createTabixIndex(SEXP arg_tabixFile,
                           SEXP arg_seqnameColumn,
                           SEXP arg_startColumn,
                           SEXP arg_endColumn,
                           SEXP arg_commentChar,
                           SEXP arg_skipLine) {
  SEXP ret = R_NilValue;
  std::string fn = CHAR(STRING_ELT(arg_tabixFile, 0));
  int chrom = INTEGER(arg_seqnameColumn)[0];
  int startPos = INTEGER(arg_startColumn)[0];
  int endPos = INTEGER(arg_endColumn)[0];
  char metaChar = CHAR(STRING_ELT(arg_commentChar, 0))[0];
  int skip = INTEGER(arg_skipLine)[0];
  ti_conf_t meta_conf = {0, chrom, startPos, endPos, metaChar, skip};
  if (ti_index_build(fn.c_str(), &meta_conf)) {
    REprintf("Create tabix index failed for [ %s ]!\n", fn.c_str());
  }
  return ret;
}
