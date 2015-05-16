#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#include <R.h>
#include <Rinternals.h>

extern "C" SEXP impl_readTabixByRange(SEXP arg_tabixFile, SEXP arg_range);
extern "C" SEXP impl_readTabixHeader(SEXP arg_tabixFile);
extern "C" SEXP impl_readTabixSkippedLine(SEXP arg_tabixFile);
extern "C" SEXP impl_createTabixIndex(SEXP arg_tabixFile,
                                      SEXP arg_seqnameColumn,
                                      SEXP arg_startColumn,
                                      SEXP arg_endColumn,
                                      SEXP arg_commentChar,
                                      SEXP arg_skipLine);

