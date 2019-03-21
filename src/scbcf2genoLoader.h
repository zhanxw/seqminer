#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#include <R.h>
#include <Rinternals.h>

extern "C" SEXP impl_readSingleChromosomeBCFToMatrixByRange(SEXP fileName, SEXP indexFile, SEXP range);
extern "C" SEXP impl_createSingleChromosomeBCFIndex(SEXP fileName, SEXP indexFileName);

