#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#include <R.h>
#include <Rinternals.h>

extern "C" SEXP impl_readBGENToMatrixByRange(SEXP arg_fileName, SEXP arg_range);
extern "C" SEXP impl_readBGENToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName);

extern "C" SEXP impl_readBGENToListByRange(SEXP arg_fileName, SEXP arg_range);
extern "C" SEXP impl_readBGENToListByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName);

