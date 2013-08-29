#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#include <R.h>
#include <Rinternals.h>

extern "C" SEXP impl_readVCFToMatrixByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType);
extern "C" SEXP impl_readVCFToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType);

extern "C" SEXP impl_readVCFToListByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag);
extern "C" SEXP impl_readVCFToListByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag);

