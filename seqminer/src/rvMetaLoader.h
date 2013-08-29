#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#include <R.h>
#include <Rinternals.h>

extern "C" SEXP impl_rvMetaReadDataByGene(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_geneFile, SEXP arg_gene);
extern "C" SEXP impl_rvMetaReadDataByRange(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_range);

extern "C" SEXP impl_readCovByRange(SEXP arg_covFile, SEXP arg_range);
extern "C" SEXP impl_readScoreByRange(SEXP arg_scoreFile, SEXP arg_range);

