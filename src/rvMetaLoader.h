#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#include <R.h>
#include <Rinternals.h>

typedef enum GroupingUnit {
  GROUP_BY_CHROM_POS = 1,
  GROUP_BY_CHROM_POS_REF_ALT = 2} GroupingUnit;

extern "C" SEXP impl_rvMetaReadDataByGene(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_geneFile, SEXP arg_gene, SEXP arg_multiAllelic);
extern "C" SEXP impl_rvMetaReadDataByRange(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_range, SEXP arg_multiAllelic);

extern "C" SEXP impl_readCovByRange(SEXP arg_covFile, SEXP arg_range);
extern "C" SEXP impl_readScoreByRange(SEXP arg_scoreFile, SEXP arg_range);
extern "C" SEXP impl_readSkewByRange(SEXP arg_scoreFile, SEXP arg_range);

extern "C" SEXP impl_rvMetaWriteScoreData(SEXP arg_data, SEXP arg_outPrefix);
extern "C" SEXP impl_rvMetaWriteCovData(SEXP arg_data, SEXP arg_outPrefix);

extern "C" SEXP impl_isInRange(SEXP arg_position, SEXP arg_range);
