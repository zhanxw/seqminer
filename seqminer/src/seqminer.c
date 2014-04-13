#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>

extern SEXP impl_readVCFToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType);
SEXP readVCFToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType) {
  return impl_readVCFToMatrixByGene(arg_fileName, arg_geneFile, arg_geneName, arg_annoType);
}

extern SEXP impl_readVCFToMatrixByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType);
SEXP readVCFToMatrixByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType) {
  return impl_readVCFToMatrixByRange(arg_fileName, arg_range, arg_annoType);
}


extern SEXP impl_readVCFToListByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag);
SEXP readVCFToListByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag) {
  return impl_readVCFToListByGene(arg_fileName, arg_geneFile, arg_geneName, arg_annoType, arg_columns, arg_infoTag, arg_indvTag);
}

extern SEXP impl_readVCFToListByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag);
SEXP readVCFToListByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag) {
  return impl_readVCFToListByRange(arg_fileName, arg_range, arg_annoType, arg_columns, arg_infoTag, arg_indvTag);
}

extern SEXP impl_rvMetaReadDataByGene(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_geneFile, SEXP arg_gene);
SEXP rvMetaReadDataByGene(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_geneFile, SEXP arg_gene) {
  return impl_rvMetaReadDataByGene(arg_pvalFile, arg_covFile, arg_geneFile, arg_gene);
}

extern SEXP impl_rvMetaReadDataByRange(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_range);
SEXP rvMetaReadDataByRange(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_range) {
  return impl_rvMetaReadDataByRange(arg_pvalFile, arg_covFile, arg_range);
}

extern SEXP impl_readCovByRange(SEXP arg_covFile, SEXP arg_range);
SEXP readCovByRange(SEXP arg_covFile, SEXP arg_range) {
  return impl_readCovByRange(arg_covFile, arg_range);
}

extern SEXP impl_readScoreByRange(SEXP arg_covFile, SEXP arg_range);
SEXP readScoreByRange(SEXP arg_covFile, SEXP arg_range) {
  return impl_readScoreByRange(arg_covFile, arg_range);
}

extern SEXP impl_readSkewByRange(SEXP arg_covFile, SEXP arg_range);
SEXP readSkewByRange(SEXP arg_covFile, SEXP arg_range) {
  return impl_readSkewByRange(arg_covFile, arg_range);
}

extern SEXP impl_readTabixByRange(SEXP arg_tabixFile, SEXP arg_range);
SEXP readTabixByRange(SEXP arg_tabixFile, SEXP arg_range) {
  return impl_readTabixByRange(arg_tabixFile, arg_range);
}

extern SEXP impl_readTabixSkippedLine(SEXP arg_tabixFile);
SEXP readTabixSkippedLine(SEXP arg_tabixFile) {
  return impl_readTabixSkippedLine(arg_tabixFile);
}

extern SEXP impl_readTabixHeader(SEXP arg_tabixFile);
SEXP readTabixHeader(SEXP arg_tabixFile) {
  return impl_readTabixHeader(arg_tabixFile);
}

extern SEXP impl_rvMetaWriteScoreData(SEXP arg_data, SEXP arg_outPrefix);
SEXP rvMetaWriteScoreData(SEXP arg_data, SEXP arg_outPrefix) {
  return impl_rvMetaWriteScoreData(arg_data, arg_outPrefix);
}

extern SEXP impl_rvMetaWriteCovData(SEXP arg_data, SEXP arg_outPrefix);
SEXP rvMetaWriteCovData(SEXP arg_data, SEXP arg_outPrefix) {
  return impl_rvMetaWriteCovData(arg_data, arg_outPrefix);
}

extern SEXP impl_createTabixIndex(SEXP arg_tabixFile,
                                      SEXP arg_seqnameColumn,
                                      SEXP arg_startColumn,
                                      SEXP arg_endColumn,
                                      SEXP arg_commentChar,
                                      SEXP arg_skipLine);
SEXP createTabixIndex(SEXP arg_tabixFile,
                      SEXP arg_seqnameColumn,
                      SEXP arg_startColumn,
                      SEXP arg_endColumn,
                      SEXP arg_commentChar,
                      SEXP arg_skipLine){
  return impl_createTabixIndex(arg_tabixFile,
                               arg_seqnameColumn,
                               arg_startColumn,
                               arg_endColumn,
                               arg_commentChar,
                               arg_skipLine);
}
