#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP anno(SEXP, SEXP, SEXP);
extern SEXP annotateGene(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP createTabixIndex(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getRefBase(SEXP, SEXP, SEXP, SEXP);
extern SEXP isInRange(SEXP, SEXP);
extern SEXP readCovByRange(SEXP, SEXP);
extern SEXP readScoreByRange(SEXP, SEXP);
extern SEXP readSkewByRange(SEXP, SEXP);
extern SEXP readTabixByRange(SEXP, SEXP);
extern SEXP readTabixHeader(SEXP);
extern SEXP readTabixSkippedLine(SEXP);
extern SEXP readVCFToListByGene(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP readVCFToListByRange(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP readVCFToMatrixByGene(SEXP, SEXP, SEXP, SEXP);
extern SEXP readVCFToMatrixByRange(SEXP, SEXP, SEXP);
extern SEXP rvMetaReadDataByGene(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rvMetaReadDataByRange(SEXP, SEXP, SEXP, SEXP);
extern SEXP rvMetaWriteCovData(SEXP, SEXP);
extern SEXP rvMetaWriteScoreData(SEXP, SEXP);
extern SEXP readBGENToListByGene(SEXP, SEXP, SEXP);
extern SEXP readBGENToListByRange(SEXP, SEXP);
extern SEXP readBGENToMatrixByGene(SEXP, SEXP, SEXP);
extern SEXP readBGENToMatrixByRange(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"anno",                   (DL_FUNC) &anno,                   3},
    {"annotateGene",           (DL_FUNC) &annotateGene,           5},
    {"createTabixIndex",       (DL_FUNC) &createTabixIndex,       6},
    {"getRefBase",             (DL_FUNC) &getRefBase,             4},
    {"isInRange",              (DL_FUNC) &isInRange,              2},
    {"readCovByRange",         (DL_FUNC) &readCovByRange,         2},
    {"readScoreByRange",       (DL_FUNC) &readScoreByRange,       2},
    {"readSkewByRange",        (DL_FUNC) &readSkewByRange,        2},
    {"readTabixByRange",       (DL_FUNC) &readTabixByRange,       2},
    {"readTabixHeader",        (DL_FUNC) &readTabixHeader,        1},
    {"readTabixSkippedLine",   (DL_FUNC) &readTabixSkippedLine,   1},
    {"readVCFToListByGene",    (DL_FUNC) &readVCFToListByGene,    7},
    {"readVCFToListByRange",   (DL_FUNC) &readVCFToListByRange,   6},
    {"readVCFToMatrixByGene",  (DL_FUNC) &readVCFToMatrixByGene,  4},
    {"readVCFToMatrixByRange", (DL_FUNC) &readVCFToMatrixByRange, 3},
    {"rvMetaReadDataByGene",   (DL_FUNC) &rvMetaReadDataByGene,   5},
    {"rvMetaReadDataByRange",  (DL_FUNC) &rvMetaReadDataByRange,  4},
    {"rvMetaWriteCovData",     (DL_FUNC) &rvMetaWriteCovData,     2},
    {"rvMetaWriteScoreData",   (DL_FUNC) &rvMetaWriteScoreData,   2},
    {"readBGENToListByGene",    (DL_FUNC) &readBGENToListByGene,    3},
    {"readBGENToListByRange",   (DL_FUNC) &readBGENToListByRange,   2},
    {"readBGENToMatrixByGene",  (DL_FUNC) &readBGENToMatrixByGene,  3},
    {"readBGENToMatrixByRange", (DL_FUNC) &readBGENToMatrixByRange, 2},
    {NULL, NULL, 0}
};

void R_init_seqminer(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
