#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#if 0
// 
#define R_NO_REMAP
#endif

#include <R.h>
#include <Rinternals.h>

extern "C" SEXP impl_annotateGene(SEXP param,
                                  SEXP chrom,
                                  SEXP pos,
                                  SEXP ref,
                                  SEXP alt);

extern "C" SEXP impl_anno(SEXP s_inFile,
               SEXP s_outFile,
               SEXP s_param);

extern "C" SEXP impl_getRefBase(SEXP param,
                                SEXP chrom,
                                SEXP pos,
                                SEXP length);
