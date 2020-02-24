#include <R.h>
#include <Rinternals.h>

extern SEXP impl_annotateGene(SEXP arg_param,
                              SEXP arg_chrom,
                              SEXP arg_pos,
                              SEXP arg_ref,
                              SEXP arg_alt);
SEXP annotateGene(SEXP arg_param,
                  SEXP arg_chrom,
                  SEXP arg_pos,
                  SEXP arg_ref,
                  SEXP arg_alt) {
  return impl_annotateGene(arg_param, arg_chrom, arg_pos, arg_ref, arg_alt);
}

extern SEXP impl_anno(SEXP arg_inFile,
                      SEXP arg_outFile,
                      SEXP arg_param);
SEXP anno(SEXP arg_inFile,
          SEXP arg_outFile,
          SEXP arg_param) {
  return impl_anno(arg_inFile, arg_outFile, arg_param);
}

extern SEXP impl_getRefBase(SEXP arg_reference,
                            SEXP arg_chrom,
                            SEXP arg_pos,
                            SEXP arg_len);
SEXP getRefBase(SEXP arg_reference,
                SEXP arg_chrom,
                SEXP arg_pos,
                SEXP arg_len) {
  return impl_getRefBase(arg_reference, arg_chrom, arg_pos, arg_len);
}
