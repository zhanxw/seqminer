#ifndef _PLINK2GENOLOADER_H_
#define _PLINK2GENOLOADER_H_

#include <iostream>  // have to include this to avoid R screw up "length" in iostream/sstream...

#include <R.h>
#include <Rinternals.h>


extern "C" SEXP impl_readPlinkToMatrixByIndex(SEXP arg_fileName,
                                              SEXP arg_sampleIdx,
                                              SEXP arg_markerIdx);

extern "C" SEXP impl_readBedToMatrixByIndex(SEXP arg_fileName,
                                            SEXP arg_numSample,
                                            SEXP arg_numMarker,
                                            SEXP arg_sampleIdx,
                                            SEXP arg_markerIdx);

#endif /* _PLINK2GENOLOADER_H_ */
