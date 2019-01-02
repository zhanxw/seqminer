#include "plink2genoLoader.h"

#include <vector>

#include "R_CPP_interface.h"

#include <R.h>

#define UNUSED(x) ((void)(x))

/**
 * @param arg_fileName: a string character
 * @param arg_range: which range to extract. NOTE: only use first element
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g.
 * Synonymous|Nonsynonymous)
 * @param arg_columns: a list of which columns to extract (e.g. CHROM, POS ...)
 * @param arg_infoTag: a list of which tag under INFO tag will be extracted
 * (e.g. ANNO, ANNO_FULL, AC ...)
 * @param arg_indvTag: a list of which tag given in individual's column (e.g.
 * GT, GD, GQ ...)
 */
SEXP impl_readPlinkToMatrixByIndex(SEXP arg_fileName, SEXP arg_sampleIdx,
                                   SEXP arg_markerIdx) {
  // begin
  REprintf("start\n");
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  REprintf("file = %s\n", FLAG_fileName.c_str());
  
  std::vector<int> FLAG_indvIndex;
  std::vector<int> FLAG_markerIndex;
  extractIntArray(arg_sampleIdx, &FLAG_indvIndex);
  extractIntArray(arg_markerIdx, &FLAG_markerIndex);

  SEXP ret;
  int allocated = 0;
  int numMarker = FLAG_markerIndex.size();
  int numSample = FLAG_indvIndex.size();
  // usually, PLINK are snp-major, and R is column-based
  // to optimze for a sequential I/O, read 1 snp from PLINK and store 1 column
  // in an R matrix
  // therefore, R matrix is better sample (row) by marker (column)
  allocated += createDoubleArray(numSample * numMarker, &ret);

  // check BED to make sure it is snp-major
  FILE* fpBed = fopen((FLAG_fileName + ".bed").c_str(), "rb");

  char c;
  // magic number
  const char magic1 = 0x6c;  // 0b01101100;
  int numRead = fread(&c, sizeof(char), 1, fpBed);
  UNUSED(numRead);
  // assert(numRead == 1);
  if (c != magic1) {
    REprintf("Magic number of binary PLINK file does not match!\n");
    REprintf("Critical error happening!\n");  // abort();
  }
  int magic2 = 0x1b;  // 0b00011011;
  numRead = fread(&c, sizeof(char), 1, fpBed);
  // assert(numRead == 1);
  if (c != magic2) {
    REprintf("Magic number of binary PLINK file does not match!\n");
    REprintf("Critical error happening!\n");  // abort();
  }

  // snp major mode
  const int SNP_MAJOR_MODE = 0x01;  // 0b00000001;
  const int INDV_MAJOR_MODE = 0x00;
  bool snpMajorMode;
  numRead = fread(&c, sizeof(char), 1, fpBed);
  // assert(numRead == 1);
  if (c == SNP_MAJOR_MODE) {
    snpMajorMode = true;
    REprintf("binary PLINK BED file is ready to read\n");
  } else {
    snpMajorMode = false;
    REprintf("individual-major mode PLINK is not supported yet!\n");
  }

  const static unsigned char HOM_REF = 0x0;     //0b00;
  const static unsigned char HET = 0x2;         //0b10;
  const static unsigned char HOM_ALT = 0x3;     //0b11;
  const static unsigned char MISSING = 0x1;     //0b01;
  const unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100, 0b110000, 0b11000000
  int retIdx = 0;
  if (snpMajorMode) {
    for (int p = 0; p < numSample; p++) {
      for (int m = 0; m < numMarker; m++) {
        // get file position
        int pos_ = 3 + (numSample / 4 + 1) * FLAG_markerIndex[m] +
                   FLAG_indvIndex[p] / 4;
        int offset = FLAG_indvIndex[p] % 4;
        unsigned char c;
        fseek(fpBed, pos_, SEEK_SET);
        int numRead = fread(&c, sizeof(unsigned char), 1, fpBed);
        UNUSED(numRead);
        unsigned char geno = (c & mask[offset]) >> (offset << 1);
        switch (geno) {
          case HOM_REF:
            REAL(ret)[retIdx++] = 0;
            break;
          case HET:
            REAL(ret)[retIdx++] = 1;
            break;
          case HOM_ALT:
            REAL(ret)[retIdx++] = 2;
            break;
          case MISSING:
            REAL(ret)[retIdx++] = -9;
            break;
          default:
            REprintf("Read PLINK genotype error!\n");
            break;
        };
      }
    }
  } else {
    REprintf("individual-major mode PLINK is not supported yet!");
  }

  fclose(fpBed);
  allocated += setDim(numMarker, numSample, &ret);
  UNPROTECT(allocated);
  
  return ret;
}  // impl_readVCFToListByRange
