#include "plink2genoLoader.h"

#include <algorithm>
#include <vector>

#include "FileIO.h"
#include "R_CPP_interface.h"
#include "TimeUtil.h"

#include <R.h>

#define UNUSED(x) ((void)(x))

void printTime(const char* s) {
  REprintf("%s - %s\n", currentTime().c_str(), s);
}

int filterIndex(int numElement, std::vector<int>* in) {
  if (!in) return 0;
  std::vector<int>& vec = *in;
  size_t inIndex = 0;
  size_t outIndex = 0;
  int droppedElement = 0;
  for (; inIndex < vec.size(); ++inIndex) {
    if (vec[inIndex] < 0 || vec[inIndex] >= numElement) {
      droppedElement++;
      continue;
    }
    if (inIndex == outIndex) {
      continue;
    }
    vec[outIndex] = vec[inIndex];
    outIndex++;
  }
  vec.resize(outIndex);
  return droppedElement;
}

std::vector<std::string> keepByIndex(const std::vector<std::string>& vec,
                                     const std::vector<int>& idx) {
  // REprintf("vec len = %d, idx len = %d \n", (int)vec.size(),
  // (int)idx.size());
  std::vector<std::string> ret;
  ret.resize(idx.size());
  for (size_t i = 0; i != idx.size(); ++i) {
    ret[i] = vec[idx[i]];
  }
  return ret;
}
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

  // read bim
  printTime("read bim");
  std::vector<std::string> markerNames;  // 1:100_A/C_rs1234
  std::string s;
  LineReader* lr = new LineReader((FLAG_fileName + ".bim").c_str());
  std::vector<std::string> fd;
  while (lr->readLineBySep(&fd, " \t")) {
    if (fd.size() != 6) {
      REprintf("Wrong format in bim file.\n");
      continue;
    }
    s.clear();
    s.append(fd[0]);
    s.append(":");
    s.append(fd[3]);
    // s.append("_");
    // s.append(fd[4]);
    // s.append("/");
    // s.append(fd[5]);
    // s.append("_");
    // s.append(fd[1]);
    markerNames.push_back(s);
  }
  delete lr;

  // read fam
  printTime("read fam");
  std::vector<std::string> sampleNames;  // fid->iid
  lr = new LineReader((FLAG_fileName + ".fam").c_str());
  while (lr->readLineBySep(&fd, " \t")) {
    if (fd.size() != 6) {
      REprintf("Wrong format in fam file.\n");
      continue;
    }
    s.clear();
    // s.append(fd[0]);
    // s.append("->");
    s.append(fd[1]);
    sampleNames.push_back(s);
  }
  delete lr;

  // prepare outputs
  const int numMarker = markerNames.size();
  const int numSample = sampleNames.size();
  REprintf("extract %d marker and %d sample out of %d marker and %d sample\n",
           (int)FLAG_markerIndex.size(), (int)FLAG_indvIndex.size(), numMarker,
           numSample);

  // remove out-of-bound index
  if (filterIndex(numMarker, &FLAG_markerIndex) ||
      filterIndex(numSample, &FLAG_indvIndex)) {
    REprintf(
        "Some indice are invalid, now extract %d marker and %d sample out of "
        "%d marker and %d sample\n",
        (int)FLAG_markerIndex.size(), (int)FLAG_indvIndex.size(), numMarker,
        numSample);
  }

  // usually, PLINK are snp-major, and R is column-based
  // to optimze for a sequential I/O, read 1 snp from PLINK and store 1 column
  // in an R matrix
  // therefore, R matrix is better sample (row) by marker (column)
  allocated +=
      createDoubleArray(FLAG_indvIndex.size() * FLAG_markerIndex.size(), &ret);

  printTime("read bed");
  // check BED to make sure it is snp-major
  FILE* fpBed = fopen((FLAG_fileName + ".bed").c_str(), "rb");

  unsigned char c;
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

  const static unsigned char HOM_REF = 0x0;  // 0b00;
  const static unsigned char HET = 0x2;      // 0b10;
  const static unsigned char HOM_ALT = 0x3;  // 0b11;
  const static unsigned char MISSING = 0x1;  // 0b01;

  // construct look-up table
  REprintf("build look-up table\n");
  double table[256][4];
  for (int i = 0; i < 256; ++i) {
    for (int j = 0; j < 4; ++j) {
      unsigned char geno = (i >> (j << 1)) & 3;
      switch (geno) {
        case HOM_REF:
          table[i][j] = 0;
          break;
        case HET:
          table[i][j] = 1;
          break;
        case HOM_ALT:
          table[i][j] = 2;
          break;
        case MISSING:
          table[i][j] = -9;
          break;
      }
    }
  }

  const int sampleStride = (numSample + 3) / 4;  // every 4 samples take 1 bytes
  int pos_;
  std::vector<unsigned char> geno;
  geno.resize(sampleStride);
  std::vector<double> buffer;
  buffer.resize(sampleStride * 4 + 4);
  if (snpMajorMode) {
    double* pRet = REAL(ret);
    for (int m = 0; m < (int)FLAG_markerIndex.size(); m++) {
      pos_ = 3 + sampleStride * FLAG_markerIndex[m];
      // // extract genotype for all individuals
      // REprintf("read marker %d\n", m);
      fseek(fpBed, pos_, SEEK_SET);
      int numRead =
          fread(geno.data(), sizeof(unsigned char), sampleStride, fpBed);
      UNUSED(numRead);
      // // expand genotypes
      // REprintf("expand genotypes for marker %d\n", m);
      for (int p = 0; p < sampleStride; ++p) {
        memcpy(buffer.data() + 4 * p, &(table[geno[p]][0]), sizeof(double) * 4);
      }
      // // assign results
      // REprintf("assign genotypes for marker %d\n", m);
      for (int p = 0; p < (int)FLAG_indvIndex.size(); ++p) {
        *pRet++ = buffer[FLAG_indvIndex[p]];
      }
    }  // loop marker
    REprintf("assigned %d values \n", pRet - REAL(ret));
  } else {
    REprintf("individual-major mode PLINK is not supported yet!");
  }

  fclose(fpBed);

  REprintf("allocate dim and dimnames\n");
  allocated +=
      setDim((int)FLAG_indvIndex.size(), (int)FLAG_markerIndex.size(), &ret);
  std::vector<std::string> v1 = keepByIndex(sampleNames, FLAG_indvIndex);
  std::vector<std::string> v2 = keepByIndex(markerNames, FLAG_markerIndex);
  allocated += setDimNames(v1, v2, &ret);
  allocated += setDimNames(keepByIndex(sampleNames, FLAG_indvIndex),
                           keepByIndex(markerNames, FLAG_markerIndex), &ret);

  printTime("end");
  UNPROTECT(allocated);

  return ret;
}  // impl_readVCFToListByRange