#include "bgen2genoLoader.h"

#include <map>
#include <set>
#include <string>
#include <vector>

#include "BGenFile.h"
#include "R_CPP_interface.h"

#include <R.h>

#include "GeneLoader.h"

/**
 * Read from @param vin and return a matrix of marker by people
 */
SEXP readBGEN2Matrix(BGenFile* bin) {
  std::vector<double> genoVec;
  std::vector<std::string> posVec;
  std::vector<std::string> idVec;
  std::string posString;

  // print header
  const int N = bin->getNumSample();
  // const int M = bin->getNumMarker();
  std::vector<std::string> sm = bin->getSampleIdentifier();  // all sample names
  std::vector<std::string>& names = idVec;
  if (!sm.size()) {
    char buf[1024];
    for (int i = 0; i < N; ++i) {
      sprintf(buf, "sample_%d", i);
      sm.push_back(buf);
    }
  }

  const size_t sampleSize = bin->getNumEffectiveSample();
  for (size_t i = 0; i != sampleSize; ++i) {
    names.push_back(sm[bin->getEffectiveIndex(i)]);
  }

  while (bin->readRecord()) {
    // REprintf("read a record\n");
    const BGenVariant& var = bin->getVariant();

    // store all results here
    posString = var.chrom;
    posString += ':';
    posString += toString(var.pos);
    posVec.push_back(posString);

    for (size_t i = 0; i < sampleSize; i++) {
      const BGenVariant& var = bin->getVariant();
      genoVec.push_back(var.computeDosage(i));
      // Rprintf( "\t%d", g);
    }
    // Rprintf( "\n");
  }  // end while

  //  REprintf("posVec = %zu, idVec = %zu, genoVec = %zu\n", posVec.size(),
  //  idVec.size(), genoVec.size());

  // pass value back to R (see Manual Chapter 5)

  int nx = (int)posVec.size();
  int ny = (int)idVec.size();

  SEXP ans = R_NilValue;

  PROTECT(ans = allocMatrix(REALSXP, nx, ny));
  double* rans = REAL(ans);
  int idx = 0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      // Rprintf("idx = %d, i = %d, j=%d, geno = %g\n", idx, i, j,
      // genoVec[idx]);
      if (genoVec[idx] < 0) {
        rans[i + nx * j] = NA_REAL;
      } else {
        rans[i + nx * j] = genoVec[idx];
      }
      ++idx;
    }
  }

  // set row and col names
  SEXP dim;
  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = nx;
  INTEGER(dim)[1] = ny;
  setAttrib(ans, R_DimSymbol, dim);

  SEXP rowName;
  PROTECT(rowName = allocVector(STRSXP, nx));
  for (int i = 0; i < nx; i++)
    SET_STRING_ELT(rowName, i, mkChar(posVec[i].c_str()));
  SEXP colName;
  PROTECT(colName = allocVector(STRSXP, ny));
  for (int i = 0; i < ny; i++)
    SET_STRING_ELT(colName, i, mkChar(idVec[i].c_str()));

  SEXP dimnames;
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, rowName);
  SET_VECTOR_ELT(dimnames, 1, colName);
  setAttrib(ans, R_DimNamesSymbol, dimnames);

  // finish up
  UNPROTECT(5);
  return (ans);
}  // end readBGEN2Matrix

/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (just allow One gene
 * name).
 */
SEXP impl_readBGENToMatrixByRange(SEXP arg_fileName, SEXP arg_range) {
  SEXP ans = R_NilValue;

  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  std::vector<std::string> FLAG_range;
  extractStringArray(arg_range, &FLAG_range);

  if (FLAG_fileName.size() == 0) {
    error("Please provide BGEN file name");
    return ans;
  }
  if (FLAG_range.empty()) {
    error("Please provide a given range, e.g. '1:100-200'");
    return ans;
  }

  int nGene = FLAG_range.size();
  Rprintf("%d region to be extracted.\n", nGene);
  int numAllocated = 0;

  // allocate return value
  PROTECT(ans = allocVector(VECSXP, nGene));
  numAllocated++;
  setListNames(FLAG_range, &ans);

  for (int i = 0; i < nGene; ++i) {
    // REprintf("range = %s\n", FLAG_range[i].c_str());
    BGenFile bin(FLAG_fileName);
    bin.setRangeList(FLAG_range[i].c_str());

    // real working part
    SET_VECTOR_ELT(ans, i, readBGEN2Matrix(&bin));
  }
  UNPROTECT(numAllocated);
  return ans;
}  // end impl_readBGENToMatrixByRange

/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (just allow One gene
 * name).
 */
SEXP impl_readBGENToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile,
                                 SEXP arg_geneName) {
  SEXP ans = R_NilValue;

  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  std::string FLAG_geneFile = CHAR(STRING_ELT(arg_geneFile, 0));
  std::vector<std::string> FLAG_geneName;
  extractStringArray(arg_geneName, &FLAG_geneName);

  if (FLAG_fileName.size() == 0) {
    error("Please provide BGEN file name");
  }
  if (FLAG_geneName.size() && FLAG_geneFile.size() == 0) {
    error("Please provide gene file name when extract genotype by gene");
  }

  int nGene = FLAG_geneName.size();
  Rprintf("%d region to be extracted.\n", nGene);
  int numAllocated = 0;

  // allocate return value
  PROTECT(ans = allocVector(VECSXP, nGene));
  numAllocated++;
  setListNames(FLAG_geneName, &ans);

  OrderedMap<std::string, std::string> geneRange;
  loadGeneFile(FLAG_geneFile, FLAG_geneName, &geneRange);
  for (int i = 0; i < nGene; ++i) {
    // REprintf("range = %s\n", FLAG_geneName[i].c_str());
    const std::string& range = geneRange[FLAG_geneName[i]];

    // Rprintf( "range = %s\n", range.c_str());
    BGenFile bin(FLAG_fileName);
    if (range.size())
      bin.setRangeList(range.c_str());
    else {
      warning("Gene name [ %s ] does not exists in provided gene file",
              FLAG_geneName[i].c_str());
      UNPROTECT(numAllocated);
      return (ans);
    };

    // real working part
    SET_VECTOR_ELT(ans, i, readBGEN2Matrix(&bin));
  }
  UNPROTECT(numAllocated);
  return ans;
}

void mypause() {
  REprintf("--------------------------------------------------\n");
}

SEXP readBGEN2List(BGenFile* bin) {
  // Rprintf("vcfColumn.size() = %u\n", FLAG_vcfColumn.size());
  // Rprintf("vcfInfo.size() = %u\n", FLAG_infoTag.size());
  // Rprintf("vcfIndv.size() = %u\n", FLAG_indvTag.size());
  // also append sample names at the end
  // 7: chrom, pos, varId, rsId, alleles, isPhased, prob, sampleId
  int retListLen = 8;
  if (retListLen == 0) {
    return R_NilValue;
  }

  int numAllocated =
      0;  // record how many times we allocate (using PROTECT in R);
  SEXP ret;
  PROTECT(ret = allocVector(VECSXP, retListLen));
  numAllocated++;

  //  store results
  std::vector<std::string> idVec;
  std::vector<std::string> chrom;
  std::vector<int> pos;
  std::vector<std::string> varId;
  std::vector<std::string> rsId;
  std::vector<std::string> alleles;
  // std::vector<std::vector<bool> > missing;
  std::vector<bool> isPhased;
  std::vector<std::vector<double> >
      prob;  // prob[variant][each_sample * (prob1, prob2, ...)]

  // std::map<std::string, std::vector<std::string> > infoMap;

  // std::map<std::string, std::vector<std::string> > indvMap;
  /// int nRow = 0;  // # of positions that will be outputed

  // get effective sample names
  const int N = bin->getNumSample();
  std::vector<std::string> sm = bin->getSampleIdentifier();  // all sample names
  std::vector<std::string>& names = idVec;
  if (!sm.size()) {
    char buf[1024];
    for (int i = 0; i < N; ++i) {
      sprintf(buf, "sample_%d", i);
      sm.push_back(buf);
    }
  }

  const size_t sampleSize = bin->getNumEffectiveSample();
  for (size_t i = 0; i != sampleSize; ++i) {
    names.push_back(sm[bin->getEffectiveIndex(i)]);
  }

  // real working part
  int nRecord = 0;
  const int numProbValues =
      3;  // if multi-allelic/multi-haploid, this value can be different
  int maxProbValues = -1;
  while (bin->readRecord()) {
    // REprintf("read a record\n");
    const BGenVariant& var = bin->getVariant();
    const size_t sampleSize = bin->getNumEffectiveSample();

    // store results here
    nRecord++;
    chrom.push_back(var.chrom);
    pos.push_back(var.pos);
    varId.push_back(var.varid);
    rsId.push_back(var.rsid);
    alleles.push_back(toString(var.alleles, ","));
    isPhased.push_back(var.isPhased);
    prob.resize(nRecord);

    std::vector<double>& p = prob[nRecord - 1];
    p.reserve(sampleSize * numProbValues);

    for (size_t i = 0; i != sampleSize; ++i) {
      int beg = var.index[bin->getEffectiveIndex(i)];
      int end = var.index[bin->getEffectiveIndex(i) + 1];
      if (end - beg > maxProbValues) {
        maxProbValues = end - beg;
      }
      for (int j = 0; j < numProbValues; ++j) {
        if (j < numProbValues) {
          p.push_back(var.prob[beg + j]);
        } else {
          p.push_back(-9);
        }
      }
      // REprintf("beg = %d, end = %d, prob[%d][%d] len = %d\n", beg,end,
      // nRecord - 1, i, p[i].size());
    }

    // Rprintf("Done add indv\n");
  }  // end while
  if (maxProbValues > numProbValues) {
    REprintf("some sample has more than %d > %d probabilities per variant!\n",
             maxProbValues, numProbValues);
  }

  // pass value back to R (see Manual Chapter 5)
  std::vector<std::string> listNames;
  int retListIdx = 0;
  storeResult(chrom, ret, retListIdx++);
  storeResult(pos, ret, retListIdx++);
  storeResult(varId, ret, retListIdx++);
  storeResult(rsId, ret, retListIdx++);
  storeResult(alleles, ret, retListIdx++);
  storeResult(isPhased, ret, retListIdx++);
  storeResult(prob, ret, retListIdx);
  for (size_t i = 0; i != prob.size(); ++i) {
    SEXP s = VECTOR_ELT(VECTOR_ELT(ret, retListIdx), i);
    setDim(numProbValues, sampleSize, s);
  }

  retListIdx++;
  listNames.push_back("chrom");
  listNames.push_back("pos");
  listNames.push_back("varid");
  listNames.push_back("rsid");
  listNames.push_back("alleles");
  listNames.push_back("isPhased");
  listNames.push_back("probability");

  // store sample ids
  // Rprintf("set sample id");
  listNames.push_back("sampleId");
  storeResult(idVec, ret, retListIdx++);

  // Rprintf("set list names\n");
  SEXP sListNames;
  PROTECT(sListNames = allocVector(STRSXP, listNames.size()));
  numAllocated++;
  for (unsigned int i = 0; i != listNames.size(); ++i) {
    SET_STRING_ELT(sListNames, i, mkChar(listNames[i].c_str()));
  }
  setAttrib(ret, R_NamesSymbol, sListNames);

  // finish up
  UNPROTECT(numAllocated);
  // Rprintf("Unprotected: %d\n", (retListLen + 1));
  return (ret);
}

/**
 * @param arg_fileName: a string character
 * @param arg_range: which range to extract. NOTE: only use first element
 */
SEXP impl_readBGENToListByRange(SEXP arg_fileName, SEXP arg_range) {
  // begin
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  std::string FLAG_range = CHAR(STRING_ELT(arg_range, 0));

  // Rprintf( "range = %s\n", range.c_str());
  BGenFile bin(FLAG_fileName.c_str());
  if (FLAG_range.size())
    bin.setRangeList(FLAG_range.c_str());
  else {
    error("Please provide a range before we can continue.\n");
  }
  return readBGEN2List(&bin);
}  // impl_readBGENToListByRange

/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (NOTE: only first one gene
 * is used).
 */
SEXP impl_readBGENToListByGene(SEXP arg_fileName, SEXP arg_geneFile,
                               SEXP arg_geneName) {
  // begin
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  std::string FLAG_geneFile = CHAR(STRING_ELT(arg_geneFile, 0));
  std::string FLAG_geneName = CHAR(STRING_ELT(arg_geneName, 0));

  OrderedMap<std::string, std::string> geneRange;
  loadGeneFile(FLAG_geneFile, FLAG_geneName, &geneRange);
  std::string range;
  int n = geneRange.size();
  for (int i = 0; i < n; ++i) {
    if (range.size() > 0) {
      range += ",";
    }
    range += geneRange.valueAt(i);
  }

  REprintf("range = %s\n", range.c_str());
  BGenFile bin(FLAG_fileName.c_str());
  if (range.size())
    bin.setRangeList(range.c_str());
  else {
    error("Please provide a valid gene name before we can continue.\n");
  };

  return readBGEN2List(&bin);
}  // end readBGEN2List
