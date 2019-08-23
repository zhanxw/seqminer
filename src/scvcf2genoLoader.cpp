#include "scvcf2genoLoader.h"

#include "R_CPP_interface.h"
#include "RangeList.h"
#include "SingleChromosomeVCFIndex.h"

SEXP impl_readSingleChromosomeVCFToMatrixByRange(SEXP arg_fileName,
                                                 SEXP arg_indexFileName,
                                                 SEXP arg_range) {
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  std::string FLAG_indexFileName = CHAR(STRING_ELT(arg_indexFileName, 0));
  std::vector<std::string> FLAG_range;
  extractStringArray(arg_range, &FLAG_range);
  SingleChromosomeVCFIndex sc(FLAG_fileName, FLAG_indexFileName);
  if (sc.openIndex()) {
    REprintf("failed to open index!\n");
  }
  SEXP ans;
  // read header to get sample names
  int64_t offset;
  if (sc.query(0, &offset) <= 0) {
    REprintf("Cannot find the header block!\n");
  }

  std::string line;
  std::vector<std::string> fd;

  if (sc.readLine(offset, &line) < 0) {
    REprintf("Cannot readline()!\n");
  }
  stringTokenize(line, '\t', &fd);
  std::vector<std::string> sampleNames(fd.begin() + 9, fd.end());
  const int numSample = sampleNames.size();
  REprintf("Inferred %d samples from header\n", numSample);

  int nGene = FLAG_range.size();
  Rprintf("%d region to be extracted.\n", nGene);
  PROTECT(ans = allocVector(VECSXP, nGene));
  setListNames(FLAG_range, &ans);

  std::string chromName;
  unsigned int chromPosBeg;
  unsigned int chromPosEnd;
  std::vector<std::string> ranges;
  for (int i = 0; i < nGene; ++i) {
    stringTokenize(FLAG_range[i], ',', &ranges);

    // RangeList rl;
    // rl.addRangeList(FLAG_range[i]);
    // int chromPosBeg = rl.begin().getBegin();
    // int chromPosEnd = rl.begin().getEnd();
    // REprintf("query position %d-%d\n", chromPosBeg, chromPosEnd);
    std::vector<double> buf;
    int cumNumVariant = 0;
    std::vector<std::string> markerNames;

    for (size_t j = 0; j != ranges.size(); ++j) {
      parseRangeFormat(ranges[j], &chromName, &chromPosBeg, &chromPosEnd);
      int nVariant = sc.query(chromPosBeg, chromPosEnd, &offset);
      if (nVariant <= 0) {
        REprintf("Cannot find the variant!\n");
        continue;
      }

      // create double array
      // SEXP val;
      // PROTECT(val = allocVector(REALSXP, numSample * nVariant));
      // double* pVal = REAL(val);

      if (sc.readLine(offset, &line) < 0) {
        REprintf("Cannot readline()!\n");
      }
      stringTokenize(line, '\t', &fd);
      markerNames.push_back(fd[0]);
      markerNames.back() += ':';
      markerNames.back() += fd[1];
      markerNames.back() += '_';
      markerNames.back() += fd[3];
      markerNames.back() += '/';
      markerNames.back() += fd[4];      
      for (size_t j = 9; j != fd.size(); ++j) {
        // TODO: marker rigorous check on genotype format is needed
        buf.push_back((fd[j][0] - '0') + (fd[j][2] - '0'));
      }
      for (int remainVariant = nVariant - 1; remainVariant > 0;
           --remainVariant) {
        if (sc.nextLine(&line) < 0) {
          REprintf("Cannot readline()!\n");
        }
        stringTokenize(line, '\t', &fd);
        markerNames.push_back(fd[0]);
        markerNames.back() += ':';
        markerNames.back() += fd[1];
        markerNames.back() += '_';
        markerNames.back() += fd[3];
        markerNames.back() += '/';
        markerNames.back() += fd[4];      
        for (size_t j = 9; j != fd.size(); ++j) {
          // TODO: marker rigorous check on genotype format is needed
          buf.push_back((fd[j][0] - '0') + (fd[j][2] - '0'));
        }
      }
      cumNumVariant += nVariant;
    }
    SEXP val;
    PROTECT(val = allocVector(REALSXP, numSample * cumNumVariant));
    memcpy(REAL(val), buf.data(), sizeof(double) * numSample * cumNumVariant);
    setDim(numSample, cumNumVariant, val);
    setDimNames(sampleNames, markerNames, val);
    UNPROTECT(1);
    SET_VECTOR_ELT(ans, i, val);
  }
  UNPROTECT(1);
  return ans;
}

SEXP impl_createSingleChromosomeVCFIndex(SEXP arg_fileName,
                                         SEXP arg_indexFileName) {
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  std::string FLAG_indexFileName = CHAR(STRING_ELT(arg_indexFileName, 0));

  SingleChromosomeVCFIndex sc(FLAG_fileName, FLAG_indexFileName);
  if (sc.createIndex()) {
    REprintf("create index file successfully!\n");
  }
  REprintf("created index file [ %s ]\n", FLAG_indexFileName.c_str());
  return arg_indexFileName;
}
