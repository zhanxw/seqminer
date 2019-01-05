#include "vcf2genoLoader.h"

#include <string>
#include <map>
#include <vector>
#include <set>

#include "VCFUtil.h"
#include "tabix.h"
#include "TabixReader.h"
#include "R_CPP_interface.h"

#include <R.h>

#include "GeneLoader.h"

/**
 * Read from @param vin and return a matrix of marker by people
 */
SEXP readVCF2Matrix(VCFExtractor* vin) {
  std::vector<double> genoVec;
  std::vector<std::string> posVec;
  std::vector<std::string> idVec;
  std::string posString;

  // print header
  std::vector<std::string>& names = idVec;
  vin->getVCFHeader()->getPeopleName(&names);

  while (vin->readRecord()){
    // REprintf("read a record\n");
    VCFRecord& r = vin->getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    // store all results here
    posString = r.getChrom();
    posString += ':';
    posString += r.getPosStr();
    posVec.push_back(posString);

    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];
      int g = indv->justGet(0).getGenotype();
      //Rprintf( "\t%d", g);
      genoVec.push_back(g);
    }
    //Rprintf( "\n");
  }; // end while

  //  REprintf("posVec = %zu, idVec = %zu, genoVec = %zu\n", posVec.size(), idVec.size(), genoVec.size());

  // pass value back to R (see Manual Chapter 5)


  int nx = (int) posVec.size();
  int ny = (int) idVec.size();

  SEXP ans = R_NilValue;

  PROTECT(ans = allocMatrix(REALSXP, nx, ny));
  double* rans = REAL(ans);
  int idx = 0;
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j <ny ; j++ ) {
      // Rprintf("idx = %d, i = %d, j=%d, geno = %g\n", idx, i, j, genoVec[idx]);
      rans[i + nx * j] = genoVec[idx];
      ++idx;
    }
  }

  // set row and col names
  SEXP dim;
  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = nx; INTEGER(dim)[1] = ny;
  setAttrib(ans, R_DimSymbol, dim);

  SEXP rowName;
  PROTECT(rowName=allocVector(STRSXP, nx));
  for (int i = 0; i < nx; i++ )
    SET_STRING_ELT(rowName, i, mkChar(posVec[i].c_str()));
  SEXP colName;
  PROTECT(colName=allocVector(STRSXP, ny));
  for (int i = 0; i < ny; i++ )
    SET_STRING_ELT(colName, i, mkChar(idVec[i].c_str()));

  SEXP dimnames;
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, rowName);
  SET_VECTOR_ELT(dimnames, 1, colName);
  setAttrib(ans, R_DimNamesSymbol, dimnames);

  // finish up
  UNPROTECT(5);
  return(ans);
} // end readVCF2Matrix

/**
 * @param file name
 * @return check whether VCF files has ANNO tag
 */
bool vcfHasAnnotation(const std::string& fn) {
  //Rprintf( "range = %s\n", range.c_str());
  VCFInputFile vin(fn);
  while (vin.readRecord()) {
    VCFRecord& r = vin.getVCFRecord();
    VCFInfo& info = r.getVCFInfo();
    bool tagMissing;
    info.getTag("ANNO", &tagMissing);
    if (tagMissing) {
      return false;
    }
    return true;
  }
  return false;
}

/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (just allow One gene name).
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g. Synonymous|Nonsynonymous)
 */
SEXP impl_readVCFToMatrixByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType) {
  SEXP ans = R_NilValue;

  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName,0));
  std::vector<std::string> FLAG_range;
  extractStringArray(arg_range, &FLAG_range);
  std::string FLAG_annoType = CHAR(STRING_ELT(arg_annoType,0));

  if (FLAG_fileName.size() == 0) {
    error("Please provide VCF file name");
    return ans;
  }
  if (FLAG_range.empty()) {
    error("Please provide a given range, e.g. '1:100-200'");
    return ans;
  }

  if (!FLAG_annoType.empty() && !vcfHasAnnotation(FLAG_fileName)) {
    REprintf("Please use annotated VCF as input (cannot find ANNO in the INFO field);\n");
    REprintf("Prefer using ANNO from https://github.com/zhanxw/anno  \n");
    return ans;
  }

  int nGene = FLAG_range.size();
  Rprintf("%d region to be extracted.\n", nGene);
  int numAllocated = 0;
  
  // allocate return value
  PROTECT(ans = allocVector(VECSXP, nGene));
  numAllocated++;
  numAllocated += setListNames(FLAG_range, &ans);
  
  for (int i = 0; i < nGene; ++i) {
    // REprintf("range = %s\n", FLAG_range[i].c_str());
    VCFExtractor vin(FLAG_fileName.c_str());
    vin.setRangeList(FLAG_range[i].c_str());

    if (FLAG_annoType.size()) {
      vin.setAnnoType(FLAG_annoType.c_str());
    }
    // real working part
    SET_VECTOR_ELT(ans, i, readVCF2Matrix(&vin));
  }
  UNPROTECT(numAllocated);
  return ans;
} //end impl_readVCFToMatrixByRange

/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (just allow One gene name).
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g. Synonymous|Nonsynonymous)
 */
SEXP impl_readVCFToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType) {
  SEXP ans = R_NilValue;

  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName,0));
  std::string FLAG_geneFile = CHAR(STRING_ELT(arg_geneFile,0));
  std::vector<std::string> FLAG_geneName;
  extractStringArray(arg_geneName, &FLAG_geneName);
  std::string FLAG_annoType = CHAR(STRING_ELT(arg_annoType,0));

  if (FLAG_fileName.size() == 0) {
    error("Please provide VCF file name");
  }
  if (FLAG_geneName.size() && FLAG_geneFile.size() == 0) {
    error("Please provide gene file name when extract genotype by gene");
  }
  if (!FLAG_annoType.empty() && !vcfHasAnnotation(FLAG_fileName)) {
    REprintf("Please use annotated VCF as input (cannot find ANNO in the INFO field);\n");
    REprintf("Prefer using ANNO from https://github.com/zhanxw/anno  \n");
    return ans;
  }

  int nGene = FLAG_geneName.size();
  Rprintf("%d region to be extracted.\n", nGene);
  int numAllocated = 0;
  
  // allocate return value
  PROTECT(ans = allocVector(VECSXP, nGene));
  numAllocated++;
  numAllocated += setListNames(FLAG_geneName, &ans);

  OrderedMap< std::string, std::string> geneRange;
  loadGeneFile(FLAG_geneFile, FLAG_geneName, &geneRange);  
  for (int i = 0; i < nGene; ++i) {
    // REprintf("range = %s\n", FLAG_geneName[i].c_str());
    const std::string& range = geneRange[FLAG_geneName[i]];

    //Rprintf( "range = %s\n", range.c_str());
    VCFExtractor vin(FLAG_fileName.c_str());
    if (range.size())
      vin.setRangeList(range.c_str());
    else {
      warning("Gene name [ %s ] does not exists in provided gene file", FLAG_geneName[i].c_str());
      UNPROTECT(numAllocated);
      return (ans);
    };

    if (FLAG_annoType.size()) {
      vin.setAnnoType(FLAG_annoType.c_str());
    }

    // real working part
    SET_VECTOR_ELT(ans, i, readVCF2Matrix(&vin));
  }
  UNPROTECT(numAllocated);
  return ans;
}

SEXP readVCF2List(VCFInputFile* vin,
                  const std::set<std::string>& FLAG_vcfColumn,
                  const std::vector<std::string>& FLAG_infoTag,
                  const std::vector<std::string>& FLAG_indvTag) {
  // Rprintf("vcfColumn.size() = %u\n", FLAG_vcfColumn.size());
  // Rprintf("vcfInfo.size() = %u\n", FLAG_infoTag.size());
  // Rprintf("vcfIndv.size() = %u\n", FLAG_indvTag.size());
  // also append sample names at the end
  int retListLen = FLAG_vcfColumn.size() + FLAG_infoTag.size() + FLAG_indvTag.size() + 1;
  if (retListLen == 0) {
    return R_NilValue;
  }

  int numAllocated = 0; // record how many times we allocate (using PROTECT in R);
  SEXP ret;
  PROTECT(ret = allocVector(VECSXP, retListLen));
  numAllocated ++;

  //  store results
  std::vector<std::string> idVec;
  std::vector<std::string> chrom;
  std::vector<int> pos;
  std::vector<std::string> rsId;
  std::vector<std::string> ref;
  std::vector<std::string> alt;
  std::vector<std::string> qual;
  std::vector<std::string> filt;
  std::vector<std::string> info;
  std::vector<std::string> format;

  std::map<std::string, std::vector<std::string> > infoMap;

  // std::vector<int> gtVec;
  // std::vector<int> gdVec;
  // std::vector<int> gqVec;

  std::map<std::string, std::vector<std::string> > indvMap;
  int nRow = 0; // # of positions that will be outputed

  // print header
  std::vector<std::string>& names = idVec;
  vin->getVCFHeader()->getPeopleName(&names);


  bool FLAG_variantOnly = false;
  // real working part
  int nonVariantSite = 0;
  while (vin->readRecord()){
    // REprintf("read a record\n");
    VCFRecord& r = vin->getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    if (FLAG_variantOnly) {
      // REprintf("filter by var\n");
      bool hasVariant = false;
      int geno;
      int GTidx = r.getFormatIndex("GT");
      for (size_t i = 0; i < people.size() ;i ++) {
        indv = people[i];
        geno = indv->justGet(GTidx).getGenotype();
        if (geno != 0 && geno != MISSING_GENOTYPE)
          hasVariant = true;
      }
      if (!hasVariant) {
        nonVariantSite++;
        continue;
      }
    }

    // store results here
    nRow++;
    if (FLAG_vcfColumn.count("CHROM")){
      chrom.push_back(r.getChrom());
    }
    if (FLAG_vcfColumn.count("POS")){
      pos.push_back(r.getPos());
    }
    if (FLAG_vcfColumn.count("ID")){
      rsId.push_back(r.getID());
    }
    if (FLAG_vcfColumn.count("REF")){
      ref.push_back(r.getRef());
    }
    if (FLAG_vcfColumn.count("ALT")){
      alt.push_back(r.getAlt());
    }
    if (FLAG_vcfColumn.count("QUAL")){
      qual.push_back(r.getQual());
    }
    if (FLAG_vcfColumn.count("FILTER")){
      filt.push_back(r.getFilt());
    }
    if (FLAG_vcfColumn.count("INFO")){
      info.push_back(r.getInfo());
    }
    if (FLAG_vcfColumn.count("FORMAT")){
      format.push_back(r.getFormat());
    }

    // store INFO field
    for (std::vector<std::string>::const_iterator it = FLAG_infoTag.begin(); it != FLAG_infoTag.end(); ++it) {
      bool missing;
      VCFValue v = r.getInfoTag(it->c_str(), &missing);
      if (missing) {
        infoMap[ *it ].push_back("");
      } else {
        infoMap[ *it ].push_back(v.toStr());
        // Rprintf("add info field [ %s ] = %s\n", it->c_str(), v.toStr());
      }
    };
    // Rprintf("Done add info\n");

    // store indv values
    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];

      for (std::vector<std::string>::const_iterator it = FLAG_indvTag.begin(); it != FLAG_indvTag.end(); ++it) {
        int idx = r.getFormatIndex(it->c_str());
        if (idx < 0) {
          indvMap[ *it ].push_back("");
        } else {
          bool missing;
          VCFValue v = indv->get(idx, &missing);
          if (missing) {
            indvMap[ *it ].push_back("");
          } else{
            indvMap[ *it ].push_back(v.toStr());
            // Rprintf("add indv field [ %s ] = %s\n", it->c_str(), v.toStr());
          }
        }
      };
    }
    // Rprintf("Done add indv\n");
  }; // end while
  //   Rprintf("indvMap.size() = %zu\n", indvMap.size());
  // REprintf("posVec = %zu, idVec = %zu, genoVec = %zu\n", posVec.size(), idVec.size(), genoVec.size());

  // pass value back to R (see Manual Chapter 5)
  std::vector<std::string> listNames;
  int retListIdx = 0;
  if (FLAG_vcfColumn.count("CHROM")) {
    numAllocated += storeResult(chrom, ret, retListIdx++);
    listNames.push_back("CHROM");
  }
  if (FLAG_vcfColumn.count("POS")) {
    numAllocated += storeResult(pos, ret, retListIdx++);
    listNames.push_back("POS");
  }
  if (FLAG_vcfColumn.count("ID"))  {
    numAllocated += storeResult(rsId, ret, retListIdx++);
    listNames.push_back("ID");
  }
  if (FLAG_vcfColumn.count("REF")) {
    numAllocated += storeResult(ref, ret, retListIdx++);
    listNames.push_back("REF");
  }
  if (FLAG_vcfColumn.count("ALT")) {
    numAllocated += storeResult(alt, ret, retListIdx++);
    listNames.push_back("ALT");
  }
  if (FLAG_vcfColumn.count("QUAL")) {
    numAllocated += storeResult(qual, ret, retListIdx++);
    listNames.push_back("QUAL");
  }
  if (FLAG_vcfColumn.count("FILTER")) {
    numAllocated += storeResult(filt, ret, retListIdx++);
    listNames.push_back("FILTER");
  }
  if (FLAG_vcfColumn.count("INFO")) {
    numAllocated += storeResult(info, ret, retListIdx++);
    listNames.push_back("INFO");
  }
  if (FLAG_vcfColumn.count("FORMAT")) {
    numAllocated += storeResult(format, ret, retListIdx++);
    listNames.push_back("FORMAT");
  }
  // pass info values to R
  for ( std::map<std::string, std::vector<std::string> >::iterator it = infoMap.begin();
        it != infoMap.end();
        ++it) {
    numAllocated += storeResult(it->first, it->second, ret, retListIdx++);
    listNames.push_back(it->first);
  }
  // pass indv tags to R
  // Rprintf("pass idnv tags\n");
  for ( std::map<std::string, std::vector<std::string> >::iterator it = indvMap.begin();
        it != indvMap.end();
        ++it) {

    // dump(it->second);
    numAllocated += storeResult(it->first, it->second, ret, retListIdx);
    // Rprintf("results done\n");
    // NOTE: R internally store values into matrix by column first!
    // thus the matrix is people by marker
    numAllocated += setDim(idVec.size(), nRow, ret, retListIdx);
    retListIdx ++;
    listNames.push_back(it->first);
  }
  // Rprintf("pass idnv tags done.\n");

  // store sample ids
  // Rprintf("set sample id");
  listNames.push_back("sampleId");
  numAllocated += storeResult(idVec, ret, retListIdx++);

  // Rprintf("set list names\n");
  SEXP sListNames;
  PROTECT(sListNames = allocVector(STRSXP, listNames.size()));
  numAllocated ++;
  for (unsigned int i = 0; i != listNames.size(); ++i){
    SET_STRING_ELT(sListNames, i, mkChar(listNames[i].c_str()));
  }
  setAttrib(ret, R_NamesSymbol, sListNames);

  // finish up
  UNPROTECT(numAllocated);
  // Rprintf("Unprotected: %d\n", (retListLen + 1));
  return(ret);
}

/**
 * @param arg_fileName: a string character
 * @param arg_range: which range to extract. NOTE: only use first element
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g. Synonymous|Nonsynonymous)
 * @param arg_columns: a list of which columns to extract (e.g. CHROM, POS ...)
 * @param arg_infoTag: a list of which tag under INFO tag will be extracted (e.g. ANNO, ANNO_FULL, AC ...)
 * @param arg_indvTag: a list of which tag given in individual's column (e.g. GT, GD, GQ ...)
 */
SEXP impl_readVCFToListByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag){
  // begin
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName,0));
  std::string FLAG_range = CHAR(STRING_ELT(arg_range,0));
  std::string FLAG_annoType = CHAR(STRING_ELT(arg_annoType,0));

  std::set<std::string> FLAG_vcfColumn;
  std::vector<std::string> FLAG_infoTag, FLAG_indvTag;

  extractStringSet(arg_columns, &FLAG_vcfColumn);
  extractStringArray(arg_infoTag, &FLAG_infoTag);
  extractStringArray(arg_indvTag, &FLAG_indvTag);

  //Rprintf( "range = %s\n", range.c_str());
  VCFExtractor vin(FLAG_fileName.c_str());
  if (FLAG_range.size())
    vin.setRangeList(FLAG_range.c_str());
  else {
    error("Please provide a range before we can continue.\n");
  };
  if (!FLAG_annoType.empty() && !vcfHasAnnotation(FLAG_fileName)) {
    REprintf("Please use annotated VCF as input (cannot find ANNO in the INFO field);\n");
    REprintf("Prefer using ANNO from https://github.com/zhanxw/anno  \n");
    SEXP ans = R_NilValue;
    return ans;
  }
  
  if (FLAG_annoType.size()) {
    vin.setAnnoType(FLAG_annoType.c_str());
  }
  return readVCF2List(&vin, FLAG_vcfColumn, FLAG_infoTag, FLAG_indvTag);
} // impl_readVCFToListByRange

/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (NOTE: only first one gene is used).
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g. Synonymous|Nonsynonymous)
 * @param arg_columns: a list of which columns to extract (e.g. CHROM, POS ...)
 * @param arg_infoTag: a list of which tag under INFO tag will be extracted (e.g. ANNO, ANNO_FULL, AC ...)
 * @param arg_indvTag: a list of which tag given in individual's column (e.g. GT, GD, GQ ...)
 */
SEXP impl_readVCFToListByGene(SEXP arg_fileName,
                              SEXP arg_geneFile,
                              SEXP arg_geneName,
                              SEXP arg_annoType,
                              SEXP arg_columns,
                              SEXP arg_infoTag,
                              SEXP arg_indvTag){
  // begin
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName,0));
  std::string FLAG_geneFile = CHAR(STRING_ELT(arg_geneFile,0));
  std::string FLAG_geneName = CHAR(STRING_ELT(arg_geneName,0));
  std::string FLAG_annoType = CHAR(STRING_ELT(arg_annoType,0));

  std::set<std::string> FLAG_vcfColumn;
  std::vector<std::string> FLAG_infoTag, FLAG_indvTag;

  extractStringSet(arg_columns, &FLAG_vcfColumn);
  extractStringArray(arg_infoTag, &FLAG_infoTag);
  extractStringArray(arg_indvTag, &FLAG_indvTag);

  if (!FLAG_annoType.empty() && !vcfHasAnnotation(FLAG_fileName)) {
    REprintf("Please use annotated VCF as input (cannot find ANNO in the INFO field);\n");
    REprintf("Prefer using ANNO from https://github.com/zhanxw/anno  \n");
    SEXP ans = R_NilValue;
    return ans;
  }
  
  // Rprintf("vcfColumn.size() = %u\n", FLAG_vcfColumn.size());
  // Rprintf("vcfInfo.size() = %u\n", FLAG_infoTag.size());
  // Rprintf("vcfIndv.size() = %u\n", FLAG_indvTag.size());
  // also append sample names at the end
  int retListLen = FLAG_vcfColumn.size() + FLAG_infoTag.size() + FLAG_indvTag.size() + 1;
  if (retListLen == 0) {
    return R_NilValue;
  }

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

  REprintf( "range = %s\n", range.c_str());
  VCFExtractor vin(FLAG_fileName.c_str());
  if (range.size())
    vin.setRangeList(range.c_str());
  else {
    error("Please provide a valid gene name before we can continue.\n");
  };

  if (FLAG_annoType.size()) {
    vin.setAnnoType(FLAG_annoType.c_str());
  }

  return readVCF2List(&vin, FLAG_vcfColumn, FLAG_infoTag, FLAG_indvTag);
} // end readVCF2List

