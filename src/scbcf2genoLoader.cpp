#include "scbcf2genoLoader.h"

#include "BCFHeader.h"
#include "R_CPP_interface.h"
#include "RangeList.h"
#include "SingleChromosomeBCFIndex.h"

int parseBCFVariant(const BCFHeader& bcfHeader,
                    const uint32_t l_shared, const uint32_t l_indiv,
                    const std::vector<char> line,
                    std::vector<double>* buf,
                    std::vector<std::string>* markerNames);

SEXP impl_readSingleChromosomeBCFToMatrixByRange(SEXP arg_fileName,
                                                 SEXP arg_indexFileName,
                                                 SEXP arg_range) {
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  std::string FLAG_indexFileName = CHAR(STRING_ELT(arg_indexFileName, 0));
  std::vector<std::string> FLAG_range;
  extractStringArray(arg_range, &FLAG_range);
  SingleChromosomeBCFIndex sc(FLAG_fileName, FLAG_indexFileName);
  if (sc.openIndex()) {
    REprintf("failed to open index!\n");
  }
  SEXP ans = R_NilValue;

  BGZF* fp = bgzf_open(FLAG_fileName.c_str(), "rb");
  if (!fp) {
    REprintf("Cannot open file: %s\n", FLAG_fileName.c_str());
    return ans;
  }
  bgzf_seek(fp, 0, SEEK_SET);

  // check magic number
  char magic[5];
  if (5 != bgzf_read(fp, magic, 5)) {
    REprintf("Encounted fatal error!\n"); return ans; // exit(1);
  }
  if (!(magic[0] == 'B' && magic[1] == 'C' && magic[2] == 'F' &&
        magic[3] == 2 && (magic[4] == 1 || magic[4] == 2))) {
    REprintf("Encounted fatal error!\n"); return ans; // exit(1);
  }

  // read header to get sample names
  uint32_t l_text;
  if (4 != bgzf_read(fp, &l_text, 4)) {
    REprintf("Encounted fatal error!\n"); return ans; // exit(1);
  }
  Rprintf("l_text = %d\n", l_text);

  std::string s;
  // int64_t bgzf_offset_before_header = bgzf_tell(fp); // the beginning of header block
  s.resize(l_text);
  if (bgzf_read(fp, (void*)s.data(), l_text) != l_text) {
    REprintf( "Read failed!\n");
  }
  BCFHeader bcfHeader;
  if (bcfHeader.parseHeader(s,
                  &bcfHeader.header_contig_id,
                  &bcfHeader.header_id,
                  &bcfHeader.header_number,
                  &bcfHeader.header_type,
                  &bcfHeader.header_description)) {
    REprintf( "Parse header failed!\n"); REprintf("Encounted fatal error!\n"); return ans; // exit(1);
  }

  // locate #CHROM line
  size_t ptr_chrom_line = s.find("#CHROM"); // the index of "#CHROM", also the size between beginning of header to '#CHROM'
  if (ptr_chrom_line == std::string::npos) {
    REprintf( "Cannot find the \"#CHROM\" line!\n");
    REprintf("Encounted fatal error!\n"); return ans; // exit(1);
  }
  s = s.substr(ptr_chrom_line, s.size() - ptr_chrom_line);
  // load sample names
  while (s.back() == '\n' || s.back() == '\0') {
    s.resize(s.size() - 1);
  } 
  stringTokenize(s, "\t", &bcfHeader.sample_names);
  bcfHeader.sample_names.erase(bcfHeader.sample_names.begin(), bcfHeader.sample_names.begin() + 9);
  
  const int numSample = (int)bcfHeader.sample_names.size() ; // vcf header has 9 columns CHROM...FORMAT before actual sample names
  const std::vector<std::string>& sampleNames = bcfHeader.sample_names;
  REprintf("Inferred %d samples from header\n", numSample);

  int nGene = FLAG_range.size();
  Rprintf("%d region to be extracted.\n", nGene);
  PROTECT(ans = allocVector(VECSXP, nGene));
  setListNames(FLAG_range, &ans);

  std::string chromName;
  unsigned int chromPosBeg;
  unsigned int chromPosEnd;
  std::vector<std::string> ranges;
  int64_t offset;
  std::vector<char> line;
  for (int i = 0; i < nGene; ++i) {
    stringTokenize(FLAG_range[i], ',', &ranges);

    std::vector<double> buf;
    int cumNumVariant = 0;
    std::vector<std::string> markerNames;
    uint32_t l_shared;
    uint32_t l_indiv;

    for (size_t j = 0; j != ranges.size(); ++j) {
      parseRangeFormat(ranges[j], &chromName, &chromPosBeg, &chromPosEnd);
      // note, bcf index is 0-based, but query is 1-based
      --chromPosBeg;
      --chromPosEnd;
      int nVariant = sc.query(chromPosBeg, chromPosEnd, &offset);
      if (nVariant <= 0) {
        REprintf("Cannot find the variant!\n");
        continue;
      }

      // read and process a bcf block for a variant
      if (sc.readLine(offset, &l_shared, &l_indiv, &line) < 0) {
        REprintf("Cannot readline()!\n");
      }
      parseBCFVariant(bcfHeader, l_shared, l_indiv, line, &buf, &markerNames);
      for (int remainVariant = nVariant - 1; remainVariant > 0;
           --remainVariant) {
        if (sc.nextLine(&l_shared, &l_indiv, &line) < 0) {
          REprintf("Cannot readline()!\n");
        }
        parseBCFVariant(bcfHeader, l_shared, l_indiv, line, &buf, &markerNames);
      }
      cumNumVariant += nVariant;
    }
    SEXP val;
    PROTECT(val = allocVector(REALSXP, numSample * cumNumVariant));
    memcpy(REAL(val), buf.data(), sizeof(double) * numSample * cumNumVariant);
    setDim(numSample, cumNumVariant, val);
    Rprintf("sampleNames.size() = %d, markerNames.size() = %d\n", (int) sampleNames.size(), (int) markerNames.size());
    // for(int i = 0; i < sampleNames.size(); ++i) {
    //   Rprintf("%d = %s\n", i, sampleNames[i].c_str());
    // }
    setDimNames(sampleNames, markerNames, val);
    UNPROTECT(1);
    SET_VECTOR_ELT(ans, i, val);
  }
  if (fp) {bgzf_close(fp);}
  UNPROTECT(1);
  return ans;
}

SEXP impl_createSingleChromosomeBCFIndex(SEXP arg_fileName,
                                         SEXP arg_indexFileName) {
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName, 0));
  std::string FLAG_indexFileName = CHAR(STRING_ELT(arg_indexFileName, 0));

  SingleChromosomeBCFIndex sc(FLAG_fileName, FLAG_indexFileName);
  if (sc.createIndex()) {
    REprintf("create index file successfully!\n");
  }
  REprintf("created index file [ %s ]\n", FLAG_indexFileName.c_str());
  return arg_indexFileName;
}

// read in one integer (int8_t, int16_t, or int32_t)
// @return -1 if error happens, or the number of actual bytes used
int readOneInteger(const char* fp, int* len) {
  int retVal = 0;
  int nRead = 0;
  uint8_t val_type = *fp;
  nRead = 1;
  retVal += nRead;
  fp += nRead;
  // read the next integer
  // if (1 != bgzf_read(fp, &val_type, 1)) {
  //   fprintf(stderr, "Wrong read!\n"); REprintf("Encounted fatal error!\n"); return ans; // exit(1);
  // }
  switch (val_type & 0x0F) {
    case 1:  // 8 bit int
      int8_t tmp8;
      memcpy(&tmp8, fp, sizeof(int8_t));
      nRead = sizeof(int8_t);
      *len = tmp8;
      // if (1 != bgzf_read(fp,  &tmp8, 1)) {
      //   *len = tmp8;
      // }
      break;
    case 2: // 16 bit int
      int16_t tmp16;
      memcpy(&tmp16, fp, sizeof(int16_t));
      nRead = sizeof(int16_t);
      *len = tmp16;
      // if (1 != bgzf_read(fp,  &tmp16, 1)) {
      //   *len = tmp16;
      // }
      break;
    case 3: // 32 bit int
      int32_t tmp32;
      memcpy(&tmp32, fp, sizeof(int32_t));
      nRead = sizeof(int32_t);
      *len = tmp32;
      // if (1 != bgzf_read(fp,  &tmp32, 1)) {
      //   *len = tmp32;
      // }
      break;
    default:
      REprintf("Wrong type!\n"); REprintf("Encounted fatal error!\n"); return retVal; // exit(1);
  }
  retVal += nRead;
  fp += nRead;
  if (val_type >> 4  != 1) {
    REprintf("Wrong array dimension!\n"); REprintf("Encounted fatal error!\n"); return retVal; // exit(1);
  }
  return retVal;
}

// read one or two bytes of given @param type from @param fp, and report array $param len
// @return -1 if error happens, or the number of actual bytes used
int readArray(const char* fp, const int type, int* len) {
  int retVal = 0;
  int nRead = 0;
  uint8_t val_type = *fp;
  nRead = 1;
  retVal += nRead;
  fp += nRead;
  // if (1 != bgzf_read(fp, &val_type, 1)) {
  //   fprintf(stderr, "Wrong read!\n"); REprintf("Encounted fatal error!\n"); return ans; // exit(1);
  // }
  if ( (val_type & 0x0F) != type) {
    REprintf("Wrong type %d != %d!\n", val_type & 0x0F, type); REprintf("Encounted fatal error!\n"); return retVal; // exit(1);
  }
  uint8_t val_len = (val_type >> 4);
  if (val_len == 0) { // missing
    *len = 0;
  } else if (val_len == 15) { // overflowed
    nRead = readOneInteger(fp, len);
    retVal += nRead;
    fp += nRead;
  } else {
    *len = val_len;
  }
  return retVal;
}

// @return -1 if error happens, or the number of actual bytes used
int readInt(const char* fp, std::vector<int8_t>* ret) {
  int retVal = 0;
  int nRead;
  int len;
  nRead = readArray(fp, 1, &len); // 1 means 8bit integer
  if (nRead < 0) { 
    REprintf("Wrong read array!\n"); REprintf("Encounted fatal error!\n"); return retVal; // exit(1);
  }
  retVal += nRead;
  fp += nRead;
  // Rprintf("len of int = %d\n", len);
  ret->resize(len);
  memcpy((void*)ret->data(), fp, len*sizeof(int8_t));
  nRead = len * sizeof(int8_t);
  fp += nRead;
  retVal += nRead;
  // if ((ssize_t)(len*sizeof(int8_t)) != bgzf_read(fp, ) ) {
  //   fprintf(stderr, "Wrong read string!\n"); REprintf("Encounted fatal error!\n"); return ans; // exit(1);
  // }
  return retVal;
}

int parseBCFVariant(const BCFHeader& bcfHeader,
                    const uint32_t l_shared, const uint32_t l_indiv,
                    const std::vector<char> line,
                    std::vector<double>* buf,
                    std::vector<std::string>* markerNames) {
  const size_t sampleSize = bcfHeader.sample_names.size();
  const char* p = line.data();
  const char* pIndv = p + l_shared;

  // parse from pIndv and expected to get GT as the first fiedl
  // 1. check GT tag
  std::vector<int8_t> format_key;
  pIndv += readInt(pIndv, &format_key);
  if (bcfHeader.header_id[format_key[0]] != "GT") {
    REprintf("The first element in FORMAT is not GT!\n");
    return -1;
  }
  // Rprintf("format key = %d [%s], number = [%s], type = [%s], description = [%s] \n",
  //        format_key[0],
  //        bcfHeader.header_id[format_key[0]].c_str(),
  //        bcfHeader.header_number[format_key[0]].c_str(),
  //        bcfHeader.header_type[format_key[0]].c_str(),
  //        bcfHeader.header_description[format_key[0]].c_str());
  // 2. read GT type
  int8_t format_type = *pIndv;
  pIndv ++;
  // Rprintf("format type = 0x%0x ", format_type);
  // int format_len_per_indv = (format_type >> 4) * // (num of types per indv)
  //     ( format_type & ( (1<<4) - 1) );                // (bytes per type, e.g. 1 for int8_t, which is 1 byte
  // // Rprintf("format len per indv = %d\n", format_len_per_indv);

  // 3. read genotypes
  for (size_t i = 0; i != sampleSize; ++i) {
    buf->push_back( ((*pIndv) >> 1) - 1);
    pIndv++;
    buf->back() += ( ((*pIndv) >> 1) - 1);
    pIndv++;
    if (buf->back() < 0) {
      buf->back() = -9;
    }
  }
  
  int32_t pos; // 0-based position
  memcpy(&pos, p+4, 4);
  markerNames->push_back(toString(pos + 1)); // report 1-based position

  return 0;
}

