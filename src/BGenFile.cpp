#include "BGenFile.h"

#include <cassert>
#include "BitReader.h"
#include "CommonFunction.h"
#include "FileIO.h"
#include "RangeList.h"
#include "TypeConversion.h"

#define UNUSED(x) ((void)(x))

BGenFile::BGenFile(const std::string& fn) : var(this->N), autoMergeRange(false) {
  this->bgenFileName = fn;
  this->mode = BGEN_LINE_MODE;

  int nRead;
  fp = fopen(fn.c_str(), "rb");
  if (fp == NULL) {
    REprintf("Cannot open file [ %s ]\n", fn.c_str());
    // exit(1);
    return;
  }
  // first 4 bytes
  parseUint32(fp, &offset);
#ifdef DEBUG
  printf("nRead = %d, offset = %d\n", nRead, (int)offset);
#endif

  // header block
  parseUint32(fp, &LH);
  parseUint32(fp, &M);
  parseUint32(fp, &N);
#ifdef DEBUG
  Rprintf("nRead = %d, LH = %d\n", nRead, (int)LH);
  Rprintf("nRead = %d, M = %d\n", nRead, (int)M);
  Rprintf("nRead = %d, N = %d\n", nRead, (int)N);
#endif
  sampleMask.resize(N);
  std::fill(sampleMask.begin(), sampleMask.end(), false);

  nRead = fread(&magic, sizeof(magic[0]), 4, fp);
  UNUSED(nRead);
#ifdef DEBUG
  Rprintf("nRead = %d, magic = %c%c%c%c\n", nRead, magic[0], magic[1], magic[2],
          magic[3]);
#endif
  assert(nRead == 4);
  if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' &&
        magic[3] == 'n')) {
    // assert(false);
    REprintf("bgen magic number [ %c%c%c%c ] does not match ['bgen']! \n",
             magic[0], magic[1], magic[2], magic[3]);
  }

  int freeDataLen = LH - 20;
  freeData.resize(freeDataLen);
  nRead = fread(freeData.data(), sizeof(freeData[0]), freeDataLen, fp);
#ifdef DEBUG
  Rprintf("nRead = %d, freeData = ", nRead);
  for (int i = 0; i < freeDataLen; ++i) {
    Rprintf("%c", freeData[i]);
  }
  Rprintf("\n");
#endif

  parseUint32(fp, &flag);
#ifdef DEBUG
  Rprintf("nRead = %d, flag = %zu\n", nRead, flag);
#endif

  snpCompression = static_cast<SNP_COMPRESSION>(flag & 3);
  layout = static_cast<LAYOUT>((flag >> 2) & 0xf);
  flagSampleIdentifier = static_cast<SAMPLE_IDENTIFIER>(flag >> 31);
#ifdef DEBUG
  Rprintf("snpCompression = %d\t", (int)snpCompression);
  switch (snpCompression) {
    case NO_COMPRESSION:
      Rprintf("No compression\n");
      break;
    case GZIP:
      Rprintf("GZIP\n");
      break;
    case ZSTD:
      Rprintf("ZSTD\n");
      break;
    default:
      Rprintf("Wrong value!\n");
  }

  Rprintf("layout= %d\n", (int)layout);
  Rprintf("flagSampleIdentifier = %d %s\n", (int)flagSampleIdentifier,
          (int)flagSampleIdentifier == 1 ? "Have sample id"
                                         : "Do not have sample id");
#endif
  // sample identifier block
  if (flagSampleIdentifier == HAS_SAMPLE_IDENTIFIER) {
    uint32_t LSI;
    parseUint32(fp, &LSI);
#ifdef DEBUG
    Rprintf("nRead = %d, LSI = %d\n", nRead, (int)LSI);
#endif
    if (!(LSI + LH <= offset)) {
      assert(false);
    }

    uint32_t N2;
    parseUint32(fp, &N2);
#ifdef DEBUG
    Rprintf("nRead = %d, N2 = %d\n", nRead, (int)N2);
#endif
    assert(N2 == N);
    sampleIdentifier.resize(N);
    uint16_t siLen;
    for (uint32_t i = 0; i < N; ++i) {
      nRead = fread(&siLen, sizeof(siLen), 1, fp);
      assert(nRead == 1);

      sampleIdentifier[i].resize(siLen);
      nRead = fread(&sampleIdentifier[i][0], sizeof(sampleIdentifier[i][0]),
                    siLen, fp);
#ifdef DEBUG
      Rprintf("nRead = %d, si[%d] = %d, %s\n", nRead, i, (int)siLen,
              sampleIdentifier[i].c_str());
#endif
    }
  }
  buildEffectiveIndex();
  
  // advance fp to variant data blocks
  // NOTE: the assertion below may not work
  // sometimes, unspecified bytes may exists between header block and variant
  // data blocks.
  // long variantOffset = ftell(fp);
  // if (variantOffset != (long)offset + 4) {
  // }
  fseek(fp, offset + 4, SEEK_SET);

  fileSize = getFileSize(fn);
}

bool BGenFile::readRecord() {
  if (mode == BGEN_RANGE_MODE) {
    long int file_pos, bytes;
    if (index.next(&file_pos, &bytes)) {
      fseek(fp, file_pos, SEEK_SET);
    } else {
      return false;
    }
  }

  switch (layout) {
    case LAYOUT1:
      return parseLayout1();
    case LAYOUT2:
      return parseLayout2();
    default:
      assert(false);
  }
  return false;
}

bool BGenFile::parseLayout1() {
  if (isFileEnd(fp)) {
    return false;
  }
  // variant identifying data
  uint32_t NinRow;
  int nRead = fread(&NinRow, sizeof(NinRow), 1, fp);
  UNUSED(nRead);
#ifdef DEBUG
  Rprintf("nRead = %d, NinRow = %d\n", nRead, (int)NinRow);
#endif
  assert(nRead == 1);

  parseString(fp, 2, &var.varid);
  parseString(fp, 2, &var.rsid);
  parseString(fp, 2, &var.chrom);
  parseUint32(fp, &var.pos);

  //      const uint16_t K = 2;  // omitted in layout 1 by specification
  var.K = 2;
  var.alleles.resize(var.K);
  for (int i = 0; i < var.K; ++i) {
    parseString(fp, 4, &var.alleles[i]);
  }

  // genotype data block
  assert(snpCompression ==
         GZIP);  // do not deal with no-compression case for now
  nRead = fread(&C, sizeof(C), 1, fp);
  assert(nRead == 1);
#ifdef DEBUG
  Rprintf("C = %zu\n", C);
#endif

  // std::vector<uint8_t> buf(NinRow * 6);
  D = NinRow * 6;
  buf.resize(D);
  // std::vector<uint8_t> bufCompress(C);
  compressedBuf.resize(C);
  nRead = fread(compressedBuf.data(), sizeof(uint8_t), C, fp);
  assert(nRead == (int)C);
  // unsigned long bufLen = NinRow * 6;
  unsigned long decompressedByte = NinRow * 6;
  int zlibStatus =
      uncompress(buf.data(), &decompressedByte, compressedBuf.data(), C);
  if(zlibStatus != Z_OK) {
    REprintf("decompress zlib failed!\n");
  }

  // parse probility
  var.missing.resize(N);
  var.ploidy.resize(N);
  var.isPhased = false;
  var.index.resize(N + 1);
  var.prob.resize(N * 3);
  float p[3];
  for (size_t i = 0; i < NinRow; ++i) {
    var.ploidy[i] = 2;
    var.index[i] = 3 * i;
    uint16_t* probPtr = (uint16_t*)(buf.data() + i * 6);
    p[0] = (float)(*probPtr) / 32768;
    ++probPtr;
    p[1] = (float)(*probPtr) / 32768;
    ++probPtr;
    p[2] = (float)(*probPtr) / 32768;
#ifdef DEBUG
    Rprintf("prob = (%g, %g, %g)\n", i, p[0], p[1], p[2]);
#endif

    if (p[0] == 0 && p[1] == 0 && p[2] == 0) {
      var.missing[i] = true;
    } else {
      var.missing[i] = false;
    }
    var.prob[i * 3] = p[0];
    var.prob[i * 3 + 1] = p[1];
    var.prob[i * 3 + 2] = p[2];
  }
  var.index.push_back(3 * N);
#ifdef DEBUG
  Rprintf("feof = %d\n", feof(fp));
#endif

  return true;
}

bool BGenFile::parseLayout2() {
  if (isFileEnd(fp)) {
    return false;
  }
  // variant identifying data
  parseString(fp, 2, &var.varid);
  parseString(fp, 2, &var.rsid);
  parseString(fp, 2, &var.chrom);
  parseUint32(fp, &var.pos);

  parseUint16(fp, &var.K);
#ifdef DEBUG
  Rprintf("K = %zu\n", var.K);
#endif

  var.alleles.resize(var.K);
  for (int i = 0; i < var.K; ++i) {
    parseString(fp, 4, &var.alleles[i]);
  }

  uint32_t C;
  parseUint32(fp, &C);
#ifdef DEBUG
  Rprintf("C = %zu\n", C);
#endif

  uint32_t D;
  if (snpCompression == NO_COMPRESSION) {
    D = C;
  } else {
    parseUint32(fp, &D);
#ifdef DEBUG
    Rprintf("D = %zu\n", D);
#endif
  }
  // genotype data block
  std::vector<uint8_t> compressedBuf;
  if (snpCompression == NO_COMPRESSION) {
    compressedBuf.resize(C);
  } else {
    compressedBuf.resize(C - 4);
  }
  std::vector<uint8_t> buf(D);
  size_t nRead;
  UNUSED(nRead);
  if (snpCompression == GZIP) {
    nRead = fread(compressedBuf.data(), sizeof(uint8_t), C - 4, fp);
    assert(nRead == C - 4);
    unsigned long bufLen = D;
    nRead = uncompress(buf.data(), &bufLen, compressedBuf.data(), C);
    assert(nRead == Z_OK);
  } else if (snpCompression == ZSTD) {
    nRead = fread(compressedBuf.data(), sizeof(uint8_t), C - 4, fp);
    assert(nRead == C - 4);
    unsigned long bufLen = D;
    // TODO: create ZSTD context to save time
    size_t ret =
        ZSTD_decompress(buf.data(), bufLen, compressedBuf.data(), C - 4);
    if (ret > bufLen) {
      if (ZSTD_isError(ret)) {
#ifdef DEBUG
        Rprintf("ZSTD error: %s\n", ZSTD_getErrorName(ret));
#endif
      }
    }
    assert(ret == bufLen);
  } else if (snpCompression == NO_COMPRESSION) {
    // nRead = fread(compressedBuf.data(), sizeof(uint8_t), C , fp);
    // assert(nRead == C);
    nRead = fread(buf.data(), sizeof(uint8_t), D, fp);
    assert(nRead == D);
  }

  const uint32_t nIndv = *(uint32_t*)buf.data();
  assert(nIndv == N);
  const uint16_t nAllele = *(uint16_t*)(buf.data() + 4);
  UNUSED(nAllele);
  const uint8_t nMinAllele = *(uint8_t*)(buf.data() + 6);
  UNUSED(nMinAllele);  
  const uint8_t nMaxAllele = *(uint8_t*)(buf.data() + 7);
  UNUSED(nMaxAllele);
  assert(0 <= nMinAllele && nMinAllele <= 63);
  assert(0 <= nMaxAllele && nMaxAllele <= 63);
  const uint8_t* ploidityAndMissing = (uint8_t*)(buf.data() + 8);
  const uint8_t isPhased = *(uint8_t*)(buf.data() + 8 + N);
  var.isPhased = isPhased != 0;
  var.B = *(uint8_t*)(buf.data() + 8 + N + 1);  // bits
  assert(1 <= var.B && var.B <= 32);
#ifdef DEBUG
  int totalBit = 0;
  int cumBit = 0;
#endif

  var.missing.resize(N);
  var.ploidy.resize(N);
  var.index.reserve(N + 1);
  var.index.resize(0);
  var.prob.resize(0);
  BitReader br(buf.data() + 8 + N + 2, (D - 8 - N - 2), var.B);
  for (uint32_t i = 0; i < nIndv; ++i) {
    var.index.push_back(var.prob.size());
    const uint8_t ploidy = ploidityAndMissing[i] & 0x3f;
    const int Z = ploidy;
    const bool missing = (ploidityAndMissing[i] & 0x80) != 0;
    var.ploidy[i] = ploidy;
    var.missing[i] = missing;
#ifdef DEBUG
    Rprintf("ploidy = %d, missing = %s\n", Z, missing ? "true" : "false");
#endif

    // ploidy = Z, allele = K
    if (var.isPhased) {
// total Z * (K- 1) * B bits bits used
#ifdef DEBUG
      totalBit = Z * (var.K - 1) * B;
      Rprintf("Total phased bits = %d\n", totalBit);
#endif

      for (int j = 0; j < ploidy; ++j) {
        float remainProb = 1.0;
        for (int k = 0; k < var.K - 1; ++k) {
          float p = br.next();
          var.prob.push_back(p);
          remainProb -= p;
        }
        var.prob.push_back(remainProb);
      }
    } else {
      // total ( ( (Z+K-1) choose (K-1) ) -1 ) * B bits used
      const int nCombination = choose((Z + var.K - 1), (var.K - 1));
#ifdef DEBUG
      totalBit = (nCombination - 1) * B;
      Rprintf("Total unphased bits = %d\n", totalBit);
#endif
      float remainProb = 1.0;
      for (int j = 0; j < nCombination - 1; ++j) {
        float p = br.next();
        var.prob.push_back(p);
        remainProb -= p;
      }
      var.prob.push_back(remainProb);
    }
#ifdef DEBUG
    cumBit += totalBit;
#endif
  }
  var.index.push_back(var.prob.size());

#ifdef DEBUG
  Rprintf("CumBit = %d\n", cumBit);
  Rprintf("Total chunk = %d\n", 10 + N + cumBit / 8);
  Rprintf("feof = %d\n", feof(fp));
#endif

  return true;
}

void BGenFile::parseString(FILE* fp, int lenByte, std::string* out) {
  if (lenByte == 2) {
    uint16_t siLen;
    int nRead = fread(&siLen, sizeof(siLen), 1, fp);
    UNUSED(nRead);    
    assert(nRead == 1);

    (*out).resize(siLen);
    // (*out)[siLen] = '\0';
    nRead = fread(&(*out)[0], sizeof((*out)[0]), siLen, fp);
#ifdef DEBUG
    Rprintf("nRead = %d, read = %s\n", nRead, (*out).c_str());
#endif

  } else if (lenByte == 4) {
    uint32_t siLen;
    int nRead = fread(&siLen, sizeof(siLen), 1, fp);
    UNUSED(nRead);
    assert(nRead == 1);

    (*out).resize(siLen);
    nRead = fread(&(*out)[0], sizeof((*out)[0]), siLen, fp);
#ifdef DEBUG
    Rprintf("nRead = %d, read = %s\n", nRead, (*out).c_str());
#endif

  } else {
    assert(false);
  }
}

void BGenFile::parseUint32(FILE* fp, uint32_t* value) {
  int nRead = fread(value, sizeof(value[0]), 1, fp);
  UNUSED(nRead);
  assert(nRead == 1);
}
void BGenFile::parseUint16(FILE* fp, uint16_t* value) {
  int nRead = fread(value, sizeof(value[0]), 1, fp);
  UNUSED(nRead);  
  assert(nRead == 1);
}
// @return choose m out of n elements
int BGenFile::choose(int n, int m) {
  if (m == 1) {
    return n;
  }
  if (n == 1) {
    return 1;
  }
  int ret = 1;
  for (int i = 0; i < m; ++i) {
    ret *= (n - i);
  }
  for (int i = 0; i < m; ++i) {
    ret /= (i + 1);
  }
  return ret;
}

bool BGenFile::isFileEnd(FILE* fp) { return (fileSize - ftell(fp)) == 0; }
long BGenFile::getFileSize(const std::string& fn) {
  FILE* fp = fopen(fn.c_str(), "rb");
  fseek(fp, 0, SEEK_END);
  long ret = ftell(fp);
  fclose(fp);
  return ret;
}

//////////////////////////////////////////////////
// Sample inclusion/exclusion
void BGenFile::setPeopleMask(const std::string& s, bool b) {
  for (size_t i = 0; i != sampleIdentifier.size(); ++i) {
    if (sampleIdentifier[i] == s) {
      sampleMask[i] = b;
    }
  }
  buildEffectiveIndex();
}
void BGenFile::includePeople(const std::string& s) { setPeopleMask(s, false); }

void BGenFile::includePeople(const std::vector<std::string>& v) {
  for (size_t i = 0; i != v.size(); ++i) {
    includePeople(v[i].c_str());
  }
}
void BGenFile::setPeopleMaskFromFile(const char* fn, bool b) {
  if (!fn || strlen(fn) == 0) {
    return;
  }
  LineReader lr(fn);
  std::vector<std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")) {
    for (unsigned int i = 0; i < fd.size(); i++) {
      setPeopleMask(fd[i].c_str(), b);
    }
  }
  buildEffectiveIndex();
}
void BGenFile::includePeopleFromFile(const char* fn) {
  setPeopleMaskFromFile(fn, false);
}
void BGenFile::includeAllPeople() {
  std::fill(sampleMask.begin(), sampleMask.end(), false);
  buildEffectiveIndex();
}

void BGenFile::excludePeople(const std::string& s) { setPeopleMask(s, true); }
void BGenFile::excludePeople(const std::vector<std::string>& v) {
  for (size_t i = 0; i != v.size(); ++i) {
    excludePeople(v[i]);
  }
}
void BGenFile::excludePeopleFromFile(const char* fn) {
  setPeopleMaskFromFile(fn, true);
}
void BGenFile::excludeAllPeople() {
  std::fill(sampleMask.begin(), sampleMask.end(), true);
  buildEffectiveIndex();
}

//////////////////////////////////////////////////
// Adjust range collections
void BGenFile::enableAutoMerge() { this->autoMergeRange = true; }
void BGenFile::disableAutoMerge() { this->autoMergeRange = false; }
// void clearRange();
void BGenFile::setRangeFile(const char* fn) {
  if (!fn || strlen(fn) == 0) return;
  RangeList r;
  r.addRangeFile(fn);
  this->setRange(r);
}
// @param l is a string of range(s)
void BGenFile::setRange(const char* chrom, int begin, int end) {
  RangeList r;
  r.addRange(chrom, begin, end);
  this->setRange(r);
}
void BGenFile::setRange(const RangeList& rl) { this->setRangeList(rl); }
void BGenFile::setRangeList(const std::string& l) {
  if (l.empty()) return;

  RangeList r;
  r.addRangeList(l);
  this->setRange(r);
}
// this function the entry point for all function add/change region list
void BGenFile::setRangeList(const RangeList& rl) {
  if (rl.size() == 0) return;

  this->setRangeMode();

  RangeList l;
  l.setRange(rl);
  if (this->autoMergeRange) {
    l.sort();
  }

  if (mode == BGEN_RANGE_MODE) {
    index.setRange(rl);
    // this->tabixReader->setRange(l);
  } else {
    REprintf("[ERROR] invalid reading mode, quitting...\n");
    // abort();
    return;
  }
}

void BGenFile::setRangeMode() {
  if (this->index.init(bgenFileName + ".bgi")) {
    REprintf("Cannot open BGEN index file [ %s ]!\n",
             (bgenFileName + ".bgi").c_str());
    // abort();
    return;
  }
  this->mode = BGEN_RANGE_MODE;
}

int BGenFile::setSiteFile(const std::string& fn) {
  if (fn.empty()) return 0;

  std::vector<std::string> fd;
  LineReader lr(fn);
  int pos;
  std::string chromPos;
  while (lr.readLineBySep(&fd, "\t ")) {
    if (fd.empty()) continue;
    if (fd[0].find(':') != std::string::npos) {
      this->allowedSite.insert(fd[0]);
      continue;
    }
    if (fd.size() >= 2 && str2int(fd[1], &pos) && pos > 0) {
      chromPos = fd[0];
      chromPos += ':';
      chromPos += fd[1];
      this->allowedSite.insert(chromPos);
      continue;
    }
  }
  return 0;
}

int BGenFile::loadSampleFile(const std::string& fn) {
  if (fn.empty()) return 0;

  LineReader lr(fn);
  std::vector<std::string> fd;
  std::vector<std::string> sn;
  int lineNo = 0;
  while (lr.readLineBySep(&fd, " \t")) {
    ++lineNo;
    if (fd.size() < 3) {
      REprintf(
          "Error: Line [ %d ] of the input file [ %s ] has less than 3 "
          "columns!\n",
          lineNo, fn.c_str());
      return -1;
    }
    if (lineNo == 1) {
      if (fd[0] != "ID_1" || fd[1] != "ID_2" || fd[2] != "missing") {
        REprintf(
            "ERROR: The header line does not start with ID_1, ID_2 and "
            "missing!\n");
        return -1;
      }
    } else if (lineNo == 2) {
      if (fd[0] != "0" || fd[1] != "0" || fd[2] != "0") {
        REprintf(
            "ERROR: The second line (variable type line) does not start "
            "with 0, 0 and 0!\n");
        return -1;
      }
    } else {
      if (fd[0] != fd[1]) {
        REprintf(
            "WARN: Line [ %d ], ID_1 [ %s ] is different than ID_2 [ %s ], "
            "by default ID_2 column is used.\n",
            lineNo, fd[0].c_str(), fd[1].c_str());
      }
      sn.push_back(fd[1]);
    }
  }
  REprintf("INFO: Loaded [ %d ] samples from .sample file [ %s ]\n",
           (int)sn.size(), fn.c_str());
  if (sn.size() != N) {
    REprintf(
        "ERROR: Sample file has [ %d ] samples, but BGEN file have [ %d ] "
        "samples, skipped loading this sample file\n",
        (int)sn.size(), N);
    return -1;
  } else {
    this->sampleIdentifier = sn;
    return 0;
  }
}

int BGenFile::getNumEffectiveSample() const {
  size_t ret = 0;
  for (size_t i = 0; i != sampleMask.size(); ++i) {
    if (sampleMask[i]) continue;
    ret++;
  }
  return ret;
}

void BGenFile::getIncludedSampleName(std::vector<std::string>* p) const {
  if (!p) return;
  p->clear();
  for (size_t i = 0; i != sampleMask.size(); ++i) {
    if (sampleMask[i]) continue;
    p->push_back(sampleIdentifier[i]);
  }
}

void BGenFile::buildEffectiveIndex() {
  effectiveIndex.resize(0);
  for (size_t i = 0; i != N; ++i) {
    if (sampleMask[i]) continue;
    effectiveIndex.push_back(i);
  }
  effectiveIndex.push_back(N);
}

int BGenFile::getEffectiveIndex(int idx) const {
  return this->effectiveIndex[idx];
}

void BGenFile::printInfo() {
  Rprintf("--First 4 bytes--\n");
  Rprintf("\toffset = %d\n", (int)offset);

  // header block
  Rprintf("--Header block--\n");
  Rprintf("\tLH = %d\n", (int)LH);
  Rprintf("\tM = %d\n", (int)M);
  Rprintf("\tN = %d\n", (int)N);

  if (freeData.empty()) {
    Rprintf("\tfreeData = []\n");
  } else {
    Rprintf("\tfreeData = [");
    for (size_t i = 0; i < freeData.size(); ++i) {
      Rprintf("%c", freeData[i]);
    }
    Rprintf("]\n");
  }
  Rprintf("\tflag = %x\n", flag);
  Rprintf("\tsnpCompression = %d ", (int)snpCompression);
  switch (snpCompression) {
    case NO_COMPRESSION:
      Rprintf("(No compression)\n");
      break;
    case GZIP:
      Rprintf("(GZIP)\n");
      break;
    case ZSTD:
      Rprintf("(ZSTD)\n");
      break;
    default:
      Rprintf("Wrong value!\n");
  }

  Rprintf("\tlayout= %d\n", (int)layout);
  Rprintf("\tflagSampleIdentifier = %d %s\n", (int)flagSampleIdentifier,
          (int)flagSampleIdentifier == 1 ? "(Have sample id)"
                                         : "(Do not have sample id)");

  // sample identifier block
  if (flagSampleIdentifier == HAS_SAMPLE_IDENTIFIER) {
    Rprintf("--Sample identifier block--\n");
    // Rprintf("LSI = %d\n", (int)LSI);
    // Rprintf("N2 = %d\n", (int)N2);
    // assert(N2 == N);
    for (uint32_t i = 0; i < N; ++i) {
      Rprintf("\tsi[%d] = %s\n", i, sampleIdentifier[i].c_str());
    }
  }

  // variant data blocks
  Rprintf("--Variant data block--\n");
  for (size_t i = 0; i < M; ++i) {
    if (!readRecord()) {
      Rprintf("\tNo variants presented, file truncated?\n");
      break;
    }

    Rprintf("\n\t[Variant %d]\n", (int)i);

    Rprintf("\tChromosomal position: %s:%d\n", var.chrom.c_str(), var.pos);
    Rprintf("\tRSID = %s\n", var.rsid.c_str());
    Rprintf("\tVarID = %s\n", var.varid.c_str());
    Rprintf("\tAlleles = %s ", var.alleles[0].c_str());
    for (size_t j = 1; j != var.alleles.size(); ++j) {
      Rprintf("/ %s ", var.alleles[j].c_str());
    }
    Rprintf("\n");

    if (layout == LAYOUT1) {
      Rprintf("\tPolidy = 2\n");
      Rprintf("\tUnphased\n");
      Rprintf("\tAlleles = 2\n");
      Rprintf("\tBitsPerGenotypeProbability = 16\n");  // 2 bytes per genotype
                                                       // probability
      int nMissing = 0;
      for (size_t j = 0; j != N; ++j) {
        if (var.prob[j * 3] == 0.0 && var.prob[j * 3 + 1] == 0.0 &&
            var.prob[j * 3 + 2] == 0.0)
          ++nMissing;
      }
      Rprintf("Missing = %d\t", nMissing);
    } else if (layout == LAYOUT2) {
      int nMissing = 0;
      for (size_t i = 0; i < var.missing.size(); ++i) {
        if (var.missing[i]) ++nMissing;
      }
      std::set<uint8_t> s = makeSet(var.ploidy);
      std::string ss = toString(s, ",");
      Rprintf("\tPolidy = %s\n", ss.c_str());
      Rprintf("\t%s\n", var.isPhased ? "Phased" : "Unphased");
      Rprintf("\tAlleles = %d\n", var.K);
      Rprintf("\tBitsPerGenotypeProbability = %d\n", (int)var.B);
      Rprintf("\tMissing = %d\n", nMissing);
    } else {
      assert(false);
    }
  }
}
