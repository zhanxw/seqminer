#include "SingleChromosomeVCFIndex.h"

#define STRICT_R_HEADERS
#include <R.h>

// #include <stdlib.h>
#include <algorithm>
#include <string>
#include <vector>

#include "MmapFile.h"
#include "StringUtil.h"
#include "bgzf.h"

struct Record {
  int64_t pos;
  int64_t offset;
  bool operator<(Record& o) { return (pos < o.pos); }
};

static bool comparator(const Record& a, const Record& b) { return a.pos < b.pos; }

SingleChromosomeVCFIndex::SingleChromosomeVCFIndex(
    const std::string& vcfFile, const std::string& indexFile) {
  vcfFile_ = vcfFile;
  indexFile_ = indexFile;
  fVcfFile_ = bgzf_open(vcfFile_.c_str(), "rb");
  data_ = NULL;
  mmapFile_ = NULL;
  str = (kstring_t*)calloc(1, sizeof(kstring_t));
}

SingleChromosomeVCFIndex::~SingleChromosomeVCFIndex() { this->close(); }

void SingleChromosomeVCFIndex::close() {
  if (str) {
    free(str);
    str = NULL;
  }
  if (fVcfFile_) {
    bgzf_close(fVcfFile_);
    fVcfFile_ = NULL;
  }
  closeIndex();
}

void SingleChromosomeVCFIndex::closeIndex() {
  if (mmapFile_) {
    delete mmapFile_;
    data_ = NULL;
  }
  if (data_) {
    delete[](uint8_t*) data_;
    data_ = NULL;
  }
}

/**
 * Create single chromosome index file
 * the file content is a 2-column matrix of int64_t type
 * line1:  num_sample  num_marker
 * line2:  0           bgzf_offset_for_#CHROM_line
 * line3:  var_1_pos   bgzf_offset_for_var_1
 * ...
 */
int SingleChromosomeVCFIndex::createIndex() {
  const char* fn = vcfFile_.c_str();
  BGZF* fp = fVcfFile_;  // bgzf_open(fn, "rb");
  bgzf_seek(fp, 0, SEEK_SET);
  // kstring_t* str;
  // str = (kstring_t*)calloc(1, sizeof(kstring_t));
  kstring_t& s = *str;
  FILE* fIndex = fopen(indexFile_.c_str(), "wb");
  int ret;

  int64_t numSample = 0;
  int64_t numMarker = 0;
  int64_t pos = -1;
  int64_t offset = -1;

  std::string line;
  std::vector<std::string> fd;

  // will rewrite with N (sample size) and M (number of marker) later
  fwrite(&numSample, sizeof(int64_t), 1, fIndex);
  fwrite(&numMarker, sizeof(int64_t), 1, fIndex);

  do {
    offset = bgzf_tell(fp);
    // REprintf("offset = %ld\n", offset);
    ret = bgzf_getline(fp, '\n', &s);
    if (ret <= 0) {  // file end or error
      break;
    }
    // process header
    if (s.s[0] == '#') {
      if (s.s[1] == '#') {  // begin with '##', header line, skip
        continue;
      } else if (s.s[1] == 'C') {  // likely '#CHROM' line
        line = s.s;
        stringTokenize(line, '\t', &fd);
        numSample =
            fd.size() - 9;  // 9 is CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
        REprintf("header line has %d samples\n", numSample);
        pos = 0;
        fwrite(&pos, sizeof(int64_t), 1, fIndex);
        fwrite(&offset, sizeof(int64_t), 1, fIndex);
        REprintf("offset = %ld\n", offset);
        continue;
      } else {
        REprintf("Strange header line!\n");
      }
    }
    // find POS
    size_t beg = 0;
    for (; beg < s.l; ++beg) {
      if (s.s[beg] == '\t') {
        pos = strtol(s.s + beg + 1, NULL, 0);
        break;
      }
    }

    // printf("%ld %ld\n", pos, offset);
    fwrite(&pos, sizeof(int64_t), 1, fIndex);
    fwrite(&offset, sizeof(int64_t), 1, fIndex);
    numMarker++;

  } while (1);

  if (fseek(fIndex, 0, SEEK_SET)) {
    REprintf("fseek failed\n!");
  }
  fwrite(&numSample, sizeof(int64_t), 1, fIndex);
  fwrite(&numMarker, sizeof(int64_t), 1, fIndex);

  // bgzf_close(fp);
  fclose(fIndex);

  REprintf("Indexing finished with %d samples and %d markers\n", numSample,
           numMarker);
  return 0;
}

int SingleChromosomeVCFIndex::openIndex() {
  closeIndex();
  // read everything
  size_t fsize = getFileSize(indexFile_.c_str());
  REprintf("fsize = %ld\n", (long int)fsize);
  data_ = new uint8_t[fsize];
  FILE* fp = fopen(indexFile_.c_str(), "rb");
  if (fread(data_, sizeof(uint8_t), fsize, fp) != fsize) {
    REprintf("Read incomplete index\n");
    return -1;
  }

  // verify file integrity
  int64_t* d = (int64_t*) data_;
  if (fsize != sizeof(Record)  * (2L + d[1])) {
    REprintf("Check file integrity!\n");
    REprintf("d = %ld %ld fsize = %ld\n", d[0], d[1], (int)fsize);
    return -1;
  }
  return 0;
}

int SingleChromosomeVCFIndex::mapIndex() {
  closeIndex();
  // read everything
  mmapFile_ = new MmapFile;
  MmapFile& mmapFile = *mmapFile_;
  if (mmapFile.open(indexFile_.c_str())) {
    return -1;
  }
  size_t Nrecord = mmapFile.getFileSize() / 16 -
                   1;  // -1: skip first block, 16: two bytes for uint64_t
  data_ = mmapFile.data;

  // verify file integrity
  int64_t* d = (int64_t*)data_;
  if (Nrecord != d[1]) {
    REprintf("Check file integrity!\n");
    return -1;
  }
  return 0;
}


int SingleChromosomeVCFIndex::query(int chromPos, int64_t* pVirtualOffset) {
  return this->query(chromPos, chromPos, pVirtualOffset);
}

int SingleChromosomeVCFIndex::query(int chromPosBeg, int chromPosEnd,
                                    int64_t* voffset) {
  if (!data_) {
    REprintf("open index first!\n");
    return -1;
  }

  if (!voffset) {
    return -1;
  }
  REprintf("query [%d, %d]\n", chromPosBeg, chromPosEnd);

  Record* r = (Record*)data_;
  const int64_t Nrecord = r[0].offset;
  
  ++r;  // skip the first block, as first block is (#sample, #marker)

  // binary search for file position
  *voffset = -1;
  Record query;
  query.pos = chromPosBeg;
  // Comparator comparator;
  Record* lb =
      std::lower_bound(r, r + Nrecord + 1, query,
                       comparator);  // r[lb].pos >= query.pos = chromPosBeg
  query.pos = chromPosEnd;
  Record* ub =
      std::upper_bound(lb, r + Nrecord + 1, query,
                       comparator);  // r[ub].pos > query.pos = chromPosEnd
  REprintf("Found %d results\n", (int)(ub - lb));
  for (Record* pi = lb; pi != ub; ++pi) {
    // REprintf("%ld %ld\n", pi->pos, pi->offset);
    *voffset = lb->offset;
    break;
  }

  if (*voffset < 0) {
    REprintf("Cannot find position!\n");
    return -1;
  } else {
    REprintf("found %d position, e.g. %ld %ld\n", (int)(ub - lb), (*lb).pos,
             (*lb).offset);
    return ub - lb;
  }
}

int SingleChromosomeVCFIndex::readLine(int64_t offset, std::string* line) {
  if (bgzf_seek(fVcfFile_, offset, SEEK_SET)) {
    REprintf("seek error!\n");
  }
  kstring_t& s = *str;
  int ret = bgzf_getline(fVcfFile_, '\n', &s);
  if (ret <= 0) {
    REprintf("getline error, ret = %d!\n", ret);
  }
  // for (size_t i = 0; i < s.l; ++i) {
  //   if (i >= 50) break;
  //   REprintf("%c", s.s[i]);
  // }
  // REprintf("\n");

  *line = s.s;

  return s.l;
}

int SingleChromosomeVCFIndex::nextLine(std::string* line) {
  kstring_t& s = *str;
  int ret = bgzf_getline(fVcfFile_, '\n', &s);
  if (ret <= 0) {
    REprintf("getline error, ret = %d!\n", ret);
  }
  // for (size_t i = 0; i < s.l; ++i) {
  //   if (i >= 50) break;
  //   REprintf("%c", s.s[i]);
  // }
  // REprintf("\n");

  *line = s.s;
  return s.l;
}
