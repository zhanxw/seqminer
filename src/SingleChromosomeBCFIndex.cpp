#include "SingleChromosomeBCFIndex.h"

#define STRICT_R_HEADERS
#include <R.h>

// #include <stdlib.h>
#include <algorithm>
#include <string>
#include <vector>

#include "BCFHeader.h"
#include "MmapFile.h"
#include "StringUtil.h"
#include "bgzf.h"

struct Record {
  int64_t pos;
  int64_t offset;
  bool operator<(Record& o) { return (pos < o.pos); }
};

static bool comparator(const Record& a, const Record& b) {
  return a.pos < b.pos;
}

SingleChromosomeBCFIndex::SingleChromosomeBCFIndex(
    const std::string& bcfFile, const std::string& indexFile) {
  bcfFile_ = bcfFile;
  indexFile_ = indexFile;
  fBcfFile_ = bgzf_open(bcfFile_.c_str(), "rb");
  data_ = NULL;
}

SingleChromosomeBCFIndex::~SingleChromosomeBCFIndex() { this->close(); }

void SingleChromosomeBCFIndex::close() {
  if (fBcfFile_) {
    bgzf_close(fBcfFile_);
    fBcfFile_ = NULL;
  }
  closeIndex();
}

void SingleChromosomeBCFIndex::closeIndex() {
  if (data_) {
    delete[] (uint8_t*)data_;
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
int SingleChromosomeBCFIndex::createIndex() {
  // const char* fn = bcfFile_.c_str();
  BGZF* fp = fBcfFile_;  // bgzf_open(fn, "rb");
  bgzf_seek(fp, 0, SEEK_SET);

  // check magic number
  char magic[5];
  if (5 != bgzf_read(fp, magic, 5)) {
    return -1;  // exit(1);
  }
  if (!(magic[0] == 'B' && magic[1] == 'C' && magic[2] == 'F' &&
        magic[3] == 2 && (magic[4] == 1 || magic[4] == 2))) {
    return -1;  // exit(1);
  }

  // read header
  uint32_t l_text;
  if (4 != bgzf_read(fp, &l_text, 4)) {
    return -1;  // exit(1);
  }
  Rprintf("l_text = %d\n", l_text);

  std::string s;
  int64_t bgzf_offset_before_header =
      bgzf_tell(fp);  // the beginning of header block
  s.resize(l_text);
  if (bgzf_read(fp, (void*)s.data(), l_text) != l_text) {
    REprintf("Read failed!\n");
  }
  BCFHeader bcfHeader;
  if (bcfHeader.parseHeader(s, &bcfHeader.header_contig_id,
                            &bcfHeader.header_id, &bcfHeader.header_number,
                            &bcfHeader.header_type,
                            &bcfHeader.header_description)) {
    REprintf("Parse header failed!\n");
    return -1;  // exit(1);
  }

  // locate #CHROM line
  int64_t bgzf_offset_after_header = bgzf_tell(fp);  // the end of header block
  size_t ptr_chrom_line =
      s.find("#CHROM");  // the index of "#CHROM", also the size between
                         // beginning of header to '#CHROM'
  if (ptr_chrom_line == std::string::npos) {
    REprintf("Cannot find the \"#CHROM\" line!\n");
    return -1;  // exit(1);
  }
  Rprintf("offset_header = %d\n", (int)ptr_chrom_line);

  bgzf_seek(fp, bgzf_offset_before_header,
            SEEK_SET);  // rewind fp to the beginning of header
  s.resize(ptr_chrom_line);
  int64_t before_chrom_size = bgzf_read(fp, (void*)s.data(), ptr_chrom_line);
  int64_t bgzf_offset_before_chrom = bgzf_tell(fp);  // the offset to #CHROM
  s.resize(l_text - before_chrom_size);
  int64_t after_chrom_size =
      bgzf_read(fp, (void*)s.data(), l_text - before_chrom_size);
  int32_t last_character = s[after_chrom_size - 1];
  // load sample names
  while (s.back() == '\n' || s.back() == '\0') {
    s.resize(s.size() - 1);
  }
  stringTokenize(s, "\t", &bcfHeader.sample_names);
  const int64_t num_sample =
      (int)bcfHeader.sample_names.size() -
      9;  // vcf header has 9 columns CHROM...FORMAT before actual sample names
#if defined(__MINGW64__) || ( defined(__APPLE__) && defined(__arm64__))
  Rprintf("sample size = %g\n", (double) num_sample);
#else
  Rprintf("sample size = %ld\n", (unsigned long int)num_sample);
#endif
  Rprintf("last character is s[after_chrom_size-1] = %d\n",
          last_character);  // should be 0, the null terminator character
  // quality check
  if (bgzf_offset_after_header != bgzf_tell(fp)) {
    REprintf("Messed up bgzf header\n");
    return -1;  // exit(1);
  }

  // create index file
  FILE* fIndex = fopen(indexFile_.c_str(), "wb");
  int64_t num_marker = 0;
  int64_t pos = 0;
  fwrite(&num_sample, sizeof(int64_t), 1, fIndex);
  fwrite(&num_marker, sizeof(int64_t), 1, fIndex);
  fwrite(&pos, sizeof(int64_t), 1, fIndex);
  fwrite(&bgzf_offset_before_chrom, sizeof(int64_t), 1, fIndex);

  uint32_t l_shared;
  uint32_t l_indiv;
  std::vector<char> data;
  int64_t offset;
  do {
    offset = bgzf_tell(fp);
    if (4 != bgzf_read(fp, &l_shared, sizeof(uint32_t))) {
      break;  // REprintf( "Wrong read!\n"); exit(1);
    }
    if (4 != bgzf_read(fp, &l_indiv, sizeof(uint32_t))) {
      break;  // REprintf( "Wrong read!\n"); exit(1);
    }
    data.resize(l_shared + l_indiv);
    if (l_shared + l_indiv !=
        bgzf_read(fp, data.data(), (l_shared + l_indiv) * sizeof(char))) {
      break;  // REprintf( "Wrong read!\n"); exit(1);
    }
    memcpy(&pos, data.data() + 4, 4);
    fwrite(&pos, sizeof(int64_t), 1, fIndex);
    fwrite(&offset, sizeof(int64_t), 1, fIndex);

    num_marker++;
    if (num_marker % 10000 == 0) {
#if defined(__MINGW64__) || ( defined(__APPLE__) && defined(__arm64__))
      Rprintf("\rprocessed %g markers", (double) num_marker);
#else
      Rprintf("\rprocessed %ld markers", (unsigned long int)num_marker);
#endif
    }
  } while (true);

  if (fseek(fIndex, 0, SEEK_SET)) {
    REprintf("fseek failed\n!");
  }
  fwrite(&num_sample, sizeof(int64_t), 1, fIndex);
  fwrite(&num_marker, sizeof(int64_t), 1, fIndex);
  fclose(fIndex);
#if defined(__MINGW64__) || ( defined(__APPLE__) && defined(__arm64__))
  Rprintf("Indexing finished with %g samples and %g markers\n", 
  (double) num_sample,
  (double) num_marker);
#else
  Rprintf("Indexing finished with %ld samples and %ld markers\n",
          (unsigned long int)num_sample, (unsigned long int)num_marker);
#endif
  return 0;
}

int SingleChromosomeBCFIndex::openIndex() {
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
  int64_t* d = (int64_t*)data_;
  if (fsize !=
      sizeof(Record) *
          (2L + d[1])) {  // d[0, 1]: number of sample; number of marker
    REprintf("Check file integrity!\n");
#if defined(__MINGW64__) || ( defined(__APPLE__) && defined(__arm64__))
    REprintf("d = %g %g fsize = %g\n", (double) d[0], (double) d[1], (double) fsize);
#else
    REprintf("d = %ld %ld fsize = %ld\n", (unsigned long int)d[0],
             (unsigned long int)d[1], (long int)fsize);
#endif
    return -1;
  }
  return 0;
}

int SingleChromosomeBCFIndex::query(int chromPos, int64_t* pVirtualOffset) {
  return this->query(chromPos, chromPos, pVirtualOffset);
}

int SingleChromosomeBCFIndex::query(int chromPosBeg, int chromPosEnd,
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
#if defined(__MINGW64__) || ( defined(__APPLE__) && defined(__arm64__))
    REprintf("%g %g\n", (double) pi->pos, (double) pi->offset);
#else
    REprintf("%ld %ld\n", (unsigned long int)pi->pos,
             (unsigned long int)pi->offset);
#endif
    // REprintf("%ld %ld\n", ub->pos, ub->offset);
    *voffset = lb->offset;
    break;
  }
  if (*voffset < 0) {
    REprintf("Cannot find position!\n");
    return -1;
  } else {
#if defined(__MINGW64__) || ( defined(__APPLE__) && defined(__arm64__))
    REprintf("found %d position, e.g. %g %g\n", (int)(ub - lb), 
    (double) (*lb).pos,
    (double) (*lb).offset);
#else
    REprintf("found %d position, e.g. %ld %ld\n", (int)(ub - lb),
             (unsigned long int)(*lb).pos, (unsigned long int)(*lb).offset);
#endif
    return ub - lb;
  }
}

int SingleChromosomeBCFIndex::readLine(int64_t offset, uint32_t* l_shared,
                                       uint32_t* l_indiv,
                                       std::vector<char>* line) {
  if (bgzf_seek(fBcfFile_, offset, SEEK_SET)) {
    REprintf("seek error!\n");
  }

  if (4 != bgzf_read(fBcfFile_, l_shared, sizeof(uint32_t)) ||
      4 != bgzf_read(fBcfFile_, l_indiv, sizeof(uint32_t))) {
    REprintf("readLine error!\n");
  }
  uint32_t totalLen = *l_shared + *l_indiv;
  line->resize(totalLen);
  if (totalLen != bgzf_read(fBcfFile_, line->data(), totalLen)) {
    REprintf("readLine bgzf_read error!\n");
  }

  return totalLen;
}

int SingleChromosomeBCFIndex::nextLine(uint32_t* l_shared, uint32_t* l_indiv,
                                       std::vector<char>* line) {
  if (4 != bgzf_read(fBcfFile_, l_shared, sizeof(uint32_t)) ||
      4 != bgzf_read(fBcfFile_, l_indiv, sizeof(uint32_t))) {
    REprintf("readLine error!\n");
  }
  uint32_t totalLen = *l_shared + *l_indiv;
  line->resize(totalLen);
  if (totalLen != bgzf_read(fBcfFile_, line->data(), totalLen)) {
    REprintf("readLine bgzf_read error!\n");
  }

  return totalLen;
}
