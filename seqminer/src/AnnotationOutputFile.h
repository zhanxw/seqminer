#ifndef _ANNOTATIONOUTPUTFILE_H_
#define _ANNOTATIONOUTPUTFILE_H_

#ifdef HAVE_UNISTD_H
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include "IO.h"
#include "AnnotationInputFile.h"
#include "OrderedMap.h"

#ifdef HAVE_UNISTD_H
bool fileExists(std::string& s) {
  struct stat res;
  int ret = stat(s.c_str(), &res);
  return (ret == 0);
}

time_t getFileMtime(std::string& s) {
  struct stat res;
  int ret = stat(s.c_str(), &res);
  if (ret != 0) {
    // perror("stat");
    return -1;
  }
  return res.st_mtime;
}
#endif

class AnnotationOutputFile{
 public:
AnnotationOutputFile(const std::string& out):headerOutputted(false), totalVariants(0),outputFileName(out) {
    /* if (strncmp (out, "stdout", 6) ) { */
    /*   this->fout = fopen(out, "wt"); */
    /*   if (!this->fout) { */
    /*     REprintf( "Cannot open otuput file %s for write.\n", out); */
    /*   } */
    /* } else { */
    /*   this->fout = stdout; */
    /* } */
    if (hasSuffix(out, ".gz")) {
      this->fout = new FileWriter(out.c_str(), BGZIP);
    } else {
      this->fout = new FileWriter(out.c_str());
    }
  };
  ~AnnotationOutputFile() {
    this->close();
  }
  void close() {
    if (this->fout) {
      REprintf( "DONE: %d varaints are annotated.\n", totalVariants);
      REprintf( "DONE: Generated annotation output in [ %s ].\n", outputFileName.c_str());
      delete this->fout;
      this->fout = NULL;
    }
  };
  void linkToInput(const AnnotationInputFile& in) {
    this->aif =  &in;
  };
  void writeHeader() {
    writeHeader(this->aif->getHeader());
  };
  void writeHeader(const std::vector< std::string> & h) {
    for (size_t i = 0; i < h.size(); ++i) {
      if (i != h.size() - 1) {
        this->fout->write(h[i].c_str());
        this->fout->write("\n");
      } else { // last line
        if ( (aif->getFormat() == IN_FMT_PLAIN || aif->getFormat() == IN_FMT_PLINK) &&
             (h[i].substr(0, 5) == "CHROM" || h[i].substr(0, 6) == "#CHROM") ) {
          this->fout->write(h[i].c_str());
          this->fout->write("\tANNO\tANNO_FULL\n");
        } else {
          this->fout->write(h[i].c_str());
          this->fout->write("\n");
        }
      }
    }
  }
  void writeResult(const OrderedMap<std::string, std::string>& res) {
    //REprintf( "Need to link output to input!\n");
    assert(aif);

    if (!this->headerOutputted) {
      this->writeHeader();
      this->headerOutputted = true;
    }

    const std::vector< std::string> & field = aif->getFields();
    std::string key;
    std::string val;
    switch (aif->getFormat()) {
      case IN_FMT_VCF:
        {
          for (size_t i = 0; i < field.size(); ++i) {
            if (i) this->fout->write("\t");
            this->fout->write(field[i].c_str()) ;
            if (i == 7) { // 7: the 8th column in 0-based index
              if (!field[i].empty())
                this->fout->write(VCF_INFO_SEPARATOR);
              for (size_t j = 0; j < res.size(); ++j) {
                if (j)
                  this->fout->write(VCF_INFO_SEPARATOR);
                res.at(j, &key, &val);
                this->fout->write(key.c_str());
                this->fout->write("=");
                this->fout->write(val.c_str());
              }
            }
          };
        }
        break;
      case IN_FMT_PLAIN:
      case IN_FMT_PLINK:
        {
          for (size_t i = 0; i < field.size(); ++i) {
            if (i) this->fout->write("\t");
            this->fout->write(field[i].c_str());
          };
          for (size_t j = 0; j < res.size(); ++j) {
            this->fout->write("\t");
            res.at(j, &key, &val);
            // this->fout->write(key.c_str());
            // this->fout->write("=");
            this->fout->write(val.c_str());
          }
        }
        break;
      case IN_FMT_EPACTS:
        {
          //this->fout->write(field[0].substr(0, sep).c_str());
          this->fout->write(aif->getEpactsPrefix(field[0]).c_str());
          this->fout->write("_");
          std::string s;
          res.value("ANNO", &s);
          this->fout->write(s.c_str());

          for (unsigned int i = 1; i < field.size(); i++ ){
            this->fout->write("\t");
            this->fout->write(field[i].c_str()) ;
          }
        }
        break;
      default:
        break;
    };
    this->fout->write("\n");
    ++this->totalVariants;
  };
  int indexOutput() {
    // back up existing/older index file
    size_t len = outputFileName.size();
    std::string indexFileName = outputFileName.substr(0, len - 3); // strip ".gz"
#ifdef HAVE_UNISTD_H
    if (fileExists(indexFileName) &&
        getFileMtime(indexFileName) <= getFileMtime(outputFileName)) { // index too old
      remove(indexFileName.c_str());
      REprintf( "Done: Removed old index file ...\n");
    }
#endif
    if ( bgzf_is_bgzf(outputFileName.c_str())!=1 )
    {
      REprintf("[tabix] was bgzip used to compress this file? %s\n", outputFileName.c_str());
      return 1;
    }

    // set up format
    ti_conf_t conf = ti_conf_vcf;
    switch (aif->getFormat()) {
      case IN_FMT_VCF:
        conf = ti_conf_vcf;
        break;
      case IN_FMT_PLAIN:
      case IN_FMT_PLINK:
        conf.sc = 1;
        conf.bc = 2;
        conf.ec = 2;
        conf.meta_char = '#';
        conf.line_skip = 0;
        for (size_t i = 0; i < aif->getHeader().size(); ++i) {
          if (aif->getHeader()[i][0] != '#')
            conf.line_skip ++;
        }
        break;
      case IN_FMT_EPACTS:
        REprintf( "EPACTS outputs are not tab-delimited and not indextable for now.\n");
        return -1;
    }

    // index
    return ti_index_build(outputFileName.c_str(), &conf);    
  };
  
private:
  // don't copy
  AnnotationOutputFile(const AnnotationOutputFile& fw);
  AnnotationOutputFile& operator=(const AnnotationOutputFile& fw);
 private:
  bool headerOutputted;
  const AnnotationInputFile* aif;
  FileWriter* fout;
  int totalVariants;
  std::string outputFileName;
}; // class AnnotationOutputFile



#endif /* _ANNOTATIONOUTPUTFILE_H_ */
