#ifndef _ANNOTATIONINPUTFILE_H_
#define _ANNOTATIONINPUTFILE_H_

typedef enum {IN_FMT_VCF = 0 , IN_FMT_PLAIN, IN_FMT_PLINK, IN_FMT_EPACTS} InputFileFormat;

// hold input file
// and return (chrom, pos, ref, alt) iteratively
class AnnotationInputFile{
 public:
  AnnotationInputFile(): checkReference(true), failedReferenceSite(0), lr(NULL){}
  AnnotationInputFile(const std::string& inputFileName, const std::string& inputFormatStr){
    if (inputFileName.empty() || inputFormatStr.empty()) {
      return;
    }
    // check inputFormat
    std::string inputFormat = tolower(inputFormatStr);
    if (!inputFormat.empty() &&
        inputFormat != "vcf" &&
        inputFormat != "plain" &&
        inputFormat != "plink" &&
        inputFormat != "epacts") {
      REprintf( "Unsupported input format [ %s ], we support VCF, plain, plink and EPACTS formats.\n", inputFormatStr.c_str());
      // LOG << "Unsupported input format [ " << inputFormatStr << " ], we support VCF, plain, plink and EPACTS formats.\n";
      return; // abort();
    };

    if (inputFormat == "vcf" || inputFormat.empty()) {
      // ga.annotateVCF(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = IN_FMT_VCF;
    } else if (inputFormat == "plain") {
      // ga.annotatePlain(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = IN_FMT_PLAIN;
    } else if (inputFormat == "plink") {
      // ga.annotatePlink(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = IN_FMT_PLINK;
    } else if (inputFormat == "epacts") {
      // ga.annotateEpacts(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = IN_FMT_EPACTS;
    } else{
      REprintf( "Cannot recognize input file format: %s \n", inputFileName.c_str());
      return; // abort();
    };

    // open input files
    this->lr = new LineReader(inputFileName);

    // init values
    checkReference = true;
    failedReferenceSite = 0;
  };
  ~AnnotationInputFile() {
    this->close();
  }
  void close() {
    if (this->lr) {
      // REprintf("Delete lr\n");
      delete lr;
      this->lr = NULL;

      if (checkReference && failedReferenceSite > 0) {
        REprintf( "ERROR: Total [ %d ] sites have unmatched reference alleles\n", failedReferenceSite);
        // LOG << "ERROR: Total [ " << failedReferenceSite << " ] sites have unmatched reference alleles\n";
      }
    }
  };


  int openReferenceGenome(const std::string& referenceGenomeFileName) {
    return this->gs.open(referenceGenomeFileName);
  }
  // // extract one header line
  // bool extractHeader(std::string* l) {
  //   this->line.clear();
  //   bool ret = this->lr->readByLine(this->line);
  //   *l = this->line;
  //   return ret;
  // };
  // // extract one header line
  // bool extractHeader(std::string* l) {
  //   this->line.clear();
  //   bool ret = this->lr->readByLine(this->line);
  //   *l = this->line;
  //   return ret;
  // };

  // check ref and alt alleles (may switch ref and alt)
  //@return true if (1) ref match reference; (2) after switch ref and alt, match reference genome.
  bool forceReferenceStrand(const std::string& chrom,
                            const int& pos,
                            std::string* ref,
                            std::string* alt) {
    // determine ref base from reference
    bool refMatchRef = true;
    for (size_t i = 0; i < ref->size(); i++ ) {
      if ((*ref)[i] != gs[chrom][pos - 1 + i]) {
        refMatchRef = false;
        break;
      }
    }
    if (!refMatchRef) {
      bool altMatchRef = true;
      for (size_t i = 0; i < alt->size(); i++ ) {
        if ( (*alt)[i] != gs[chrom][pos - 1 + i]) {
          altMatchRef = false;
          break;
        }
      }
      if (!altMatchRef) {
        REprintf( "Ref [%s] and alt [%s] does not match reference: %s:%d\n", ref->c_str(), alt->c_str(),
                chrom.c_str(), pos);
        return false;
      } else {
        std::swap(*ref, *alt);
      }
    }
    return true;
  }
  void setCheckReference(bool b) {
    this->checkReference = b;
  };
  // if reach end or experience something wrong, will @return false
  bool extract(std::string* chrom,
               int* pos,
               std::string* ref,
               std::string* alt) {
    bool ret;
    do {
      ret = this->lr->readLine(&this->line);
      if (ret == false)
        return ret;
    } while (this->line.empty());

    // for any line beginning with '#', store headers
    // " CHR " is for PLINK header, "CHROM" is for other headers    
    while (line[0] == '#' || line.substr(0,5) == "CHROM" || line.substr(0, 5) == " CHR ") { 
      this->header.push_back(line);
      do {
        ret = this->lr->readLine(&this->line);
        if (ret == false)
          return ret;
      } while (this->line.empty());
    }



    switch (this->format){
      case IN_FMT_VCF:
        stringTokenize(line, "\t ", &fd);
        if (fd.size() < 5) return false;
        *chrom = fd[0];
        *pos = toInt(fd[1]);
        *ref = fd[3];
        *alt = fd[4];
        break;
      case IN_FMT_PLAIN:
        stringNaturalTokenize(line, "\t ", &fd);
        if (fd.size() < 4) return false;
        *chrom = fd[0];
        *pos = toInt(fd[1]);
        *ref = fd[2];
        *alt = fd[3];
        break;
      case IN_FMT_PLINK:
        stringNaturalTokenize(line, "\t ", &fd);
        // will access fd[6], so at least 7 elements here
        if (fd.size() <= 7) return false;
        *chrom = fd[0];
        *pos = toInt(fd[2]);
        *ref = fd[3];
        *alt = fd[6];

        if (!forceReferenceStrand(*chrom, *pos, ref, alt))
          return false;
        break;
      case IN_FMT_EPACTS:
        {
          stringNaturalTokenize(line, "\t ", &fd);
          // e.g.
          // 20:139681_G/A   266     1       1       0.0018797       NA      NA
          // find
          int beg = 0;
          int sep = fd[0].find(':', beg);
          *chrom = fd[0].substr(beg, sep - beg);

          beg = sep + 1;
          sep = fd[0].find('_', beg);
          *pos = toInt(fd[0].substr(beg, sep - beg));

          beg = sep + 1;
          sep = fd[0].find('/', beg);
          *ref = toupper(fd[0].substr(beg, sep - beg));

          beg = sep + 1;
          sep = fd[0].find_first_of(" _\t", beg);
          *alt = toupper(fd[0].substr(beg, sep - beg));

          epactsPrefixLength = sep;
          if ( chrom->empty() || *pos <= 0 || ref->empty() || alt->empty()) {
            REprintf( "Skip line: %s ...." , fd[0].c_str());
            // LOG << "Skip: " << fd[0];
            return false;
          }
        }
        break;
      default:
        REprintf( "Unknown format, quitting!\n");
        return false; // abort();
        break;
    }// end switch

    // verify reference
    if (this->checkReference) {
      // getBase() use 0-based index
      std::string refFromGenome = this->gs.getBase(*chrom, *pos - 1, *pos + ref->size() - 1);
      if ((*ref) != refFromGenome) {
        ++ failedReferenceSite;
        if (failedReferenceSite <= 10) { // output at most 10 warnings
          REprintf( "ERROR: Reference allele [ %s ] does not match reference genome [ %s ] at [ %s:%d ]\n", ref->c_str(), refFromGenome.c_str(), chrom->c_str(), *pos);
        }
        // LOG << "ERRROR: Reference allele [" << (*ref) <<   "]  does not match reference genome [" << refFromGenome << "] at " << *chrom << ":" << *pos << "\n";
      };
    }
    return true;
  };
  InputFileFormat getFormat() const {
    return this->format;
  };
  std::string getEpactsPrefix(const std::string& s) const{
    return s.substr(0, this->epactsPrefixLength);
  };
  const std::vector< std::string>& getHeader() const {
    return this->header;
  };
  const std::vector< std::string>& getFields() const {
    return this->fd;
  };
 private:
  bool checkReference;
  int failedReferenceSite;
  InputFileFormat format;
  LineReader* lr;
  std::vector< std::string> fd;
  std::string line;
  std::vector< std::string> header;
  GenomeSequence gs; // check if ref alleles matched reference genome
  size_t epactsPrefixLength;
}; // end class AnnotationInputFile



#endif /* _ANNOTATIONINPUTFILE_H_ */
