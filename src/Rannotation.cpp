#include "Rannotation.h"

#include "AnnotationController.h"
#include "R_CPP_interface.h"
#include "Rinternals.h"
#include "BedReader.h"
#include "GenomeScore.h"
#include "TabixReader.h"
#include "AnnotationOutputFile.h"

OutputAnnotationString AnnotationString; // global variable

/**
 * @param param is VALIDATED annotation parameter
 * @param result is parsing result
 * @return true if parsing is successful
 */
std::string parseParameter(SEXP param,
                           const std::string& key,
                           const std::string& def) {
  std::string ret;
  SEXP val = getListElement(param, key.c_str());
  if (val == R_NilValue) {
    return def;
  } else {
    ret.assign(CHAR(STRING_ELT(val, 0)));
    return ret;
  }
}

std::string parseParameter(SEXP param,
                           const char* arg_key,
                           const char* arg_def) {
  std::string key = arg_key;
  std::string def = arg_def;
  return parseParameter(param, key, def);
}

int parseParameter(SEXP param, const std::string& key, const int def) {
  int ret;
  SEXP val = getListElement(param, key.c_str());
  if (val == R_NilValue) {
    return def;
  } else {
    if (Rf_isInteger(val)) {
      ret = INTEGER(val)[0];
      return ret;
    } else if (Rf_isReal(val)){
      ret = (int) (REAL(val)[0]);
      return ret;
    }
    return def;
  }
}

int parseParameter(SEXP param, const char* key, const int def) {
  std::string k = key;
  return parseParameter(param, k, def);
}

bool parseParameter(SEXP param, const std::string& key, const bool def) {
  bool ret = false;
  SEXP val = getListElement(param, key.c_str());
  if (val == R_NilValue) {
    return def;
  } else {
    if (Rf_isLogical(val)) {
      ret = LOGICAL(val)[0];
      return ret;
    } else {
      return ret;
    }
    return def;
  }
}

bool parseParameter(SEXP param, const char* key, const bool def) {
  std::string k  =key;
  return parseParameter(param, k, def);
}

SEXP impl_annotateGene(SEXP s_param,
                       SEXP s_chrom,
                       SEXP s_pos,
                       SEXP s_ref,
                       SEXP s_alt) {
  SEXP ret = R_NilValue;

  const std::string reference    = parseParameter(s_param, "reference", "/net/fantasia/home/zhanxw/anno/resources/human.g1k.v37.fa");
  const std::string codonFile    = parseParameter(s_param, "codonFile", "/net/fantasia/home/zhanxw/anno/codon.txt");
  const std::string priorityFile = parseParameter(s_param, "priorityFile", "/net/fantasia/home/zhanxw/anno/priority.txt");
  const std::string geneFile     = parseParameter(s_param, "geneFile", "/net/fantasia/home/zhanxw/mycode/seqForStat/resources/refFlat_hg19.txt.gz");

  GeneAnnotationParam param;
  param.upstreamRange    = parseParameter(s_param, "upstreamRange", 50);
  param.downstreamRange  = parseParameter(s_param, "downstreamRange", 50);
  param.spliceIntoExon   = parseParameter(s_param, "spliceIntoExon", 3);
  param.spliceIntoIntron = parseParameter(s_param, "spliceIntoIntron", 8);

  // REprintf("param parsed\n");

  AnnotationInputFile aif;
  AnnotationController controller(aif);
  controller.geneAnnotation.setAnnotationParameter(param);
  controller.geneAnnotation.openReferenceGenome(reference);
  controller.geneAnnotation.openCodonFile(codonFile);
  controller.geneAnnotation.openPriorityFile(priorityFile);
  controller.geneAnnotation.openGeneFile(geneFile, "refFlat");

  // REprintf("ready to annotate\n");

  // loop input
  int n = LENGTH(s_chrom);
  std::vector<std::string> anno(n);
  std::vector<std::string> annoFull(n);

  std::string chrom;
  int pos;
  std::string ref;
  std::string alt;
  for (int i = 0; i < n ; ++i) {
    // REprintf("i = %d\n", i);
    chrom = CHAR(STRING_ELT(s_chrom, i));
    pos = INTEGER(s_pos)[i];
    ref = CHAR(STRING_ELT(s_ref, i));
#ifdef __sun
    // Solaris will treat the last argument as pair list
    // very puzzling why this will happen...
    if (TYPEOF(s_alt) == LISTSXP) {
      alt = (const char*)(CAR(s_alt));
      s_alt = CDR(s_alt);
    }
#else
    alt = CHAR(STRING_ELT(s_alt, i));
#endif

    AnnotationString.setFormat("default");
    std::string choppedChr = chopChr(chrom);

    controller.annotate(choppedChr, pos, ref, alt);
    const OrderedMap<std::string, std::string>& result = controller.getResult();
    result.value("ANNO", &anno[i]);
    result.value("ANNOFULL", &annoFull[i]);
  }

  // store results
  REprintf("store results\n");
  int numAllocated = 0;
  numAllocated += createList(2, &ret);
  numAllocated += storeResult(anno, ret, 0);
  numAllocated += storeResult(annoFull, ret, 1);

  std::vector<std::string> dimNames;
  dimNames.push_back("ANNO");
  dimNames.push_back("ANNOFULL");
  numAllocated += setListNames(dimNames, &ret);

  UNPROTECT(numAllocated);
  return ret;
}

SEXP impl_anno(SEXP s_inFile,
               SEXP s_outFile,
               SEXP s_param) {
  SEXP ret = R_NilValue;
  std::string inputFile;
  extractString(s_inFile, &inputFile);
  if (inputFile.empty()) {
    REprintf("No input file provided\n");
    return ret;
  }

  std::string outputFile;
  extractString(s_outFile, &outputFile);
  if (outputFile.empty()) {
    REprintf("No output file provided\n");
    return ret;
  }

  // parse all parameters
  const std::string reference      = parseParameter(s_param, "reference", "/net/fantasia/home/zhanxw/anno/resources/human.g1k.v37.fa");
  const std::string codonFile      = parseParameter(s_param, "codonFile", "/net/fantasia/home/zhanxw/anno/codon.txt");
  const std::string priorityFile   = parseParameter(s_param, "priorityFile", "/net/fantasia/home/zhanxw/anno/priority.txt");
  const std::string geneFile       = parseParameter(s_param, "geneFile", "/net/fantasia/home/zhanxw/mycode/seqForStat/resources/refFlat_hg19.txt.gz");
  const std::string geneFileFormat = parseParameter(s_param, "geneFileFormat", "refFlat");
  const std::string FLAG_bedFile = parseParameter(s_param, "bed", "");
  const std::string FLAG_genomeScore = parseParameter(s_param, "genomeScore", "");
  const std::string FLAG_tabixFile = parseParameter(s_param, "tabix", "");
  const bool FLAG_indexOutput = parseParameter(s_param, "indexOutput", false);
  const std::string FLAG_inputFormat = parseParameter(s_param, "inputFormat", "vcf");
  const bool FLAG_checkReference = parseParameter(s_param, "checkReference", TRUE);

  GeneAnnotationParam param;
  param.upstreamRange    = parseParameter(s_param, "upstreamRange", 50);
  param.downstreamRange  = parseParameter(s_param, "downstreamRange", 50);
  param.spliceIntoExon   = parseParameter(s_param, "spliceIntoExon", 3);
  param.spliceIntoIntron = parseParameter(s_param, "spliceIntoIntron", 8);

  AnnotationInputFile aif(inputFile, FLAG_inputFormat.c_str());
  aif.openReferenceGenome(reference);
  aif.setCheckReference(FLAG_checkReference);

  AnnotationController controller(aif);
  controller.geneAnnotation.setAnnotationParameter(param);
  controller.geneAnnotation.openReferenceGenome(reference);
  controller.geneAnnotation.openCodonFile(codonFile);
  controller.geneAnnotation.openPriorityFile(priorityFile);
  controller.geneAnnotation.openGeneFile(geneFile, geneFileFormat);

  // load resources
  if (!FLAG_bedFile.empty()) {
    REprintf("Use bed option: %s\n", FLAG_bedFile.c_str() );
    std::vector< std::string> fd;
    std::vector< std::string> bed;
    stringTokenize(FLAG_bedFile, ",", &fd);
    for (size_t i = 0; i < fd.size(); ++i ){
      stringTokenize(fd[i], "=", &bed);
      if (bed.size() == 2) {
        controller.openBedFile(bed[0].c_str(), bed[1].c_str());
      } else {
        REprintf("ERROR: Cannot recognized format [ %s ].\n", fd[i].c_str());
        // exit 1;
        return ret;
      };
    };
  }
  if (!FLAG_genomeScore.empty()){
    // REprintf("Use binary GERP score: %s\n", FLAG_genomeScore.c_str());
    // ga.addGenomeScore("GERP", FLAG_genomeScore.c_str());
    REprintf("Use binary score file: %s\n", FLAG_genomeScore.c_str() );
    std::vector< std::string> fd;
    std::vector< std::string> gs;
    stringTokenize(FLAG_genomeScore, ",", &fd);
    for (size_t i = 0; i < fd.size(); ++i ){
      stringTokenize(fd[i], "=", &gs);
      if (gs.size() == 2) {
        controller.openGenomeScoreFile(gs[0].c_str(), gs[1].c_str());
      } else {
        REprintf("ERROR: Cannot recognized format [ %s ].\n", fd[i].c_str());
        //exit(1);
        return ret;
      };
    };
  }
  // parse something like:
  // abc.txt.gz(chrom=1,pos=7,ref=3,alt=4,SIFT=7,PolyPhen=10)
  if(!FLAG_tabixFile.empty()){
    REprintf("Use tabix option: %s\n", FLAG_tabixFile.c_str() );
    std::vector< std::string> fd;
    stringTokenize(FLAG_bedFile, ")", &fd); /// use ')' to split is for convenience
    for (size_t i = 0; i < fd.size(); ++i) {
      fd[i].push_back(')');
      ModelParser mp;
      mp.parse(FLAG_tabixFile);
      std::string fn = mp.getName();
      int chrom, pos, ref, alt;
      mp.assign("chrom", &chrom, 1).assign("pos", &pos, 2).assign("ref", &ref, 3).assign("alt", &alt, 4);
      REprintf("Column %d, %d, %d and %d in tabix file will be matched to chromosome, position, reference allele, alternative allele respectively.\n", chrom, pos, ref, alt);
      TabixData* tabix = new TabixData(fn.c_str(), chrom, pos, ref, alt);

      for (size_t i = 0; i < mp.getParam().size(); ++i) {
        if ( tolower(mp.getParam()[i]) == "chrom" ||
             tolower(mp.getParam()[i]) == "pos" ||
             tolower(mp.getParam()[i]) == "ref" ||
             tolower(mp.getParam()[i]) == "alt") {
          continue;
        }
        int intValue;
        if (str2int(mp.getValue(i), &intValue)) {
          tabix->addTag(mp.getParam()[i], intValue);
          REprintf("Tag %s will be from column %d in tabix file\n", mp.getParam()[i].c_str(), intValue);
        } else {
          tabix->addTag(mp.getParam()[i], mp.getValue(i));
          REprintf("Tag %s will be from column %s (from header) in tabix file\n", mp.getParam()[i].c_str(), mp.getValue(i).c_str());
        }
      }
      controller.addTabixData(tabix);
    }
  } // end halding tabix database

  // annotate, store results
  std::string chrom;
  int pos;
  std::string ref;
  std::string alt;
  AnnotationOutputFile aof(outputFile);
  aof.linkToInput(aif);
  std::string choppedChr; // chop leading 'chr'
  while (aif.extract(&chrom, &pos, &ref, &alt)) {
    choppedChr = chopChr(chrom);
    controller.annotate(choppedChr, pos, ref, alt);
    aof.writeResult(controller.getResult());
  };
  // aof.writeResult(controller.getResult()); // TODO: will add this to handle when input only have comment lines

  // output stats
  controller.geneAnnotation.outputAnnotationStats(outputFile);

  aof.close();
  aif.close();

  if (FLAG_indexOutput) {
    if (aof.indexOutput() == 0) {
      REprintf("DONE: Indexing succeed!\n");
    } else {
      REprintf("ERROR: Indexing failed!\n");
    }
  };

  // return results
  Rprintf("Output file are written to [ %s ].\n", CHAR(STRING_ELT(s_outFile, 0)));
  return ret;
}

SEXP impl_getRefBase(SEXP reference,
                     SEXP chrom,
                     SEXP pos,
                     SEXP length) {
  SEXP ret = R_NilValue;
  GenomeSequence gs;
  std::string fn;
  extractString(reference, &fn);
  REprintf("to open %s\n", fn.c_str());
  if (!gs.open(fn)) {
    return ret;
  }

  int n = LENGTH(chrom);
  std::vector<std::string> seq(n);
  for (int i = 0; i < n; ++i) {
    const std::string _chrom = CHAR(STRING_ELT(chrom, i));
    const int _pos = INTEGER(pos)[i];
    const int _len = INTEGER(length)[i];
    seq[i] = gs.getBase(_chrom, _pos, _pos + _len);
  }
  int numAllocated = 0;
  numAllocated += createStringArray(n, &ret);
  initStringArray(ret);
  numAllocated += storeResult(seq, &ret);
  UNPROTECT(numAllocated);

  return ret;
}
