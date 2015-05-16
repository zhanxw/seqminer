#ifndef _ANNOTATIONCONTROLLER_H_
#define _ANNOTATIONCONTROLLER_H_

#include "Type.h"
#include "Codon.h"
#include "GenomeSequence.h"
#include "Gene.h"
#include "GeneFormat.h"
#include "SequenceUtil.h"
#include "FreqTable.h"
#include "StringTemplate.h"
#include "GeneAnnotation.h"
#include "BedReader.h"
#include "GenomeScore.h"
#include "ModelParser.h"
#include "TabixData.h"
#include "AnnotationInputFile.h"

// run annotations (gene based, bed file, genomeScores, tabix database)
// and store results
class AnnotationController{
 public:
  AnnotationController(AnnotationInputFile& in):aif(in) {
  };
  virtual ~AnnotationController() {
    for (size_t i = 0; i < bedReader.size() ; ++i) {
      delete bedReader[i];
    }
    for (size_t i = 0; i < genomeScore.size(); ++i ) {
      delete genomeScore[i];
    };
    for (size_t i = 0; i < tabixData.size(); ++i ) {
      delete tabixData[i];
    };
  };
  void openBedFile(const char* tag, const char* fn) {
    // check duplication
    for (size_t i = 0; i < this->bedTag.size(); ++i ) {
      if (this->bedTag[i] == tag) {
        REprintf( "ERROR: Duplicated tag [ %s ] \n", tag);
        return;
      }
    }

    // add bedFile
    BedReader* p = new BedReader;
    int ret = p->open(fn);
    if (ret < 0) {
      REprintf( "Cannot open BED file: [ %s ]\n", fn);
      delete p;
      return ;
    } else {
      REprintf( "DONE: Load %d regions from BED file\n", ret);
    };

    this->bedTag.push_back(tag);
    this->bedReader.push_back(p);
  };
  void openGenomeScoreFile(const char* tag, const char* fn){
    // check duplication
    for (size_t i = 0; i < this->genomeScoreTag.size(); ++i) {
      if (this->genomeScoreTag[i] == tag ) {
        REprintf( "ERROR: Duplicated tag [ %s ] \n", tag);
        return;
      }
    }

    // add genome score file
    GenomeScore* p = new GenomeScore(fn);
    this->genomeScoreTag.push_back(tag);
    this->genomeScore.push_back(p);
  };

  void addTabixData(TabixData* t) {
    this->tabixData.push_back(t);
  };

  void annotate(std::string& chrom,
                int& pos,
                std::string& ref,
                std::string& alt) {
    this->geneAnnotation.annotate(chrom, pos, ref, alt);

    this->result.clear();
    this->result["ANNO"] = this->geneAnnotation.getTopPriorityAnnotation();
    this->result["ANNOFULL"] = this->geneAnnotation.getFullAnnotation();

    std::vector<std::string> bedString;
    for (size_t br = 0; br < this->bedReader.size(); ++ br) {
      if (this->bedReader[br]->find(chrom.c_str(), pos, &bedString)){
        if (!bedString.empty())  {
          this->result[this->bedTag[br]] = stringJoin(bedString, ',');
        } else {
          this->result[this->bedTag[br]] = "";
        }
      }
    }
    for (size_t gs = 0; gs < this->genomeScore.size(); ++ gs) {
      this->result[this->genomeScoreTag[gs]] = toStr(this->genomeScore[gs]->baseScore(chrom.c_str(), pos));
    }
    for (size_t tb = 0; tb < this->tabixData.size(); ++ tb) {
      TabixData& tabix = *this->tabixData[tb];
      tabix.addAnnotation(chrom, pos, ref, alt);
      size_t s = tabix.getTag().size();
      for (size_t i = 0; i < s; ++i) {
        this->result[tabix.getTag()[i]] = tabix.getAnnotation()[i];
      }
    }
  };

  const OrderedMap<std::string, std::string>& getResult() const {
    return this->result;
  }
 public:
  AnnotationInputFile& aif;
  // various annotation types
  GeneAnnotation geneAnnotation;
 private:
  // various annotation types
  std::vector<BedReader*> bedReader;
  std::vector<std::string> bedTag;
  std::vector<GenomeScore*> genomeScore;
  std::vector<std::string> genomeScoreTag;
  std::vector<TabixData*> tabixData;
  OrderedMap<std::string, std::string> result; // store all types of annotation results
};

#endif /* _ANNOTATIONCONTROLLER_H_ */
