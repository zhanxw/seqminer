#ifndef _FILEFORMAT_H_
#define _FILEFORMAT_H_

class MetaFileFormat{
 public:
  /// check if all the necessary headers are set up
  virtual bool isComplete() const = 0 ;
  void setHeader(const std::map<std::string, int>& header) {
    std::map<std::string, int>::const_iterator iter = header.begin();
    for ( ;
          iter != header.end();
          ++iter) {
      if (setHeader(iter->first, iter->second) < 0) {
        REprintf("Problem when using the header [ %s ] for column [ %d ]\n", iter->first.c_str(), iter->second);
      }
    }
  }
  int setHeader(const std::string& key, const int& val) {
    if (val < 0) return -1;
    if (get(key) >= 0) // duplicate
      return -1;
    data[key] = val;
    return 0;
  }
  int get(const std::string& key) {
    std::string k = toupper(key);
    std::map<std::string, int>::const_iterator it = data.find(k);

    if (it != data.end()) {
      return it->second;
    }

    // check synonym
    if (synonym.count(k) != 0) {
      const std::set<std::string> & s = synonym.find(k)->second;
      std::set<std::string>::const_iterator i;
      for ( i = s.begin();
            i != s.end();
            ++i) {
        if (data.find(*i) != data.end()) {
          return data.find(*i)->second;
        }
      }
    }

    missingKey.insert(key);
    return -1;
  }

  int peakHeader(const std::string& fn, std::map<std::string, int>* header) {
    header->clear();
    LineReader lr(fn);
    std::vector<std::string> fd;
    std::string line;
    while(lr.readLine(&line)) {
      stringNaturalTokenize(line, "\t ", &fd);
      if (fd.empty()) continue;
      if (fd[0][0] == '#' && fd[0] != "#CHROM") continue;
      // strip out the '#' prefix
      if (fd[0][0] == '#') {
        fd[0] = fd[0].substr(1, fd[0].size() - 1);
      }
      for (size_t i = 0; i < fd.size(); ++i ) {
        if (header->count(fd[i]) ) {
          REprintf("Duplicatd header [ %s ]\n", fd[i].c_str());
          continue;
        }
        (*header) [toupper(fd[i])] = i;
      }
      return 0;
    };
    return -1;
  }
  int open(const std::string& fn) {
    return peakHeader(fn, &data);
  }
  void dump() {
    REprintf("Missing header:\n");
    for (std::set<std::string>::const_iterator it = missingKey.begin();
         it != missingKey.end();
         ++it) {
      REprintf("[ %s ] \n", it->c_str());
    }
    REprintf("Known header:\n");
    for (std::map<std::string, int>::const_iterator it = data.begin();
         it != data.end();
         ++it) {

      REprintf("[ %s ] => [ %d ]\n", it->first.c_str(), it->second);
    }
    REprintf("Synonym headers:\n");
    for (std::map<std::string, std::set<std::string> >::const_iterator it = synonym.begin();
         it != synonym.end();
         ++it) {
      REprintf("[ %s ] => ", it->first.c_str());

      const std::set<std::string> & s = synonym.find(it->first)->second;
      std::set<std::string>::const_iterator i;
      for ( i = s.begin();
            i != s.end();
            ++i) {
        REprintf("[ %s ] ", i->c_str());
      }
      REprintf("\n");
    }
  };
  /// add synonym of s1 and s2
  int addSynonym(const std::string& key1, const std::string& key2) {
    std::string s1 = toupper(key1);
    std::string s2 = toupper(key2);
    if (s1 == s2) return 0;
    if (synonym.find(s1) != synonym.end() ) { // already synonym
      if (synonym[s1].count(s2) != 0) {
        return 0;
      }
    }
    if (synonym.find(s2) != synonym.end() ) {
      if (synonym[s2].count(s1) != 0) {
        return 0;
      }
    }

    // add synonym
    {
      synonym[s1].insert(s2);
      const std::set<std::string> & s = synonym.find(s1)->second;
      std::set<std::string>::const_iterator i;
      for ( i = s.begin();
            i != s.end();
            ++i) {
        synonym[*i].insert(s2);
      }
    }
    {
      synonym[s2].insert(s1);
      const std::set<std::string> & s = synonym.find(s2)->second;
      std::set<std::string>::const_iterator i;
      for ( i = s.begin();
            i != s.end();
            ++i) {
        synonym[*i].insert(s1);
      }
    }
    return 0;
  }

 private:
  std::map<std::string, int> data;
  std::set<std::string> missingKey;
  std::map<std::string, std::set <std::string> > synonym; // key: word val: set of synonym of the word
};


/**
   CHROM   POS     REF     ALT     N_INFORMATIVE   AF      INFORMATIVE_ALT_AC      CALL_RATE       HWE_PVALUE      N_REF   N_HET   N_ALT   U_STAT  SQRT_V_STAT     ALT_EFFSIZE     PVALUE
*/
class PvalFileFormat: public MetaFileFormat {
 public:
  PvalFileFormat() {
    addSynonym("AF", "ALL_AF");
  }
  bool isComplete() const {
    return true;
  }
};
/**
   CHROM   START_POS       END_POS NUM_MARKER      MARKER_POS      COV
*/
class CovFileFormat: public MetaFileFormat {
 public:
  CovFileFormat() {
    addSynonym("CURRENT_POS", "START_POS");
    addSynonym("MARKERS_IN_WINDOW", "MARKER_POS");
    addSynonym("COV_MATRICES", "COV");
    addSynonym("CURRENT_POS", "END_POS");
  }
  bool isComplete() const {
    return true;
  }
};

#endif /* _FILEFORMAT_H_ */
