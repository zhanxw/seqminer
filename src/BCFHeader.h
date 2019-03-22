#ifndef _BCFHEADER_H_
#define _BCFHEADER_H_

#include "StringUtil.h"
#include "TypeConversion.h"

class BCFHeader {
 public:
  // 1-st level keys: contig, FILTER, INFO, FORMAT
  // 2-st level keys:
  //   contig: ID, length
  //   FILTER: ID, description, (will manually add Number=1, Type=String)
  //   INFO: ID, Number, Type, Description
  //   FORMAT: ID, Number, Type, Description
  // we use a vector to store contigs (header_contig_id)
  // and four vectors to store FILTER, INFO, FORMATs (header_id, header_number, header_type, header_description)
  std::vector<std::string> header_contig_id;
  std::vector<std::string> header_id;
  std::vector<std::string> header_number;
  std::vector<std::string> header_type;
  std::vector<std::string> header_description;
  std::vector<std::string> sample_names; // sample names

  class BCFHeaderParser {
   public:
    std::string parseValue(const std::string& s, const std::string& key) {
      size_t begin = s.find("<");
      size_t end = s.rfind(">");
      if (begin == std::string::npos ||
          end == std::string::npos) {
        REprintf( "Wrong intput string during parsing!\n");
      }
      std::string ss = s.substr(begin + 1, end - begin - 1);
      // printf("ss = %s\n", ss.c_str());
      begin = ss.find(key);
      ss = ss.substr(begin, ss.size() - begin);
      if (ss.substr(0, key.size()) != key ||
          ss[key.size()] != '=')  {
        REprintf( "Cannot find the key\n");
      }
      ss = ss.substr(key.size() + 1, ss.size() - key.size() - 1);
      // printf("ss = %s\n", ss.c_str());
      if (ss[0] == '"') {
        end = ss.find("\"", 1);
        ss = ss.substr(1, end - 1);
      } else {
        end = ss.find_first_of(",>");
        ss = ss.substr(0, end);
      }
    
      // printf("ss = %s\n", ss.c_str());
      return ss;
    }
    int parse(const std::string& s) {
      key = id = number = type = desc = "";
      if (s.find("##contig=<") == 0) {
        key = "contig";
        id = parseValue(s, "ID");
      } else if (s.find("##FILTER=<") == 0) {
        key = "filter";
        id = parseValue(s, "ID");
        number = "1";
        type = "String";
        desc = parseValue(s, "Description");
        idx = atoi(parseValue(s, "IDX"));      
      } else if (s.find("##INFO=<") == 0) {
        key = "info";
        id = parseValue(s, "ID");
        number = parseValue(s, "Number");
        type = parseValue(s, "Type");
        desc = parseValue(s, "Description");
        idx = atoi(parseValue(s, "IDX"));      
      } else if (s.find("##FORMAT=<") == 0) {
        key = "format";
        id = parseValue(s, "ID");
        number = parseValue(s, "Number");
        type = parseValue(s, "Type");
        desc = parseValue(s, "Description");
        idx = atoi(parseValue(s, "IDX"));
      }
      return 0;
    }
    std::string key;
    std::string id;
    std::string number;
    std::string type;
    std::string desc;
    int idx;
  };

  int parseHeader(const std::string& header,
                  std::vector<std::string>* header_contig_id,
                  std::vector<std::string>* header_id,
                  std::vector<std::string>* header_number,
                  std::vector<std::string>* header_type,
                  std::vector<std::string>* header_desc) {
    std::vector<std::string> lines;
    stringTokenize(header, "\n", &lines);
    BCFHeaderParser parser;
    for (size_t i = 0; i != lines.size(); ++i) {
      if (parser.parse(lines[i]) < 0) {
        REprintf( "Parser encountered error!\n");
        return -1;
      }
      if (parser.key == "contig") {
        header_contig_id -> push_back(parser.id);
      } else if (parser.key == "filter" ||
                 parser.key == "info" ||
                 parser.key == "format") {
        // need to check the hidden IDX= value
        // e.g. "ID=PASS,Description="All filters passed",IDX=0"
        // in BCF implementation, the default and hidden IDX=0 specified the dictionary index
        // this is important as FILTEr and INFO can share the same ID, for example:
        // ##FILTER=<ID=SVM,Description="Variant failed SVM filter",IDX=6>
        // ##FILTER=<ID=DISC,Description="Mendelian or duplicate genotype discordance is high (3/5% or more)",IDX=7>
        // ...
        // ##INFO=<ID=SVM,Number=1,Type=Float,Description="Milk-SVM score for variant quality, passing -0.5 or greater,IDX=6">

        if (parser.idx == (int) header_id->size()) {
          header_id->push_back(parser.id);
          header_number->push_back(parser.number);
          header_type->push_back(parser.type);
          header_desc->push_back(parser.desc);
        } else if (parser.idx < (int) header_id->size()) {
          Rprintf("BCF index (IDX=) is reused for [%s] with IDX=%d\n", parser.id.c_str(), parser.idx);
          (*header_id)[parser.idx] = parser.id;
          (*header_number)[parser.idx] = parser.number;
          (*header_type)[parser.idx]= parser.type;
          (*header_desc)[parser.idx]= parser.desc;
        } else {
          Rprintf("BCF index is invalid for [%s] with IDX=%d, skipped!\n", parser.id.c_str(), parser.idx);
        }
      } else {
        // do nothing
      }
    }
    Rprintf("Total contig parse = %d, total header index used = %d\n",
           (int) header_contig_id->size(),
           (int) header_id->size());
    return 0;
  }

};

#endif /* _BCFHEADER_H_ */
