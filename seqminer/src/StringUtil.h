#ifndef _STRINGUTIL_H_
#define _STRINGUTIL_H_

#include "R.h"
#include <string.h> // for strlen
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <set>
#include <string>
#include <cassert>
#include <algorithm>

/** tokenize the string
 * @return number of tokens we obtained
 * Special case:
 * For empty input string, we will return 1, and @param result will have only 1 element (the empty string)
 * When delim is empty, we will give warning, return 1, and @param result will have the whole input string
 */
inline int stringTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){
  assert(result);
  result->clear();
  if (!delim.size()) {
    REprintf( "stringTokenize() using an empty delim");
    result->push_back(str);
    return -1;
  }

  std::string s;
  unsigned int l = str.size();
  unsigned int i = 0;
  while (i < l) {
    if (delim.find(str[i]) != std::string::npos) { // it's a delimeter
      result->push_back(s);
      s.clear();
    } else {
      s.push_back(str[i]);
    }
    ++i;
  };
  result->push_back(s);
  return result->size();
}

inline int stringTokenize(const std::string& str, const char delim, std::vector<std::string>* result){
    std::string d;
    d.push_back(delim);
    return (stringTokenize(str, d, result));
}

inline bool isEmptyString(const std::string& s){
    return s.size() == 0;
}
/* int stringNaturalTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){ */
/*     int ret = stringTokenize(str, delim, result); */
/*     ret = remove_if(result->begin(), result->end(), isEmptyString) - result->begin(); */
/*     result->resize(ret); */
/*     return ret; */
/* }; */
/* int stringNaturalTokenize(const std::string& str, const char delim, std::vector<std::string>* result){ */
/*     std::string d; */
/*     return (stringNaturalTokenize(str, d, result)); */
/* }; */

// pretty much like stringTokenize, but @param result will not contain empty string
inline int stringNaturalTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){
  assert(result);
  result->clear();
  if (!delim.size()) {
    REprintf( "stringTokenize() using an empty delim");
    result->push_back(str);
    return -1;
  }
  std::string s;
  unsigned int l = str.size();
  unsigned int i = 0;
  while (i < l) {
    if (delim.find(str[i]) != std::string::npos) { // it's a delimeter
      if (s.size()>0){
        result->push_back(s);
        s.clear();
      }
    } else {
      s.push_back(str[i]);
    }
    ++i;
  }
  if (s.size() > 0)
    result->push_back(s);
  return result->size();
}

inline int stringNaturalTokenize(const std::string& str, const char delim, std::vector<std::string>* result){
  std::string d(1, delim);
  return (stringNaturalTokenize(str, d, result));
}

//remove leading and trailing characters
inline void stringStrip(std::string* input, const char* characters = " ") {
  if (!input || !input->size()) return;
    size_t beg = input->find_first_not_of(characters);
    size_t end = input->find_last_not_of(characters);
    input->assign( input->substr(beg, end - beg + 1) );
}

// remove the leading and trailing white spaces
inline std::string stringStrip(const std::string& s){
  unsigned int beg = s.find_first_not_of(' ');
  unsigned int end = s.find_last_not_of(' ');
  return s.substr(beg, end-beg);
}

//extract piece of string from @param beg(inclusive) to @param end(exclusive)
//NOTE: beg and end can be negative meaning count from right hand side. 
inline void stringSlice(std::string* input, int beg, int end) {
    assert(input);
    unsigned int len = input->size();
    if (beg < 0) beg += len;
    if (end < 0) end += len;
    assert (beg >= 0 && end >= 0);
    input -> assign ( input->substr(beg, end- beg)) ;
}

template <class T>
std::string stringJoin(const std::vector<std::string>& input, const T delim) {
    std::string s;
    if (input.size() == 0) {
        return s;
    }
    s = input[0];
    for (unsigned int i = 1; i < input.size(); i++) {
        s+= delim;
        s+= input[i];
    }
    return s;
}
/**
 * for std::string type, we use reference to save memory.
 */
template <>
inline std::string stringJoin<const std::string&>(const std::vector<std::string>& input, const std::string& delim) {
    std::string s;
    if (input.size() == 0) {
        return s;
    }
    s = input[0];
    for (unsigned int i = 1; i < input.size(); i++) {
        s+= delim;
        s+= input[i];
    }
    return s;
}

inline void tolower(std::string* s) {
  for (std::string::iterator i = s->begin();
       i != s->end();
       ++i)
    (*i) = tolower(*i);
}

inline std::string tolower(const std::string& s) {
  std::string ret(s);
  tolower(&ret);
  return ret;
}

inline void toupper(std::string* s) {
  for (std::string::iterator i = s->begin();
       i != s->end();
       ++i)
    (*i) = toupper(*i);
}

inline std::string toupper(const std::string& s) {
  std::string ret(s);
  toupper(&ret);
  return ret;
}
/* std::string toUpper(const std::string& s) { */
/*     std::string r; */
/*     for (unsigned int i = 0; i < s.size(); i++) { */
/*         r.push_back(toupper(s[i])); */
/*     } */
/*     return r; */
/* }; */

/* std::string toLower(const std::string& s) { */
/*     std::string r; */
/*     for (unsigned int i = 0; i < s.size(); i++) { */
/*         r.push_back(tolower(s[i])); */
/*     } */
/*     return r; */
/* }; */

// remove the leading 'chr' if any
inline std::string chopChr(const std::string& s) {
    if (s.size() > 3 && 
        (s[0] == 'c' || s[0] == 'C') &&
        (s[1] == 'h' || s[1] == 'H') &&
        (s[2] == 'r' || s[2] == 'R')){
        return s.substr(3);
    }
    return s;
}

static char _bufferStr[128];
// convert number to char*
inline const char* toStr(const int i) {
    sprintf(_bufferStr, "%d", i);
    return _bufferStr;
}

inline const char* toStr(const double d) {
    sprintf(_bufferStr, "%lf", d);
    return _bufferStr;
}

/* void tolower(std::string* s) { */
/*   for (std::string::iterator i = s->begin(); */
/*        i != s->end(); */
/*        ++i) */
/*     (*i) = tolower(*i); */
/* }; */

/* std::string tolower(const std::string& s) { */
/*   std::string ret(s); */
/*   tolower(&ret); */
/*   return ret; */
/* }; */
/**
 * print out the content for debug only
 */
inline void dumpStringVector(const std::vector<std::string> s) {
  for (unsigned int i = 0; i < s.size(); i++) {
    Rprintf( "%u: %s\n", i, s[i].c_str());
  }
}

/**
 * @return true if @param s ends with @param tail
 */
inline bool endsWith(const std::string& s, const std::string& tail) {
  if (s.size() < tail.size()) return false;
  size_t l = tail.size();
  size_t idx = s.size() - l;
  for (size_t i = 0; i != l; ++i) {
    if (s[idx + i] != tail[i]) {
      return false;
    }
  }
  return true;
}

/**
 * Remove duplicated element from @param input and store it to @param output
 * @return number of duplicated elements
 */
template<class T>
inline int dedup(const std::vector<T>& input,
          std::vector<T>* output) {
  int ret = 0;
  assert(output);
  output->clear();
  std::set<T> checked;
  for (typename std::vector<T>::const_iterator it = input.begin();
       it != input.end();
       ++it) {
    if (checked.count(*it)) {
      ++ret;
      continue;
    }
    checked.insert(*it);
    output->push_back(*it);
  }
  return ret;
}

/**
 * Remove duplicated element from @param input
 */
template<class T>
int dedup(std::vector<T>* input) {
  assert(input);
  int ret = 0;
  if (input->empty()) return ret;
  
  std::set<T> checked;
  size_t n = input->size();
  size_t i = 0;
  size_t p = 0;
  while (p < n) {
    if (checked.count(input->at(p))) {
      ++p;
      continue;
    }
    if (i != p) {
      (*input)[i] = (*input)[p];
    }
    ++i;
    ++p;
  }
  if (i != n) {
    input->resize(i);
  }
  return (n - i);
}

// @return true: if @param s has leading "chr", "CHR", "Chr"...
inline bool hasLeadingChr(const std::string& s) {
  if (s.size() > 3 &&
      (s[0] == 'c' || s[0] == 'C') &&
      (s[1] == 'h' || s[1] == 'H') &&
      (s[2] == 'r' || s[2] == 'R')){
    return true;
  }
  return false;
}

/**
 * split " ",
 *   for "a b", split to "a", "b"
 *   for "a", split to "a"
 *   for "a ", split to "a", ""
 *   for "", split to ""
 */
class StringTokenizer{
 public:
  StringTokenizer(const std::string& i, char token):
            data(i) {
    this->token = token;
    reset();
  }
  StringTokenizer(const std::string& i, const std::string& token):
            data(i) {
    this->token = token;
    reset();
  }
  void reset() {
    this->begin = 0;
    this->end = data.size();
  }
  /**
   * @param piece parsed a piece of string
   * @return true if there are more parsed results
   */
  bool next(std::string* piece) {
    std::string& s = *piece;
    s.clear();
    while (begin <= end) {
          if (begin == end) {
            ++ begin;
            return true;
          }


          const char& c = data[begin];
          if (token.find(c) == std::string::npos) {
            // not a token
            s.push_back(c);
            ++begin;
          } else {
            ++begin;
            return begin < end;
          }
    }
    return begin <= end;
  }
 private:
  const std::string& data;
  std::string token;
  size_t begin;
  size_t end;
}; // StringTokenizer

#endif /* _STRINGUTIL_H_ */
