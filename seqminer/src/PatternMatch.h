#ifndef _PATTERNMATCH_H_
#define _PATTERNMATCH_H_

#define ERROR_BUF_LEN 64

#include "config.h"

#ifndef HAVE_REG_STARTEND
// e.g. in Solaris, it doesnot implement REG_STARTEND,
// we have to skip using regex related functions in lib
#undef HAVE_POSIX_REGEX
// #pragma message "undef posix (1)"
#endif

#ifdef HAVE_POSIX_REGEX
// We use PCRE here, use 'man pcreposix' for more information
// accordig to http://lh3lh3.users.sourceforge.net/reb.shtml
// PCRE-posix is fast
// #include <pcreposix.h> --> this header is only for PCRE
#include <regex.h>
#include <R.h>
#include <string>

// #pragma message "undef posix (2)"
class Regex {
public:
    /**
     * read pattern like "=Synonymous,=Indel"
     */
    int readPattern(const char* argRegex) {
        // REprintf("pattern = %s\n", argRegex);
        if (this->initialized){
            regfree(&this->pattern);
            this->initialized = false;
        }
        // int cflags = 0;
        int ret = regcomp(& this->pattern, argRegex, 0);
        if (ret) {
            regerror(ret, & this->pattern, error_buf, ERROR_BUF_LEN);
            // REprintf("[ERROR] %s\n", error_buf);
            // fputs(error_buf, stderr);
            return -1;
        }
        // Rprintf("Regcomp OK\n");
        this->initialized = true;
        return 0;
    };
    void readPattern(std::string& argRegex) {
        readPattern(argRegex.c_str());
    };
    /**
     * @return true if matches.
     */
    bool match(const char* text) {
/*         REprintf("To match: %s\n", text); */
/*         if (!this->initialized) { */
/*             REprintf("Uninitialized regex!\n"); */
/*             return false; */
/*         } */
        if (text[0] == '\0') 
            return false;
        int ret = regexec(&this->pattern, text, 1, &this->matchResult, 0);
/*         Rprintf("ret = %d, REG_NOMATCH = %d\n", ret, REG_NOMATCH); */
        if (ret == 0) {
/*             REprintf("Match: %s\n", text); */
            return true;
        } else if (ret == REG_NOMATCH){
/*             REprintf("Nomatch: %s\n", text); */
            return false;
        } else {
            regerror(ret, & this->pattern, error_buf, ERROR_BUF_LEN);
            REprintf("[ERROR] %s\n", error_buf);
            //fputs(error_buf, stderr);
            return false;
        }
        return false;
    };

    /**
     * @return if any pattern matches the @param text[begin...end], will return true
     * @param begin: inclusive
     * @param end: exclusive
     */
    bool match(const char* text, int begin, int end) {
/*         for (int i = begin; i < end; ++i) { */
/*             REprintf("%c", text[i]); */
/*         } */
/*         REprintf("\n"); */
/*         Rprintf("match [ %d - %d ] %s\n", begin, end, text); */
        if (!this->initialized) {
            REprintf("Uninitialized regex!\n");
            return false;
        }

        if (begin == end){
            //empty string
            return false;
        };
        /* size_t nmatch = 1; */
        /* regmatch_t pmatch[nmatch]; */
        this->matchResult.rm_so = begin;
        this->matchResult.rm_eo = end;
        int eflags = REG_STARTEND;
        int ret = regexec(&this->pattern, text + begin, 1, &this->matchResult, eflags);
        if (ret == 0) {
            // check range
            if (this->matchResult.rm_eo <= end) {
                return true;
            } else {
                return false;
            }
        } else if (ret == REG_NOMATCH){
/*             REprintf("Nomatch: %s\n", text); */
            return false;
        } else {
            regerror(ret, & this->pattern, error_buf, ERROR_BUF_LEN);
/*             REprintf("[ERROR] %s\n", error_buf); */
            //fputs(error_buf, stderr);
            return false;
        }
        return false;
    };
    Regex() {
        this->initialized = false;
    }
    ~Regex(){
        if (this->initialized)
            regfree(&pattern);
        this->initialized = false;
    }
    bool isInitialized() const {
      return this->initialized;
    }
private:
    bool initialized;
    regex_t pattern;
    char error_buf[ERROR_BUF_LEN];
    regmatch_t matchResult;
};

#else

#include <string>
#include <vector>
#include "Utils.h"
class Regex {
public:
    /**
     * parse patterns by '|', e.g. "=Synonymous|=Indel"
     * @return 0 if success
     */
    int readPattern(const char* argRegex) {
      if (!argRegex ||  (*argRegex == '\0')) return -1;
      this->initialized = false;
      stringTokenize(argRegex, "|", &this->pattern);
      this->initialized = true;
      return 0;
    };
    void readPattern(std::string& argRegex) {
        readPattern(argRegex.c_str());
    };
    /**
     * @return true if matches.
     */
    bool match(const char* text) const{
        if (!this->initialized) {
            REprintf("Uninitialized regex!\n");
            return false;
        }
        if (text[0] == '\0') 
            return false;
        std::string t;
        for (size_t i = 0 ; i < this->pattern.size(); ++i) {
          if (t.find(this->pattern[i]) != std::string::npos) { //found
            return true;
          }
        }
        return false;
    };

    /**
     * @return if any pattern matches the @param text[begin...end], will return true
     * @param begin: inclusive
     * @param end: exclusive
     */
    bool match(const char* text, int begin, int end) {
        if (!this->initialized) {
            REprintf("Uninitialized regex!\n");
            return false;
        }
        if (begin == end){
            //empty string
            return false;
        };
        std::string t(text + begin, end - begin);
        for (size_t i = 0 ; i < this->pattern.size(); ++i) {
          if (t.find(this->pattern[i]) != std::string::npos) { //found
            return true;
          }
        }
        return false;
    };
    Regex() {
        this->initialized = false;
    }
    ~Regex(){
        this->initialized = false;
    }
    bool isInitialized() const {
      return this->initialized;
    }
private:
    bool initialized;
    std::vector<std::string> pattern;
    char error_buf[ERROR_BUF_LEN];
};

#endif


#endif /* _PATTERNMATCH_H_ */
