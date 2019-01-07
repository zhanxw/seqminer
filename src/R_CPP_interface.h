#ifndef _R_CPP_INTERFACE_H_
#define _R_CPP_INTERFACE_H_

// include <iostream> to avoid ''length'' macro problem on MacOS machines
// see https://stat.ethz.ch/pipermail/r-help/2005-January/065117.html
#include <iostream>

#include <R.h>
#include <Rdefines.h>  // define SEXP

#include <set>
#include <string>
#include <vector>

#include <cstring>

void extractString(SEXP in, std::string* out);

/**
 * NOTE
 */
void extractStringArray(SEXP in, std::vector<std::string>* out);
void extractIntArray(SEXP in, std::vector<int>* out);

void extractStringSet(SEXP in, std::set<std::string>* out);

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char* str);

void dump(std::vector<std::string>& s);

// store @param in to @param ret at index @param idx

void storeResult(const std::vector<bool>& in, SEXP& ret, int idx);
void storeResult(const std::vector<int>& in, SEXP& ret, int idx);
void storeResult(const std::vector<double>& in, SEXP& ret, int idx);
void storeResult(const std::vector<std::string>& in, SEXP ret, int idx);

void storeResult(const std::vector<std::vector<double> >& in, SEXP& ret,
                 int idx);
void storeResult(const std::vector<std::vector<std::vector<double> > >& in,
                 SEXP& ret, int idx);

// store result @param in to a SEXP @param ret
void storeResult(const std::vector<std::string>& in, SEXP* ret);
void storeResult(const std::vector<double>& in, SEXP* ret);
void storeResult(const std::vector<std::vector<double> >& in, SEXP* ret);

// store @param key and @param val to @param ret at index @param idx
void storeResult(const std::string& key, const std::vector<std::string>& val,
                 SEXP ret, int idx);
void storeResult(const std::string& key, const std::vector<int>& val, SEXP& ret,
                 int idx);

// convert and store results
void storeDoubleResult(const std::vector<std::string>& in, SEXP& ret, int idx);
void storeIntResult(const std::vector<std::string>& in, SEXP& ret, int idx);

void setDim(int nrow, int ncol, SEXP s);
void setDim(int i, int j, int k, SEXP s);
/**
 * Set dim attributes for ret[idx]
 */
void setDim(int nrow, int ncol, SEXP ret, int idx);
void setDimNames(const std::vector<std::string>& nrow,
                 const std::vector<std::string>& ncol, SEXP s);
void setDimNames(const std::vector<std::string>& ni,
                 const std::vector<std::string>& nj,
                 const std::vector<std::string>& nk, SEXP s);

#if 0
int createList(int n, SEXP* s); 

int createStringArray(int n, SEXP* s);

int createDoubleArray(int n, SEXP* s);

int createIntArray(int n, SEXP* s);
#endif

void setListNames(std::vector<std::string>& names, SEXP* s);

void initDoubleArray(SEXP s);

void initIntArray(SEXP s);

void initStringArray(SEXP s);

int getDim(SEXP s, std::vector<int>* d);

/**
 * Print the type of @param x
 */
void printType(SEXP x);

#endif /* _R_CPP_INTERFACE_H_ */
