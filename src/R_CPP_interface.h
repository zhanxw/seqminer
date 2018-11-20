#ifndef _R_CPP_INTERFACE_H_
#define _R_CPP_INTERFACE_H_

// include <iostream> to avoid ''length'' macro problem on MacOS machines
#include <iostream>   
#include <R.h>
#include <Rdefines.h> // define SEXP

#include <string>
#include <vector>
#include <set>

#include <cstring>

void extractString(SEXP in, std::string* out);

/**
 * NOTE
 */
void extractStringArray(SEXP in, std::vector<std::string>* out);

void extractStringSet(SEXP in, std::set<std::string>* out);

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str);

void dump(std::vector<std::string> & s);

int storeResult(const std::vector<std::string>& in , SEXP ret, int idx) ;

int storeResult(const std::vector<int>& in , SEXP& ret, int idx);

int storeIntResult(const std::vector<std::string>& in , SEXP& ret, int idx);

int storeResult(const std::vector<double>& in , SEXP& ret, int idx);

int storeDoubleResult(const std::vector<std::string>& in , SEXP& ret, int idx);

int storeResult(const std::vector<bool>& in ,  SEXP& ret, int idx) ;

int storeResult(const std::vector<std::vector<double> >& in ,  SEXP& ret, int idx) ;

int storeResult(const std::vector<std::vector<std::vector<double> > >& in ,  SEXP& ret, int idx);

int storeResult(const std::string& key, const std::vector<std::string>& val , SEXP ret, int idx); 

int storeResult(const std::string& key, const std::vector<int>& val , SEXP& ret, int idx);


int setDim(int nrow, int ncol, SEXP* s);
int setDim(int i, int j, int k, SEXP* s);
int setDimNames(const std::vector<std::string>& nrow,
                const std::vector<std::string>& ncol,
                SEXP* s);
int setDimNames(const std::vector<std::string>& ni,
                const std::vector<std::string>& nj,
                const std::vector<std::string>& nk,                
                SEXP* s);

/**
 * Set dim attributes for ret[idx]
 */
int setDim(int nrow, int ncol, SEXP ret, int idx);

int createList(int n, SEXP* s); 

int createStringArray(int n, SEXP* s);

int createDoubleArray(int n, SEXP* s);

int createIntArray(int n, SEXP* s);

int setListNames(std::vector<std::string>& names, SEXP* s);

void initDoubleArray(SEXP s); 

void initIntArray(SEXP s); 

void initStringArray(SEXP s); 

/**
 * Another set of utility function
 */
int storeResult(const std::vector<std::string>& in , SEXP* ret) ;
int storeResult(const std::vector<double>& in, SEXP* ret);
int storeResult(const std::vector<std::vector<double> >& in, SEXP* ret);

int getDim(SEXP s, std::vector<int>* d);

/**
 * Print the type of @param x
 */
void printType(SEXP x);
  
#endif /* _R_CPP_INTERFACE_H_ */
