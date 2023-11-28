#include "R_CPP_interface.h"

#include "TypeConversion.h"

void extractString(SEXP in, std::string* out) {
  *out = CHAR(STRING_ELT(in, 0));
}

/**
 * NOTE
 */
void extractStringArray(SEXP in, std::vector<std::string>* out) {
  out->clear();
  std::string s;
  for (R_len_t i = 0; i < length(in); i++) {
    s = CHAR(STRING_ELT(in, i));
    if (s.size()) {
      out->push_back(s);
      // Rprintf("extractStringArray: [%d] %s\n", out->size(), s.c_str());
    }
  }
}

void extractIntArray(SEXP in, std::vector<int>* out) {
  out->clear();
  for (R_len_t i = 0; i < length(in); i++) {
    out->push_back(INTEGER(in)[i]);
  }
}

void extractStringSet(SEXP in, std::set<std::string>* out) {
  std::string s;
  for (R_len_t i = 0; i < length(in); i++) {
    s = CHAR(STRING_ELT(in, i));
    out->insert(s);
    // Rprintf("extractStringArray: [%d] %s\n", out->size(), s.c_str());
  }
}

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char* str) {
  SEXP elmt = R_NilValue;
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (R_len_t i = 0; i < length(list); i++)
    if (std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

void dump(std::vector<std::string>& s) {
#ifndef __MINGW64__
  Rprintf("Dump %zu elements:\n", s.size());
#else
  Rprintf("Dump %lu elements:\n", (unsigned long int)s.size());
#endif
  for (unsigned int i = 0; i != s.size(); i++) {
    Rprintf("s[%u] = \"%s\"\n", i, s[i].c_str());
  }
}

void storeResult(const std::vector<bool>& in, SEXP& ret, int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = in.size();
  PROTECT(s = allocVector(LGLSXP, n));
  for (int i = 0; i < n; ++i) {
    LOGICAL(s)[i] = in[i];
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(1);
}

void storeResult(const std::vector<int>& in, SEXP& ret, int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = in.size();
  PROTECT(s = allocVector(INTSXP, n));
  for (int i = 0; i < n; ++i) {
    INTEGER(s)[i] = in[i];
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(1);
}

void storeResult(const std::vector<std::string>& in, SEXP ret, int idx) {
  SEXP s;  //  = VECTOR_ELT(ret, i);
  int n = in.size();
  PROTECT(s = allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) {
    SET_STRING_ELT(s, i, mkChar(in[i].c_str()));
    // Rprintf("storeResult [%d] = %s\n", i, in[i].c_str());
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(1);
}

void storeResult(const std::vector<double>& in, SEXP& ret, int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = in.size();
  PROTECT(s = allocVector(REALSXP, n));
  for (int i = 0; i < n; ++i) {
    REAL(s)[i] = in[i];
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(1);
}

void storeIntResult(const std::vector<std::string>& in, SEXP& ret, int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = in.size();
  int tmp;
  PROTECT(s = allocVector(INTSXP, n));
  for (int i = 0; i < n; ++i) {
    if (str2int(in[i], &tmp))
      INTEGER(s)[i] = tmp;
    else
      INTEGER(s)[i] = NA_INTEGER;
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(1);
}

void storeDoubleResult(const std::vector<std::string>& in, SEXP& ret, int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = in.size();
  double tmp;
  PROTECT(s = allocVector(REALSXP, n));
  for (int i = 0; i < n; ++i) {
    if (str2double(in[i], &tmp))
      REAL(s)[i] = tmp;
    else
      REAL(s)[i] = NA_REAL;
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(1);
}

void storeResult(const std::vector<std::vector<double> >& in, SEXP& ret,
                 int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = in.size();
  int numAllocated = 0;
  PROTECT(s = allocVector(VECSXP, n));
  numAllocated++;
  for (int i = 0; i < n; ++i) {
    SEXP si;
    storeResult(in[i], &si);
    SET_VECTOR_ELT(s, i, si);
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(numAllocated);
}

void storeResult(const std::vector<std::vector<std::vector<double> > >& in,
                 SEXP& ret, int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = in.size();
  int numAllocated = 0;
  PROTECT(s = allocVector(VECSXP, n));
  numAllocated++;
  for (int i = 0; i < n; ++i) {
    SEXP si;
    storeResult(in[i], &si);
    SET_VECTOR_ELT(s, i, si);
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(numAllocated);
}

void storeResult(const std::string& key, const std::vector<std::string>& val,
                 SEXP ret, int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = val.size();
  PROTECT(s = allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) {
    SET_STRING_ELT(s, i, mkChar(val[i].c_str()));
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(1);
}

void storeResult(const std::string& key, const std::vector<int>& val, SEXP& ret,
                 int idx) {
  SEXP s;  // = VECTOR_ELT(ret, i);
  int n = val.size();
  PROTECT(s = allocVector(INTSXP, n));
  for (int i = 0; i < n; ++i) {
    INTEGER(s)[i] = val[i];
  }
  SET_VECTOR_ELT(ret, idx, s);
  UNPROTECT(1);
}

void setDim(int nrow, int ncol, SEXP s) {
  SEXP dim;
  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = nrow;
  INTEGER(dim)[1] = ncol;
  setAttrib(s, R_DimSymbol, dim);
  UNPROTECT(1);
}

void setDim(int i, int j, int k, SEXP s) {
  SEXP dim;
  PROTECT(dim = allocVector(INTSXP, 3));
  INTEGER(dim)[0] = i;
  INTEGER(dim)[1] = j;
  INTEGER(dim)[2] = k;
  setAttrib(s, R_DimSymbol, dim);
  UNPROTECT(1);
}

/**
 * Set dim attributes for ret[idx]
 */
void setDim(int nrow, int ncol, SEXP ret, int idx) {
  SEXP s = VECTOR_ELT(ret, idx);
  setDim(nrow, ncol, s);
  SET_VECTOR_ELT(ret, idx, s);
}

void setDimNames(const std::vector<std::string>& nrow,
                 const std::vector<std::string>& ncol, SEXP s) {
  int numAllocated = 0;
  SEXP dimnames;
  PROTECT(dimnames = allocVector(VECSXP, 2));
  ++numAllocated;

  storeResult(nrow, dimnames, 0);
  storeResult(ncol, dimnames, 1);
  setAttrib(s, R_DimNamesSymbol, dimnames);
  UNPROTECT(numAllocated);
}

void setDimNames(const std::vector<std::string>& ni,
                 const std::vector<std::string>& nj,
                 const std::vector<std::string>& nk, SEXP s) {
  int numAllocated = 0;
  SEXP dimnames;
  PROTECT(dimnames = allocVector(VECSXP, 3));
  ++numAllocated;

  storeResult(ni, dimnames, 0);
  storeResult(nj, dimnames, 1);
  storeResult(nk, dimnames, 2);
  setAttrib(s, R_DimNamesSymbol, dimnames);

  UNPROTECT(numAllocated);
}

#if 0
int createList(int n, SEXP* s) {
  PROTECT((*s) = allocVector(VECSXP, n));
  return 1;
}

int createStringArray(int n, SEXP* s) {
  PROTECT((*s) = allocVector(STRSXP, n));
  return 1;
}

int createDoubleArray(int n, SEXP* s) {
  PROTECT((*s) = allocVector(REALSXP, n));
  return 1;
}

int createIntArray(int n, SEXP* s) {
  PROTECT((*s) = allocVector(INTSXP, n));
  return 1;
}
#endif

void setListNames(std::vector<std::string>& names, SEXP* s) {
  SEXP sListNames;
  PROTECT(sListNames = allocVector(STRSXP, names.size()));
  for (unsigned int i = 0; i != names.size(); ++i) {
    SET_STRING_ELT(sListNames, i, mkChar(names[i].c_str()));
  }
  setAttrib((*s), R_NamesSymbol, sListNames);
  UNPROTECT(1);
}

void initDoubleArray(SEXP s) {
  double* r = REAL(s);
  for (int i = 0; i < length(s); i++) {
    r[i] = NA_REAL;
  }
}

void initIntArray(SEXP s) {
  int* r = INTEGER(s);
  for (int i = 0; i < length(s); i++) {
    r[i] = NA_INTEGER;
  }
}

void initStringArray(SEXP s) {
  for (int i = 0; i < length(s); i++) {
    SET_STRING_ELT(s, i, NA_STRING);
  }
}

void storeResult(const std::vector<std::string>& in, SEXP* ret) {
  int alloc = 0;
  PROTECT((*ret) = allocVector(STRSXP, in.size()));
  alloc++;
  for (size_t i = 0; i < in.size(); ++i) {
    SET_STRING_ELT((*ret), i, mkChar(in[i].c_str()));
  }
  UNPROTECT(alloc);
}

void storeResult(const std::vector<double>& in, SEXP* ret) {
  int alloc = 0;
  PROTECT((*ret) = allocVector(REALSXP, in.size()));
  alloc++;
  for (size_t i = 0; i < in.size(); ++i) {
    REAL(*ret)[i] = in[i];
  }
  UNPROTECT(alloc);
}

void storeResult(const std::vector<std::vector<double> >& in, SEXP* ret) {
  int alloc = 0;
  PROTECT((*ret) = allocVector(VECSXP, in.size()));
  alloc++;
  for (size_t i = 0; i < in.size(); ++i) {
    SEXP si;
    storeResult(in[i], &si);
    SET_VECTOR_ELT((*ret), i, si);
  }
  UNPROTECT(alloc);
}

/**
 * Store dimensions of @param s in @param d
 * @return 0 if success
 */
int getDim(SEXP s, std::vector<int>* d) {
  SEXP r = getAttrib(s, R_DimSymbol);
  if (isNull(r)) return -1;
  int n = length(r);
  d->resize(n);
  for (int i = 0; i < n; ++i) {
    (*d)[i] = INTEGER(r)[i];
  }
  return 0;
}

/**
 * Print the type of @param x
 */
void printType(SEXP x) {
  switch (TYPEOF(x)) {
    case NILSXP:
      REprintf("NILSXP");
      break;
    case SYMSXP:
      REprintf("SYMSXP");
      break;
    case LISTSXP:
      REprintf("LISTSXP");
      break;
    case CLOSXP:
      REprintf("CLOSXP");
      break;
    case ENVSXP:
      REprintf("ENVSXP");
      break;
    case PROMSXP:
      REprintf("PROMSXP");
      break;
    case LANGSXP:
      REprintf("LANGSXP");
      break;
    case SPECIALSXP:
      REprintf("SPECIALSXP");
      break;
    case BUILTINSXP:
      REprintf("BUILTINSXP");
      break;
    case CHARSXP:
      REprintf("CHARSXP");
      break;
    case LGLSXP:
      REprintf("LGLSXP");
      break;
    case INTSXP:
      REprintf("INTSXP");
      break;
    case REALSXP:
      REprintf("REALSXP");
      break;
    case CPLXSXP:
      REprintf("CPLXSXP");
      break;
    case STRSXP:
      REprintf("STRSXP");
      break;
    case DOTSXP:
      REprintf("DOTSXP");
      break;
    case ANYSXP:
      REprintf("ANYSXP");
      break;
    case VECSXP:
      REprintf("VECSXP");
      break;
    case EXPRSXP:
      REprintf("EXPRSXP");
      break;
    case BCODESXP:
      REprintf("BCODESXP");
      break;
    case EXTPTRSXP:
      REprintf("EXTPTRSXP");
      break;
    case WEAKREFSXP:
      REprintf("WEAKREFSXP");
      break;
    case S4SXP:
      REprintf("S4SXP");
      break;
    case RAWSXP:
      REprintf("RAWSXP");
      break;
    default:
      REprintf("<unknown>");
      break;
  }
  REprintf("\n");
}
