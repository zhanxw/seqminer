#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "R.h"

#if 0
/**
 * The following is taken from:
 * http://c.learncodethehardway.org/book/learn-c-the-hard-waych21.html#x26-10700021.2
 */
#ifdef NDEBUG
#define debug(M, ...)
#else 
#define debug(M, ...) REprintf("DEBUG %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__) 
#endif

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
 
#define log_error(M, ...) REprintf("[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__) 
 
#define log_warn(M, ...) REprintf("[WARN] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__) 
 
#define log_info(M, ...) REprintf("[INFO] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__) 

#define check(A, M, ...) if(!(A)) { log_err(M, ##__VA_ARGS__); errno=0; goto error; } 

#define sentinel(M, ...)  { log_err(M, ##__VA_ARGS__); errno=0; goto error; } 

#define check_mem(A) check((A), "Out of memory.")
 
#define check_debug(A, M, ...) if(!(A)) { debug(M, ##__VA_ARGS__); errno=0; goto error; } 
#else
inline void log_error(const char* s) {
  REprintf("[ERROR] %s\n", s);
}

#endif

inline void REPORT(const char* x) { 
    REprintf("Report '%s' to zhanxw@umich.edu\n", x ); 
}

inline void FATAL(const char* x) {
    REPORT(x);
    REprintf("Critical error happening!\n"); //abort();
}

#endif /* _EXCEPTION_H_ */
