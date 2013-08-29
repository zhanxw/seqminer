#include <R.h>
void bz_internal_error ( int errcode ){
  REprintf("[ERROR] Bzip2 error %d", errcode);
}
