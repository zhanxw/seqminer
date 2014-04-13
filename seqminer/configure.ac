AC_INIT(seqminer, 2.0, zhanxw@gmail.com)
AC_CONFIG_HEADERS(src/config.h)
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

dnl skip checking C compilers as we also need to check C++ compiler
dnl AC_LANG(C)
dnl AC_PROG_CC
AC_LANG(C++)
AC_PROG_CXX

LIB_BZIP2=""
AC_CHECK_LIB([bz2], [BZ2_bzDecompressInit], [
  AC_CHECK_HEADER([bzlib.h], [
    HAVE_BZIP2=1
    AC_DEFINE([HAVE_BZIP2], [1], [Define to 1 if you have bzip2 library])
    LIB_BZIP2="-lbz2"
  ]) 
])

AC_CHECK_LIB([c], [regcomp], [
  AC_CHECK_HEADER([regex.h], [
    HAVE_POSIX_REGEX=1
    AC_DEFINE([HAVE_POSIX_REGEX], [1], [Define to 1 if you have POSIX Regular Expressionlibrary])
  ])
])

dnl copy from http://ftp.netbsd.org/pub/pkgsrc/current/pkgsrc/textproc/nbsed/files/configure.ac
AC_MSG_CHECKING([for working REG_STARTEND])
AC_TRY_RUN([
#undef NDEBUG
#include <assert.h>
#include <regex.h>
int main(void) {
	regex_t re;
	regmatch_t rm;

	assert(regcomp(&re, "lo", 0) == 0);
	rm.rm_so = 0, rm.rm_eo = 7;
	assert(regexec(&re, "hel\0lo\n", 1, &rm, REG_STARTEND) == 0);
	assert(rm.rm_so == 4 && rm.rm_eo == 6);
	return 0;
}], [AC_DEFINE([HAVE_REG_STARTEND], [1], [Do you have a working REG_STARTEND?])
AC_MSG_RESULT([yes])], [AC_MSG_RESULT([no])])


dnl AC_SUBST(LIB_BZIP2)
dnl -I. @CFLAGS@ -DNDEBUG 
dnl -I. @CXXFLAGS@ -DNDEBUG
dnl -lz -lm @LIB_BZIP2@

AC_SUBST([PKG_CFLAGS], ["-I. -DNDEBUG ${CFLAGS}"])
AC_SUBST([PKG_CPPFLAGS], ["-I. -DNDEBUG ${CPPFLAGS} ${CXXFLAGS}"])
AC_SUBST([PKG_LIBS], ["-lz -lm ${LIB_BZIP2}"])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
