define([pkgversion], esyscmd([sh -c "grep Version: DESCRIPTION | cut -d' ' -f2 | tr -d '\n'"]))dnl
AC_INIT(seqminer, [pkgversion], zhanxw@gmail.com)
AC_CONFIG_HEADERS(src/config.h)
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi

dnl adopted from: https://github.com/cran/dfmta/blob/master/configure.ac
AC_DEFUN([AX_GET_R_CONF_VAR],[$1=`"${R_HOME}/bin/R" CMD config $1`])
AX_GET_R_CONF_VAR([CC])
AX_GET_R_CONF_VAR([CFLAGS])
AX_GET_R_CONF_VAR([CXX11])
AX_GET_R_CONF_VAR([CXX11STD])
AX_GET_R_CONF_VAR([CXX11FLAGS])
CXX="${CXX11} ${CXX11STD}"

dnl print out flags to ease debugging
AC_MSG_NOTICE([CC = $CC])
AC_MSG_NOTICE([CFLAGS = $CFLAGS])
AC_MSG_NOTICE([CXX = $CXX11 $CXX11STD])
AC_MSG_NOTICE([CXXFLAGS = $CXX11FLAGS])

dnl skip checking C compilers as we also need to check C++ compiler
dnl AC_LANG(C)
dnl AC_PROG_CC
AC_LANG(C++)
AC_PROG_CXX

dnl check libraries
dnl refer:
dnl   https://www.gnu.org/savannah-checkouts/gnu/autoconf/manual/autoconf-2.69/html_node/Libraries.html
dnl   AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
LIB_BZIP2=""
AC_CHECK_LIB([bz2], [BZ2_bzDecompressInit], [
  AC_CHECK_HEADER([bzlib.h], [
    HAVE_BZIP2=1
    AC_DEFINE([HAVE_BZIP2], [1], [Define to 1 if you have bzip2 library])
    LIB_BZIP2="-lbz2"
    AC_MSG_NOTICE([Will use system-wide bzip2 library])    
  ]) 
])

LIB_ZSTD=""
AC_CHECK_LIB([zstd], [ZSTD_compress], [
  AC_CHECK_HEADERS([zstd.h], [
    HAVE_ZSTD=1
    AC_DEFINE([HAVE_ZSTD], [1], [Define to 1 if you have zstd library])
    LIB_ZSTD="-lzstd"
    AC_MSG_NOTICE([Will use system-wide zstd library])    
  ])                     
], [AC_MSG_NOTICE([Will build zstd library from the package])])
             

LIB_SQLITE=""
AC_CHECK_LIB([sqlite], [sqlite3_open], [
  AC_CHECK_HEADERS([sqlite3.h], [
    HAVE_SQLITE=1
    AC_DEFINE([HAVE_SQLITE], [1], [Define to 1 if you have sqlite library])
    LIB_SQLITE="-lsqlite"
    AC_MSG_NOTICE([Will use system-wide sqlite library])
  ])                     
], [AC_MSG_NOTICE([Will build sqlite library from the package])])

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

dnl compile sqlite-amalgamation-3250300
dnl specify sqlite.o in Makevars.in

dnl dnl compile zstd-1.3.7
dnl AC_MSG_NOTICE([Finished building zstd library])
dnl base_dir=$(pwd)
dnl ZSTD_CFLAGS=$("${R_HOME}/bin/R" CMD config CFLAGS)
dnl ZSTD_LIBS="$base_dir/src/zstd/lib/libzstd.a"
dnl AC_CONFIG_COMMANDS([zstd compiling],
dnl                    [make -C src/zstd/lib lib CFLAGS="-fPIC ${ZSTD_CFLAGS}"])


AC_SUBST([PKG_CFLAGS], ["-D_USE_KNETFILE -I. ${CFLAGS}"])
AC_SUBST([PKG_CPPFLAGS], ["-D_USE_KNETFILE -I. ${CXX11FLAGS}"])
AC_SUBST([PKG_LIBS], ["-lz -lm ${LIB_BZIP2} ${LIB_ZSTD} ${LIB_SQLITE}"])
AC_SUBST([HAVE_BZIP2])
AC_SUBST([HAVE_ZSTD])
AC_SUBST([HAVE_SQLITE])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
