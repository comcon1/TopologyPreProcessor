#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Topology PreProcessor], [0.2], [comconadin@gmail.com])
AM_INIT_AUTOMAKE([])
AC_CONFIG_SRCDIR([src/logger.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET

AC_LANG([C++])

AX_BOOST_BASE([1.60],, [AC_MSG_ERROR([TppMKTOP needs Boost, but it was not found in your system])])
AX_BOOST_THREAD
AX_BOOST_LOG
AX_BOOST_LOG_SETUP
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_SYSTEM
AX_BOOST_CHRONO
AX_BOOST_REGEX
AX_BOOST_DATE_TIME
#AX_BOOST_ATOMIC
AC_CHECK_LIB([boost_atomic], [main],[],
             [AC_MSG_ERROR([Can not link to libboost_atomic!])] )

echo
echo '>> After BOOST checks'
echo

echo $CPPFLAGS
echo $CXXFLAGS
echo $LIBS
echo $LDFLAGS


echo
echo '>> Checking for common requirements of TPPMKTOP'
echo


# Check for libs required for BOOST logger

#AC_MSG_CHECKING([pthread])

# This line works, but AC_PTHREAD creates empty $PTHREAD_CFLAGs and $PTHREAD_LIBS
AC_CHECK_LIB(pthread, pthread_create, [PTHREAD_CFLAGS+=-pthread])
#AC_PTHREAD()


echo 'Pthread flags:'
echo $CXXFLAGS
echo $PTHREAD_CFLAGS
echo $PTHREAD_CXXFLAGS
echo $PTHREAD_CC
echo $PTHREAD_LIBS
echo 'end of flags'

CPPFLAGS="$BOOST_CPPFLAGS $CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
CXXFLAGS="-DBOOST_LOG_DYN_LINK $PTHREAD_CFLAGS $CXXFLAGS"
LIBS="$PTHREAD_LIBS $LIBS $BOOST_LOG_LIB $BOOST_LOG_SETUP_LIB \
 $BOOST_PROGRAM_OPTIONS_LIB $BOOST_FILESYSTEM_LIB $BOOST_SYSTEM_LIB \
 $BOOST_THREAD_LIB $BOOST_ATOMIC_LIB $BOOST_CHRONO_LIB $BOOST_REGEX_LIB \
 $BOOST_DATE_TIME_LIB"

# Not every required header is included into check but I think that's enough

AC_CHECK_HEADERS_ONCE([boost/log/trivial.hpp boost/log/sources/severity_logger.hpp boost/log/utility/setup/file.hpp],,
                [AC_MSG_ERROR([BOOST::LOG headers are not present!])])

AC_MSG_CHECKING([for BOOST::LOG correct behavior])

AC_TRY_LINK([#include <boost/log/trivial.hpp>],
            [
    BOOST_LOG_TRIVIAL(trace) << "A trace severity message";
    BOOST_LOG_TRIVIAL(debug) << "A debug severity message";
    BOOST_LOG_TRIVIAL(info) << "An informational severity message";
            ],
            [AC_MSG_RESULT([works])],
            [AC_MSG_ERROR([Error compiling example BOOST::LOG code! May be wrong version.])]
)

AC_CHECK_HEADERS_ONCE([boost/program_options/options_description.hpp boost/program_options/parsers.hpp boost/program_options/variables_map.hpp boost/program_options/errors.hpp],,
                [AC_MSG_ERROR([BOOST headers are required! Install libboost-dev!])] )

AC_MSG_CHECKING([checking for options_description interface])
AC_TRY_LINK([#include <boost/program_options/options_description.hpp>],
              [boost::program_options::options_description desc("testprg");],
              [AC_MSG_RESULT([works])],
              [AC_MSG_ERROR([Incorrect boost::program_options library!])])


# Check for MySQL++ features
# FIXME: option for redefining mysql include path

AC_CHECK_HEADER([/usr/include/mysql/mysql.h],
                [CPPFLAGS="-I/usr/include/mysql $CPPFLAGS"],
                [AC_MSG_ERROR([Mysql headers are required! Install libmysql-dev!])]
)

# Check for MySQL++
# do not link mysql for all project targets
SAVELIBS="$LIBS"
AC_CHECK_LIB([mysqlpp], [main],[],
             [AC_MSG_ERROR([Can not link to libmysqlpp!])] )
AC_CHECK_HEADER([mysql++/mysql++.h],,
                [AC_MSG_ERROR([Mysql++ headers are required! Install mysqlpp-dev!])],
                [#include <mysql.h>] )

AC_MSG_CHECKING([checking for MySQL++ correct behavior])
AC_TRY_LINK([#include <mysql.h>
#include <mysql++/mysql++.h>],
            [
mysqlpp::Connection con (false);
            ],
            [AC_MSG_RESULT([works])],
            [AC_MSG_ERROR([Incorrect MySQL++ usage! May be wrong version.])]
)
LIBS="$SAVELIBS"

# checking for openbabel functions
PKG_CHECK_MODULES([OPENBABEL], [openbabel-2.0],
		  [CPPFLAGS="$CPPFLAGS $OPENBABEL_CFLAGS"])

AC_CHECK_LIB([openbabel], [main],,
             [AC_MSG_ERROR([Can not link to libopenbabel!])]
)
AC_CHECK_HEADERS_ONCE([openbabel/parsmart.h openbabel/obiter.h openbabel/obconversion.h openbabel/mol.h],,
                [AC_MSG_ERROR([OpenBabel headers are required! Install libopenbabel-dev!])]
)

AC_MSG_CHECKING([checking for OpenBabel correct behavior])
AC_TRY_LINK([#include <openbabel/parsmart.h>],
            [
OpenBabel::OBSmartsPattern pat;
pat.Init("CC(=O)-[NH2]");
            ],
            [AC_MSG_RESULT([works])],
            [AC_MSG_ERROR([Incorrect OpenBabel SMARTS usage! May be wrong version.])]
)


# Checks for STL header
AC_CHECK_HEADERS_ONCE([iostream string sstring exception])

if test x$HAVE_IOSTREAM = x1; then
  AC_MSG_ERROR([Error with STL headers!])
fi

# Check for BOOST header-only libs
AC_CHECK_HEADERS_ONCE([boost/lexical_cast.hpp boost/cast.hpp boost/date_time/posix_time/posix_time.hpp boost/format.hpp],,
                 [AC_MSG_ERROR([BOOST headers are incomplete!])])

AC_CHECK_HEADERS_ONCE([boost/spirit/core.hpp boost/spirit/home/classic/core.hpp])

if test x$HAVE_BOOST_SPIRIT_CORE_HPP = x0; then
  if test x$HAVE_BOOST_SPIRIT_HOME_CLASSIC_CORE_HPP = x0; then
      AC_MSG_ERROR([Unknown version of BOOST::SPIRIT!])
  fi
fi

# Check for BOOST ublas headers
AC_CHECK_HEADERS_ONCE([boost/numeric/ublas/vector.hpp boost/numeric/ublas/matrix_proxy.hpp boost/numeric/ublas/matrix.hpp boost/numeric/ublas/io.hpp],,
                 [AC_MSG_ERROR([BOOST numeric ublas headers are required!])])

# activate GAMESS features (as HESSIAN load)
AC_ARG_ENABLE([gamess-features],
     [  --enable-gamess-features    Turn on GAMESS features <experimental>],
     [case "${enableval}" in
       yes) gamessfe=1 ;;
       no)  gamessfe=0 ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-gamess-features]) ;;
     esac],[gamessfe=0])

AC_DEFINE_UNQUOTED([ENABLE_GAMESS_FEATURES],[$gamessfe],[Turn on GAMESS features.])
AM_CONDITIONAL([GAMESSFEATURES], [test x$gamessfe = x1])

if test x$gamessfe = x1; then

  echo
  echo '>> Checking for requirements for GAMESS FEATURES '
  echo
  AC_CHECK_HEADERS([boost/numeric/bindings/lapack/gesv.hpp boost/numeric/bindings/lapack/gesvd.hpp boost/numeric/bindings/traits/ublas_matrix.hpp],,
                  [AC_MSG_ERROR([BOOST numeric bindings are required for GAMESS FEATURES!])])

  SAVELIBS="$LIBS"
  AC_CHECK_LIB([lapack], [sgesvd_],[true],
              [AC_MSG_ERROR([Lapack SVD should work for GAMESS FEATURES!])]
  )
  LIBS="$SAVELIBS"
fi

# activate TEST features
AC_ARG_ENABLE([selftest-features],
     [  --enable-selftest-features  Turn on self test features (write access to DB is required)],
     [case "${enableval}" in
       yes) selftestfe=1 ;;
       no)  selftestfe=0 ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-selftest-features]) ;;
     esac],[selftestfe=0])

AC_DEFINE_UNQUOTED([ENABLE_SELFTEST_FEATURES],[$selftestfe],[Turn on self test features.])
AM_CONDITIONAL([SELFTESTFEATURES], [test x$selftestfe = x1])

if test x$selftestfe = x1; then

  echo
  echo '>> Checking for requirements for SELFTEST FEATURES '
  echo

  AC_MSG_RESULT([>> No specific tests are required yet.])

fi

#TODO: make compiler flags conditional
CXXFLAGS="-g $CXXFLAGS"

# FIXME: unit test with ENABLE option

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

tmpvar=`LANG="C" date`
AC_DEFINE_UNQUOTED([CONFIGURE_CDATE],["$tmpvar"],[Date of last confiugration.])
#CONFIGURE_CDATE=`LANG="C" date`

# Database version format:
#
# NNNNx (N - digit, x - letter)
#
# NNNN - fatal database update
# x - non-fatal update

DB_MN="0002a"
AC_DEFINE_UNQUOTED([DB_MAGICNUMBER],["$DB_MN"],
                   [Number of Database version the code is designed for.])
echo
AC_MSG_RESULT([>> DataBase must have magic number (version) $DB_MN])
echo

AC_CONFIG_FILES([
  Makefile
  src/Makefile])
AC_OUTPUT
