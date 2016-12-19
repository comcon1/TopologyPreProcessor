#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Topology PreProcessor], [0.1], [comconadin@gmail.com])
AM_INIT_AUTOMAKE([])
AC_CONFIG_SRCDIR([src/runtime.cpp])
AC_CONFIG_HEADERS([config.h])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET

AC_LANG([C++])

echo
echo '>> Checking for common requirements of TPPMKTOP'
echo

# Checks for BOOST::PROGRAM_OPTIONS
AC_CHECK_LIB([boost_program_options], [main],,
             [AC_MSG_ERROR([Can not link to libboost_program_options!])]
)
AC_CHECK_HEADERS([boost/program_options/options_description.hpp boost/program_options/parsers.hpp boost/program_options/variables_map.hpp boost/program_options/errors.hpp],,
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


# checking for openbabel functions
PKG_CHECK_MODULES([OPENBABEL], [openbabel-2.0],
		  [CPPFLAGS="$CPPFLAGS $OPENBABEL_CFLAGS"])

AC_CHECK_LIB([openbabel], [main],,
             [AC_MSG_ERROR([Can not link to libopenbabel!])]
)
AC_CHECK_HEADERS([openbabel/parsmart.h openbabel/obiter.h openbabel/obconversion.h openbabel/mol.h],,
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
AC_CHECK_HEADERS([iostream])

if test x$HAVE_IOSTREAM = x1; then
  AC_MSG_ERROR([Error with STL headers!])
fi

# Check for BOOST header-only libs
AC_CHECK_HEADERS([boost/lexical_cast.hpp boost/cast.hpp boost/date_time/posix_time/posix_time.hpp boost/format.hpp],,
                 [AC_MSG_ERROR([BOOST headers are incomplete!])])

AC_CHECK_HEADERS([boost/spirit/core.hpp])
AC_CHECK_HEADERS([boost/spirit/home/classic/core.hpp])

if test x$HAVE_BOOST_SPIRIT_CORE_HPP = x0; then
  if test x$HAVE_BOOST_SPIRIT_HOME_CLASSIC_CORE_HPP = x0; then
      AC_MSG_ERROR([Unknown version of BOOST::SPIRIT!])
  fi
fi

# Check for BOOST ublas headers
AC_CHECK_HEADERS([boost/numeric/ublas/vector.hpp boost/numeric/ublas/matrix_proxy.hpp boost/numeric/ublas/matrix.hpp boost/numeric/ublas/io.hpp],,
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

  AC_CHECK_LIB([lapack], [sgesvd_],[true],
              [AC_MSG_ERROR([Lapack SVD should work for GAMESS FEATURES!])]
  )
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


# Override LIBS to make program-specific linking in AM
LIBS="-lpthread -lboost_system -lboost_thread  -lboost_log_setup  -lboost_log -lopenbabel -lboost_program_options"

#TODO: make compiler flags conditional
CXXFLAGS="-g"

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

DB_MN="0001a"
AC_DEFINE_UNQUOTED([DB_MAGICNUMBER],["$DB_MN"],
                   [Number of Database version the code is designed for.])
echo
AC_MSG_RESULT([>> DataBase must have magic number (version) $DB_MN])
echo

AC_CONFIG_FILES([
  Makefile
  src/Makefile])
AC_OUTPUT
