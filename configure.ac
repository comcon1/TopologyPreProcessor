#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Topology PreProcessor], [0.1], [comconadin@gmail.com])
AM_INIT_AUTOMAKE([])
AC_CONFIG_SRCDIR([src/global.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET

AC_LANG([C++])

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
AC_CHECK_LIB([mysqlpp], [main],,
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
AC_CHECK_LIB([openbabel], [main],,
             [AC_MSG_ERROR([Can not link to libopenbabel!])]
)
AC_CHECK_HEADERS([openbabel/parsmart.h openbabel/obiter.h openbabel/obconversion.h],,
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


#if test ${}

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

AC_CONFIG_FILES([
  Makefile
  src/Makefile])
AC_OUTPUT
