#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# This autoconf based configure is replaced by more portable R package 'configure' 
AC_PREREQ([2.65])
AC_INIT([cytolib], [2.0], [mike@ozette.ai])
AC_CONFIG_SRCDIR([cytolibConfig.h.in])
AC_PROG_SED

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo "could not determine R_HOME"
    exit 1
fi

CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX11FLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

	   
AC_MSG_NOTICE([Parsing cytolib version info ...])
CYTOLIB_VERSION=`bash getverion.sh`

AC_MSG_NOTICE([CYTOLIB_VERSION=${CYTOLIB_VERSION}])
AC_SUBST(CYTOLIB_VERSION)
AC_CONFIG_FILES([inst/include/cytolib/cytolibConfig.h:cytolibConfig.h.in])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
