#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([cytolib], [2.0], [wjiang2@fhcrc.org])
AC_CONFIG_SRCDIR([cytolibConfig.h.in])
AC_PROG_SED

AC_MSG_NOTICE([Parsing cytolib version info ...])
CYTOLIB_VERSION=`bash getverion.sh`

AC_MSG_NOTICE([CYTOLIB_VERSION=${CYTOLIB_VERSION}])
AC_SUBST(CYTOLIB_VERSION)
AC_CONFIG_FILES([inst/include/cytolib/cytolibConfig.h:cytolibConfig.h.in])
AC_OUTPUT
