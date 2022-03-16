#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([SpherePack], [3.2], [https://github.com/NCAR/NCAR-Classic-Libraries-for-Geophysics])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
LT_INIT()

# Checks for programs.
AC_PROG_F77

# Checks for libraries.
AC_CHECK_LIB([ncarg], [set])
AC_CHECK_LIB([fftpack], [costi])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AM_CONDITIONAL([HAVE_FFTPACK], [test x$HAVE_LIBFFTPACK != x])
AM_COND_IF([HAVE_FFTPACK], [AC_CHECK_FUNCS([sinti costi cosqi sint cost cosqb sinqb])])

AM_CONDITIONAL([HAVE_NCAR_GRAPHICS], [test x$HAVE_LIBNCARG != x])
AM_COND_IF([HAVE_NCAR_GRAPHICS], [AC_CHECK_FUNCS([set line frame])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
