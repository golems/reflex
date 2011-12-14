#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([reflex], [0.20110914], [ntd@gatech.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([lqg.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_CC_C99

# Checks for libraries.
AC_SEARCH_LIBS([cblas_dgemm],[blas], [true], [echo "Failed: need BLAS library" && exit 1])
AC_SEARCH_LIBS([dgesvd_],[lapack], [true], [echo "Failed: need LAPACK library" && exit 1])
AC_SEARCH_LIBS([aa_la_dlsnp],[amino])


# Checks for header files.
AC_CHECK_HEADERS([inttypes.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T

DX_HTML_FEATURE(ON)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(reflex)


# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([memset])

# Enable maximum warnings
AX_CFLAGS_WARN_ALL

dnl APPEND_FLAG( [FLAG] )
dnl Check if FLAG is supported, and if so, add it
AC_DEFUN([APPEND_FLAG],
         [AX_CHECK_COMPILER_FLAGS([$1], [CFLAGS="$1 $CFLAGS"])])

APPEND_FLAG([-Wextra])
APPEND_FLAG([-Wconversion])
APPEND_FLAG([-Wpointer-arith])
APPEND_FLAG([-Wfloat-equal])
APPEND_FLAG([-Wshadow])
APPEND_FLAG([-Wwrite-strings])
APPEND_FLAG([-Wc++-compat])


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
