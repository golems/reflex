#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([reflex], [0.20120504], [ntd@gatech.edu])
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([include/reflex.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_PROG_FC


# Enable C++11
dnl AC_LANG([C++])
dnl AC_LANG_CPLUSPLUS
m4_ifdef([AX_CXX_HEADER_STDCXX_11],
         [AX_CXX_HEADER_STDCXX_11],
         [AX_CXX_COMPILE_STDCXX_0X
         if test x"$ax_cv_cxx_compile_cxx0x_gxx" = xyes; then CXXFLAGS="$CXXFLAGS -std=gnu++0x"; fi])

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

m4_ifdef([AX_CHECK_COMPILE_FLAG],
         [AC_DEFUN([APPEND_CFLAG],
                   [AX_CHECK_COMPILE_FLAG([$1], [CFLAGS="$1 $CFLAGS"])])
          AC_DEFUN([APPEND_FCFLAG],
                   [AX_CHECK_COMPILE_FLAG([$1], [FCFLAGS="$1 $FCFLAGS"])])])

m4_ifdef([AX_CHECK_COMPILER_FLAGS],
         [AC_DEFUN([APPEND_CFLAG],
                   [AX_CHECK_COMPILER_FLAGS([$1], [CFLAGS="$1 $CFLAGS"])])
          AC_DEFUN([APPEND_FCFLAG],
                   [AX_CHECK_COMPILER_FLAGS([$1], [FCFLAGS="$1 $FCFLAGS"])])])


APPEND_CFLAG([-Wextra])
APPEND_CFLAG([-Wconversion])
APPEND_CFLAG([-Wpointer-arith])
APPEND_CFLAG([-Wfloat-equal])
APPEND_CFLAG([-Wshadow])
APPEND_CFLAG([-Wwrite-strings])
APPEND_CFLAG([-Wc++-compat])


AC_LANG([Fortran])
APPEND_FCFLAG([-fimplicit-none])
APPEND_FCFLAG([-ffree-line-length-none])
APPEND_FCFLAG([-Wall])
APPEND_FCFLAG([-pedantic])
APPEND_FCFLAG([-Wextra])
APPEND_FCFLAG([-Wconversion])
APPEND_FCFLAG([-Wconversion-extra])
APPEND_FCFLAG([-Wshadow])
APPEND_FCFLAG([-Wimplicit-interface])
APPEND_FCFLAG([-Wimplicit-procedure])
APPEND_FCFLAG([-Wunderflow])
APPEND_FCFLAG([-Warray-temporaries])
APPEND_FCFLAG([-Wtabs])
APPEND_FCFLAG([-Wunused-parameter])
APPEND_FCFLAG([-Wintrinsic-shadow])


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
