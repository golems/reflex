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

m4_ifdef([AX_CXX_HEADER_STDCXX_11],
         [AX_CXX_HEADER_STDCXX_11],
         [AX_CXX_COMPILE_STDCXX_0X
         if test x"$ax_cv_cxx_compile_cxx0x_gxx" = xyes; then CXXFLAGS="$CXXFLAGS -std=gnu++0x"; fi])


# Enable C++11
dnl AC_LANG([C++])
dnl AC_LANG_CPLUSPLUS
m4_ifdef([AX_CXX_HEADER_STDCXX_11],
         [AX_CXX_HEADER_STDCXX_11],
         [AX_CXX_COMPILE_STDCXX_0X
         if test x"$ax_cv_cxx_compile_cxx0x_gxx" = xyes; then CXXFLAGS="$CXXFLAGS -std=gnu++0x"; fi])

# Checks for libraries.
#AC_SEARCH_LIBS([cblas_dgemm],[blas], [true], [echo "Failed: need BLAS library" && exit 1])
#AC_SEARCH_LIBS([dgesvd_],[lapack], [true], [echo "Failed: need LAPACK library" && exit 1])
#AC_SEARCH_LIBS([aa_la_dlsnp],[amino])


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


AC_SEARCH_LIBS([dgemm_],[openblas blas])
AC_SEARCH_LIBS([dgesvd_],[lapack])

# Default Install Prefix
AS_IF( [test "x$prefix" = xNONE], [prefix=/usr/local] )

# Find amino fortran modules
FCFLAGS="$FCFLAGS -I$prefix/include"

# Set LP solve shared library location
AC_SEARCH_LIBS([solve], [lpsolve55],
	       [],
	       [LDFLAGS="-L/usr/lib/lp_solve $LDFLAGS"
		$as_unset ac_cv_search_solve
		AC_SEARCH_LIBS([solve], [lpsolve55],
			       [],
			       [AC_MSG_ERROR([unable to find lpsolve])])
	       ])

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


AC_CHECK_PROGS([HELP2MAN], [help2man])
AC_CHECK_PROGS([MAN2HTML], [man2html])
AM_CONDITIONAL([HAVE_HELP2MAN], [test x$HELP2MAN != x])
AM_CONDITIONAL([HAVE_MAN2HTML], [test x$MAN2HTML != x])

AC_CHECK_FILE( [$srcdir/doc/man/rfx-trajgen.1],           [FOUND_MANPAGES=yes])
AC_CHECK_FILE( [$srcdir/doc/man/rfx-trajgen.html],        [FOUND_MANHTML=yes])
AM_CONDITIONAL([HAVE_MANPAGES], [test x$HELP2MAN != x -o x$FOUND_MANPAGES = xyes])
AM_CONDITIONAL([HAVE_MANHTML],  [test '(' x$HELP2MAN != x -a x$MAN2HTML != x ')' -o x$FOUND_MANHTML = xyes])


dnl JAVA
AS_IF( [test "x$JAVAPREFIX" = x],
       [
        AS_IF( [test "x$JAVA" = x],  [AC_CHECK_PROGS(JAVA, java kaffe)],  )
        AS_IF( [test "x$JAVAC" = x], [AC_CHECK_PROGS(JAVAC, javac "gcj -C" jikes guavac)] )
        AS_IF( [test "x$JAVAH" = x], [AC_CHECK_PROGS(JAVAH, javah)] )
        AS_IF( [test "x$JAR" = x],   [AC_CHECK_PROGS(JAR, jar)] )
       ],
       [
        AS_IF( [test "x$JAVA" = x],  [AC_CHECK_PROGS(JAVA, java kaffe)],
               [$JAVAPREFIX] )
        AS_IF( [test "x$JAVAC" = x], [AC_CHECK_PROGS(JAVAC, javac "gcj -C" jikes guavac)],
               [$JAVAPREFIX] )
        AS_IF( [test "x$JAVAH" = x], [AC_CHECK_PROGS(JAVAH, javah)],
               [$JAVAPREFIX] )
        AS_IF( [test "x$JAR" = x],   [AC_CHECK_PROGS(JAR, jar)],
               [$JAVAPREFIX] )
       ]
)

AS_IF([test "x$JAVAC" = "x"],
      [BUILD_JAVA=no],
      [
       BUILD_JAVA=yes
       AX_JNI_INCLUDE_DIR
       for JNI_INCLUDE_DIR in $JNI_INCLUDE_DIRS
       do
           JNICPPFLAGS="$JNICPPFLAGS -I$JNI_INCLUDE_DIR"
       done
       dnl Guess the classpath
       AS_IF([test "x$CLASSPATH" = x],
             [CLASSPATH=.:$prefix/share/java/\*])
       AC_SUBST([CLASSPATH])
      ])

AC_SUBST([JNICPPFLAGS])
AM_CONDITIONAL([HAVE_JAVAC],   [test "x$JAVAC" != x])

dnl TODO: better (any) handling of CLASSPATH variable

AC_CONFIG_FILES([Makefile doc/Makefile])
AC_OUTPUT

AC_MSG_NOTICE([CONFIGURATION SUMMARY])
AC_MSG_NOTICE([=====================])
AC_MSG_NOTICE([PREFIX:          $prefix])
AC_MSG_NOTICE([BUILD JAVA LIB:  $BUILD_JAVA])
AC_MSG_NOTICE([CLASSPATH:       $CLASSPATH])
AC_MSG_NOTICE([CFLAGS:          $CFLAGS])
AC_MSG_NOTICE([FCFLAGS:         $FCFLAGS])
