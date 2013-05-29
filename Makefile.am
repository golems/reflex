INCLUDES = -I$(top_srcdir)/include

EXTRA_DIST = cartpend.gnuplot

include_HEADERS = include/reflex.h include/reflex.hpp

# pkginclude_HEADERS =

#TESTS =

lib_LTLIBRARIES = libreflex.la

noinst_PROGRAMS = lqgtest cartpend

cartpend_SOURCES = src/cartpend.c
cartpend_LDADD = libreflex.la

lqgtest_SOURCES = src/lqgtest.c
lqgtest_LDADD = libreflex.la

noinst_PROGRAMS += plot_test
plot_test_SOURCES = src/plot_test.c
plot_test_LDADD = libreflex.la

libreflex_la_SOURCES =              \
	src/control.c               \
	src/trajectory.cpp          \
	src/trajectory_plot.cpp     \
	src/lqg.c                   \
	src/plot.c                  \
	src/traj.c

clean-local:
	-rm -rf debian/tmp debian/files debian/libreflex? debian/libreflex-dev debian/*.log debian/*.substvars texput.log debian/*.debhelper

distclean-local:
	-rm -rf doxygen-doc

cp.eps: cartpend cartpend.gnuplot
	./cartpend
	gnuplot ./cartpend.gnuplot

include aminclude.am
