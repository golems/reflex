INCLUDES = -I$(top_srcdir)/include

EXTRA_DIST = cartpend.gnuplot

include_HEADERS = include/reflex.h include/reflex.hpp

pkginclude_HEADERS = include/reflex/kinematics.h include/reflex/trajectory.h include/reflex/trajx.h

nodist_include_HEADERS = reflex_kin.mod

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
	src/lqg.c                   \
	src/lqgf.f90                \
	src/plot.c                  \
	src/trajq.c                 \
	src/kin.c                   \
	src/trajx.c                 \
	src/traj.c

libreflex_la_LIBADD =               \
	libreflex_kin.la

noinst_LTLIBRARIES =                \
	libreflex_kin.la

libreflex_kin_la_SOURCES = src/kinematics.f90

reflex_kin.mod: $(libreflex_kin_la_OBJECTS)

clean-local:
	-rm -rf debian/tmp debian/files debian/libreflex? debian/libreflex-dev debian/*.log debian/*.substvars texput.log debian/*.debhelper

distclean-local:
	-rm -rf doxygen-doc

cp.eps: cartpend cartpend.gnuplot
	./cartpend
	gnuplot ./cartpend.gnuplot


## JAVA ##
JSOURCES =                                        \
	$(top_srcdir)/java/Lib.java               \
	$(top_srcdir)/java/Test.java

EXTRA_DIST += $(JSOURCES) doc/javadoc

if HAVE_JAVAC

BUILT_SOURCES = org_golems_reflex_Lib.h

JCLASSES =                                        \
	org/golems/reflex/Lib.class               \
	org/golems/reflex/Test.class

org/golems/reflex/Lib.class: $(JSOURCES)
	$(JAVAC) -d $(top_builddir) $(JSOURCES)

org_golems_reflex_Lib.h: $(top_builddir)/org/golems/reflex/Lib.class
	$(JAVAH) org.golems.reflex.Lib
	touch $@

$(top_builddir)/org/golems/reflex/Lib.class: org.golems.reflex.jar

org.golems.reflex.jar: $(JSOURCES)
	$(JAVAC) -d $(top_builddir) $(JSOURCES)
	$(JAR) cf $@ $(JCLASSES)


# install the jar file somewhere reasonable
javadir = $(datadir)/java
java_DATA = org.golems.reflex.jar

lib_LTLIBRARIES += libreflexj.la
libreflexj_la_CPPFLAGS = $(CPPFLAGS) $(JNICPPFLAGS)
libreflexj_la_SOURCES = java/reflexj.c
libreflexj_la_LIBADD = libreflex.la

# Javadocs
doc/.javadoc.stamp: $(JSOURCES)
	javadoc -d doc/javadoc $(JSOURCES)
	touch $@

doc/javadoc: doc/.javadoc.stamp

javadoc: doc/.javadoc.stamp


endif


include aminclude.am
