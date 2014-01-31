INCLUDES = -I$(top_srcdir)/include

SUBDIRS = . doc

EXTRA_DIST = cartpend.gnuplot

include_HEADERS = include/reflex.h include/reflex.hpp

pkginclude_HEADERS = \
	include/reflex/kinematics.h \
	include/reflex/trajq.h \
	include/reflex/trajx.h \
	include/reflex/tf.h \
	include/reflex/lqg.h \
	include/reflex/body.h

nodist_include_HEADERS = reflex.mod

# pkginclude_HEADERS =

#TESTS =

lib_LTLIBRARIES = libreflex.la

libreflex_la_SOURCES =              \
	src/control.c               \
	src/lqg.c                   \
	src/lqg_tf.c                \
	src/rfx_tf.c                \
	src/tf/dud.c                \
	src/plot.c                  \
	src/trajq.c                 \
	src/kin.c                   \
	src/trajx.cpp               \
	src/traj.c                  \
	src/body/body.cpp

libreflex_la_LIBADD = libreflex_mod.la

noinst_PROGRAMS = cartpend
cartpend_SOURCES = src/cartpend.c
cartpend_LDADD = libreflex.la -lamino -llapack -lblas

noinst_PROGRAMS += lqgtest
lqgtest_SOURCES = src/lqgtest.c
lqgtest_LDADD = libreflex.la -lamino -llapack -lblas

noinst_PROGRAMS += plot_test
plot_test_SOURCES = src/plot_test.c
plot_test_LDADD = libreflex.la -lamino -llapack -lblas


noinst_PROGRAMS += test-tf-filter
test_tf_filter_SOURCES = src/test/test-tf-filter.c
test_tf_filter_LDADD = libreflex.la -lamino -llapack -lblas


bin_PROGRAMS = rfx-trajgen
rfx_trajgen_SOURCES = src/demo/rfx-trajgen.c
rfx_trajgen_LDADD = libreflex.la -lamino -lblas -llapack


bin_PROGRAMS += rfx-camcal
rfx_camcal_SOURCES = src/demo/rfx-camcal.c
rfx_camcal_LDADD = libreflex.la -lamino -lblas -llapack

reflex_mod.f90: src/reflex_mod.F90 src/kinematics.f90 src/lqg_tf_f.f90 src/tf/rfx_tff.f90
	$(CPP) $(INCLUDES) -x c  $< -o $@

noinst_LTLIBRARIES = libreflex_mod.la
libreflex_mod_la_SOURCES = reflex_mod.f90
reflex.mod: $(libreflex_mod_la_OBJECTS)


clean-local:
	-rm -rf debian/tmp debian/files debian/libreflex? debian/libreflex-dev debian/*.log debian/*.substvars texput.log debian/*.debhelper *.jar org

distclean-local:
	-rm -rf doxygen-doc

cp.eps: cartpend cartpend.gnuplot
	./cartpend
	gnuplot ./cartpend.gnuplot



tf-filter.eps: test-tf-filter src/test/plot-tf-filter
	./test-tf-filter
	gnuplot ./src/test/plot-tf-filter

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
	$(JAVAC) -classpath $(CLASSPATH)  -d $(top_builddir) $(JSOURCES)

org_golems_reflex_Lib.h: $(top_builddir)/org/golems/reflex/Lib.class
	$(JAVAH) org.golems.reflex.Lib
	touch $@

$(top_builddir)/org/golems/reflex/Lib.class: org.golems.reflex.jar

org.golems.reflex.jar: $(JSOURCES)
	$(JAVAC) -classpath $(CLASSPATH) -d $(top_builddir) $(JSOURCES)
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
