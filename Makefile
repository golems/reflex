# Project Name
PROJECT := reflex

# Project Version
VERSION := 20100801

# Binary Files
#BINFILES :=

# Library files
SHAREDLIBS := reflex

all: default

include /usr/share/make-common/common.1.mk

#CFLAGS += -O0 -Wno-conversion
CFLAGS += --std=gnu99 -O2
FFLAGS += -I/usr/include -O2 -fPIC

default: $(LIBFILES) $(BINFILES) test lqgtest

test: rfx_test lqgtest cartpend cp.eps
	./lqgtest

OBJS := control.o trajectory.o trajectory_plot.o lqg.o

$(call LINKLIB, reflex, $(OBJS))
$(call LINKBIN, rfx_test, rfx_test.o $(OBJS), amino stdc++ blas lapack rt)
$(call LINKBIN, lqgtest, lqgtest.o $(OBJS), amino stdc++ blas lapack rt)
$(call LINKBIN, cartpend, cartpend.o $(OBJS), amino stdc++ blas lapack rt)

cp.eps: cartpend cartpend.gnuplot
	./cartpend
	gnuplot ./cartpend.gnuplot


.PHONY: default clean doc

doc:
	doxygen

clean:
	rm -fr *.o $(BINFILES) $(LIBFILES) $(BUILDDIR)/*.o .dep debian *.deb *.lzma rfx_test


