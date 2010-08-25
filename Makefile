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
FFLAGS += -I/usr/include

default: $(LIBFILES) $(BINFILES) test

test: rfx_test


$(call LINKLIB, reflex, control.o trajectory.o)
$(call LINKBIN, rfx_test, rfx_test.o trajectory.o control.o, amino stdc++ )



.PHONY: default clean doc

doc:
	doxygen

clean:
	rm -fr *.o $(BINFILES) $(LIBFILES) $(BUILDDIR)/*.o .dep debian *.deb *.lzma rfx_test


