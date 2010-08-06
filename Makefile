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

default: $(LIBFILES) $(BINFILES)


$(call LINKLIB, reflex, control.o trajectory.o)

.PHONY: default clean doc

doc:
	doxygen

clean:
	rm -fr *.o $(BINFILES) $(LIBFILES) $(BUILDDIR)/*.o .dep debian *.deb *.lzma


