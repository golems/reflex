
if HAVE_HELP2MAN

man/rfx-trajgen.1: $(top_srcdir)/src/demo/rfx-trajgen.c $(top_srcdir)/configure.ac $(top_srcdir)/doc/rfx-trajgen.inc
	$(HELP2MAN) -i rfx-trajgen.inc -h -h -v -V --no-info -n "Generate trajectory" $(top_builddir)/rfx-trajgen$(EXEEXT) -o $@

if HAVE_MAN2HTML
.1.html:
	$(MAN2HTML) -r $< | sed -e '1d;2d;8d;s/@\|&amp;/{at}/g' > $@
endif
endif



if HAVE_MANPAGES
dist_man_MANS = man/rfx-trajgen.1
endif

maintainer-clean-local:
	-rm -rf *.1 *.html man/*.html man/*.1
