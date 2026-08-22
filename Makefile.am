ACLOCAL_AMFLAGS = -I m4 -I config
SUBDIRS=src

dep depclean:
	$(MAKE) -C src $@

# Run the standalone hot-kernel regression tests (see test/hotkernels).
test-hotkernels:
	$(MAKE) -C test/hotkernels run

.PHONY: test-hotkernels
