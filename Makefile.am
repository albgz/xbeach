ACLOCAL_AMFLAGS = -I m4 -I config
SUBDIRS = src test/hotkernels

dep depclean:
	$(MAKE) -C src $@

# Run the hot-kernel regression tests through Automake's test driver.
test-hotkernels:
	$(MAKE) -C test/hotkernels check

.PHONY: test-hotkernels
