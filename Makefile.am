ACLOCAL_AMFLAGS = -I m4 -I config
SUBDIRS = src test/hotkernels

EXTRA_DIST = \
	makedepo.py \
	scripts/generate.py \
	tools/benchmark_delilah.py \
	case-delilah/README.md \
	case-delilah/oracle.json \
	case-delilah/params.txt \
	case-delilah/params_original.txt \
	case-delilah/RF_table.txt \
	case-delilah/tide.txt \
	case-delilah/vardens.txt \
	case-delilah/x.grd \
	case-delilah/y.grd \
	case-delilah/zneumann.grd

dep depclean:
	$(MAKE) -C src $@

# Run the hot-kernel regression tests through Automake's test driver.
test-hotkernels:
	$(MAKE) -C test/hotkernels check

.PHONY: test-hotkernels
