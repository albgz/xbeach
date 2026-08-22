AM_FCFLAGS = -I$(top_builddir)/src/xbeachlibrary

check_PROGRAMS = hotkernels_test hotkernels2_test
TESTS = $(check_PROGRAMS)

hotkernels_test_SOURCES = hotkernels_test.f90
hotkernels2_test_SOURCES = hotkernels2_test.f90

LDADD = $(top_builddir)/src/xbeachlibrary/libxbeach.la
