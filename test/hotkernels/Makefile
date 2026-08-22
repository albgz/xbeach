# Standalone regression tests for the hot wave kernels.
#
#   hotkernels_test   - the shoaling advection kernels (advecxho/advecyho/
#                       advecthetaho), pure-dummy interfaces
#   hotkernels2_test  - the state-carrying kernels (iteratedispersion,
#                       compute_wave_direction_velocities, wave_dispersion)
#
# `make` (or `make run`) builds and runs both against the already-built
# libxbeach. They need the library and its .mod files, so run a full
# `make` at the tree root first.

LIBDIR = $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../../src/xbeachlibrary)
LIBS   = -L$(LIBDIR)/.libs -lxbeach -Wl,-rpath,"$(LIBDIR)/.libs"
FC     = gfortran
FCFLAGS = -g -O2
CPPFLAGS = -I$(LIBDIR)
AM_FCFLAGS = -ffree-form

hotkernels_test: hotkernels_test.f90
	$(FC) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) hotkernels_test.f90 $(LDFLAGS) $(LIBS) -o $@

hotkernels2_test: hotkernels2_test.f90
	$(FC) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) hotkernels2_test.f90 $(LDFLAGS) $(LIBS) -o $@

run: hotkernels_test hotkernels2_test
	DYLD_LIBRARY_PATH=$(LIBDIR)/.libs ./hotkernels_test
	DYLD_LIBRARY_PATH=$(LIBDIR)/.libs ./hotkernels2_test

all: run

clean:
	rm -f hotkernels_test hotkernels2_test

.PHONY: run all clean test
test: run
