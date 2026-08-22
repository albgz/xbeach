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
AM_FCFLAGS = -ffree-form -I$(LIBDIR)

hotkernels_test: hotkernels_test.f90
	$(FC) $(AM_FCFLAGS) -g -O2 hotkernels_test.f90 $(LIBS) -o $@

hotkernels2_test: hotkernels2_test.f90
	$(FC) $(AM_FCFLAGS) -g -O2 hotkernels2_test.f90 $(LIBS) -o $@

run: hotkernels_test hotkernels2_test
	DYLD_LIBRARY_PATH=$(LIBDIR)/.libs ./hotkernels_test
	DYLD_LIBRARY_PATH=$(LIBDIR)/.libs ./hotkernels2_test

all: run

clean:
	rm -f hotkernels_test hotkernels2_test

.PHONY: run all clean test
test: run
