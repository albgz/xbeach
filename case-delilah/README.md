# Deterministic DELILAH modernization benchmark

This directory freezes a deterministic derivative of the 177 x 70,
nine-direction DELILAH surfbeat case used to validate the Fortran
modernization. It is not the original case unchanged: output scheduling and
the boundary-phase seed were deliberately made deterministic.

- `params_original.txt` preserves the imported parameter values and original
  NetCDF/point-output request (with line endings normalized).
- `params.txt` is the benchmark definition. It disables randomized boundary
  phases (`random = 0`), uses the built-in Fortran output backend, and writes
  `H`, `E`, `Fx`, and `Fy` once at `t = 3800 s`.

From a configured and built source tree, run:

```sh
python3 tools/benchmark_delilah.py --warmups 1 --runs 3 \
  --work-root /tmp/xbeach-delilah-$(date +%s) --output result.json
```

For a quick harness check without changing the committed workload:

```sh
python3 tools/benchmark_delilah.py --tstop 20 --runs 2 \
  --work-root /tmp/xbeach-delilah-smoke-$(date +%s)
```

The runner creates a clean directory per run, measures wall time, records the
model-reported duration and timestep count, hashes every retained field, and
fails if retained runs do not agree. It validates `dims.dat`, field dimensions,
FP64 value counts, finiteness, and normal termination. Its JSON also records
the source and runner trees, working-diff identity, input/executable/library
hashes, dynamic dependencies, configured compiler/flags, CPU affinity, load,
warnings, and host. Existing work roots are rejected rather than overwritten.
Promotion evidence should come from a clean tree; dirty-tree reports are for
candidate experiments. Compare performance only when the workload, compiler,
host, timestep count, output sizes, and all four hashes match.
