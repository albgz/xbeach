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
  --work-root /tmp/xbeach-delilah-$(date +%s) \
  --output /tmp/xbeach-delilah-result-$(date +%s).json
```

For a quick harness check without changing the committed workload:

```sh
python3 tools/benchmark_delilah.py --tstop 20 --runs 2 \
  --work-root /tmp/xbeach-delilah-smoke-$(date +%s)
```

The runner creates a clean directory per run and emits a schema-4 JSON record.
It fails closed when a result root is reused; the frozen `177 x 70 x 9`,
`random=0`, Fortran-output contract or requested output time is not observed;
outputs are missing, non-finite, wrongly sized, or hash-inconsistent; XBeach
terminates abnormally; or an unclassified warning appears. It records both the
source fixture and effective per-run input hashes, so `--tstop` smoke overrides
cannot be confused with the full workload. Before and after every sample it
binds the clean/dirty source and runner trees, runner, compiler, ELF executable,
actually loaded `libxbeach.so`, and dynamic dependency closure. Every retained
science field is hashed with full SHA-256. For the full `tstop=3800`
workload, the runner also fails unless timesteps, dimensions, sizes, and all four
field hashes match `oracle.json`, frozen from the clean first current-compiler
build-compatible commit (`8ccf245`). Short smoke overrides record that the
oracle was not applied and cannot be mistaken for full regression evidence.
Compiler provenance is mandatory and includes the configured command, wrappers,
underlying driver, GNU Fortran frontend where available, flags, version output,
and hashes. Result publication is rejected inside either reviewed repository;
an external temporary result is written, identities are sealed again, and the
file is then atomically published.

Promotion evidence should come from a clean tree; dirty-tree reports are for
candidate experiments. Compare performance only when the workload, compiler,
host, timestep count, output sizes, exact runtime artifacts, and all four hashes
match.
