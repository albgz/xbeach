#!/usr/bin/env python3
"""Run the frozen DELILAH case and record comparable timing and output hashes."""

from __future__ import annotations

import argparse
import array
import hashlib
import json
import math
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import cast

OUTPUTS = ("H.dat", "E.dat", "Fx.dat", "Fy.dat")
TIMESTEPS_RE = re.compile(r"Timesteps\s*:\s*(\d+)")
DURATION_RE = re.compile(r"Duration\s*:\s*([0-9.Ee+-]+)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_provenance(repo: Path) -> dict[str, object]:
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repo, text=True
    ).strip()
    status = subprocess.check_output(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        cwd=repo,
        text=True,
    ).splitlines()
    unstaged = subprocess.check_output(
        ["git", "diff", "--binary", "HEAD"], cwd=repo
    )
    return {
        "repository": str(repo),
        "commit": commit,
        "dirty": bool(status),
        "status": status,
        "working_diff_sha256": sha256_bytes(unstaged),
    }


def discover_git_root(path: Path) -> Path:
    return Path(
        subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "--show-toplevel"], text=True
        ).strip()
    ).resolve()


def make_variables(path: Path, names: tuple[str, ...]) -> dict[str, str]:
    if not path.is_file():
        return {}
    wanted = set(names)
    values: dict[str, str] = {}
    for line in path.read_text(errors="replace").splitlines():
        match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$", line)
        if match and match.group(1) in wanted:
            values[match.group(1)] = match.group(2)
    return values


def compiler_provenance(source_repo: Path) -> dict[str, object]:
    variables = make_variables(
        source_repo / "src" / "xbeach" / "Makefile",
        ("FC", "FCFLAGS", "AM_FCFLAGS", "LDFLAGS", "LIBS"),
    )
    command = variables.get("FC", "gfortran").split()[0]
    compiler = shutil.which(command)
    if compiler is None:
        return {"configured_variables": variables, "path": None}
    compiler_path = Path(compiler).resolve()
    version = subprocess.check_output(
        [str(compiler_path), "--version"], text=True, errors="replace"
    ).splitlines()[0]
    return {
        "configured_variables": variables,
        "path": str(compiler_path),
        "sha256": sha256(compiler_path),
        "version": version,
    }


def runtime_executable(launcher: Path) -> Path:
    candidate = launcher.parent / ".libs" / launcher.name
    if launcher.read_bytes()[:4] == b"\x7fELF":
        return launcher
    if candidate.is_file() and candidate.read_bytes()[:4] == b"\x7fELF":
        return candidate.resolve()
    raise RuntimeError(f"cannot resolve an ELF XBeach executable from {launcher}")


def runtime_environment(executable: Path) -> dict[str, str]:
    env = os.environ.copy()
    library_dir = executable.parents[2] / "xbeachlibrary" / ".libs"
    for key in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
        env[key] = os.pathsep.join(
            item for item in (str(library_dir), env.get(key, "")) if item
        )
    return env


def dynamic_dependencies(
    executable: Path, env: dict[str, str]
) -> list[dict[str, object]]:
    if shutil.which("ldd") is None:
        return []
    output = subprocess.check_output(["ldd", str(executable)], text=True, env=env)
    paths: set[Path] = set()
    for line in output.splitlines():
        match = re.search(r"=>\s+(/\S+)", line)
        if match is None:
            match = re.match(r"\s*(/\S+)", line)
        if match:
            path = Path(match.group(1)).resolve()
            if path.is_file():
                paths.add(path)
    return [
        {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(paths)
    ]


def runtime_artifacts(executable: Path, env: dict[str, str]) -> dict[str, object]:
    dependencies = dynamic_dependencies(executable, env)
    loaded = [
        item for item in dependencies if Path(str(item["path"])).name.startswith("libxbeach.so")
    ]
    if len(loaded) != 1:
        raise RuntimeError(f"expected one loaded libxbeach shared object, found {len(loaded)}")
    return {
        "executable": str(executable),
        "executable_size": executable.stat().st_size,
        "executable_sha256": sha256(executable),
        "dynamic_dependencies": dependencies,
        "loaded_library": loaded[0],
    }


def cpu_model() -> str:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(errors="replace").splitlines():
            if line.lower().startswith("model name"):
                return line.partition(":")[2].strip()
    return platform.processor()


def load_average() -> list[float] | None:
    try:
        return list(os.getloadavg())
    except OSError:
        return None


def replace_scalar(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?im)^(\s*{re.escape(key)}\s*=\s*)[^\r\n]+")
    updated, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise ValueError(f"expected one {key}= entry, found {count}")
    return updated


def parameter_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(errors="replace").splitlines():
        match = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(\S+)", line)
        if match:
            values[match.group(1).lower()] = match.group(2)
    return values


def validate_effective_parameters(case_dir: Path, expected_tstop: float) -> dict[str, object]:
    values = parameter_values(case_dir / "params.txt")
    required = ("nx", "ny", "thetamin", "thetamax", "dtheta", "random", "tstart", "tstop", "outputformat")
    missing = [key for key in required if key not in values]
    if missing:
        raise RuntimeError(f"effective params.txt lacks: {', '.join(missing)}")
    nx, ny = int(values["nx"]), int(values["ny"])
    dtheta = float(values["dtheta"])
    directions = int(round((float(values["thetamax"]) - float(values["thetamin"])) / dtheta))
    observed = {
        "nx": nx,
        "ny": ny,
        "directions": directions,
        "random": int(values["random"]),
        "tstart": float(values["tstart"]),
        "tstop": float(values["tstop"]),
        "outputformat": values["outputformat"].lower(),
    }
    expected = {
        "nx": 177,
        "ny": 70,
        "directions": 9,
        "random": 0,
        "tstart": expected_tstop,
        "tstop": expected_tstop,
        "outputformat": "fortran",
    }
    if observed != expected:
        raise RuntimeError(f"effective workload differs: observed={observed}, expected={expected}")
    return observed


def validate_warnings(path: Path) -> list[str]:
    lines = [line.strip() for line in path.read_text(errors="replace").splitlines() if line.strip()]
    allowed = (
        re.compile(r"^Warning: Specification of instat using parameter 'wbctype'$"),
        re.compile(r"^(morstart|tintg|tintp|tintm)\s*=.*Warning: value > recommended value of"),
        re.compile(r"^Warning: Lax Wendroff .* scheme is not supported, changed$"),
        re.compile(r"^to Warming and Beam \[scheme=SCHEME_WARMBEAM\]$"),
        re.compile(r"^Unknown, unused or multiple statements of parameter[A-Z0-9_]+ inparams\.txt$"),
        re.compile(r"^Warning: shallow water so long wave variance is reduced using par%nmax$"),
    )
    unknown = [line for line in lines if not any(pattern.match(line) for pattern in allowed)]
    if unknown:
        raise RuntimeError(f"unclassified XBeach warnings: {unknown}")
    return lines


def prepare_case(source: Path, destination: Path, tstop: float | None) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)
    if tstop is not None:
        params = destination / "params.txt"
        text = params.read_text()
        formatted = f"{tstop:.12g}"
        text = replace_scalar(text, "tstop", formatted)
        text = replace_scalar(text, "tstart", formatted)
        params.write_text(text)


def parse_run_log(text: str) -> tuple[int, float]:
    steps = TIMESTEPS_RE.findall(text)
    durations = DURATION_RE.findall(text)
    if not steps or not durations:
        raise RuntimeError("XBeach log did not contain Duration and Timesteps")
    if "End of program xbeach" not in text:
        raise RuntimeError("XBeach log did not contain the normal termination marker")
    return int(steps[-1]), float(durations[-1])


def read_native_doubles(path: Path) -> array.array[float]:
    data = path.read_bytes()
    values = array.array("d")
    if len(data) % values.itemsize != 0:
        raise RuntimeError(f"{path.name} size is not a multiple of FP64 storage")
    values.frombytes(data)
    return values


def validate_outputs(case_dir: Path, expected_tstop: float) -> dict[str, object]:
    effective = validate_effective_parameters(case_dir, expected_tstop)
    dims_path = case_dir / "dims.dat"
    if not dims_path.is_file():
        raise RuntimeError("missing output metadata: dims.dat")
    dims = read_native_doubles(dims_path)
    if len(dims) != 11:
        raise RuntimeError(f"dims.dat contains {len(dims)} values instead of 11")
    nx, ny, directions = (int(dims[1]), int(dims[2]), int(dims[3]))
    if any(float(int(value)) != value for value in dims[:10]):
        raise RuntimeError("dims.dat integer metadata contains non-integral values")
    observed = (nx, ny, directions, dims[-1])
    expected = (177, 70, 9, expected_tstop)
    if observed != expected:
        raise RuntimeError(f"output metadata differs: observed={observed}, expected={expected}")
    expected_values = (nx + 1) * (ny + 1)
    for name in OUTPUTS:
        values = read_native_doubles(case_dir / name)
        if len(values) != expected_values:
            raise RuntimeError(
                f"{name} contains {len(values)} values instead of {expected_values}"
            )
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(f"{name} contains a non-finite value")
    warning_path = case_dir / "XBwarning.txt"
    if not warning_path.is_file():
        raise RuntimeError("missing XBwarning.txt")
    warnings = validate_warnings(warning_path)
    return {
        "effective_parameters": effective,
        "dims_sha256": sha256(dims_path),
        "dims_values": list(dims),
        "nx": nx,
        "ny": ny,
        "directions": directions,
        "output_time": dims[-1],
        "values_per_field": expected_values,
        "warning_sha256": sha256(warning_path),
        "warning_lines": warnings,
    }


def input_hashes(case_dir: Path) -> dict[str, str]:
    return {
        str(path.relative_to(case_dir)): sha256(path)
        for path in sorted(case_dir.iterdir())
        if path.is_file()
    }


def validate_frozen_oracle(
    case_source: Path, expected_tstop: float, runs: list[dict[str, object]]
) -> dict[str, object]:
    path = case_source / "oracle.json"
    oracle = cast(dict[str, object], json.loads(path.read_text()))
    workload = cast(dict[str, object], oracle["workload"])
    applied = expected_tstop == float(cast(float, workload["tstop"]))
    result = {
        "path": str(path),
        "sha256": sha256(path),
        "source_commit": oracle["source_commit"],
        "applied": applied,
    }
    if not applied:
        return result

    first = runs[0]
    validation = cast(dict[str, object], first["validation"])
    observed = {
        "timesteps": first["timesteps"],
        "values_per_field": validation["values_per_field"],
        "output_hashes": first["hashes"],
        "output_sizes": first["sizes"],
        "dims_sha256": validation["dims_sha256"],
    }
    expected = {
        "timesteps": workload["timesteps"],
        "values_per_field": workload["values_per_field"],
        "output_hashes": oracle["output_hashes"],
        "output_sizes": oracle["output_sizes"],
        "dims_sha256": oracle["dims_sha256"],
    }
    if observed != expected:
        raise RuntimeError(
            f"full DELILAH result differs from frozen oracle: "
            f"observed={observed}, expected={expected}"
        )
    return result


def run_once(
    executable: Path,
    case_dir: Path,
    env: dict[str, str],
    expected_tstop: float,
    expected_artifacts: dict[str, object],
) -> dict[str, object]:
    effective_inputs = input_hashes(case_dir)
    artifacts_before = runtime_artifacts(executable, env)
    if artifacts_before != expected_artifacts:
        raise RuntimeError("runtime artifacts changed before a sample")

    started_at = datetime.now(timezone.utc).isoformat()
    load_before = load_average()
    started = time.perf_counter()
    completed = subprocess.run(
        [str(executable)],
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    wall_seconds = time.perf_counter() - started
    (case_dir / "console.log").write_text(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError(
            f"XBeach exited with {completed.returncode}; see {case_dir / 'console.log'}"
        )

    log_text = completed.stdout
    xblog = case_dir / "XBlog.txt"
    if xblog.exists():
        log_text += "\n" + xblog.read_text(errors="replace")
    timesteps, reported_seconds = parse_run_log(log_text)

    missing = [name for name in OUTPUTS if not (case_dir / name).is_file()]
    if missing:
        raise RuntimeError(f"missing output files: {', '.join(missing)}")
    validation = validate_outputs(case_dir, expected_tstop)
    artifacts_after = runtime_artifacts(executable, env)
    if artifacts_after != artifacts_before:
        raise RuntimeError("runtime artifacts changed during a sample")

    return {
        "started_at": started_at,
        "effective_input_hashes": effective_inputs,
        "runtime_artifacts": artifacts_before,
        "load_average_before": load_before,
        "load_average_after": load_average(),
        "wall_seconds": wall_seconds,
        "reported_seconds": reported_seconds,
        "timesteps": timesteps,
        "validation": validation,
        "hashes": {name: sha256(case_dir / name) for name in OUTPUTS},
        "sizes": {name: (case_dir / name).stat().st_size for name in OUTPUTS},
    }


def ensure_consistent(runs: list[dict[str, object]]) -> None:
    reference = (
        runs[0]["timesteps"],
        runs[0]["hashes"],
        runs[0]["sizes"],
        runs[0]["validation"],
        runs[0]["effective_input_hashes"],
        runs[0]["runtime_artifacts"],
    )
    for index, run in enumerate(runs[1:], start=2):
        current = (
            run["timesteps"],
            run["hashes"],
            run["sizes"],
            run["validation"],
            run["effective_input_hashes"],
            run["runtime_artifacts"],
        )
        if current != reference:
            raise RuntimeError(f"retained run {index} differs from retained run 1")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    repo = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--executable",
        type=Path,
        default=repo / "src" / "xbeach" / ".libs" / "xbeach",
    )
    parser.add_argument("--runs", type=int, default=1)
    parser.add_argument("--warmups", type=int, default=0)
    parser.add_argument(
        "--work-root", type=Path, default=Path("/tmp/xbeach-delilah-benchmark")
    )
    parser.add_argument(
        "--tstop",
        type=float,
        help="Override tstop and move the single global output to that time (smoke runs)",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON result to this path")
    args = parser.parse_args()

    if args.runs < 1 or args.warmups < 0:
        parser.error("--runs must be at least 1 and --warmups cannot be negative")

    launcher = args.executable.resolve()
    if not launcher.is_file():
        parser.error(f"executable does not exist: {launcher}")
    executable = runtime_executable(launcher)
    source_repo = discover_git_root(executable.parent)
    case_source = repo / "case-delilah"
    if not (case_source / "params.txt").is_file():
        parser.error(f"benchmark input is missing: {case_source / 'params.txt'}")

    expected_tstop = args.tstop if args.tstop is not None else 3800.0
    env = runtime_environment(executable)
    artifacts = runtime_artifacts(executable, env)
    provenance = git_provenance(source_repo)
    runner_provenance = git_provenance(repo)
    runner_hash = sha256(Path(__file__).resolve())
    compiler = compiler_provenance(source_repo)
    case_hashes = input_hashes(case_source)

    work_root = args.work_root.resolve()
    if work_root.exists():
        parser.error(f"work root already exists; choose a fresh path: {work_root}")
    work_root.mkdir(parents=True)
    retained: list[dict[str, object]] = []
    total = args.warmups + args.runs
    for ordinal in range(total):
        label = "warmup" if ordinal < args.warmups else "run"
        number = ordinal + 1 if label == "warmup" else ordinal - args.warmups + 1
        run_dir = work_root / f"{label}-{number}"
        prepare_case(case_source, run_dir, args.tstop)
        result = run_once(executable, run_dir, env, expected_tstop, artifacts)
        print(
            f"{label} {number}: {result['wall_seconds']:.6f} s, "
            f"{result['timesteps']} timesteps",
            flush=True,
        )
        if label == "run":
            retained.append(result)

    ensure_consistent(retained)
    oracle = validate_frozen_oracle(case_source, expected_tstop, retained)
    if runtime_artifacts(executable, env) != artifacts:
        raise RuntimeError("runtime artifacts changed across the benchmark")
    if git_provenance(source_repo) != provenance:
        raise RuntimeError("source tree changed across the benchmark")
    if git_provenance(repo) != runner_provenance:
        raise RuntimeError("runner tree changed across the benchmark")
    if sha256(Path(__file__).resolve()) != runner_hash:
        raise RuntimeError("benchmark runner changed across the benchmark")
    if compiler_provenance(source_repo) != compiler:
        raise RuntimeError("compiler identity changed across the benchmark")
    if input_hashes(case_source) != case_hashes:
        raise RuntimeError("source case changed across the benchmark")
    walls = [cast(float, run["wall_seconds"]) for run in retained]
    relevant_environment = {
        key: value
        for key, value in sorted(os.environ.items())
        if key.startswith(("OMP_", "GOMP_", "KMP_"))
        or key in {"LANG", "LC_ALL", "LD_PRELOAD"}
    }
    output_validation = cast(dict[str, object], retained[0]["validation"])
    report = {
        "schema": 4,
        "measured_at": datetime.now(timezone.utc).isoformat(),
        "source_commit": provenance["commit"],
        "source_tree": provenance,
        "runner_tree": runner_provenance,
        "launcher": str(launcher),
        "runtime_artifacts": artifacts,
        "runner_sha256": runner_hash,
        "command": [str(Path(__file__).resolve()), *sys.argv[1:]],
        "compiler": compiler,
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": cpu_model(),
            "cpu_count": os.cpu_count(),
            "cpu_affinity": sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else None,
        },
        "environment": relevant_environment,
        "case_input_hashes": case_hashes,
        "effective_input_hashes": retained[0]["effective_input_hashes"],
        "frozen_oracle": oracle,
        "workload": {
            "case": "DELILAH",
            "nx": output_validation["nx"],
            "ny": output_validation["ny"],
            "directions": output_validation["directions"],
            "random": 0,
            "tstop": expected_tstop,
            "output_time": output_validation["output_time"],
        },
        "warmups": args.warmups,
        "runs": retained,
        "median_wall_seconds": statistics.median(walls),
        "minimum_wall_seconds": min(walls),
        "maximum_wall_seconds": max(walls),
        "output_hashes": retained[0]["hashes"],
        "output_sizes": retained[0]["sizes"],
        "timesteps": retained[0]["timesteps"],
    }
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"benchmark failed: {error}", file=sys.stderr)
        raise SystemExit(1)
