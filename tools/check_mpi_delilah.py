#!/usr/bin/env python3
"""Compare serial and MPI DELILAH runs across explicit processor grids.

XBeach reserves one MPI rank for output, so an MxN compute grid is launched
with M*N+1 ranks.  The command fails if a run crashes, takes a different
number of timesteps, produces non-finite data, or exceeds the declared
fieldwise physical tolerances.
"""

from __future__ import annotations

import argparse
import array
import hashlib
import json
import math
import platform
import re
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from benchmark_delilah import OUTPUTS, parse_run_log, prepare_case, replace_scalar, sha256

DISTRIBUTION_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$"
)
DEFAULT_TOLERANCES = {
    "H.dat": (1.0e-8, 1.0e-7),
    "E.dat": (1.0e-6, 1.0e-7),
    "Fx.dat": (1.0e-6, 1.0e-7),
    "Fy.dat": (1.0e-6, 1.0e-7),
}


@dataclass(frozen=True)
class Layout:
    m: int
    n: int

    @property
    def label(self) -> str:
        return f"{self.m}x{self.n}"

    @property
    def ranks(self) -> int:
        return self.m * self.n + 1


@dataclass
class FieldComparison:
    field: str
    values: int
    max_abs: float
    max_rel: float
    max_abs_index_zero_based: int
    max_abs_ij_one_based: tuple[int, int]
    violating_values: int
    seam_max_abs: float | None
    passed: bool


def parse_layout(text: str) -> Layout:
    match = re.fullmatch(r"([1-9]\d*)x([1-9]\d*)", text.strip().lower())
    if match is None:
        raise argparse.ArgumentTypeError(f"invalid processor grid: {text!r}")
    return Layout(int(match.group(1)), int(match.group(2)))


def add_mpi_layout(params_path: Path, layout: Layout) -> None:
    text = params_path.read_text()
    for key, value in (
        ("mpiboundary", "man"),
        ("mmpi", str(layout.m)),
        ("nmpi", str(layout.n)),
    ):
        pattern = re.compile(rf"(?im)^\s*{key}\s*=")
        if pattern.search(text):
            text = replace_scalar(text, key, value)
        else:
            if not text.endswith("\n"):
                text += "\n"
            text += f"{key} = {value}\n"
    params_path.write_text(text)


def read_doubles(path: Path) -> array.array[float]:
    values = array.array("d")
    raw = path.read_bytes()
    if len(raw) % values.itemsize:
        raise RuntimeError(f"{path} is not native FP64 storage")
    values.frombytes(raw)
    if not values:
        raise RuntimeError(f"{path} is empty")
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError(f"{path} contains non-finite values")
    return values


def read_dims(path: Path) -> tuple[int, int]:
    values = read_doubles(path)
    if len(values) < 3:
        raise RuntimeError(f"{path} lacks nx/ny metadata")
    nx, ny = values[1:3]
    if nx < 1 or ny < 0 or nx != int(nx) or ny != int(ny):
        raise RuntimeError(f"{path} has invalid nx/ny metadata: {nx}, {ny}")
    return int(nx) + 1, int(ny) + 1


def source_identity(repo_root: Path) -> dict[str, object]:
    def git(*args: str) -> str:
        completed = subprocess.run(
            ["git", *args],
            cwd=repo_root,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return completed.stdout.strip()

    status = git("status", "--porcelain=v1")
    diff = subprocess.run(
        ["git", "diff", "--binary", "HEAD"],
        cwd=repo_root,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    return {
        "repo_root": str(repo_root),
        "commit": git("rev-parse", "HEAD"),
        "dirty": bool(status),
        "status": status.splitlines(),
        "working_diff_sha256": hashlib.sha256(diff).hexdigest(),
    }


def tree_identity(root: Path) -> dict[str, object]:
    files: dict[str, dict[str, object]] = {}
    digest = hashlib.sha256()
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        relative = path.relative_to(root).as_posix()
        payload = path.read_bytes()
        file_hash = hashlib.sha256(payload).hexdigest()
        files[relative] = {"bytes": len(payload), "sha256": file_hash}
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(file_hash))
    return {"tree_sha256": digest.hexdigest(), "files": files}


def runtime_identity(mpirun: str) -> dict[str, str]:
    completed = subprocess.run(
        [mpirun, "--version"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return {
        "platform": platform.platform(),
        "python": sys.version,
        "mpirun": completed.stdout.strip(),
    }


def parse_subdomains(log: str) -> list[tuple[int, int, int, int, int]]:
    rows: list[tuple[int, int, int, int, int]] = []
    in_table = False
    for line in log.splitlines():
        if "proc   is   lm   js   ln" in line:
            in_table = True
            continue
        if not in_table:
            continue
        match = DISTRIBUTION_ROW_RE.match(line)
        if match is not None:
            rank, is_, lm, js, ln = map(int, match.groups())
            rows.append((rank, is_, lm, js, ln))
            continue
        if rows:
            break
    if not rows:
        raise RuntimeError("MPI log lacks the processor distribution table")
    return rows


def validate_subdomains(
    rows: list[tuple[int, int, int, int, int]], layout: Layout
) -> None:
    expected_ranks = layout.m * layout.n
    ranks = sorted(row[0] for row in rows)
    if ranks != list(range(expected_ranks)):
        raise RuntimeError(
            f"processor table ranks {ranks} do not match 0..{expected_ranks - 1}"
        )
    actual_m = len({row[1] for row in rows})
    actual_n = len({row[3] for row in rows})
    if (actual_m, actual_n) != (layout.m, layout.n):
        raise RuntimeError(
            f"requested {layout.label}, log reports {actual_m}x{actual_n}"
        )


def internal_seams(
    rows: list[tuple[int, int, int, int, int]], width: int, height: int
) -> tuple[set[int], set[int]]:
    # `is` and `js` are one-based starts including two overlap cells. Centre
    # the diagnostic band on the first independently owned cell.
    m_starts = {row[1] + 2 for row in rows if row[1] > 1}
    n_starts = {row[3] + 2 for row in rows if row[3] > 1}
    m_band = {
        index
        for start in m_starts
        for index in range(max(1, start - 2), min(width, start + 2) + 1)
    }
    n_band = {
        index
        for start in n_starts
        for index in range(max(1, start - 2), min(height, start + 2) + 1)
    }
    return m_band, n_band


def compare_field(
    serial: array.array[float],
    candidate: array.array[float],
    field: str,
    width: int,
    atol: float,
    rtol: float,
    m_band: set[int],
    n_band: set[int],
) -> FieldComparison:
    if len(serial) != len(candidate):
        raise RuntimeError(
            f"{field} value-count mismatch: serial={len(serial)}, MPI={len(candidate)}"
        )
    max_abs = -1.0
    max_rel = -1.0
    max_index = 0
    violating = 0
    seam_max: float | None = None
    for index, (reference, observed) in enumerate(zip(serial, candidate)):
        abs_error = abs(observed - reference)
        rel_error = abs_error / max(abs(reference), abs(observed), 1.0e-300)
        if abs_error > max_abs:
            max_abs = abs_error
            max_index = index
        max_rel = max(max_rel, rel_error)
        if abs_error > atol and rel_error > rtol:
            violating += 1
        i = index % width + 1
        j = index // width + 1
        if i in m_band or j in n_band:
            seam_max = abs_error if seam_max is None else max(seam_max, abs_error)
    return FieldComparison(
        field=field,
        values=len(serial),
        max_abs=max_abs,
        max_rel=max_rel,
        max_abs_index_zero_based=max_index,
        max_abs_ij_one_based=(max_index % width + 1, max_index // width + 1),
        violating_values=violating,
        seam_max_abs=seam_max,
        passed=violating == 0,
    )


def run_case(command: list[str], case_dir: Path) -> dict[str, object]:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    elapsed = time.perf_counter() - started
    (case_dir / "run.log").write_text(completed.stdout)
    (case_dir / "run.err").write_text(completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(
            f"run failed with exit {completed.returncode}: {completed.stderr[-2000:]}"
        )
    timesteps, reported_duration = parse_run_log(completed.stdout)
    missing = [name for name in OUTPUTS if not (case_dir / name).is_file()]
    if missing:
        raise RuntimeError(f"run did not produce: {', '.join(missing)}")
    return {
        "command": command,
        "returncode": completed.returncode,
        "wall_seconds": elapsed,
        "reported_duration": reported_duration,
        "timesteps": timesteps,
        "stdout_sha256": sha256(case_dir / "run.log"),
        "stderr_sha256": sha256(case_dir / "run.err"),
        "output_hashes": {name: sha256(case_dir / name) for name in OUTPUTS},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial-exe", type=Path, required=True)
    parser.add_argument("--mpi-exe", type=Path, required=True)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tstop", type=float, default=20.0)
    parser.add_argument(
        "--layouts",
        type=parse_layout,
        nargs="+",
        default=[Layout(1, 1), Layout(1, 7), Layout(7, 1)],
    )
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--oversubscribe", action="store_true")
    args = parser.parse_args()

    serial_exe = args.serial_exe.resolve()
    mpi_exe = args.mpi_exe.resolve()
    case_source = args.case_dir.resolve()
    output_dir = args.output_dir.resolve()
    for path in (serial_exe, mpi_exe, case_source / "params.txt"):
        if not path.exists():
            parser.error(f"missing required path: {path}")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    serial_dir = output_dir / "serial"
    prepare_case(case_source, serial_dir, args.tstop)
    serial_result = run_case([str(serial_exe)], serial_dir)
    serial_values = {name: read_doubles(serial_dir / name) for name in OUTPUTS}
    width, height = read_dims(serial_dir / "dims.dat")
    serial_dims = read_doubles(serial_dir / "dims.dat")
    expected_values = width * height
    if any(len(values) != expected_values for values in serial_values.values()):
        raise RuntimeError("DELILAH outputs do not have the expected 178x71 shape")

    report: dict[str, object] = {
        "schema": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "case_source": str(case_source),
        "case_inputs": tree_identity(case_source),
        "tstop": args.tstop,
        "serial_executable": {
            "path": str(serial_exe),
            "sha256": sha256(serial_exe),
        },
        "mpi_executable": {"path": str(mpi_exe), "sha256": sha256(mpi_exe)},
        "source": source_identity(Path(__file__).resolve().parents[1]),
        "runtime": runtime_identity(args.mpirun),
        "tolerances": {
            name: {"atol": values[0], "rtol": values[1]}
            for name, values in DEFAULT_TOLERANCES.items()
        },
        "serial": serial_result,
        "layouts": [],
    }
    all_passed = True
    layout_reports: list[dict[str, object]] = []
    for layout in args.layouts:
        run_dir = output_dir / f"mpi-{layout.label}"
        prepare_case(case_source, run_dir, args.tstop)
        add_mpi_layout(run_dir / "params.txt", layout)
        command = [args.mpirun]
        if args.oversubscribe:
            command.append("--oversubscribe")
        command.extend(["-np", str(layout.ranks), str(mpi_exe)])
        try:
            result = run_case(command, run_dir)
            rows = parse_subdomains((run_dir / "run.log").read_text(errors="replace"))
            validate_subdomains(rows, layout)
            if read_doubles(run_dir / "dims.dat") != serial_dims:
                raise RuntimeError("MPI dims.dat differs from serial")
            m_band, n_band = internal_seams(rows, width, height)
            comparisons = []
            for name in OUTPUTS:
                atol, rtol = DEFAULT_TOLERANCES[name]
                comparison = compare_field(
                    serial_values[name],
                    read_doubles(run_dir / name),
                    name,
                    width,
                    atol,
                    rtol,
                    m_band,
                    n_band,
                )
                comparisons.append(asdict(comparison))
            same_steps = result["timesteps"] == serial_result["timesteps"]
            passed = same_steps and all(item["passed"] for item in comparisons)
            layout_report = {
                "layout": layout.label,
                "compute_ranks": layout.m * layout.n,
                "launch_ranks": layout.ranks,
                "subdomains": rows,
                "same_timesteps": same_steps,
                "comparisons": comparisons,
                "run": result,
                "passed": passed,
            }
        except Exception as exc:
            layout_report = {
                "layout": layout.label,
                "compute_ranks": layout.m * layout.n,
                "launch_ranks": layout.ranks,
                "passed": False,
                "error": str(exc),
            }
        layout_reports.append(layout_report)
        all_passed = all_passed and bool(layout_report["passed"])
        state = "PASS" if layout_report["passed"] else "FAIL"
        print(f"{layout.label}: {state}")
        for item in layout_report.get("comparisons", []):
            print(
                f"  {item['field']}: maxabs={item['max_abs']:.6e} "
                f"maxrel={item['max_rel']:.6e} violations={item['violating_values']} "
                f"at={tuple(item['max_abs_ij_one_based'])}"
            )
        if "error" in layout_report:
            print(f"  {layout_report['error']}")
    report["layouts"] = layout_reports
    report["passed"] = all_passed
    report_path = output_dir / "mpi-correctness.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"report: {report_path}")
    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
