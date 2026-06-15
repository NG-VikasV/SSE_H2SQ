from pathlib import Path
import subprocess
import os
from typing import Dict
import pandas as pd

# --------------------------------------------------
# BASE DIRECTORY
# --------------------------------------------------
# Set project root explicitly (one level above 'files')
project_root = Path(__file__).resolve().parent.parent.parent  # adjust if qmc_runner.py is deeper
base_dir = project_root / "files"  # points to .../IT_QMC_H2SQ_Nov_2026_JAN/files

# --------------------------------------------------
# PATH BUILDER
# --------------------------------------------------
def build_run_path(params: Dict) -> Path:
    """
    Build directory path like:
    betaVp_100/fangleVp_0/L4/J2_0.00/J3_0.00/M1
    """
    required_keys = ["betaVp", "fangleVp", "L", "J2", "J3", "M"]
    missing = [k for k in required_keys if k not in params]
    if missing:
        raise ValueError(f"Missing parameters: {missing}")

    path = base_dir / f"betaVp_{params['betaVp']}" \
                   / f"fangleVp_{params['fangleVp']}" \
                   / f"L{params['L']}" \
                   / f"J2_{params['J2']}" \
                   / f"J3_{params['J3']}" \
                   / f"M{params['M']}"
    return path.resolve()  # absolute path

# --------------------------------------------------
# INTERNAL EXECUTOR
# --------------------------------------------------
def _run_command(cmd, cwd: Path):
    print("\n" + "=" * 80)
    print(f"📂 Directory : {cwd}")
    print("▶ Command   :", " ".join(cmd))
    print("=" * 80)

    process = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    for line in process.stdout:
        print(line, end="")

    process.wait()

    if process.returncode != 0:
        raise RuntimeError(f"Command failed with return code {process.returncode}")

# --------------------------------------------------
# MAIN FUNCTION
# --------------------------------------------------
def run_job(params: Dict, mode: str = "all", executable: str = "main"):
    run_dir = build_run_path(params)

    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory does not exist:\n{run_dir}")

    if mode not in {"all", "build", "run", "clean"}:
        raise ValueError("mode must be one of: all, build, run, clean")

    exe_name = f"{executable}.exe" if os.name == "nt" else f"./{executable}"

    if mode in ("all", "clean"):
        _run_command(["make", "clean"], cwd=run_dir)

    if mode in ("all", "build"):
        _run_command(["make"], cwd=run_dir)

    if mode in ("all", "run"):
        exe_path = run_dir / exe_name
        if not exe_path.exists():
            raise FileNotFoundError(f"Executable not found: {exe_name}\nDid you compile?")
        _run_command([str(exe_path)], cwd=run_dir)


def show_data(params: Dict):
    dir = build_run_path(params)

    if not dir.is_dir():
        raise FileNotFoundError(f"The directory does not exist:\n{dir}")


    file_str = "IT_data_avg.txt"

    file = dir / file_str

    if not file.exists():
        raise FileNotFoundError(f"The data file does not exist: \n{file}")


    return pd.read_csv(file, sep=r'\s+')


# --------------------------------------------------
# OPTIONAL: QUICK CHECK
# --------------------------------------------------
def print_path(params: Dict):
    """Utility to just print resolved directory path"""
    print(build_run_path(params))
