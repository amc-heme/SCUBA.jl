#!/usr/bin/env python
"""Run a command with a timeout (macOS has no `timeout`), passing --start-time through to harnesses.

  python -I tools/run_with_timeout.py --timeout 7200 -- julia --project=. julia/run.jl --scale 1M
"""
import subprocess
import sys
import time


def main():
    args = sys.argv[1:]
    timeout = None
    if args and args[0] == "--timeout":
        timeout = float(args[1]); args = args[2:]
    if args and args[0] == "--":
        args = args[1:]
    cmd = list(args) + ["--start-time", f"{time.time():.3f}"]
    try:
        r = subprocess.run(cmd, timeout=timeout)
        sys.exit(r.returncode)
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT after {timeout}s: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(124)


if __name__ == "__main__":
    main()
