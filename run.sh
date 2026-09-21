#!/bin/sh
# Runs the Demo. The library has no service; its door is the worked example.
set -eu

cd -- "$(dirname -- "$0")"

PYTHON="${PYTHON:-python3}"

verify_interpreter() {
	command -v "$PYTHON" >/dev/null 2>&1 || {
		printf '\342\235\214 wanted a Python interpreter named %s.\n' "$PYTHON" >&2
		printf '   PYTHON=/path/to/python ./run.sh\n' >&2
		exit 1
	}
}

run_demo() {
	"$PYTHON" dissociation.py
}

verify_interpreter
run_demo
