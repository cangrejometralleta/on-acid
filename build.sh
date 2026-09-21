#!/bin/sh
# Builds the Package. Format, then types, then tests, then the artefact.
set -eu

cd -- "$(dirname -- "$0")"

PYTHON="${PYTHON:-python3}"

refuse() {
	printf '\342\235\214 %s\n' "$1" >&2
	printf '   %s\n' "$2" >&2
	exit 1
}

skip() {
	printf '\342\200\224 %s not installed, skipped: %s\n' "$1" "$2"
}

verify_interpreter() {
	command -v "$PYTHON" >/dev/null 2>&1 ||
		refuse "wanted a Python interpreter named $PYTHON." \
			"PYTHON=/path/to/python ./build.sh"
}

verify_format() {
	command -v ruff >/dev/null 2>&1 || {
		skip ruff "pip install ruff"
		return 0
	}
	ruff format --check . ||
		refuse "wanted the formatting to be settled." "ruff format ."
}

verify_types() {
	command -v mypy >/dev/null 2>&1 || {
		skip mypy "pip install mypy"
		return 0
	}
	mypy dissociation.py ||
		refuse "wanted the annotations to hold." "mypy dissociation.py"
}

verify_tests() {
	"$PYTHON" -m unittest discover -p 'test_*.py' ||
		refuse "wanted every test to pass." "$PYTHON -m unittest discover -p 'test_*.py'"
}

verify_artefact() {
	"$PYTHON" -c 'import build' >/dev/null 2>&1 || {
		skip build "pip install build"
		return 0
	}
	"$PYTHON" -m build ||
		refuse "wanted a wheel and an sdist." "$PYTHON -m build"
}

verify_interpreter
verify_format
verify_types
verify_tests
verify_artefact
